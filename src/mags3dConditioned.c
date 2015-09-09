#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mags3d.h>

int main(int argc,char *argv[]){

	char *line;
	int BUFSIZE=100;
	line=(char *)malloc(sizeof(char)*BUFSIZE + 1);

	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int nz = atoi(argv[3]);
	TYPE xlo = (TYPE)atof(argv[4]);
	TYPE ylo = (TYPE)atof(argv[5]);
	TYPE zlo = (TYPE)atof(argv[6]);
	TYPE h = (TYPE)atof(argv[7]);
	TYPE a = (TYPE)atof(argv[8]);

	int mode = atoi(argv[9]);

	double mu = 0.0;
	if(mode==4){
		mu=atof(argv[10]);
	}

	int bufferNodes=ceil(a/h);
	int neighRadius=bufferNodes;
	int neighSide=2*neighRadius+1;
	int totalNeighbours=(2*neighRadius+1);  
	float totalNeighboursinv=1.0f/((float)(2*neighRadius+1));  

	int nxExt=nx+2*bufferNodes;
	int nyExt=ny+2*bufferNodes;
	int nzExt=nz+2*bufferNodes;


	FILE *fpoctave,*fpgslib,*fpweight,*fpdistanceweight,*fpharddata;

	//fpoctave = fopen("octaveout.m","w");
	if(mode==0 || mode==3 || mode==4){
		fpgslib  = fopen("gslibout.dat","w");
		fpoctave  = fopen("octaveout.m","w");
	}
	if(mode==1){
		fpweight = fopen("weightout.dat","w");
	}
	if(mode==2){
		fpdistanceweight = fopen("targetdistanceweights.dat","w");
	}
	if(mode==3 || mode==4){
		fpdistanceweight = fopen("currentbestdistanceweights.dat","r");
	}
	if(mode==4){
	      	if( (fpharddata = fopen("harddata.dat","r"))==NULL ){
			printf("ERROR: file harddata.dat does not exist. goodbye!\n");
			exit(1);
		}
	}

	//fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
	if(mode==0 || mode==3 || mode==4) fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");

	latticeParams params;
	params.nx=nx;
	params.ny=ny;
	params.nz=nz;
	params.xlo=xlo;
	params.ylo=ylo;
	params.zlo=zlo;
	params.h=h;

	TYPE *image;
	TYPE *imageNext;
	//image = genRandomImage(params.nx+2*bufferNodes,params.ny+2*bufferNodes,params.nz+2*bufferNodes);
	image = genRandomImage(nxExt,nyExt,nzExt);

	//printf("allocated %d bytes, plus %d bytes of buffer zone\n",params.nx*params.ny*params.nz*sizeof(TYPE),0);

	//int initialAddress = 	bufferNodes*((ny+2*bufferNodes)*(nx+2*bufferNodes)) + 
	//			bufferNodes*(nx+2*bufferNodes) + 
	//			bufferNodes; 
	int initialAddress = 	bufferNodes*nyExt*nxExt + 
				bufferNodes*nxExt + 
				bufferNodes; 

	//printf("initialAddress=%d\n",initialAddress);

	TYPE *imageZero = image+initialAddress;
	TYPE *imageZeroNext;

	int i,j,k,ni,nj,nk;
	TYPE sum=0.0;
	TYPE xcoord1,ycoord1,zcoord1;
	TYPE xcoord2,ycoord2,zcoord2;
	TYPE weight;
	TYPE r, rsquared, asquared, asquaredinv, cons;

	asquared = a*a;
	asquaredinv=1.0/asquared;
	//cons = pow((4.0/((TYPE)(neighRadius*neighRadius)*PI)),0.75);
	cons = pow((4.0/(asquared*PI)),0.75);

	TYPE *kernelWeight = (TYPE *)malloc(sizeof(TYPE)*neighSide*neighSide*neighSide);
	
	TYPE *distanceValue = (TYPE *)malloc(sizeof(TYPE)*MAXDISTVALS); // 101 radius lags 
	TYPE *distanceWeight = (TYPE *)malloc(sizeof(TYPE)*MAXDISTVALS); // 101 radius lags 

	int distanceCheck[MAXDISTVALS];
	int distance;
	//int totalcells=(2*neighRadius+1); 
	//int totalcells=50; 
	//int totalcells=5; 
	int totalcells=10; 
	for(i=0;i<MAXDISTVALS;i++)distanceCheck[i]=0;

	if(mode==1) fprintf(fpweight,"%d\n",neighSide);
	getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);

	int distCounter=0;
	if(mode==3 || mode==4){
		TYPE dwval,dwwei;
		for(i=0;i<MAXDISTVALS;i++){
			distanceValue[i]=-1.0;
			distanceWeight[i]  =-1.0;
		}
		i=0;
		while((fgets(line,BUFSIZE,fpdistanceweight))!=0){
			sscanf(line,"%lf %lf",&dwval,&dwwei);
			//printf("%f %f\n",dwval,dwwei);
			if(existValue(distanceValue,MAXDISTVALS,dwval)==-1){
				distanceValue[distCounter]=dwval;
				distanceWeight[distCounter]=dwwei;
				//printf("distanceValue[%d]=%f\n",distCounter,dwval);
				//printf("distanceWeight[%d]=%f\n",distCounter,dwwei);
				distCounter++;
			}
			//printf("%d %d %d %f\n",i,j,k,val);
		}
	}

	//printf("weight=zeros(%d,%d,%d);\n",neighSide,neighSide,neighSide);
	for(nk=-neighRadius;nk<=neighRadius;nk++){
	   for(nj=-neighRadius;nj<=neighRadius;nj++){
	      for(ni=-neighRadius;ni<=neighRadius;ni++){
	         getCoordinates(ni,nj,nk,params,&xcoord2,&ycoord2,&zcoord2);
	         rsquared = 	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	  			(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	  			(zcoord1-zcoord2)*(zcoord1-zcoord2);
		 r = sqrt(rsquared);

		//printf("r=%f, rsquared=%f\n",r,rsquared);
		if(mode==3 || mode==4){
			int p;
			for(p=0;p<distCounter;p++){
		    		if(abs(r-distanceValue[p])<0.00001){
					//printf("%f\t%f\tdistanceValue[%d]=%f\tdistanceWeight[%d]=%f\n",rsquared,r,p,distanceValue[p],p,distanceWeight[p]);
	         			*(kernelWeight+
					(ni+neighRadius)+
					(nj+neighRadius)*neighSide+
					(nk+neighRadius)*neighSide*neighSide)=distanceWeight[p]; 
   				}
		 	}
		}
		else{
	         *(kernelWeight+
				(ni+neighRadius)+
				(nj+neighRadius)*neighSide+
				(nk+neighRadius)*neighSide*neighSide) 
				//= cons * exp(-2.0*rsquared*asquaredinv); 	// gaussian 3D
				//= r<=a*0.5? cons : 0.0 ; 			// circular 2D
				//= (1.0/(r+2.0)) * 1.0/sqrt(2.0*PI*a)*exp(-(r+2.0)/a);	// exponetial 3D
				//= (1.0/(r+20.0)) * 1.0/sqrt(2.0*PI*a)*exp(-(r+20.0)/a);	// exponetial 3D
				= 1.0/((TYPE)totalcells); // constant
		}

                 //printf("weight(%d,%d,%d)=%f;\n",(ni+neighRadius)+1,
		 //		(nj+neighRadius)+1,
		 //		(nk+neighRadius)+1,
                 //              //cons * exp(-2.0*rsquared*asquaredinv) 		// gaussian 3D
		 //              //r<=a*0.5? cons : 0.0 				// circular 2D
		 //              //(1.0/(r+2.0)) * (1.0/sqrt(2.0*PI*a)) *exp(-(r+2.0)/a)	// exponetial 3D
		 //              //1.0/((TYPE)totalcells)
		 //       	*(kernelWeight+(ni+neighRadius)+(nj+neighRadius)*neighSide+(nk+neighRadius)*neighSide*neighSide)
                 //);

	         if(ni<=0 && nj<=0 && nk<=0){
		    if(mode==1){
	               //fprintf(fpweight,"%d %d %d %f\n",ni+neighRadius,nj+neighRadius,nk+neighRadius,(cons * exp(-2.0*rsquared*asquaredinv))); // gaussian 3D
	               //fprintf(fpweight,"%d %d %d %f\n",ni+neighRadius,nj+neighRadius,nk+neighRadius,(cons * exp(-0.5*rsquared*asquaredinv))); // gaussian 3D wrong
	               //fprintf(fpweight,"%d %d %d %f\n",ni+neighRadius,nj+neighRadius,nk+neighRadius, r<=a*0.5? cons : 0.0 ); // circular 2D
	               //fprintf(fpweight,"%d %d %d %f\n",ni+neighRadius,nj+neighRadius,nk+neighRadius,(1.0/(r+2.0))*(1.0/sqrt(2.0*PI*a))*exp(-(r+2.0)/a));//exponential3D
	               //fprintf(fpweight,"%d %d %d %f\n",ni+neighRadius,nj+neighRadius,nk+neighRadius,1.0/((TYPE)totalcells));// constant
	               fprintf(fpweight,"%d %d %d %f\n",ni+neighRadius,nj+neighRadius,nk+neighRadius,
				*(kernelWeight+(ni+neighRadius)+(nj+neighRadius)*neighSide+(nk+neighRadius)*neighSide*neighSide));
                    }
  		    if(mode==2){
		    	fprintf(fpdistanceweight,"%f %f\n",sqrt(rsquared),
							*(kernelWeight+
	  						(ni+neighRadius)+
	  						(nj+neighRadius)*neighSide+
	  						(nk+neighRadius)*neighSide*neighSide));
		    }
	         }

	      }
	   }
	}


	if(mode!=4){

		for(k=0;k<nz;k++){
		  for(j=0;j<ny;j++){
		    for(i=0;i<nx;i++){
        	      sum=0.0;
        	      getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
		      for(nk=-neighRadius;nk<=neighRadius;nk++){
		        for(nj=-neighRadius;nj<=neighRadius;nj++){
		          for(ni=-neighRadius;ni<=neighRadius;ni++){
        	            weight= *(kernelWeight+
					(ni+neighRadius)+
					(nj+neighRadius)*neighSide+
					(nk+neighRadius)*neighSide*neighSide);
 
		            sum=sum + weight * 
        	            		*(imageZero+(
						(i+ni)+
						(j+nj)*(nxExt)+
						(k+nk)*(nxExt)*(nyExt)
						)
					);
		          }
		        }
		      }

        	      //fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
			//sum = sum/((float)totalNeighbours); 	
			sum = sum*totalNeighboursinv; 	
        	      if(mode==0 || mode==3){ 
				fprintf(fpgslib,"%f %f %f %f\n",xcoord1,ycoord1,zcoord1,sum);
        	      		fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
			}

		    }
		  }
		}

	}
	else{ //conditioning (local gradual deformation)
		int BUFSIZE=100;
		char *line=(char *)malloc(sizeof(char)*BUFSIZE + 1);
		char *line2=(char *)malloc(sizeof(char)*BUFSIZE + 1);
		

		int flag=0;
		int numcols=0;
		int maxdat=0;
		while((fgets(line,BUFSIZE,fpharddata))!=0){
			if(flag==1){
				sscanf(line,"%d",&numcols);
			}
			if(flag>=numcols+2){
				//sscanf(line2,"%f[ \t]+%f[ \t]+%f",&numcols);
				//printf("datafl: %s",line2);
				maxdat=maxdat+1;
			}
			flag++;
		}	

		printf("numcols=%d maxdat=%d\n",numcols,maxdat);

		TYPE *dataCoordinate = (TYPE *)malloc(sizeof(TYPE)*maxdat*3); // 101 radius lags 
		TYPE *dataValue = (TYPE *)malloc(sizeof(TYPE)*maxdat); // 101 radius lags 

		float *dataPair=(float *)malloc(sizeof(float)*nx*ny*nz);
		float *dataDistance=(float *)malloc(sizeof(float)*nx*ny*nz);

		for(i=0;i<nx*ny*nz;i++){
			dataPair[i]=FLT_MAX;
			dataDistance[i]=FLT_MAX;
		}

		rewind(fpharddata);
		flag=0;
		i=0;
		j=0;
		k=0;
		int out=0;
		int idx,idy,idz;
		float xtmp,ytmp,ztmp,datatmp;
		int blockx,blocky,blockz;
		float remx,remy,remz;
		while((fgets(line,BUFSIZE,fpharddata))!=0){
			if(flag==1){
				sscanf(line,"%d",&numcols);
			}
			if(flag>=numcols+2){
				//sscanf(line,"%f %f %f %f",&(dataCoordinate[i]),&(dataCoordinate[i+1]),&(dataCoordinate[i+2]),&(dataValue[i]));
				sscanf(line,"%f %f %f %f",&xtmp,&ytmp,&ztmp,&datatmp);
				dataCoordinate[i]	=(TYPE)xtmp;
				dataCoordinate[i+1]	=(TYPE)ytmp;
				dataCoordinate[i+2]	=(TYPE)ztmp;
				dataValue[i]		=(TYPE)datatmp;
				//printf("x[%d]=%f\ty[%d]=%f\tz[%d]=%f\tdata[%d]=%f\n",i,dataCoordinate[i],i,dataCoordinate[i+1],i,dataCoordinate[i+2],i,dataValue[i]);

				if(
					xtmp>=params.xlo && xtmp<params.xlo + params.nx*params.h &&
					ytmp>=params.ylo && ytmp<params.ylo + params.ny*params.h &&
					ztmp>=params.zlo && ztmp<params.zlo + params.nz*params.h 
				){
					blockx = (int)floor(xtmp/params.h);
					blocky = (int)floor(ytmp/params.h);
					blockz = (int)floor(ztmp/params.h);
					remx = xtmp - (params.xlo + blockx*params.h);
					remy = ytmp - (params.ylo + blocky*params.h);
					remz = ztmp - (params.zlo + blockz*params.h);
					//if(remx<params.h*0.5 && remy<params.h*0.5 && remz<params.h*0.5){
					//	j=blockx + blocky*params.nx + blockz*(params.nx * params.ny);
					//	printf("xd[%d]=%f\tyd[%d]=%f\tzd[%d]=%f\n",i,dataCoordinate[i],i,dataCoordinate[i+1],i,dataCoordinate[i+2]);
					//	printf("xs[%d]=%f\tys[%d]=%f\tzs[%d]=%f\n",j,params.xlo + blockx*params.h,j,params.ylo + blocky*params.h,j,params.zlo + blockz*params.h);
					//}

					idx=(blockx + (int)(remx>params.h*0.5) ); 
					idy=(blocky + (int)(remy>params.h*0.5) ); 
					idz=(blockz + (int)(remz>params.h*0.5) ); 

					j=(blockx + (int)(remx>params.h*0.5) )+ (blocky + (int)(remy>params.h*0.5))*params.nx + (blockz + (int)(remz>params.h*0.5))*(params.nx * params.ny);
					//printf("xd[%d]=%f\tyd[%d]=%f\tzd[%d]=%f\n",i,dataCoordinate[i],i,dataCoordinate[i+1],i,dataCoordinate[i+2]);
					//printf("xs[%d]=%f\tys[%d]=%f\tzs[%d]=%f\n",j,params.xlo + (blockx + (int)(remx>params.h*0.5) )*params.h,j,params.ylo + (blocky + (int)(remy>params.h*0.5))*params.h,j,params.zlo + (blockz + (int)(remz>params.h*0.5))*params.h);


	         			getCoordinates(idx,idy,idz,params,&xcoord1,&ycoord1,&zcoord1);
	         			rsquared = 	(xcoord1-xtmp)*(xcoord1-xtmp)+
	  						(ycoord1-ytmp)*(ycoord1-ytmp)+
	  						(zcoord1-ztmp)*(zcoord1-ztmp);
		 			r = sqrt(rsquared);

					if(r<dataDistance[j]){
						dataDistance[j]=r;
						dataPair[j]=dataValue[i];
						++k;
					}

					//printf("%d %d\n",dataPair[k].indexdata,dataPair[k].indexlattice);
					
				}
				else{
					//printf("x[%d]=%f\ty[%d]=%f\tz[%d]=%f\tdata[%d]=%f\n",i,dataCoordinate[i],i,dataCoordinate[i+1],i,dataCoordinate[i+2],i,dataValue[i]);
					++out;
				}

				++i;
			}
			flag++;
		}
		
	

		float condCost=0.0;
		float condCostMin=FLT_MAX;
		float condCostInitial=0.0;

		int control=0;
		int harddata=k;
		
		printf("(before) condCost=%f (control=%d, out=%d)\n",condCost,k,out);
		int p=0;
		int t=0;
		int tmin=-1;


		int sizeTable=20;
		float sinTable[sizeTable+1];
		float cosTable[sizeTable+1];
		for(t=0;t<=sizeTable;t++){
			//sinTable[t]=sin( (((float)t)/((float)sizeTable))*0.5*M_PI )*sin( (((float)t)/((float)sizeTable))*0.5*M_PI );
			//cosTable[t]=cos( (((float)t)/((float)sizeTable))*0.5*M_PI )*cos( (((float)t)/((float)sizeTable))*0.5*M_PI );
			sinTable[t]=sin( (((float)t)/((float)sizeTable))*0.25*M_PI - 0.125*M_PI );
			cosTable[t]=cos( (((float)t)/((float)sizeTable))*0.25*M_PI - 0.125*M_PI );
		}
			
		imageNext = (TYPE *)malloc(sizeof(TYPE)*nxExt*nyExt*nzExt);
		imageZeroNext = imageNext+initialAddress;


		int falseStart=0;
		time_t falseStartSeed;
		time_t falseStartMin;
		int maxFalseStart=500;

		for(p=0;p<1000;p++){
			//printf("BEGIN ITERATION\n");

			condCost=0.0;
			control=0;	


			if(p==0){


				for(falseStart=1;falseStart<=maxFalseStart;falseStart++){

					falseStartSeed=((time_t)falseStart)*time(NULL); 
					//srand(falseStart*time(NULL));
					srand(falseStartSeed);
					for(k=0;k<nzExt;k++){
						for(i=0;i<nxExt;i++){
							for(j=0;j<nyExt;j++){
								image[k*(nxExt*nyExt) + i*nyExt + j]= generateGaussianNoise(); // N(0,1) random
								//image[k*(nxExt*nyExt) + i*nyExt + j]= (TYPE)( 2.0*((TYPE)rand())/((TYPE)RAND_MAX) - 1.0 ); // uniform [-1,1] random
							}
						}
					}
			
					condCost=0.0;
					control=0;	

					for(k=0;k<nz;k++){
	  					for(j=0;j<ny;j++){
	    						for(i=0;i<nx;i++){
              							sum=0.0;
              							getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
	      							for(nk=-neighRadius;nk<=neighRadius;nk++){
	        							for(nj=-neighRadius;nj<=neighRadius;nj++){
	        			  					for(ni=-neighRadius;ni<=neighRadius;ni++){
                			  
                			    						weight= *(kernelWeight+
												(ni+neighRadius)+
												(nj+neighRadius)*neighSide+
												(nk+neighRadius)*neighSide*neighSide);
 
	        			    						sum=sum + weight * 
                			    		    					*(imageZero+(
					    								(i+ni)+
					    								(j+nj)*(nxExt)+
					    								(k+nk)*(nxExt)*(nyExt)
					    								)
					    							);
	        			  					}
	        							}
	      							}
              
								//sum = sum/((float)totalNeighbours); 	
								sum = sum*totalNeighboursinv; 	
								if(dataPair[i+j*nx+k*(nx*ny)] < FLT_MAX){
									//condCost = condCost + (sum+mu-dataPair[i+j*nx+k*(nx*ny)])*(sum+mu-dataPair[i+j*nx+k*(nx*ny)]);
									condCost = condCost + (sum-dataPair[i+j*nx+k*(nx*ny)])*(sum-dataPair[i+j*nx+k*(nx*ny)]);
									control++;
								}

	    						}
	  					}
					}

					if(falseStart==1) condCostInitial=condCost;
					if(condCost < condCostMin){
						condCostMin = condCost;
						falseStartMin = falseStartSeed;
						//fprintf(stderr,"MIN: p=%d t=%d (after) condCost=%f (control=%d)\n",p,t,condCost,control);
						fprintf(stderr,"FALSE-START %d: condCost=%f\tcondCostInitial=%f\tcondCostMin=%f\t(%3.6f)--------MIN\n",falseStart,condCost,condCostInitial,condCostMin,(condCost/condCostInitial)*100.0);
					}
					else{
						//fprintf(stderr,"---: p=%d t=%d (after) condCost=%f (control=%d)\n",p,t,condCost,control);
						fprintf(stderr,"FALSE-START %d: condCost=%f\tcondCostInitial=%f\tcondCostMin=%f\t(%3.6f)\n",falseStart,condCost,condCostInitial,condCostMin,(condCost/condCostInitial)*100.0);
					}
				}

				//condCostInitial=condCostMin;
	
				falseStartSeed=falseStartMin; 
				//srand(falseStart*time(NULL));
				srand(falseStartSeed);
				for(k=0;k<nzExt;k++){
					for(i=0;i<nxExt;i++){
						for(j=0;j<nyExt;j++){
							image[k*(nxExt*nyExt) + i*nyExt + j]= generateGaussianNoise(); // N(0,1) random
							//image[k*(nxExt*nyExt) + i*nyExt + j]= (TYPE)( 2.0*((TYPE)rand())/((TYPE)RAND_MAX) - 1.0 ); // uniform [-1,1] random
						}
					}
				}
				srand(((time_t)(maxFalseStart+1))*time(NULL));

				fprintf(stderr,"ITERATION %d: condCost=%f\tcondCostInitial=%f\tcondCostBest=%f\t(%3.6f)\n",p,condCost,condCostInitial,condCostMin,(condCost/condCostInitial)*100.0);

			}
			else{
				for(k=0;k<nz;k++){
  					for(j=0;j<ny;j++){
    						for(i=0;i<nx;i++){
      							sum=0.0;
      							getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
      							for(nk=-neighRadius;nk<=neighRadius;nk++){
        							for(nj=-neighRadius;nj<=neighRadius;nj++){
        			  					for(ni=-neighRadius;ni<=neighRadius;ni++){
        			  
        			    						weight= *(kernelWeight+
											(ni+neighRadius)+
											(nj+neighRadius)*neighSide+
											(nk+neighRadius)*neighSide*neighSide);

        			    						sum=sum + weight * 
        			    		    					*(imageZero+(
				    								(i+ni)+
				    								(j+nj)*(nxExt)+
				    								(k+nk)*(nxExt)*(nyExt)
				    								)
				    							);
        			  					}
        							}
      							}
      
							//sum = sum/((float)totalNeighbours); 	
							sum = sum*totalNeighboursinv; 	
							if(dataPair[i+j*nx+k*(nx*ny)] < FLT_MAX){
								//condCost = condCost + (sum+mu-dataPair[i+j*nx+k*(nx*ny)])*(sum+mu-dataPair[i+j*nx+k*(nx*ny)]);
								condCost = condCost + (sum-dataPair[i+j*nx+k*(nx*ny)])*(sum-dataPair[i+j*nx+k*(nx*ny)]);
								control++;
							}

    						}
  					}
				}

				fprintf(stderr,"ITERATION %d: condCost=%f\tcondCostInitial=%f\tcondCostBest=%f\t(%3.6f)\n",p,condCost,condCostInitial,condCostMin,(condCost/condCostInitial)*100.0);
			}



			//if(condCost < condCostMin){
			//	condCostMin = condCost;
			//	tmin=0;
			//}


			//srand(time(NULL) + p);
			//srand(1);
			int iinf,isup;
			int jinf,jsup;
			int kinf,ksup;
			int d=0;
			for(k=0;k<nzExt;k++){
				for(j=0;j<nyExt;j++){
					for(i=0;i<nxExt;i++){
						//if(dataPair[i+j*nx+k*(nx*ny)] < FLT_MAX){
						//	k=d%(nx*ny);
						//	j=d/nx;
						//	i=d%nx;
						//	iinf=MAX(i-neighRadius,0); 
						//	jinf=MAX(j-neighRadius,0); 
						//	kinf=MAX(k-neighRadius,0); 
						//	isup=MIN(i+neighRadius,nxExt);
						//	jsup=MIN(j+neighRadius,nyExt);
						//	ksup=MIN(k+neighRadius,nzExt);
			      			//	for(nk=kinf;nk<ksup;nk++){
			        		//		for(nj=jinf;nj<jsup;nj++){
			          		//			for(ni=iinf;ni<isup;ni++){
										//imageZeroNext[
										////k*(nxExt*nyExt) + i*nyExt + j
										//(i+ni)+	(j+nj)*(nxExt)+(k+nk)*(nxExt)*(nyExt)
										//]= generateGaussianNoise(); // N(0,1) random
										*(imageNext+(
				    							(i)+
				    							(j)*(nxExt)+
				    							(k)*(nxExt)*(nyExt)
				    							)
				    						) =generateGaussianNoise(); 


						//			}
						//		}
						//	}
						//}
					}
				}
			}
			

			TYPE tmp;
		
			tmin=-1;
				
	
			for(t=0;t<=sizeTable;t++){
				//imageMiddle = genRandomImageTrigon(nxExt,nyExt,nzExt,imageZero,imageZeroNext,t);
				//imageZeroMiddle = imageMiddle+initialAddress;

				condCost=0.0;
				control=0;	
				for(k=0;k<nz;k++){
			  		for(j=0;j<ny;j++){
			    			for(i=0;i<nx;i++){
							if(dataPair[i+j*nx+k*(nx*ny)] < FLT_MAX){
	        	      				sum=0.0;
	        	      				getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
			      				for(nk=-neighRadius;nk<=neighRadius;nk++){
			        				for(nj=-neighRadius;nj<=neighRadius;nj++){
			          					for(ni=-neighRadius;ni<=neighRadius;ni++){
	        	          
	        	            						weight= *(kernelWeight+
											(ni+neighRadius)+
											(nj+neighRadius)*neighSide+
											(nk+neighRadius)*neighSide*neighSide);
	 
			            						sum=sum + weight * 
											(
												cosTable[t]*
	        	            		    						*(imageZero+(
				    									(i+ni)+
				    									(j+nj)*(nxExt)+
				    									(k+nk)*(nxExt)*(nyExt)
				    									)
				    								)
												+
												sinTable[t]*
												*(imageZeroNext+(
				    									(i+ni)+
				    									(j+nj)*(nxExt)+
				    									(k+nk)*(nxExt)*(nyExt)
				    									)
				    								)
											)
											;
			          					}
			        				}
			      				}
	        	      
							//sum = sum/((float)totalNeighbours); 	
							sum = sum*totalNeighboursinv; 	
							//if(dataPair[i+j*nx+k*(nx*ny)] < FLT_MAX){
								//condCost = condCost + (sum+mu-dataPair[i+j*nx+k*(nx*ny)])*(sum+mu-dataPair[i+j*nx+k*(nx*ny)]);
								condCost = condCost + (sum-dataPair[i+j*nx+k*(nx*ny)])*(sum-dataPair[i+j*nx+k*(nx*ny)]);
								control++;
							}
	
			    			}
			  		}
				}
					
				//printf("     p=%d t=%d (after) condCost=%f (control=%d)\n",p,t,condCost,control);
	
				if(condCost < condCostMin){
					condCostMin = condCost;
					tmin = t;
					//fprintf(stderr,"MIN: p=%d t=%d (after) condCost=%f (control=%d)\n",p,t,condCost,control);
					fprintf(stderr,"ITERATION %d/%d: condCost=%f\tcondCostInitial=%f\tcondCostBest=%f\t(%3.6f)--------MIN\n",p,t,condCost,condCostInitial,condCostMin,(condCost/condCostInitial)*100.0);
				}
				else{
					//fprintf(stderr,"---: p=%d t=%d (after) condCost=%f (control=%d)\n",p,t,condCost,control);
					fprintf(stderr,"ITERATION %d/%d: condCost=%f\tcondCostInitial=%f\tcondCostBest=%f\t(%3.6f)\n",p,t,condCost,condCostInitial,condCostMin,(condCost/condCostInitial)*100.0);
				}
			
				//freeRandomImage(imageMiddle);

			}

			//printf("tmin=%d\n",tmin);
			if(tmin!=-1){
				for(k=0;k<nzExt;k++){
				   	for(i=0;i<nxExt;i++){
						for(j=0;j<nyExt;j++){
							if(dataPair[i+j*nx+k*(nx*ny)] < FLT_MAX){
							//printf("i=%d j=%d k=%d\n",i,j,k);
							tmp 
							=
							(
								cosTable[tmin]*
	        		            				image[
					    					(j)+
					    					(i)*(nyExt)+
					    					(k)*(nxExt)*(nyExt)
					    					]
								+
								sinTable[tmin]*
									imageNext[
					    					(j)+
					    					(i)*(nyExt)+
					    					(k)*(nxExt)*(nyExt)
					    					]
							);	
							//printf("%f %f %f\n",tmp,image[(j)+(i)*(nyExt)+(k)*(nxExt)*(nyExt)],imageNext[(j)+(i)*(nyExt)+(k)*(nxExt)*(nyExt)]);
							image[
					    			(j)+
					    			(i)*(nyExt)+
					    			(k)*(nxExt)*(nyExt)
					    			]
							= tmp;
							}
						}
					}
				}
			}
			//printf("END ITERATION\n");

			//freeRandomImage(imageNext);
			//freeRandomImage(image);
			//image = genRandomImage(nxExt,nyExt,nzExt);
			//imageZero = image+initialAddress;

		}


		fprintf(stderr,"DIFF:\n");

		for(k=0;k<nz;k++){
		  for(j=0;j<ny;j++){
		    for(i=0;i<nx;i++){
        	      sum=0.0;
        	      getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
		      for(nk=-neighRadius;nk<=neighRadius;nk++){
		        for(nj=-neighRadius;nj<=neighRadius;nj++){
		          for(ni=-neighRadius;ni<=neighRadius;ni++){
        	            weight= *(kernelWeight+
					(ni+neighRadius)+
					(nj+neighRadius)*neighSide+
					(nk+neighRadius)*neighSide*neighSide);
 
		            sum=sum + weight * 
        	            		*(imageZero+(
						(i+ni)+
						(j+nj)*(nxExt)+
						(k+nk)*(nxExt)*(nyExt)
						)
					);
		          }
		        }
		      }
			//sum = sum/((float)totalNeighbours); 	
			sum = sum*totalNeighboursinv; 	

        	      //fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
			fprintf(fpgslib,"%f %f %f %f\n",xcoord1,ycoord1,zcoord1,sum);
        	      	fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
        	      		//printf("image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);

			if(dataPair[i+j*nx+k*(nx*ny)] < FLT_MAX){
				//fprintf(stderr,"%d %d %d %d %f %f %f\n",i,j,k,i+j*nx+k*(nx*ny),(sum+mu),dataPair[i+j*nx+k*(nx*ny)],((sum+mu)-dataPair[i+j*nx+k*(nx*ny)])*((sum+mu)-dataPair[i+j*nx+k*(nx*ny)]));
				fprintf(stderr,"%d %d %d %d %f %f %f\n",i,j,k,i+j*nx+k*(nx*ny),(sum),dataPair[i+j*nx+k*(nx*ny)],((sum)-dataPair[i+j*nx+k*(nx*ny)])*((sum)-dataPair[i+j*nx+k*(nx*ny)]));
			}



		    }
		  }
		}



	
		//freeRandomImage(imageNext);
		free(imageNext);

		free(dataDistance);
		free(dataPair);
		free(dataCoordinate);
		free(dataValue);
		free(line);	
		free(line2);	
	}


	free(kernelWeight);
	free(distanceValue);
	free(distanceWeight);
	freeRandomImage(image);

	free(line);

	//fclose(fpoctave);
	if(mode==0 || mode==3 || mode==4){ fclose(fpgslib); fclose(fpoctave);}
	if(mode==1) fclose(fpweight);
	if(mode==2 || mode==3) fclose(fpdistanceweight);
	if(mode==4) fclose(fpharddata);

	return 0;
}
