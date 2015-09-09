#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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


	int bufferNodes=ceil(a/h);
	int neighRadius=bufferNodes;
	int neighSide=2*neighRadius+1;

	int nxExt=nx+2*bufferNodes;
	int nyExt=ny+2*bufferNodes;
	int nzExt=nz+2*bufferNodes;


	FILE *fpoctave,*fpgslib,*fpweight,*fpdistanceweight;

	//fpoctave = fopen("octaveout.m","w");
	if(mode==0 || mode==3){
		fpgslib  = fopen("gslibout.dat","w");
	}
	if(mode==1){
		fpweight = fopen("weightout.dat","w");
	}
	if(mode==2){
		fpdistanceweight = fopen("targetdistanceweights.dat","w");
	}
	if(mode==3){
		fpdistanceweight = fopen("currentbestdistanceweights.dat","r");
	}

	//fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
	if(mode==0 || mode==3) fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");

	latticeParams params;
	params.nx=nx;
	params.ny=ny;
	params.nz=nz;
	params.xlo=xlo;
	params.ylo=ylo;
	params.zlo=zlo;
	params.h=h;

	TYPE *image;
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
	if(mode==3){
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
		if(mode==3){
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
              if(mode==0 || mode==3){ 
			fprintf(fpgslib,"%f %f %f %f\n",xcoord1,ycoord1,zcoord1,sum);
              		//printf("image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
		}

	    }
	  }
	}

	free(kernelWeight);
	free(distanceValue);
	free(distanceWeight);
	freeRandomImage(image);

	free(line);

	//fclose(fpoctave);
	if(mode==0 || mode==3) fclose(fpgslib);
	if(mode==1) fclose(fpweight);
	if(mode==2 || mode==3) fclose(fpdistanceweight);

	return 0;
}
