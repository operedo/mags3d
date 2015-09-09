#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>

int main(int argc,char *argv[]){

	int i,j,k,p,q;

	char *line;
	int BUFSIZE=100;
	ssize_t read;

	line=(char *)malloc(sizeof(char)*BUFSIZE + 1);

	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int nz = atoi(argv[3]);
	TYPE xlo = (TYPE)atof(argv[4]);
	TYPE ylo = (TYPE)atof(argv[5]);
	TYPE zlo = (TYPE)atof(argv[6]);
	TYPE h = (TYPE)atof(argv[7]);
	TYPE a = (TYPE)atof(argv[8]);

	int bufferNodes=ceil(a/h);
	int neighRadius=bufferNodes;
	int neighSide=2*neighRadius+1;

	int nxExt=nx+2*bufferNodes;
	int nyExt=ny+2*bufferNodes;
	int nzExt=nz+2*bufferNodes;

//#ifdef OCTAVE
//	printf("image=zeros(%d,%d,%d);\n",nx,ny,nz);
//#else
//#ifdef GSLIB
//	printf("Test\n4\nX\nY\nZ\ndata\n");
//#endif
//#endif

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
	int initialAddress = 	bufferNodes*nyExt*nxExt + 
				bufferNodes*nxExt + 
				bufferNodes; 
	TYPE *imageZero = image+initialAddress;

	TYPE weightTmp[MAXDISTVALS];
	TYPE distanceValue[MAXDISTVALS];
	int distanceElemNumber[MAXDISTVALS];
	TYPE weightDiff[MAXDISTVALS];
	int distanceIndex[MAXDISTVALS];

	for(i=0;i<MAXDISTVALS;i++){
		distanceValue[i]=-1.0;
		distanceIndex[i]=-1;
	}


	int retTarget = genTargetVariogram(nx,ny,nz,xlo,ylo,zlo,h,a,imageZero,0);

	FILE *fpweight, *fptargetvariogram;
	fpweight=fopen("weightout.dat","r");
	fptargetvariogram=fopen("targetvariogram.dat","r");

	int nlags,num;
	if((fgets(line,BUFSIZE,fptargetvariogram))==0){
		fprintf(stderr,"ERROR: variogram nlags not specified in variogram file.\n");
		return 1;	
	}
	sscanf(line,"%d\n",&nlags);

	TYPE *expVariogram = (TYPE *)malloc(sizeof(TYPE)*nlags);
	int *npairs=(int *)malloc(sizeof(int)*nlags);

	TYPE vgval;
	i=0;
	while((fgets(line,BUFSIZE,fptargetvariogram))!=0){
		sscanf(line,"%lf %d",&vgval,&num);
		expVariogram[i]=vgval;	
		npairs[i]=num;
		//printf("%f %d\n",vgval,num);
		i++;	
	}

	fprintf(stderr,"experimental variogram and number of pairs loaded.\n");
	fclose(fptargetvariogram);


	int side;

	if((fgets(line,BUFSIZE,fpweight))==0){
		fprintf(stderr,"ERROR: weights cell value not specified in weight file.\n");
		return 2;	
	}
	sscanf(line,"%d\n",&side);

	//printf("%d\n",side);

	nweight = side*side*side;
	TYPE *weight = (TYPE *)malloc(sizeof(TYPE)*nweight); 
	TYPE *weightPerturbed = (TYPE *)malloc(sizeof(TYPE)*nweight); 

	int midside;
	TYPE val;

	midside=ceil(side/2);
	//midsidesqrt3=ceil((TYPE)(side/2)*(TYPE)sqrt(3.0));

	// esto funciona con D+1 weights
	//int distance;
	////for(i=0;k<=midside;k++){
	//while((fgets(line,BUFSIZE,fpweight))!=0){
	//	sscanf(line,"%d %lf",&distance,&val);	 
	//	if(distance>=MAXDISTVALS){
	//		printf("ERROR: too many distance-weight pair in weights input file.\n");
	//		return 2;
	//	}
	//	distanceValue[distance]=val;
	//        printf("distanceValue(%d)=%f;\n",distance,distanceValue[distance]);
	//}
	////}
	//int idiag;
	//for(k=0;k<=midside;k++){
	//   for(j=0;j<=midside;j++){
	//      for(i=0;i<=midside;i++){
	//         idiag = nearbyint(sqrt((TYPE)(i*i+j*j+k*k)));
	//         //if(idiag<=midside){
	//	 weight[i + j*side + k*side*side]=distanceValue[idiag];
	//         //}
	//         //else{
	//	 //   weight[i + j*side + k*side*side]=0.0;
	//         //}
	//         //printf("weight(%d,%d,%d)=%f;\n",i+1,j+1,k+1,weight[i + j*side + k*side*side]);
	//      }
	//   }
	//}


	// esto funciona con (D+1)^3 weights
	//for(k=0;k<=midside;k++){
	//   for(j=0;j<=midside;j++){
	//      for(i=0;i<=midside;i++){
	//         if((fgets(line,BUFSIZE,fpweight))==0){
	//	    printf("ERROR: weights cell value not specified in weight file.\n");
	//	    return 2;	
	//         }
	//         sscanf(line,"%lf",&val);
	//	 weight[i + j*side + k*side*side]=val;
	//         //printf("%d %d %d %f\n",i,j,k,val);
	//      }
	//   }
	//}

	TYPE xcoord1,ycoord1,zcoord1;
	TYPE xcoord2,ycoord2,zcoord2,dist;
	getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);
	int distCounter=0;
	int distindex=0;
	while((fgets(line,BUFSIZE,fpweight))!=0){
		sscanf(line,"%d %d %d %lf",&i,&j,&k,&val);
		weight[i + j*side + k*side*side]=val;
		getCoordinates(i-neighRadius,j-neighRadius,k-neighRadius,params,&xcoord2,&ycoord2,&zcoord2);
               	dist = sqrt(	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	       		(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	       		(zcoord1-zcoord2)*(zcoord1-zcoord2));

		distindex = existValue(distanceValue,MAXDISTVALS,dist);
		if(distindex==-1){
			distanceValue[distCounter]=dist;
			distanceElemNumber[distCounter]=1;
			distCounter++;
		}
		//else{
		//	distanceElemNumber[distindex]++;
		//}
		//printf("%d %d %d %f\n",i,j,k,val);
	}

	quickSort(distanceValue,0,distCounter-1);

	//for(i=0;i<distCounter;i++)
	//	printf("%f ",distanceValue[i]);
	//printf("\n");



	fprintf(stderr,"weights (octant) loaded.\n");
	fclose(fpweight);

	int ieff,jeff,keff;
	for(k=0;k<side;k++){
	   for(j=0;j<side;j++){
	      for(i=0;i<side;i++){
		 getCoordinates(i-neighRadius,j-neighRadius,k-neighRadius,params,&xcoord2,&ycoord2,&zcoord2);
               	 dist = sqrt(	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	       		(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	       		(zcoord1-zcoord2)*(zcoord1-zcoord2));
		 for(p=0;p<distCounter;p++){
		    if(dist==distanceValue[p]){
			distanceIndex[p]=i + j*side + k*side*side;
			distanceElemNumber[p]++;
		    }
		 }



	         if(!(i<=midside && j<=midside && k<=midside)){
	            ieff=i;
	            jeff=j;
	            keff=k;
	            if(i>midside){
	               ieff=side-i-1;
	            }
	            if(j>midside){
	               jeff=side-j-1;
	            }
	            if(k>midside){
	               keff=side-k-1;
	            }
		    weight[i + j*side + k*side*side] = weight[ieff + jeff*side + keff*side*side];
	         }
	         //printf("%d %d %d %f\n",i,j,k,weight[i + j*side + k*side*side]);
	         //printf("weight(%d,%d,%d)=%f;\n",i+1,j+1,k+1,weight[i + j*side + k*side*side]);
		 weightPerturbed[i + j*side + k*side*side] = weight[i + j*side + k*side*side];

		

	      }
	   }
	}


	//for(i=0;i<distCounter;i++)
	//	fprintf(stderr,"%f ",distanceValue[i]);
	//printf("\n");
	//for(i=0;i<distCounter;i++)
	//	fprintf(stderr,"%d ",distanceIndex[i]);
	//printf("\n");


	fprintf(stderr,"weights loaded.\n");

	//return 1;

	int ret = genMovingAverageImage(imageZero,weight,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,0,100);

	fprintf(stderr,"starting cost function calculation...\n");
	TYPE cost = costFunctionGSLIB(nlags,expVariogram,npairs,0,0);
	fprintf(stderr,"stoping cost function calculation...\n");

	fprintf(stderr,"initial cost=%f\n",cost);

	TYPE initialCost = cost;
	TYPE tempCost=0.0;
	TYPE tempCostStep1=0.0;
	TYPE tempCostStep2=0.0;
	TYPE tempCostStep3=0.0;
	TYPE lastCost=0.0;
	TYPE currentCost= cost;

	int counter=0;

	int iters=6000;
	int linesearchiters=20;

	TYPE alphaInitial=1.0, alpha;
	TYPE lambda=0.5;

	//int iters=2000;
	int report=100;
	int updateTemperature=500;
	int resetCounter=0;
	int averageCounter=0;
	TYPE modifiedRadius;
	TYPE tol=0.001; // percentage of decrease w/r initial cost
	TYPE newval,oldval;
	//TYPE temperatureInitial=0.000075;
	//TYPE temperatureInitial=0.0003;
	TYPE temperatureInitial=0.01;
	TYPE temperature=temperatureInitial;
	TYPE probability, ran;

	TYPE delta=1.0e-03;
	TYPE percent=0.5;
	TYPE distIn;
	TYPE norm2=0.0, invnorm2;
	TYPE maxval=-1.0e+20, invmaxval;
	TYPE minweightval=1.0e+20, invminweightval;

	//for(i=1;i<=iters;i++){
	for(i=1;i<=5;i++){
		if((currentCost/initialCost)<=tol){
			fprintf(stdout,"STOP: convergence achieved. Congratulations!\n");

			//modifiedRadius = modifyWeights(	
			//			weight,
			//			neighRadius,
			//			neighSide,
			//			distanceValue,
			//			distanceIndex,
			//			distCounter,
			//			params,
			//			2,
			//			0,
			//			0.0,
			//			&newval,
			//			&oldval,
			//			i);

			return 0;
		}
		else{
			fprintf(stdout,"ITERATION %d:\n",i);

			norm2=0.0;
			maxval=-1.0e+20;
			minweightval=1.0e+20;

			for(j=0;j<distCounter;j++){

				distIn = distanceValue[i];
				oldval = weightPerturbed[distanceIndex[j]];
				newval = oldval + delta; 

				modifiedRadius = modifyWeights(	
						weightPerturbed,
						neighRadius,
						neighSide,
						distanceValue,
						distanceIndex,
						distCounter,
						params,
						1,
						distIn,
						newval,
						NULL,
						NULL,
						i);

				ret = genMovingAverageImage(imageZero,weightPerturbed,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,i,report);
				tempCost = costFunctionGSLIB(nlags,expVariogram,npairs,i,0);

				newval = oldval - delta; 
				modifiedRadius = modifyWeights(	
						weightPerturbed,
						neighRadius,
						neighSide,
						distanceValue,
						distanceIndex,
						distCounter,
						params,
						1,
						distIn,
						newval,
						NULL,
						NULL,
						i);

				ret = genMovingAverageImage(imageZero,weightPerturbed,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,i,report);
				tempCost = tempCost - costFunctionGSLIB(nlags,expVariogram,npairs,i,0);
				//weightDiff[j] = ((TYPE)distanceElemNumber[j]) * tempCost / (2.0*delta);
				weightDiff[j] = tempCost / (2.0*delta);
				norm2 = norm2 + weightDiff[j]*weightDiff[j];   
				if(abs(weightDiff[j])>maxval){
					maxval=abs(weightDiff[j]); 
				}
				if(oldval!=0.0 && abs(oldval)<minweightval){
					minweightval = oldval; 
				}

				//printf("%d %f\n",j,weightDiff[j]);

				newval = oldval; 
				modifiedRadius = modifyWeights(	
						weightPerturbed,
						neighRadius,
						neighSide,
						distanceValue,
						distanceIndex,
						distCounter,
						params,
						1,
						distIn,
						newval,
						NULL,
						NULL,
						i);

			}

			invnorm2=1.0/norm2;
			invmaxval=1.0/maxval;
			for(j=0;j<distCounter;j++){
				weightDiff[j] = weightDiff[j] * invnorm2;
				//weightDiff[j] = weightDiff[j] * invmaxval;
				printf("%d %d %f %f %f\n",j,distanceElemNumber[j],distanceValue[j],weightPerturbed[distanceIndex[j]],weightDiff[j]);
			}

			//printf("diff=zeros(%d,%d,%d);\n",side,side,side);
			//for(k=0;k<side;k++){
	   		//for(j=0;j<side;j++){
	      		//for(i=0;i<side;i++){
		 	//	getCoordinates(i-neighRadius,j-neighRadius,k-neighRadius,params,&xcoord2,&ycoord2,&zcoord2);
               	 	//	dist = sqrt(	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	       		//		(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	       		//		(zcoord1-zcoord2)*(zcoord1-zcoord2));
			//	printf("diff(%d,%d,%d)=%f;\n",i+1,j+1,k+1,weightDiff[existValue(distanceValue,MAXDISTVALS,dist)]);
			//}
			//}
			//}
		
			/* line search */
			counter=0;
			alpha=alphaInitial;
			//while(counter<linesearchiters){
			while(counter<10){
				for(j=0;j<distCounter;j++){

					distIn = distanceValue[j];
					weightTmp[j] = weightPerturbed[distanceIndex[j]];
					//newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]/((TYPE)distanceElemNumber[j]); 
					newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]; 

					modifiedRadius = modifyWeights(	
						weightPerturbed,
						neighRadius,
						neighSide,
						distanceValue,
						distanceIndex,
						distCounter,
						params,
						1,
						distIn,
						newval,
						NULL,
						NULL,
						i+1+j);

				}

				ret = genMovingAverageImage(imageZero,weightPerturbed,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,i,report);
				tempCost = costFunctionGSLIB(nlags,expVariogram,npairs,i+1,0);

				if(counter==0){ 
					tempCostStep1 = tempCost;
				}
				if(counter==1){
					tempCostStep2 = tempCostStep1;	
					tempCostStep1 = tempCost;	
				}
				if(counter>=2){
					tempCostStep3 = tempCostStep2;
					tempCostStep2 = tempCostStep1;	
					tempCostStep1 = tempCost;	
				}

				printf("iter=%d line=%d cost=%f perc=%f alpha=%f norm2=%f\n",i,counter+1,tempCost, tempCost/initialCost,alpha,norm2);

				if( counter>=2 ){
					if(tempCostStep1>tempCostStep2 && tempCostStep2>tempCostStep3 && tempCostStep3<initialCost){
						printf("LINE-SEARCH END: cost=%f perc=%f alpha=%f norm2=%f\n",tempCostStep3,tempCostStep3/initialCost,alpha/(lambda*lambda),norm2);
						
						alpha=alpha/(lambda*lambda);
						for(j=0;j<distCounter;j++){
		
							distIn = distanceValue[j];
							weightTmp[j] = weightPerturbed[distanceIndex[j]];
							//newval = weightTmp[j] - alpha * weightDiff[j]; 
							//newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]/((TYPE)distanceElemNumber[j]); 
							newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]; 
		
							modifiedRadius = modifyWeights(	
								weightPerturbed,
								neighRadius,
								neighSide,
								distanceValue,
								distanceIndex,
								distCounter,
								params,
								1,
								distIn,
								newval,
								NULL,
								NULL,
								i+1+j);
		
						}

						break;

					}

					if(tempCostStep1>tempCostStep2 && tempCostStep2>tempCostStep3 && tempCostStep3>initialCost){
						printf("LINE-SEARCH RESTART: alpha=%f\n",(6.0*alphaInitial));
						
						alpha=6.0 * alphaInitial;
						for(j=0;j<distCounter;j++){
		
							distIn = distanceValue[j];
							weightTmp[j] = weightPerturbed[distanceIndex[j]];
							//newval = weightTmp[j] - alpha * weightDiff[j]; 
							//newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]/((TYPE)distanceElemNumber[j]); 
							newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]; 
		
							modifiedRadius = modifyWeights(	
								weightPerturbed,
								neighRadius,
								neighSide,
								distanceValue,
								distanceIndex,
								distCounter,
								params,
								1,
								distIn,
								newval,
								NULL,
								NULL,
								i+1+j);
		
						}

						break;

					}

					if(tempCostStep1>tempCostStep2 && tempCostStep2<tempCostStep3){
						printf("LINE-SEARCH END: cost=%f perc=%f alpha=%f norm2=%f\n",tempCostStep2,tempCostStep2/initialCost,alpha/lambda,norm2);

						alpha=alpha/lambda;
						for(j=0;j<distCounter;j++){
		
							distIn = distanceValue[j];
							weightTmp[j] = weightPerturbed[distanceIndex[j]];
							//newval = weightTmp[j] - alpha * weightDiff[j]; 
							//newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]/((TYPE)distanceElemNumber[j]); 
							newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]; 
		
							modifiedRadius = modifyWeights(	
								weightPerturbed,
								neighRadius,
								neighSide,
								distanceValue,
								distanceIndex,
								distCounter,
								params,
								1,
								distIn,
								newval,
								NULL,
								NULL,
								i+1+j);
		
						}

						break;

					}
				}


				if(counter==9){
					if(tempCostStep1<tempCostStep2 && tempCostStep2<tempCostStep3){
						printf("LINE-SEARCH END: cost=%f perc=%f alpha=%f\n",tempCostStep1,tempCostStep1/initialCost,alpha);
						for(j=0;j<distCounter;j++){
	
							distIn = distanceValue[j];
							weightTmp[j] = weightPerturbed[distanceIndex[j]];
							//newval = weightTmp[j] - alpha * weightDiff[j]; 
							//newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]/((TYPE)distanceElemNumber[j]); 
							newval = weightTmp[j] - alpha * percent * minweightval * weightDiff[j]; 
	
							modifiedRadius = modifyWeights(	
								weightPerturbed,
								neighRadius,
								neighSide,
								distanceValue,
								distanceIndex,
								distCounter,
								params,
								1,
								distIn,
								newval,
								NULL,
								NULL,
								i+1+j);
	
						}
						break;
					}
					else{
						printf("LINE-SEARCH END: error, not local minimum found in line-search.\n");
						//counter=0;
						//alpha=alphaInitial*0.75;
					}
				}


				for(j=0;j<distCounter;j++){

					distIn = distanceValue[j];
					newval = weightTmp[j]; 

					modifiedRadius = modifyWeights(	
						weightPerturbed,
						neighRadius,
						neighSide,
						distanceValue,
						distanceIndex,
						distCounter,
						params,
						1,
						distIn,
						newval,
						NULL,
						NULL,
						i+1+j);

				}				


				alpha=alpha*lambda;
				counter++;	
			}
		}
	}

	//return 1;

	//printf("starting cost function calculation...\n");
	//TYPE value = costFunction(nlags,expVariogram,npairs,imageZero,nx,ny,nz,bufferNodes);
	//printf("stoping cost function calculation...\n");

	//printf("%f\n",value);

	free(expVariogram);
	free(npairs);
	free(weight);
	free(weightPerturbed);
	freeRandomImage(image);

	free(line);

	return 0;
}
