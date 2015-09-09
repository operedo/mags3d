#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>
#include <sys/time.h>

int main(int argc,char *argv[]){

	int i,j,k,p,q;

	char *line;
	int BUFSIZE=100;
	ssize_t read;
	float elapsed_time;
	struct timeval start_time, end_time; 

	line=(char *)malloc(sizeof(char)*BUFSIZE + 1);

	/* gamv variables */
	int nd, irepo, maxdat, MAXVAR;
	float *x=NULL, *y=NULL, *z=NULL;
	float EPSLON;
	int nlag;
	float xlag, xltol;
	int mxdlv;
	double *np=NULL, *dis=NULL, *gam=NULL, *hm=NULL, *tm=NULL, *hv=NULL, *tv=NULL;
	int numThreads;
	double *reducedVariables;
	float dismxs, tmax, tmin;
	int ndir, nvarg;
	float *uvxazm=NULL, *uvyazm=NULL, *uvzdec=NULL, *uvhdec=NULL;
	float *csatol=NULL, *csdtol=NULL, *bandwh=NULL, *bandwd=NULL;
	float *atol=NULL;
	int *ivtype=NULL, *ivtail=NULL, *ivhead=NULL;
	float *vr=NULL;
	int isill;
	double *sills=NULL;
	char paramsFile[100];
	char variogramFile[100];




	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int nz = atoi(argv[3]);
	TYPE xlo = (TYPE)atof(argv[4]);
	TYPE ylo = (TYPE)atof(argv[5]);
	TYPE zlo = (TYPE)atof(argv[6]);
	TYPE h = (TYPE)atof(argv[7]);
	TYPE a = (TYPE)atof(argv[8]);
	int generateTargetVariogram = atoi(argv[9]);
	int useNscore = atoi(argv[10]);


	int iters=atoi(argv[11]);
	int report=atoi(argv[12]);
	//int updateTemperature=(int)(0.05*((double)iters));
	int updateTemperature=atoi(argv[13]);
	int resetCounter=0;
	int averageCounter=0;
	TYPE modifiedRadius;
	TYPE tol=atof(argv[14]); // percentage of decrease w/r initial cost
	TYPE newval,oldval;
	TYPE temperatureInitial=atof(argv[15]);
	TYPE temperature=temperatureInitial;



	//strncpy(paramsFile,argv[11],100);
	strncpy(paramsFile,"gamv.par",100);

	//strncpy(variogramFile,argv[12],100);
	//printf("paramsFile=%s\n",paramsFile);


	numThreads=1;
#ifdef _OPENMP
#pragma omp parallel default(shared) 
{
	numThreads=omp_get_num_threads();
}
#endif	


	gamvReadParams(
		&nd, &irepo, &maxdat, &MAXVAR,
		&x, &y, &z,
		&EPSLON,
		&nlag,
		&xlag, &xltol,
		&mxdlv,
		&np, &dis, &gam, &hm, &tm, &hv, &tv,
		&numThreads,
		&reducedVariables,
		&dismxs, &tmax, &tmin,
		&ndir, &nvarg,
		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
		&csatol, &csdtol, &bandwh, &bandwd,
		&atol,
		&ivtype, &ivtail, &ivhead,
		&vr,
		paramsFile,
		&isill,
		&sills
	);


	//printf("Parameters loaded ok.\n");
	//printf("nd=%d\n",nd);



//	gamv(
//		&nd, &irepo, &maxdat, &MAXVAR,
//		&x, &y, &z,
//		&EPSLON,
//		&nlag,
//		&xlag, &xltol,
//		&mxdlv,
//		&np, &dis, &gam, &hm, &tm, &hv, &tv,
//		&numThreads,
//		&reducedVariables,
//		&dismxs, &tmax, &tmin,
//		&ndir, &nvarg,
//		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
//		&csatol, &csdtol, &bandwh, &bandwd,
//		&atol,
//		&ivtype, &ivtail, &ivhead,
//		&vr	
//	);





	int bufferNodes=ceil(a/h);
	int neighRadius=bufferNodes;
	int neighSide=2*neighRadius+1;

	int nxExt=nx+2*bufferNodes;
	int nyExt=ny+2*bufferNodes;
	int nzExt=nz+2*bufferNodes;


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

	TYPE distanceValue[MAXDISTVALS];
	int distanceIndex[MAXDISTVALS];

	for(i=0;i<MAXDISTVALS;i++){
		distanceValue[i]=-1.0;
		distanceIndex[i]=-1;
	}



	FILE *fpweight, *fptargetvariogram;
	fpweight=fopen("weightout.dat","r");


	if(generateTargetVariogram){
		int retTarget = genTargetVariogram(nx,ny,nz,xlo,ylo,zlo,h,a,imageZero,useNscore);
	}
	fptargetvariogram=fopen("targetvariogram.dat","r");
	/*
	else{
		if(useNscore){
			//fptargetvariogram=fopen("gamv_Cu_vertical_nscore.out","r");
			fptargetvariogram=fopen("gamv_Cu_omnihoriz_nscore.out","r");
		}
		else{
			//fptargetvariogram=fopen("gamv_Cu_vertical.out","r");
			//fptargetvariogram=fopen("gamv_Cu_omnihoriz.out","r");
			//fptargetvariogram=fopen("salida.out.log_c","r");
			fptargetvariogram=fopen("initialimage.nscore_kernelgauss_gamv","r");
		}
	}
	*/



	int nlags,num;
	if((fgets(line,BUFSIZE,fptargetvariogram))==0){
		fprintf(stderr,"ERROR: variogram nlags not specified in variogram file.\n");
		return 1;	
	}
	sscanf(line,"%d\n",&nlags);

	TYPE *expVariogram = (TYPE *)malloc(sizeof(TYPE)*nlags);
	int *npairs=(int *)malloc(sizeof(int)*nlags);
	TYPE *lagDistance = (TYPE *)malloc(sizeof(TYPE)*nlags);

	TYPE vgval,lagval;
	i=0;
	while((fgets(line,BUFSIZE,fptargetvariogram))!=0){
		sscanf(line,"%lf %d %lf",&vgval,&num,&lagval);
		expVariogram[i]=vgval;	
		npairs[i]=num;
		lagDistance[i]=lagval;	
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
	TYPE *weightTemp = (TYPE *)malloc(sizeof(TYPE)*nweight); 
	TYPE *weightPerturbed = (TYPE *)malloc(sizeof(TYPE)*nweight); 
	TYPE *weightPerturbedTemp = (TYPE *)malloc(sizeof(TYPE)*nweight); 

	int midside;
	TYPE val;

	midside=ceil(side/2);


	TYPE xcoord1,ycoord1,zcoord1;
	TYPE xcoord2,ycoord2,zcoord2,dist;
	getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);
	int distCounter=0;
	while((fgets(line,BUFSIZE,fpweight))!=0){
		sscanf(line,"%d %d %d %lf",&i,&j,&k,&val);
		weight[i + j*side + k*side*side]=val;
		getCoordinates(i-neighRadius,j-neighRadius,k-neighRadius,params,&xcoord2,&ycoord2,&zcoord2);
               	dist = sqrt(	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	       		(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	       		(zcoord1-zcoord2)*(zcoord1-zcoord2));
		if(existValue(distanceValue,MAXDISTVALS,dist)==-1){
			distanceValue[distCounter]=dist;
			distCounter++;
		}
	}

	quickSort(distanceValue,0,distCounter-1);


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


	fprintf(stderr,"weights loaded.\n");

	//int ret = genMovingAverageImage(imageZero,weight,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,0,100);
	int ret = genMovingAverageImageFast(imageZero,weight,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,0,100,&nd,&maxdat,&MAXVAR,&x,&y,&z,&vr);

	printGSLIB(nx,ny,nz,params,0,100,&vr,"initialimage.dat");

	nscore(&nd,&maxdat,&MAXVAR,&x,&y,&z,&vr);

	printGSLIB(nx,ny,nz,params,0,100,&vr,"initialimage.nscore");

	fprintf(stderr,"Moving average image generated.\n");

	irepo=100;

//	gamv(
//		&nd, &irepo, &maxdat, &MAXVAR,
//		&x, &y, &z,
//		&EPSLON,
//		&nlag,
//		&xlag, &xltol,
//		&mxdlv,
//		&np, &dis, &gam, &hm, &tm, &hv, &tv,
//		&numThreads,
//		&reducedVariables,
//		&dismxs, &tmax, &tmin,
//		&ndir, &nvarg,
//		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
//		&csatol, &csdtol, &bandwh, &bandwd,
//		&atol,
//		&ivtype, &ivtail, &ivhead,
//		&vr	
//	);

//	gamvCUDA(
//		&nd, &irepo, &maxdat, &MAXVAR,
//		&x, &y, &z,
//		&EPSLON,
//		&nlag,
//		&xlag, &xltol,
//		&mxdlv,
//		&np, &dis, &gam, &hm, &tm, &hv, &tv,
//		&numThreads,
//		&reducedVariables,
//		&dismxs, &tmax, &tmin,
//		&ndir, &nvarg,
//		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
//		&csatol, &csdtol, &bandwh, &bandwd,
//		&atol,
//		&ivtype, &ivtail, &ivhead,
//		&vr	
//	);

//	gamv(
	gamvOMP(
//	gamvCUDA(
//	gamvCUDAoptimized(
		&nd, &irepo, &maxdat, &MAXVAR,
		&x, &y, &z,
		&EPSLON,
		&nlag,
		&xlag, &xltol,
		&mxdlv,
		&np, &dis, &gam, &hm, &tm, &hv, &tv,
		&numThreads,
		&reducedVariables,
		&dismxs, &tmax, &tmin,
		&ndir, &nvarg,
		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
		&csatol, &csdtol, &bandwh, &bandwd,
		&atol,
		&ivtype, &ivtail, &ivhead,
		&vr	
	);


	gamvAverages(&ndir, &nvarg, &nlag, &isill, 
			&ivtype, &ivtail, &ivhead, &sills,
			&EPSLON,
			&np, &dis, &gam, &hm, &tm, &hv, &tv);
	int nsiz=(ndir)*(nvarg)*(nlag+2); 
//	for(i=0;i<nsiz;i++){
//		printf("%d %f %f %d %f %f\n",(i+1),dis[i],gam[i],(int)(np[i]),tm[i],hm[i]);
//	}
	FILE *fpinitialvariogram;
	fpinitialvariogram = fopen("initialvariogram.dat","w");
	fprintf(fpinitialvariogram,"%d\n",nsiz);
	for(i=0;i<nsiz;i++){
		fprintf(fpinitialvariogram,"%f %d %f\n",gam[i],(int)(np[i]),dis[i]);
	}
	fclose(fpinitialvariogram);

	
	fprintf(stderr,"starting cost function calculation...\n");
	//TYPE cost = costFunctionGSLIB(nlags,expVariogram,npairs,0,useNscore);
	TYPE cost = costFunctionFast(nsiz,gam,expVariogram,0);
	fprintf(stderr,"stoping cost function calculation...\n");

	fprintf(stderr,"initial cost=%f\n",cost);

	TYPE initialCost = cost;
	TYPE tempCost=0.0;
	TYPE currentCost= cost;
	TYPE currentCostPrev= 0.0;

/*
	//int iters=0;
	int iters=6000;
	//int iters=100;
	//int iters=2000;
	int report=100;
	//int updateTemperature=100;
	//int updateTemperature=750;
	int updateTemperature=(int)(0.05*((double)iters));
	int resetCounter=0;
	int averageCounter=0;
	TYPE modifiedRadius;
	TYPE tol=0.000001; // percentage of decrease w/r initial cost
	TYPE newval,oldval;
	//TYPE temperatureInitial=0.000075;
	//TYPE temperatureInitial=0.0003;
	//TYPE temperatureInitial=0.01; // 2D gaussian-wrong to gaussian works
	//TYPE temperatureInitial=0.005; //  2D gaussian to circular works
	TYPE temperatureInitial=0.0001; //  2D gaussian to circular works
	//TYPE temperatureInitial=0.001; //  2D circular to gaussian works
	//TYPE temperatureInitial=0.01; //  real data works
	//TYPE temperatureInitial=0.01; // 
	//TYPE temperatureInitial=0.05; // 
	TYPE temperature=temperatureInitial;
*/

	TYPE probability, ran;

	int averageFlag=0;

	gettimeofday( &start_time, NULL ); 
	//for(i=1;i<=1;i++){
	for(i=1;i<=iters;i++){
		if((currentCost/initialCost)<=tol){
			fprintf(stdout,"STOP: convergence achieved. Congratulations!\n");

			modifiedRadius = modifyWeights(	
						weight,
						neighRadius,
						neighSide,
						distanceValue,
						distanceIndex,
						distCounter,
						params,
						2,
						0,
						0.0,
						&newval,
						&oldval,
						i);

			return 0;
		}
		else{
			//if(i%50==0)
			fprintf(stdout,"ITERATION %d:\t",i);

			// smoothing of weights
//
//			if(i%(iters/10)==0 || i==iters){
//
//				modifiedRadius = averageWeights(	
//						weightPerturbed,
//						weight,
//						neighRadius,
//						neighSide,
//						distanceValue,
//						distanceIndex,
//						distCounter,
//						params,
//						0,
//						0,
//						0.0,
//						&newval,
//						&oldval,
//						i);
//
//				ret = genMovingAverageImage(imageZero,weightPerturbed,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,i,report);
//				averageFlag=1;
//			}

			modifiedRadius = modifyWeights(	
						weightPerturbed,
						neighRadius,
						neighSide,
						distanceValue,
						distanceIndex,
						distCounter,
						params,
						0,
						0,
						0.0,
						&newval,
						&oldval,
						i);

			//ret = genMovingAverageImage(imageZero,weightPerturbed,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,i,report);
			//tempCost = costFunctionGSLIB(nlags,expVariogram,npairs,i,useNscore);


			ret = genMovingAverageImageFast(imageZero,weightPerturbed,nx,ny,nz,nxExt,nyExt,nzExt,neighRadius,neighSide,params,i,report,&nd,&maxdat,&MAXVAR,&x,&y,&z,&vr);
	
			printGSLIB(nx,ny,nz,params,i,report,&vr,"currentimage.dat");

			nscore(&nd,&maxdat,&MAXVAR,&x,&y,&z,&vr);

			printGSLIB(nx,ny,nz,params,i,report,&vr,"currentimage.nscore");


//			gamvCUDA(
//				&nd, &irepo, &maxdat, &MAXVAR,
//				&x, &y, &z,
//				&EPSLON,
//				&nlag,
//				&xlag, &xltol,
//				&mxdlv,
//				&np, &dis, &gam, &hm, &tm, &hv, &tv,
//				&numThreads,
//				&reducedVariables,
//				&dismxs, &tmax, &tmin,
//				&ndir, &nvarg,
//				&uvxazm, &uvyazm, &uvzdec, &uvhdec,
//				&csatol, &csdtol, &bandwh, &bandwd,
//				&atol,
//				&ivtype, &ivtail, &ivhead,
//				&vr	
//			);

//			gamv(
			gamvOMP(
//			gamvCUDA(
//			gamvCUDAoptimized(
				&nd, &irepo, &maxdat, &MAXVAR,
				&x, &y, &z,
				&EPSLON,
				&nlag,
				&xlag, &xltol,
				&mxdlv,
				&np, &dis, &gam, &hm, &tm, &hv, &tv,
				&numThreads,
				&reducedVariables,
				&dismxs, &tmax, &tmin,
				&ndir, &nvarg,
				&uvxazm, &uvyazm, &uvzdec, &uvhdec,
				&csatol, &csdtol, &bandwh, &bandwd,
				&atol,
				&ivtype, &ivtail, &ivhead,
				&vr	
			);

		
		
			gamvAverages(&ndir, &nvarg, &nlag, &isill, 
					&ivtype, &ivtail, &ivhead, &sills,
					&EPSLON,
					&np, &dis, &gam, &hm, &tm, &hv, &tv);



			tempCost = costFunctionFast(nsiz,gam,expVariogram,0);

			//printf("tempCost=%f\t(%1.6f)\t",tempCost,tempCost/initialCost);
			
			//if(i%200==0){
			if(i%updateTemperature==0){
				//temperature=temperature*0.1;
				temperature=temperature*0.5;
				//temperature=temperature*0.9;
				resetCounter++;
				//if(resetCounter==3){
				if(resetCounter==5){
					//temperature=temperatureInitial*0.9;
					temperature=temperatureInitial;
					resetCounter=0;
				}
			}

			probability = tempCost<currentCost?1.0:exp((currentCost-tempCost)/(currentCost*temperature)); 
			ran=((TYPE)rand()/(TYPE)RAND_MAX);

			//fprintf(stderr,"probability=%f>%f?%d (%f)\n",probability,ran,probability>ran,(currentCost-tempCost)/(temperature) );

			//if(tempCost<currentCost){
			if(probability>ran){
				//accept modification
				//weight[modifiedIndex] = weightPerturbed[modifiedIndex];
				//if(i%50==0)
				if(averageFlag==0){
					fprintf(stdout,"accepted(prob=%f,r=%f,new=%f,old=%f)\ttempCost=%f\tcurrCost=%f\t%f\n",probability,modifiedRadius,newval,oldval,tempCost,currentCost,currentCost/initialCost);
				}
				else{
					fprintf(stdout,"accepted(prob=%f,r=%f,new=%f,old=%f)\ttempCost=%f\tcurrCost=%f\t%f\t(averaged)\n",probability,modifiedRadius,newval,oldval,tempCost,currentCost,currentCost/initialCost);
					averageFlag=0;
				}
				ret = modifyWeights(
							weight,
							neighRadius,
							neighSide,
							distanceValue,
							distanceIndex,
							distCounter,
							params,
							1,
							modifiedRadius,
							newval,
							NULL,
							NULL,
							i);
				currentCost = tempCost;
				if(i%report==0){
					int sysret = system("cp currentdistanceweights.dat currentbestdistanceweights.dat");

					//gamvAverages(&ndir, &nvarg, &nlag, &isill, 
					//	&ivtype, &ivtail, &ivhead, &sills,
					//	&EPSLON,
					//	&np, &dis, &gam, &hm, &tm, &hv, &tv);

					FILE *fpcurrentvariogram;
					fpcurrentvariogram = fopen("currentvariogram.dat","w");
					nsiz=(ndir)*(nvarg)*(nlag+2); 
					fprintf(fpcurrentvariogram,"%d\n",nsiz);
					int j;
					for(j=0;j<nsiz;j++){
						fprintf(fpcurrentvariogram,"%f %d %f\n",gam[j],(int)(np[j]),dis[j]);
					}
					fclose(fpcurrentvariogram);
					sysret = system("cp currentvariogram.dat currentbestvariogram.dat");
				}
			}
			else{
				//reject modification
				//weightPerturbed[modifiedIndex] = weight[modifiedIndex];
				//if(i%50==0)
				if(averageFlag==0){
					fprintf(stdout,"rejected(prob=%f,r=%f,new=%f,old=%f)\ttempCost=%f\tcurrCost=%f\t%f\n",probability,modifiedRadius,newval,oldval,tempCost,currentCost,currentCost/initialCost);
				}
				else{
					fprintf(stdout,"rejected(prob=%f,r=%f,new=%f,old=%f)\ttempCost=%f\tcurrCost=%f\t%f\t(averaged)\n",probability,modifiedRadius,newval,oldval,tempCost,currentCost,currentCost/initialCost);
					averageFlag=0;
				}
				ret = modifyWeights(
							weightPerturbed,
							neighRadius,
							neighSide,
							distanceValue,
							distanceIndex,
							distCounter,
							params,
							1,
							modifiedRadius,
							oldval,
							NULL,
							NULL,
							i);
			}
		}
	}

	gettimeofday( &end_time, NULL );
        elapsed_time = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec ) / 1e6; 
	printf("elapsed_time=%f\n",elapsed_time/3.0f);


	


	gamvFreeMemory(
		&x, &y, &z,
		&np, &dis, &gam, &hm, &tm, &hv, &tv,
		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
		&csatol, &csdtol, &bandwh, &bandwd,
		&atol,
		&ivtype, &ivtail, &ivhead,
		&vr,
		&sills,
		&reducedVariables,
		&numThreads
		);



	free(expVariogram);
	free(lagDistance);
	free(npairs);
	free(weight);
	free(weightTemp);
	free(weightPerturbed);
	free(weightPerturbedTemp);
	freeRandomImage(image);


	free(line);

	return 0;
}
