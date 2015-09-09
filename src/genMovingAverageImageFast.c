#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>

#define LINESIZE 50
//#define NUMLINES 3000
//#define NUMLINES 256
#define NUMLINES 20
#define BUFSIZE NUMLINES*LINESIZE

int genMovingAverageImageFast(	
				TYPE *imageZero, 
				TYPE* kernelWeight, 
				int nx, 
				int ny, 
				int nz, 
				int nxExt, 
				int nyExt, 
				int nzExt,
				int neighRadius, 
				int neighSide, 
				latticeParams params,
				int iter,
				int report,
				int *nd, 
				int *maxdat,
				int *MAXVAR,
				float **x, 
				float **y, 
				float **z,
				float **vr
			)
{
	int i,j,k,ni,nj,nk;
	TYPE sum=0.0;
	TYPE xcoord1,ycoord1,zcoord1;
	TYPE xcoord2,ycoord2,zcoord2;
	TYPE weight;

	int threadId, numThreads;
	int totalNeighbours=(2*neighRadius+1)*(2*neighRadius+1)*(2*neighRadius+1);  

	char line[LINESIZE];
	char buffer[BUFSIZE];
	char filename[LINESIZE];
	char filenameOct[LINESIZE];
	int counter=0;
			
	threadId=0;
#ifdef _OPENMP
#pragma omp parallel private(threadId)
{
	threadId=omp_get_thread_num();
	numThreads=omp_get_num_threads();
}
#endif

/*
	FILE *fpgslib;
	FILE *fpoctave;
	if(iter==0){
		fpoctave = fopen("initialimage.m","w");
		fpgslib  = fopen("initialimage.dat","w");
		fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");
		fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
	}
	else if(iter%report==0){
		fpoctave = fopen("currentimage.m","w");
		fpgslib  = fopen("currentimage.dat","w");
		fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");
		fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
	}
*/

	if(iter==0){
		free(*x);
		free(*y);
		free(*z);
		free(*vr);
		*maxdat = nx*ny*nz;
		*nd= *maxdat;
		*x=malloc(*maxdat * sizeof(**x));
		*y=malloc(*maxdat * sizeof(**y));
		*z=malloc(*maxdat * sizeof(**z));
		*vr=malloc(*maxdat * *MAXVAR * sizeof(**vr));
	}


/*
	if(iter==0){
#ifndef _OPENMP
		fpoctave = fopen("initialimage.m","w");
		fpgslib  = fopen("initialimage.dat","w");
		fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");
		fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
#else
		for(i=0;i<numThreads;i++){
			sprintf(filename,"initialimage0%d.dat",i);
			fpgslibThreads[i]  = fopen(filename,"w");
			sprintf(filenameOct,"initialimage0%d.m",i);
			fpoctaveThreads[i]  = fopen(filenameOct,"w");
		}
		fprintf(fpgslibThreads[0],"Test\n4\nX\nY\nZ\ndata\n");
		fprintf(fpoctaveThreads[0],"image=zeros(%d,%d,%d);\n",nx,ny,nz);
#endif
		//fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
	}
	else{
#ifndef _OPENMP
		if(iter%report==0){
			fpoctave = fopen("currentimage.m","w");
			fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
		}
		fpgslib  = fopen("currentimage.dat","w");
		fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");
#else
		for(i=0;i<numThreads;i++){
			sprintf(filename,"currentimage0%d.dat",i);
			fpgslibThreads[i]  = fopen(filename,"w");
			if(iter%report==0){
				sprintf(filenameOct,"currentimage0%d.m",i);
				fpoctaveThreads[i]  = fopen(filenameOct,"w");
			}
		}
		fprintf(fpgslibThreads[0],"Test\n4\nX\nY\nZ\ndata\n");
		if(iter%report==0)
			fprintf(fpoctaveThreads[0],"image=zeros(%d,%d,%d);\n",nx,ny,nz);
#endif
	}
*/
#pragma omp parallel default(shared) private(threadId,k,j,i,sum,xcoord1,ycoord1,zcoord1,nk,nj,ni,weight,line,buffer,counter) 
{

#ifdef _OPENMP
	threadId=omp_get_thread_num();
#else
	threadId=0;
#endif

	//memset(line, 0, strlen(line));
	//memset(buffer, 0, strlen(buffer));

#pragma omp for 
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
              	
              //printf("%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
              (*x)[i+j*nx+k*(nx*ny)]=xcoord1;
              (*y)[i+j*nx+k*(nx*ny)]=ycoord1;
              (*z)[i+j*nx+k*(nx*ny)]=zcoord1;
              (*vr)[i+j*nx+k*(nx*ny)]=sum/((float)totalNeighbours);
	    }
	  }
	}

/*
	if(counter!=0){
#ifdef _OPENMP
        	//fprintf(fpgslibThreads[threadId],"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
		fwrite(buffer,1,strlen(buffer),fpgslibThreads[threadId]);
#else
              	//fprintf(fpgslib,"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
		fwrite(buffer,1,strlen(buffer),fpgslib);
#endif
	}
*/

}



/*
#ifdef _OPENMP
	for(i=0;i<numThreads;i++)
		fclose(fpgslibThreads[i]);
	if(iter%report==0){
		for(i=0;i<numThreads;i++)
			fclose(fpoctaveThreads[i]);
	}
	if(iter==0){
		system("rm initialimage.dat; cat initialimage*.dat > initialimage.dat");
		system("rm initialimage.m; cat initialimage*.m > initialimage.m");
	}
	else{
		system("rm currentimage.dat; cat currentimage*.dat > currentimage.dat");
		if(iter%report==0)
			system("rm currentimage.m; cat currentimage*.m > currentimage.m");
	}
#else
	fclose(fpgslib);
	if(iter%report==0){
		fclose(fpoctave);
	}
#endif
*/


/*
	if(iter==0 || iter%report==0){
		for(k=0;k<nz;k++){
	  		for(j=0;j<ny;j++){
	    			for(i=0;i<nx;i++){
              				getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
					fprintf(fpgslib,"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,(*vr)[i+j*nx+k*(nx*ny)]);
              				fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,(*vr)[i+j*nx+k*(nx*ny)]);
				}
			}
		}
		
		fclose(fpgslib);
		fclose(fpoctave);
	}
*/

	return 0;
}
