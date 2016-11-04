
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mags3d.h>

#define LINESIZE 50
//#define NUMLINES 3000
//#define NUMLINES 256
#define NUMLINES 20
#define BUFSIZE NUMLINES*LINESIZE

int genMovingAverageImage(	
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
				int report
			)
{
	int i,j,k,ni,nj,nk;
	TYPE sum=0.0;
	TYPE xcoord1,ycoord1,zcoord1;
	TYPE xcoord2,ycoord2,zcoord2;
	TYPE weight;

	int threadId, numThreads;

	char line[LINESIZE];
	char buffer[BUFSIZE];
	char filename[LINESIZE];
	char filenameOct[LINESIZE];
	int counter=0;
			
#ifdef _OPENMP
#pragma omp parallel
{
	numThreads=omp_get_num_threads();
}
#endif

#ifdef _OPENMP
	FILE *fpgslibThreads[numThreads];
	FILE *fpoctaveThreads[numThreads];
#else
	FILE *fpgslib;
	FILE *fpoctave;
#endif

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

#pragma omp parallel default(shared) private(threadId,k,j,i,sum,xcoord1,ycoord1,zcoord1,nk,nj,ni,weight,line,buffer,counter) 
{

#ifdef _OPENMP
	threadId=omp_get_thread_num();
#else
	threadId=0;
#endif

	memset(line, 0, strlen(line));
	memset(buffer, 0, strlen(buffer));

#pragma omp for 
	for(k=0;k<nz;k++){
	  for(j=0;j<ny;j++){
	    for(i=0;i<nx;i++){
              sum=0.0;
              getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
	      for(nk=-neighRadius;nk<=neighRadius;nk++){
	        for(nj=-neighRadius;nj<=neighRadius;nj++){
	          for(ni=-neighRadius;ni<=neighRadius;ni++){
                    //getCoordinates(i+ni,j+nj,k+nk,params,&xcoord2,&ycoord2,&zcoord2);
                    //rsquared = 	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
		    //		(ycoord1-ycoord2)*(ycoord1-ycoord2)+
		    //		(zcoord1-zcoord2)*(zcoord1-zcoord2);
                    ////weight = cons * exp(-2.0*rsquared/((TYPE)(neighRadius*neighRadius)));
                    //weight = cons * exp(-2.0*rsquared*asquaredinv);
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
	      	//sum= *(imageZero+(
	  	//(i)+
	  	//(j)*(nxExt)+
	  	//(k)*(nxExt)*(nyExt)
	  	//)
	  	//);

		//if(iter%report==0){
              	//	fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
		//}
              	//fprintf(fpgslib,"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
              	
		sprintf(line,"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
		strcat(buffer,line);
		if(counter==NUMLINES-1){
#ifdef _OPENMP
              		//fprintf(fpgslibThreads[threadId],"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
			fwrite(buffer,1,strlen(buffer),fpgslibThreads[threadId]);
#else
              		//fprintf(fpgslib,"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
			fwrite(buffer,1,strlen(buffer),fpgslib);
#endif
			memset(buffer, 0, strlen(buffer));
		}
		counter=(counter+1)%NUMLINES;

#ifdef _OPENMP
		if(iter%report==0){
              		fprintf(fpoctaveThreads[threadId],"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
		}
#else
		if(iter%report==0){
              		fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
		}
#endif
	    }
	  }
	}

	if(counter!=0){
#ifdef _OPENMP
        	//fprintf(fpgslibThreads[threadId],"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
		fwrite(buffer,1,strlen(buffer),fpgslibThreads[threadId]);
#else
              	//fprintf(fpgslib,"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,sum);
		fwrite(buffer,1,strlen(buffer),fpgslib);
#endif
	}

}




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
	return 0;
}
