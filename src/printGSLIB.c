#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>

#define LINESIZE 50
//#define NUMLINES 3000
//#define NUMLINES 256
#define NUMLINES 20
#define BUFSIZE NUMLINES*LINESIZE

int printGSLIB(	
				int nx, 
				int ny, 
				int nz, 
				latticeParams params,
				int iter,
				int report,
				//float **x, 
				//float **y, 
				//float **z,
				float **vr,
				char *filename
			)
{
	int i,j,k;
	TYPE xcoord1,ycoord1,zcoord1;

	FILE *fpgslib;
	//FILE *fpoctave;
	if(iter==0){
		//fpoctave = fopen("initialimage.m","w");
		fpgslib  = fopen(filename,"w");
		fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");
		//fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
	}
	else if(iter%report==0){
		//fpoctave = fopen("currentimage.m","w");
		fpgslib  = fopen(filename,"w");
		fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");
		//fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
	}

	if(iter==0 || iter%report==0){
		//printf("printing...\n");
		for(k=0;k<nz;k++){
	  		for(j=0;j<ny;j++){
	    			for(i=0;i<nx;i++){
              				getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
					fprintf(fpgslib,"%3.6f %3.6f %3.6f %3.6f\n",xcoord1,ycoord1,zcoord1,(*vr)[i+j*nx+k*(nx*ny)]);
              				//fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,(*vr)[i+j*nx+k*(nx*ny)]);
				}
			}
		}
		
		fclose(fpgslib);
		//fclose(fpoctave);
	}


	return 0;
}
