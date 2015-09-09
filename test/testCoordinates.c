#include <stdio.h>
#include <stdlib.h>
#include <mags3d.h>

int main(){

	latticeParams params;
	params.nx=100;
	params.ny=100;
	params.nz=100;

	params.xlo=0.0;
	params.ylo=0.0;
	params.zlo=0.0;
	
	params.h=1.0;

	TYPE xcoord, ycoord, zcoord;

	int retcode;
	int index;
	//for(index=0;index<params.nx*params.ny*params.nz;index++){
	for(index=0;index<20;index++){
		retcode = getCoordinates(index,params,&xcoord,&ycoord,&zcoord);
		if(retcode==1){
			printf("testCoordinates: FAILED\n");
		}
		else{
			printf("testCoordinates: PASSED\n");
			printf("\tindex=%d, xcoord=%f, ycoord=%f, zcoord=%f\n",index,xcoord,ycoord,zcoord);
		}
	}
	printf("\n");
}
