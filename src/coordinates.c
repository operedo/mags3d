#include <mags3d.h>

//int getCoordinates(int index, latticeParams params, TYPE *xcoord, TYPE *ycoord, TYPE *zcoord){
//	if(index>=params.nx*params.ny*params.nz || index<0) return 1;
//	else{
//		*xcoord = params.xlo + (TYPE)(index%params.nx)*params.h;
//		*ycoord = params.ylo + (TYPE)(((int)(index/params.nx))%params.ny)*params.h;
//		*zcoord = params.zlo + (TYPE)((int)(index/(params.nx*params.ny))%params.nz)*params.h;
//		return 0;
//	}
//}

int getCoordinates(int i, int j, int k, latticeParams params, TYPE *xcoord, TYPE *ycoord, TYPE *zcoord){
	//if(index>=params.nx*params.ny*params.nz || index<0) return 1;
	//else{
		*xcoord = params.xlo + (TYPE)(i)*params.h;
		*ycoord = params.ylo + (TYPE)(j)*params.h;
		*zcoord = params.zlo + (TYPE)(k)*params.h;
		return 0;
	//}
}
