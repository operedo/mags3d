#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>

TYPE * genRandomImage(int xsize, int ysize, int zsize){
	int i,j,k;
	int row,col,slide;

	TYPE *image;	

	image=(TYPE *)malloc(sizeof(TYPE)*xsize*ysize*zsize);

	srand(time(NULL));
	//srand(1);
	for(k=0;k<zsize;k++){
		for(i=0;i<xsize;i++){
			for(j=0;j<ysize;j++){

				image[k*(xsize*ysize) + i*ysize + j]= generateGaussianNoise(); // N(0,1) random
				//image[k*(xsize*ysize) + i*ysize + j]= (TYPE)( 2.0*((TYPE)rand())/((TYPE)RAND_MAX) - 1.0 ); // uniform [0,1] random
				//image[k*(xsize*ysize) + i*ysize + j]= (TYPE)(((TYPE)rand())/((TYPE)RAND_MAX)); // uniform [0,1] random
				//image[k*(xsize*ysize) + i*ysize + j]= (TYPE)(rand()%2); // {0,1} random
				//image[k*(xsize*ysize) + i*ysize + j]= (TYPE)((i*j+k)%2); // non random
				//printf("image[%d]=%f\n",k*(xsize*ysize) + i*ysize + j,image[k*(xsize*ysize) + i*ysize + j]);
			}
		}
	}
	return image;
}

int freeRandomImage(TYPE *image){
	free(image);
}



double generateGaussianNoise()
{
	static int haveSpare = 0;
	static double rand1, rand2;
 
	if(haveSpare)
	{
		haveSpare = 0;
		return sqrt(rand1) * sin(rand2);
	}
 
	haveSpare = 1;
 
	rand1 = rand() / ((double) RAND_MAX);
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;
 
	return sqrt(rand1) * cos(rand2);
}

