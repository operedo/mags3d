#include <stdio.h>
#include <stdlib.h>
#include <mags3d.h>

int main(){
	int sizex=100,sizey=100,sizez=100;
	TYPE *image=NULL;
	image = genRandomImage(sizex,sizey,sizez);
	if(image==NULL){
		printf("testGenRandomImage: FAILED\n");
	}
	else{
		printf("testGenRandomImage: PASSED\n");
		printf("\tsizex=%d, sizey=%d, sizez=%d\n",sizex,sizey,sizez);
		int i,sum0=0,sum1=0,total=sizex*sizey*sizez;
		for(i=0;i<total;i++){
			if(image[i]==0.0) sum0++;
			if(image[i]==1.0) sum1++;
		}
		printf("\tnumber of 0s=%d\n",sum0);	
		printf("\tnumber of 1s=%d\n",sum1);	
		printf("\n");
	}
	freeRandomImage(image);
}
