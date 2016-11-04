#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mags3d.h>

TYPE costFunction(int nlags, TYPE *expVariogram, int *npairs,TYPE *image,int nx,int ny,int nz,int bufferNodes){
	TYPE valueTotal=0.0;
	TYPE valuePerLag=0.0;
	TYPE valueWeight=0.0;
	TYPE delta;
	int i,j,k,ii,jj,kk;
	srand(time(NULL));

	printf("starting nlags loop.\n");

	for(i=0;i<nlags;i++){
	   valuePerLag=0.0;
           printf("calculating lag=%d, npairs[%d]=%d...\n",i+1,i,npairs[i]);
	   for(j=0;j<npairs[i];j++){
              //printf("calculating pairs %d...\n",j);
              valueWeight=0.0;
	      for(k=0;k<nweight;k++){
                 //printf("%f\t%f\n",weight[k],(TYPE)((rand()%3)-1));

                 /* we will use horizontal lags*/
                 ii=rand()%(nx-(i+1));
                 jj=rand()%ny;
                 kk=rand()%nz;

                 //printf("%d %d %d\n",ii,jj,kk);

                 delta = *(image+ii+jj*(nx+2*bufferNodes)+kk*(nx+2*bufferNodes)*(ny+2*bufferNodes));
                 delta = delta - *(image+ii+(i+1)+jj*(nx+2*bufferNodes)+kk*(nx+2*bufferNodes)*(ny+2*bufferNodes));

                 //printf("%d %d: %f\n",j,k,delta);

	         //valueWeight = valueWeight + weight[k]*(TYPE)((rand()%3)-1);
	         valueWeight = valueWeight + weight[k]*delta;
	      }
              valuePerLag = valuePerLag + valueWeight*valueWeight;
              //printf("valuePerLag=%f\n",valuePerLag);
	   }
           valuePerLag = valuePerLag / (2.0*(TYPE)(npairs[i]));
           printf("%f\texpVariogram=%f\n",valuePerLag,expVariogram[i]);
	   valueTotal = valueTotal + (valuePerLag - expVariogram[i])*(valuePerLag - expVariogram[i]);
	}
	return valueTotal;
}
