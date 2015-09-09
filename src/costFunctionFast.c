#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mags3d.h>

TYPE costFunctionFast(int nsiz, double *gam, TYPE *expVariogram, int iter){
	int i;
	TYPE valueTotal=0.0;

	for(i=0;i<nsiz;i++){
		valueTotal = valueTotal + ((TYPE)(gam[i]) - expVariogram[i])*((TYPE)(gam[i]) - expVariogram[i])/((expVariogram[i]+1.0)*(expVariogram[i]+1.0)); 
	}
	return valueTotal;
}

