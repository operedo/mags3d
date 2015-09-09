#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mags3d.h>

TYPE costFunctionGSLIB(int nlags, TYPE *expVariogram, int *npairs,int iter,int useNscore){

	TYPE valueTotal=0.0;
	char *line;
	int BUFSIZE=100;
	ssize_t read;
	line=(char *)malloc(sizeof(char)*BUFSIZE + 1);

	int sysret;

	if(iter==0){

		if(useNscore){
			sysret=system("../gslib90/nscore initialnscore.par > /dev/null 2>&1");
			//sysret=system("../gslib90/gam-io initialgam128x128x1_nscore.par > /dev/null 2>&1");
			//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/initialgamv_Cu_vertical_nscore.par > /dev/null 2>&1");
			sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/initialgamv_Cu_omnihoriz_nscore.par > /dev/null 2>&1");

		}
		else{
			//sysret=system("../gslib90/gam-io initialgam128x128x1.par > /dev/null 2>&1");
			//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/initialgamv_Cu_vertical.par > /dev/null 2>&1");
			sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/initialgamv_Cu_omnihoriz.par > /dev/null 2>&1");
		}

		//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/initialgamv_Cu_vertical.par > /dev/null 2>&1");
		//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/initialgamv_Cu_omnihoriz.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gamv 3Dmuestras/initialgamv_Cu_vertical.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gamv initialgamv.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io initialgam.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io initialgam64x64x1.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io initialgam128x128x1.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io initialgam256x256x1.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io initialgam128x128x32.par > /dev/null 2>&1");
		sysret=system("tail -n +2 initialvariogram.out | awk '{ print $3, $4, $2 }' > initialvariogram.dat");
	}
	else{
		if(useNscore){
			sysret=system("../gslib90/nscore currentnscore.par > /dev/null 2>&1");
			//sysret=system("../gslib90/gam-io currentgam128x128x1_nscore.par > /dev/null 2>&1");
			//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/currentgamv_Cu_vertical_nscore.par > /dev/null 2>&1");
			sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/currentgamv_Cu_omnihoriz_nscore.par > /dev/null 2>&1");
		}
		else{
			//sysret=system("../gslib90/gam-io currentgam128x128x1.par > /dev/null 2>&1");
			//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/currentgamv_Cu_vertical.par > /dev/null 2>&1");
			sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/currentgamv_Cu_omnihoriz.par > /dev/null 2>&1");
		}
		//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/currentgamv_Cu_vertical.par > /dev/null 2>&1");
		//sysret=system("export OMP_NUM_THREADS=4; export OMP_SCHEDULE=\"static,32\";../../pargamv/gamv.exe 3Dmuestras/currentgamv_Cu_omnihoriz.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gamv 3Dmuestras/currentgamv_Cu_vertical.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gamv currentgamv.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io currentgam64x64x1.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io currentgam128x128x1.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io currentgam256x256x1.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io currentgam128x128x32.par > /dev/null 2>&1");
		sysret=system("tail -n +2 currentvariogram.out | awk '{ print $3, $4, $2 }' > currentvariogram.dat");
	}

	FILE *fpcurrentvariogram;
	if(iter==0){
		fpcurrentvariogram=fopen("initialvariogram.dat","r");
	}
	else{
		fpcurrentvariogram=fopen("currentvariogram.dat","r");
	}

	int num,i;


	TYPE *currentVariogram = (TYPE *)malloc(sizeof(TYPE)*nlags);
	int *currentNpairs=(int *)malloc(sizeof(int)*nlags);
	TYPE *currentLags = (TYPE *)malloc(sizeof(TYPE)*nlags);

	TYPE vgval,lagval;
	i=0;



	while((fgets(line,BUFSIZE,fpcurrentvariogram))!=0){
		sscanf(line,"%lf %d %lf",&vgval,&num,&lagval);
		currentVariogram[i]=vgval;	
		currentNpairs[i]=num;
		currentLags[i]=lagval;
		//printf("%f %d %f\n",vgval,num,lagval);

		//valueTotal = valueTotal + (log(currentVariogram[i]+2.0) - log(expVariogram[i]+2.0))*(log(currentVariogram[i]+2.0) - log(expVariogram[i]+2.0))/(log(expVariogram[i]+2.0)*log(expVariogram[i]+2.0)); 
		valueTotal = valueTotal + (currentVariogram[i] - expVariogram[i])*(currentVariogram[i] - expVariogram[i])/((expVariogram[i]+1.0)*(expVariogram[i]+1.0)); 

		i++;	
	}

	free(currentVariogram);
	free(currentNpairs);
	free(currentLags);

	free(line);
	fclose(fpcurrentvariogram);

	return valueTotal;
}
