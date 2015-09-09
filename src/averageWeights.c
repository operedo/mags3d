#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>


TYPE averageWeights(	
				TYPE* kernelWeight, 
				TYPE* kernelWeightOrig, 
				int neighRadius, 
				int neighSide,
				TYPE* distanceValue,
				int* distanceIndex,
				int distanceSize,
				latticeParams params,
				int mode,
				TYPE distIn,
				TYPE newvalIn,
				TYPE *newvalOut,
				TYPE *oldvalOut,
				int iter
			)
{
	int i,j,k,ni,nj,nk,rani,ranj,rank,index,radi,radj,radk;
	int distIndex,ldistIndex,udistIndex;
	TYPE lvalue,uvalue;
	TYPE xcoord1,ycoord1,zcoord1;
	TYPE xcoord2,ycoord2,zcoord2;
	TYPE dist=0.0;
	TYPE newval=0.0, val=0.0;


	TYPE avgWeightFwd[distanceSize], avgWeightBwd[distanceSize];
	TYPE medianWeightFwd[distanceSize], medianWeightBwd[distanceSize];
	TYPE emaWeightFwd[distanceSize];
	//TYPE alphaExp=0.8;
	TYPE alphaExp=0.6;


	TYPE avg=0.0, maxval=-1.0e+20, minval=1.0e+20;
	TYPE count=0.0;

	int counter=0;

	//int windowSize=ceil((((TYPE)(distanceSize))*0.068)); 
	int windowSize=ceil((((TYPE)(distanceSize))*0.08)); 
	//int windowSize=ceil((((TYPE)(distanceSize))*0.1)); 
	// fast forward average pass
	//for(i=windowSize;i<distanceSize;i++){
	for(i=0;i<distanceSize;i++){

		if(i==0){
			//emaWeightFwd[i]=kernelWeight[distanceIndex[i]]; 
			emaWeightFwd[i]=kernelWeightOrig[distanceIndex[i]]; 
		}
		else{
			//emaWeightFwd[i]=alphaExp*kernelWeight[distanceIndex[i]] + (1.0-alphaExp)*emaWeightFwd[i-1]; 
			emaWeightFwd[i]=alphaExp*kernelWeightOrig[distanceIndex[i]] + (1.0-alphaExp)*emaWeightFwd[i-1]; 
		}

		if(i>=windowSize && i<distanceSize-windowSize-1){
			avg=0.0;
			maxval=-1.0e+20;
			minval=1.0e+20;
			for(j=windowSize;j>0;j--){
				avg = avg + kernelWeight[distanceIndex[i-j]];
				if(kernelWeight[distanceIndex[i-j]]>maxval)
					maxval = kernelWeight[distanceIndex[i-j]]; 
				if(kernelWeight[distanceIndex[i-j]]<minval)
					minval = kernelWeight[distanceIndex[i-j]];
			}
			//kernelWeight[distanceIndex[i]] = avg/((TYPE)windowSize); 
			avgWeightFwd[i] = avg/((TYPE)windowSize); 
			medianWeightFwd[i] = 0.5*(maxval+minval); 
		}
		else if(i<windowSize){
			if(i==0){
				avgWeightFwd[i] = kernelWeight[distanceIndex[i]]; 
				medianWeightFwd[i] = kernelWeight[distanceIndex[i]]; 
			}
			else{
				avg=0.0;
				maxval=-1.0e+20;
				minval=1.0e+20;
				for(j=i;j>=0;j--){
					avg = avg + kernelWeight[distanceIndex[i-j]];
					if(kernelWeight[distanceIndex[i-j]]>maxval)
						maxval = kernelWeight[distanceIndex[i-j]]; 
					if(kernelWeight[distanceIndex[i-j]]<minval)
						minval = kernelWeight[distanceIndex[i-j]]; 

				}
				//kernelWeight[distanceIndex[i]] = avg/((TYPE)windowSize); 
				avgWeightFwd[i] = avg/((TYPE)(i+1));
				medianWeightFwd[i] = 0.5*(maxval+minval); 
			}
		}
		else{ //i>=distaceSize-windowSize-1
			avg=0.0;
			maxval=-1.0e+20;
			minval=1.0e+20;
			counter=0;
			for(j=distanceSize-i;j>=0;j--){
				avg = avg + kernelWeight[distanceIndex[i-j]];
				if(kernelWeight[distanceIndex[i-j]]>maxval)
					maxval = kernelWeight[distanceIndex[i-j]]; 
				if(kernelWeight[distanceIndex[i-j]]<minval)
					minval = kernelWeight[distanceIndex[i-j]]; 

				counter++;
			}
			//kernelWeight[distanceIndex[i]] = avg/((TYPE)windowSize); 
			avgWeightFwd[i] = avg/((TYPE)(counter));
			medianWeightFwd[i] = 0.5*(maxval+minval);
		}
	}

	//// fast backward average pass
	////for(i=distanceSize-windowSize-1;i>=0;i--){
	//for(i=distanceSize-1;i>=0;i--){
	//	if(i<=distanceSize-windowSize-1){
	//		avg=0.0;
	//		for(j=0;j<windowSize;j++){	
	//			avg = avg + kernelWeight[distanceIndex[i+j]];
	//		}
	//		//kernelWeight[distanceIndex[i]] = avg/((TYPE)windowSize); 
	//		avgWeightBwd[i] = avg/((TYPE)windowSize); 
	//	}
	//	else{
	//		avgWeightBwd[i] = kernelWeight[distanceIndex[i]]; 
	//	}
	//}

	for(i=0;i<distanceSize;i++){
		//kernelWeight[distanceIndex[i]] = avgWeightFwd[i]; 
		//kernelWeight[distanceIndex[i]] = medianWeightFwd[i]; 
		//kernelWeight[distanceIndex[i]] = 0.5*(medianWeightFwd[i] + avgWeightFwd[i]); 

		kernelWeight[distanceIndex[i]] = emaWeightFwd[i];
		kernelWeightOrig[distanceIndex[i]] = emaWeightFwd[i];

		//if(i==0){
		//	kernelWeight[distanceIndex[i]] = avgWeightFwd[i]; 
		//}
		//else if(i<windowSize){
		//	kernelWeight[distanceIndex[i]] = avgWeightFwd[i]*0.9+avgWeightBwd[i]*0.1; 
		//}
		//else if(i>=distanceSize-windowSize){
		//	kernelWeight[distanceIndex[i]] = avgWeightFwd[i]*0.1+avgWeightBwd[i]*0.9; 
		//}
		//else{
		//	kernelWeight[distanceIndex[i]] = (avgWeightFwd[i]+avgWeightBwd[i])*0.5; 
		//}
	}


/*
	//windowSize=ceil((TYPE)(windowSize)*1.4); 
	windowSize=10; 
	// medium average pass
	for(i=windowSize;i<distanceSize;i++){
		avg=0.0;
		for(j=windowSize;j>0;j--){
			
			avg = avg + kernelWeight[distanceIndex[i-j]];
		}
		//kernelWeight[distanceIndex[i]] = avg/((TYPE)windowSize); 
		avgWeight[i] = avg/((TYPE)windowSize); 
	}
	for(i=windowSize;i<distanceSize;i++){
		kernelWeight[distanceIndex[i]] = avgWeight[i]; 
	}

	//windowSize=ceil((TYPE)(windowSize)*1.4); 
	windowSize=30; 
	// slow average pass
	for(i=windowSize;i<distanceSize;i++){
		avg=0.0;
		for(j=windowSize;j>0;j--){
			
			avg = avg + kernelWeight[distanceIndex[i-j]];
		}
		//kernelWeight[distanceIndex[i]] = avg/((TYPE)windowSize); 
		avgWeight[i] = avg/((TYPE)windowSize); 
	}
	for(i=windowSize;i<distanceSize;i++){
		kernelWeight[distanceIndex[i]] = avgWeight[i]; 
	}
*/




	// backward average pass
	//for(i=distanceSize-windowSize-1;i>=0;i--){
	//	avg=0.0;
	//	for(j=0;j<windowSize;j++){
	//		
	//		avg = avg + kernelWeight[distanceIndex[i+j]];
	//	}
	//	//kernelWeight[distanceIndex[i]] = avg/((TYPE)windowSize); 
	//	avgWeight[i] = avg/((TYPE)windowSize); 
	//}
	//for(i=0;i<distanceSize-windowSize;i++){
	//	kernelWeight[distanceIndex[i]] = avgWeight[i]; 
	//}



//	for(i=distanceSize-windowSize;i>=0;i++){
//		avg=0.0;
//		for(j=windowSize;j>0;j--)
//			avg = avg + distanceValue[i-j];
//		distanceValue[i] = avg/((TYPE)windoSize); 
//	}

	getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);
	
	for(nk=-neighRadius;nk<=0;nk++){
		for(nj=-neighRadius;nj<=0;nj++){
			for(ni=-neighRadius;ni<=0;ni++){
        			getCoordinates(ni,nj,nk,params,&xcoord2,&ycoord2,&zcoord2);
        			dist = sqrt(	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
					(ycoord1-ycoord2)*(ycoord1-ycoord2)+
					(zcoord1-zcoord2)*(zcoord1-zcoord2));
		
				distIndex = existValue(distanceValue,distanceSize,dist);

				*(kernelWeight+
	  					(ni+neighRadius)+
	  					(nj+neighRadius)*neighSide+
	  					(nk+neighRadius)*neighSide*neighSide) = kernelWeight[distanceIndex[distIndex]]; 
				*(kernelWeightOrig+
	  					(ni+neighRadius)+
	  					(nj+neighRadius)*neighSide+
	  					(nk+neighRadius)*neighSide*neighSide) = kernelWeightOrig[distanceIndex[distIndex]];
	  					//(nk+neighRadius)*neighSide*neighSide) = kernelWeight[distanceIndex[distIndex]];

			}
		}
	}
	
	return dist;
}
