#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>


TYPE modifyWeights(	
				TYPE* kernelWeight, 
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

	//TYPE maxval=0.05, minval=0.0;
	TYPE maxval=0.1, minval=0.0;

	//if(iter==1 || mode==2){
	FILE *fpcurrentweight;
	
	if(iter==1)
		fpcurrentweight=fopen("initialdistanceweights.dat","w");
	else{
		if(iter%100==0)
			fpcurrentweight=fopen("currentdistanceweights.dat","w");
	}

	if(iter==1 || iter%100==0){
	getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);
	for(nk=-neighRadius;nk<=0;nk++){
 		for(nj=-neighRadius;nj<=0;nj++){
			for(ni=-neighRadius;ni<=0;ni++){
				getCoordinates(ni,nj,nk,params,&xcoord2,&ycoord2,&zcoord2);
				dist = sqrt(	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
					(ycoord1-ycoord2)*(ycoord1-ycoord2)+
					(zcoord1-zcoord2)*(zcoord1-zcoord2));

  				//fprintf(fpcurrentweight,"weight(%d,%d,%d)=%f;\n",ni+1,nj+1,nk+1,
  				fprintf(fpcurrentweight,"%f %f\n",dist,
					*(kernelWeight+
					(ni+neighRadius)+
					(nj+neighRadius)*neighSide+
					(nk+neighRadius)*neighSide*neighSide));
			}
		}
	}

	fclose(fpcurrentweight);
	}
	//}

	if(mode==0 || mode==2){

		//if(iter%50==0 || iter==1 || mode==2){


		if(mode==0){

/*
		rani = rand()%neighSide;
		ranj = rand()%neighSide;
		rank = rand()%neighSide;

		index = (rani)+
		    (ranj)*neighSide+
		    (rank)*neighSide*neighSide;

		val = *( kernelWeight + index );

		radi = rani-neighRadius;
		radj = ranj-neighRadius;
		radk = rank-neighRadius;

        	getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);
        	getCoordinates(radi,radj,radk,params,&xcoord2,&ycoord2,&zcoord2);

               	dist = sqrt(	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	       			(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	       			(zcoord1-zcoord2)*(zcoord1-zcoord2));

		distIndex = existValue(distanceValue,distanceSize,dist);
*/


		/*only for circular kernel*/
		//i=0;
		//while(i<distanceSize){
		//	if(distanceValue[i]<5.0)
		//		i++;
		//	else
		//		break;
		//}

		int counter=0;
		distIndex=0;	
		//while(distIndex==0 && counter<10){
		//while(counter<10){
			distIndex = rand() % distanceSize;
			//distIndex = rand() % i;// only for circular kernel
		//	counter++;
		//}
		dist = distanceValue[distIndex];

		val=kernelWeight[distanceIndex[distIndex]]; 


		if(distIndex==0){
			uvalue=maxval;
			lvalue=kernelWeight[distanceIndex[1]];//kernelWeight[distanceIndex[0]];
		}
		else if(distIndex==distanceSize-1){
			lvalue=minval;
			uvalue=kernelWeight[distanceIndex[distanceSize-2]];
		}
		else{
			lvalue=kernelWeight[distanceIndex[distIndex+1]];
			uvalue=kernelWeight[distanceIndex[distIndex-1]];
		}
		

		//if(lvalue>uvalue)printf("ERROR: uvalue=%f val=%f lvalue=%f\n",uvalue,val,lvalue);


		TYPE rn=((TYPE)rand()/(TYPE)RAND_MAX);
		TYPE rn2=((TYPE)rand()/(TYPE)RAND_MAX)*(maxval-minval);

		newval=0.0;

		//newval = MAX(val + (TYPE)(2*(rand()%2) -1)*0.68*val*rn,0.0); 
		counter=0;

		if(val==0.0){
			val=0.1 * MAX( kernelWeight[distanceIndex[0]] , maxval );
		}

		while(newval==0.0 && counter<10){
			//newval = MAX(val + (TYPE)(2*(rand()%2) -1)*0.68*val*rn,0.0); 

			//newval = MAX(rn*(abs(uvalue - lvalue)+0.001) + lvalue,0.0); 
			//rn=((TYPE)rand()/(TYPE)RAND_MAX);

			//newval = MAX(rn*(abs(uvalue - lvalue))*((TYPE)(2*(rand()%2) -1)) + 1.2*lvalue,0.0); 
			//newval = MAX(rn*(abs(uvalue - lvalue) + (TYPE)(2*(rand()%2) -1)*0.0001) + lvalue,0.0); 
			//newval = MAX(rn*(abs(uvalue - lvalue) + (TYPE)(2*(rand()%2) -1)*0.0001) + lvalue,0.0); 
			//newval = MAX(rn*(abs(uvalue - lvalue) + (TYPE)(2*(rand()%2) -1)*0.001) + lvalue,0.0); // ok
			//newval = MAX(rn*(abs(uvalue - lvalue) +(TYPE)(2*(rand()%2) -1)*0.005) + lvalue,0.0); 
			//newval = MAX(rn*(abs((uvalue - lvalue)) * (TYPE)(2*(rand()%2) -1)* rn2) + val,0.0); 
			//newval = MAX(rn*(0.25*val)*((TYPE)(2*(rand()%2) -1))+val,0.0); 
			//newval = MAX(rn*((maxval-minval)*0.1)*((TYPE)(2*(rand()%2) -1))+val,0.0); 
			newval = MAX(rn*(MIN(2.0*val,2.0*maxval)),0.0); 
			//newval = MAX(rn*(2.0*val),0.0); 
			rn=((TYPE)rand()/(TYPE)RAND_MAX);
			//rn2=((TYPE)rand()/(TYPE)RAND_MAX)*(maxval-minval);
			//newval = MAX(rn*(2.0*abs(uvalue - lvalue))*(TYPE)(2*(rand()%2) -1) + val,0.0); 
			counter++;
		}
		
		//val = *( kernelWeight + index );

		if(newval==0.0){
			newval=val;
			//printf("\nZERO!\n");
			//counter=0;
			//while(newval==0.0 && counter<10){
			//	newval = MAX(rn*(2.0*(val+0.0001)),0.0); 
			//	rn=((TYPE)rand()/(TYPE)RAND_MAX);
			//	counter++;
			//}	
		}

		//newval = MAX(val + (TYPE)(2*(rand()%2) -1)*1.0*(val<=0.000005?0.005:val)*rn,0.0); 
		//newval = MAX(val + (TYPE)(2*(rand()%2) -1)*0.2*val*rn,0.0); 
		*newvalOut = newval;
		*oldvalOut = val;

		}

	}
	else{
		index=-1;
		dist = distIn;
		newval=newvalIn;
	}
	//*(kernelWeight+
	//	(rani)+
	//	(ranj)*neighSide+
	//	(rank)*neighSide*neighSide) = val + MAX((TYPE)(2*(rand()%2) -1)*0.5*val,0.0);

	//printf("new val = %f\n",*(kernelWeight+
	//	(rani)+
	//	(ranj)*neighSide+
	//	(rank)*neighSide*neighSide));


        getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);

	TYPE sum=0.0,rsquared,distsquared=dist*dist;

	for(nk=-neighRadius;nk<=neighRadius;nk++){
	  for(nj=-neighRadius;nj<=neighRadius;nj++){
	    for(ni=-neighRadius;ni<=neighRadius;ni++){
               getCoordinates(ni,nj,nk,params,&xcoord2,&ycoord2,&zcoord2);
               rsquared = 	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	       			(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	       			(zcoord1-zcoord2)*(zcoord1-zcoord2);
               ////weight = cons * exp(-2.0*rsquared/((TYPE)(neighRadius*neighRadius)));
               //weight = cons * exp(-2.0*rsquared*asquaredinv);
               if(abs(rsquared-distsquared)<0.000001){
	          *(kernelWeight+
	  		(ni+neighRadius)+
	  		(nj+neighRadius)*neighSide+
	  		(nk+neighRadius)*neighSide*neighSide)=newval;
	       }
		//sum=sum+*(kernelWeight+
	  	//	(ni+neighRadius)+
	  	//	(nj+neighRadius)*neighSide+
	  	//	(nk+neighRadius)*neighSide*neighSide);
	    }
	  }
	}

	//printf("sum=%f\n",sum);

	return dist;
}
