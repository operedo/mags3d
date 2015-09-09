#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>


struct myData{
	float value;
	int index;
};


int cmpfunc (const void * a, const void * b)
{
	if ( *(float*)a <  *(float*)b ) return -1;
	else if ( *(float*)a >  *(float*)b ) return 1;
	else return 0;
}

int cmpfuncMyData (const void * a, const void * b)
{
	if ( ((struct myData*)a)->value < ((struct myData*)b)->value ) return -1;
	else if ( ((struct myData*)a)->value > ((struct myData*)b)->value ) return 1;
	else return 0;
}



void gauinv(double *p, float *xp){
	double p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim;
	lim=1.0e-10;
	p0=-0.322232431088 , p1=-1.0, p2=-0.342242088547, p3=-0.0204231210245, p4=-0.0000453642210148;
	q0=0.0993484626060, q1=0.588581570495, q2=0.531103462366, q3=0.103537752850, q4=0.0038560700634;

	if(*p<lim){
		*xp = -1.0e10f;
		return;
	}
	if(*p>(1.0-lim)){
		*xp =  1.0e10f;
		return;
	}

	
//
// Get k for an error situation:
//
	pp   = *p;
	if(*p>0.5) pp = 1 - pp;
	*xp   = 0.0;
	if(*p==0.5) return;
//
// Approximate the function:
//
	y  = sqrt(log(1.0/(pp*pp)));
	*xp = (float)( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / ((((y*q4+q3)*y+q2)*y+q1)*y+q0) );
	if((float)(*p)==(float)(pp)) *xp = -*xp;
//
// Return with G^-1(p):
//
	return;
}

/*
void dlocate(float *xx, int n, int is, int ie, double x, int *j){
//
// Initialize lower and upper methods:
//
	int jl = is-1;
	int ju = ie;
	int jm;
//
// If we are not done then compute a midpoint:
//

	while(ju-jl>1){
		jm = (ju+jl)/2;
//
// Replace the lower or upper limit with the midpoint:
//
		if((xx[ie]>xx[is]) && (x>(double)(xx[jm]))) 
			jl = jm;
		else
			ju = jm;
	}
//
// Return with the array index:
//
	*j = jl;
	return;
}
*/

void dlocate(struct myData *xx, int n, int is, int ie, double x, int *j){
//
// Initialize lower and upper methods:
//
	int jl = is-1;
	int ju = ie;
	int jm;
//
// If we are not done then compute a midpoint:
//

	while(ju-jl>1){
		jm = (ju+jl)/2;
//
// Replace the lower or upper limit with the midpoint:
//
		//if((xx[ie]>xx[is]) && (x>(double)(xx[jm]))) 
		if((xx[ie].value > xx[is].value) && (x>(double)(xx[jm].value))) 
			jl = jm;
		else
			ju = jm;
	}
//
// Return with the array index:
//
	*j = jl;
	return;
}




void dpowint(float xlow,float xhigh,double ylow,double yhigh,double xval,double power,double *ret){
//-----------------------------------------------------------------------
//
// Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
//                 for a value of x and a power pow.
//
//-----------------------------------------------------------------------
      double EPSLON=1.0e-20;

      if((xhigh-xlow)<EPSLON)
            *ret = (yhigh+ylow)/2.0;
      else
            *ret = ylow + (yhigh-ylow)*(pow(((xval-(double)(xlow))/((double)(xhigh)-(double)(xlow))),power));

      return;
}


void nscore(
				int *nd, 
				int *MAXVAR,
				int *maxdat,
				float **x, 
				float **y, 
				float **z,
				float **vr
){

	int i;
	float *vrtmp;
	double *wt_ns;
	double twt;

	struct myData* vraux = (struct myData*)malloc(*nd*sizeof(struct myData));

	for(i=0;i<*nd;i++){
		vraux[i].value=(*vr)[i];
		vraux[i].index=i;
	}

	twt=0.0;
	wt_ns=(double *)calloc(*nd,sizeof(double));
	//vrtmp=(float *)calloc(*nd,sizeof(float));

	for(i=0;i<*nd;i++){
		wt_ns[i]=1.0;
		twt=twt+wt_ns[i];
		//vrtmp[i]=(*vr)[i];
	}

	//qsort(*vr, *nd, sizeof(float), cmpfunc);
	//qsort(vrtmp, *nd, sizeof(float), cmpfunc);
	qsort(vraux, *nd, sizeof(struct myData), &cmpfuncMyData);

	double wtfac = 1.0/twt;
	double oldcp=0.0;
	double cp=0.0;
	double w;

	float vrrg;
	double vrg,vrr;

	for(i=0;i<*nd;i++){
		w = wtfac * wt_ns[i];
		cp = cp + w;
		wt_ns[i] = (cp+oldcp)*0.5;
		gauinv(&(wt_ns[i]),&vrrg);
		vrg=(double)vrrg;
		oldcp=cp;
		wt_ns[i]=vrg;
	}

	int doubone,j;
	for(i=0;i<*nd;i++){
		//vrr = (double)( (*vr)[i] );
		vrr = (double)( vraux[i].value );
		//dlocate(*vr,*nd,0,*nd-1,vrr,&j);
		//dlocate(vrtmp,*nd,0,*nd-1,vrr,&j);
		dlocate(vraux,*nd,0,*nd-1,vrr,&j);
		j=MIN(MAX(0,j),*nd-1);
		//dpowint((*vr)[j],(*vr)[j+1],wt_ns[j],wt_ns[j+1],vrr,1.0,&vrg);
		dpowint(vraux[j].value,vraux[j+1].value,wt_ns[j],wt_ns[j+1],vrr,1.0,&vrg);
		//printf("vrg=%f\n",vrg);
		//vrtmp[i]=(float)vrg;	
		(*vr)[vraux[i].index]=(float)vrg;	
	}
	//for(i=0;i<*nd;i++){
	//	(*vr)[i]=vrtmp[i];
	//}

	free(wt_ns);
	//free(vrtmp);
	free(vraux);
}
