#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>
#include <sys/time.h>
	
#define BUFFSIZE 20
static char buffer[BUFFSIZE];

int main(int argc,char *argv[]){
//      integer nd,irepo,maxdat,MAXVAR
	int nd, irepo, maxdat, MAXVAR;
//      real x(maxdat),y(maxdat),z(maxdat)
	float *x=NULL, *y=NULL, *z=NULL;
//      real EPSLON
	float EPSLON;
//      integer nlag
	int nlag;
//      real xlag,xltol
	float xlag, xltol;
//      integer mxdlv
	int mxdlv;
//      real*8 np(mxdlv),dis(mxdlv),gam(mxdlv),hm(mxdlv),
//     + tm(mxdlv),hv(mxdlv),tv(mxdlv)
	double *np=NULL, *dis=NULL, *gam=NULL, *hm=NULL, *tm=NULL, *hv=NULL, *tv=NULL;
//      integer numThreads
	int numThreads;
//      real*8 reducedVariables(7,mxdlv,numThreads)
	double *reducedVariables;
//      real dismxs,tmax,tmin
	float dismxs, tmax, tmin;
//      integer ndir,nvarg
	int ndir, nvarg;
//      real uvxazm(100),uvyazm(100),uvzdec(100),uvhdec(100)
	float *uvxazm=NULL, *uvyazm=NULL, *uvzdec=NULL, *uvhdec=NULL;
//      real csatol(100),csdtol(100),bandwh(ndir),bandwd(ndir)
	float *csatol=NULL, *csdtol=NULL, *bandwh=NULL, *bandwd=NULL;
//      real atol(ndir)
	float *atol=NULL;
//      integer ivtype(nvarg),ivtail(nvarg),ivhead(nvarg)
	int *ivtype=NULL, *ivtail=NULL, *ivhead=NULL;
//      real vr(maxdat,MAXVAR)
	float *vr=NULL;
	int isill;
	double *sills=NULL;

	printf("begin...\n");
	strncpy(buffer,"gamv-large.par",BUFFSIZE);
	printf("end...\n");

	gamvReadParams(
		&nd, &irepo, &maxdat, &MAXVAR,
		&x, &y, &z,
		&EPSLON,
		&nlag,
		&xlag, &xltol,
		&mxdlv,
		&np, &dis, &gam, &hm, &tm, &hv, &tv,
		&numThreads,
		&reducedVariables,
		&dismxs, &tmax, &tmin,
		&ndir, &nvarg,
		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
		&csatol, &csdtol, &bandwh, &bandwd,
		&atol,
		&ivtype, &ivtail, &ivhead,
		&vr,
		buffer,
		&isill,
		&sills
	);

	printf("Parameters loaded ok.\n");
	printf("nd=%d\n",nd);

	int i=0;
	//for(i=0;i<maxdat;i++)
	//	printf("vr[%d]=%f\n",i,vr[i]);
/*
	for(i=0;i<ndir;i++)
		printf("%f %f %f %f\n",uvxazm[i],uvyazm[i],uvzdec[i],uvhdec[i]);
*/



	gamv(
		&nd, &irepo, &maxdat, &MAXVAR,
		&x, &y, &z,
		&EPSLON,
		&nlag,
		&xlag, &xltol,
		&mxdlv,
		&np, &dis, &gam, &hm, &tm, &hv, &tv,
		&numThreads,
		&reducedVariables,
		&dismxs, &tmax, &tmin,
		&ndir, &nvarg,
		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
		&csatol, &csdtol, &bandwh, &bandwd,
		&atol,
		&ivtype, &ivtail, &ivhead,
		&vr	
	);


/*
	gamvCUDA(
		&nd, &irepo, &maxdat, &MAXVAR,
		&x, &y, &z,
		&EPSLON,
		&nlag,
		&xlag, &xltol,
		&mxdlv,
		&np, &dis, &gam, &hm, &tm, &hv, &tv,
		&numThreads,
		&reducedVariables,
		&dismxs, &tmax, &tmin,
		&ndir, &nvarg,
		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
		&csatol, &csdtol, &bandwh, &bandwd,
		&atol,
		&ivtype, &ivtail, &ivhead,
		&vr	
	);
*/

	int nsiz=(ndir)*(nvarg)*(nlag+2); 
	//for(i=0;i<nsiz;i++){
	//	printf("%f\n",gam[i]);
	//}

	gamvAverages(&ndir, &nvarg, &nlag, &isill, 
			&ivtype, &ivtail, &ivhead, &sills,
			&EPSLON,
			&np, &dis, &gam, &hm, &tm, &hv, &tv);


	nsiz=(ndir)*(nvarg)*(nlag+2); 
	for(i=0;i<nsiz;i++){
		printf("%d %f %f %d %f %f\n",(i+1),dis[i],gam[i],(int)(np[i]),tm[i],hm[i]);
	}

	gamvFreeMemory(
		&x, &y, &z,
		&np, &dis, &gam, &hm, &tm, &hv, &tv,
		&uvxazm, &uvyazm, &uvzdec, &uvhdec,
		&csatol, &csdtol, &bandwh, &bandwd,
		&atol,
		&ivtype, &ivtail, &ivhead,
		&vr,
		&sills,
		&reducedVariables,
		&numThreads);

}
