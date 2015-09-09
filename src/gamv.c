
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX(x,y)  ((x) >= (y) ? (x) : (y))
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#define PI 3.14159265
#define MV 500

int gamvReadParams(
//      integer nd,irepo,maxdat,MAXVAR
	int *nd, int *irepo, int *maxdat, int *MAXVAR,
//      real x(maxdat),y(maxdat),z(maxdat)
	float **x, float **y, float **z,
//      real EPSLON
	float *EPSLON,
//      integer nlag
	int *nlag,
//      real xlag,xltol
	float *xlag, float *xltol,
//      integer mxdlv
	int *mxdlv,
//      real*8 np(mxdlv),dis(mxdlv),gam(mxdlv),hm(mxdlv),
//     + tm(mxdlv),hv(mxdlv),tv(mxdlv)
	double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv,
//      integer numThreads
	int *numThreads,
//      real*8 reducedVariables(7,mxdlv,numThreads)
	double **reducedVariables,
//      real dismxs,tmax,tmin
	float *dismxs, float *tmax, float *tmin,
//      integer ndir,nvarg
	int *ndir, int *nvarg,
//      real uvxazm(100),uvyazm(100),uvzdec(100),uvhdec(100)
	float **uvxazm, float **uvyazm, float **uvzdec, float **uvhdec,
//      real csatol(100),csdtol(100),bandwh(ndir),bandwd(ndir)
	float **csatol, float **csdtol, float **bandwh, float **bandwd,
//      real atol(ndir)
	float **atol,
//      integer ivtype(nvarg),ivtail(nvarg),ivhead(nvarg)
	int **ivtype, int **ivtail, int **ivhead,
//      real vr(maxdat,MAXVAR)
	float **vr,
	char *str,
	int *isill,
	double **sills
	){

	int i,j,flag,flag2,flag3,numcols;
	char *line,*line2;
	char datafl[50],outfl[50],dummyfl[50],dummyfl2[50];
	FILE *fpparams, *fpdata;
	int BUFSIZE=100;

	int ixl,iyl,izl,ncut,iv;
	int nvar, *ivar,ivar0,ivar1,ivar2;
	float fvar0,fvar1,fvar2;
	float xtmp,ytmp,ztmp;
	int nvari;
	float *azm;
	float *dip;
	float *dtol;
	float dummy;

	int num[MV];
	float var[MV],cut[MV];
	double avg[MV],ssq[MV];
	float *vrmin,*vrmax;

	//printf("dentro de gamvReadParams: params %s\n",str);

	line=(char *)malloc(sizeof(char)*BUFSIZE + 1);
	line2=(char *)malloc(sizeof(char)*BUFSIZE + 1);

	//EPSLON=(float *)malloc(sizeof(float));
      	*EPSLON  = 1.0e-20;

	//printf("dentro de gamvReadParams: abriendo\n");
      	//fpparams = fopen(str,"r");
      	if( (fpparams = fopen(str,"r"))==NULL ){
		printf("ERROR: file %s does not exist. goodbye!\n",str);
		exit(1);
	}
	//printf("dentro de gamvReadParams: abrio ok?\n");

//
// Find Start of Parameters:
//
	i=0;
	flag3=0;
	while((fgets(line,BUFSIZE,fpparams))!=0){
		if(i==1){
			sscanf(line,"%s[ \t]+%s",datafl,dummyfl);
			//printf("parse: %s\n",datafl);

		      	if( (fpdata = fopen(datafl,"r"))==NULL ){
				printf("ERROR: file %s does not exist. Be careful!\n",datafl);
				//printf("ERROR: file %s does not exist. goodbye!\n",datafl);
				//exit(1);
			}			
	
			//maxdat=(int *)malloc(sizeof(int));

			*maxdat=0;
			flag2=0;
			while((fgets(line2,BUFSIZE,fpdata))!=0){
				if(flag2==1){
					sscanf(line2,"%d",&numcols);
				}
				if(flag2>=numcols+2){
					//sscanf(line2,"%f[ \t]+%f[ \t]+%f",&numcols);
					//printf("datafl: %s",line2);
					*maxdat=*maxdat+1;
				}
				flag2++;
				
			}
			//fclose(fpdata);	


		}

		if(i==2){
			//sscanf(line,"%d[ \t]+%d[ \t]+%d[ \t]+%s",&ixl,&iyl,&izl,dummyfl);
			sscanf(line,"%d %d %d[ \t]+%s",&ixl,&iyl,&izl,dummyfl);
			//printf("parse: %d %d %d\n",ixl,iyl,izl);
		}

		if(i==3){
			sscanf(line,"%d[ \t]+%s",&nvar,dummyfl);
			ivar=(int *)malloc(nvar*sizeof(int));
			if(nvar==1){
				sscanf(line,"%d %d[ \t]+%s",&nvar,&(ivar[0]),dummyfl);
				//ivar[0]=ivar0;
				//printf("parse: %d %d\n",nvar,ivar[0]);
			}
			else if(nvar==2){
				sscanf(line,"%d %d %d[ \t]+%s",&nvar,&(ivar[0]),&(ivar[1]),dummyfl);
				//ivar[0]=ivar0;
				//ivar[1]=ivar1;
				//printf("parse: %d %d %d\n",nvar,ivar[0],ivar[1]);

			}
			else if(nvar==3){
				sscanf(dummyfl,"%d %d %d %d[ \t]+%s",&nvar,&(ivar[0]),&(ivar[1]),&(ivar[2]),dummyfl);
				//ivar[0]=ivar0;
				//ivar[1]=ivar1;
				//ivar[2]=ivar2;
				//printf("parse: %d %d %d %d\n",nvar,ivar[0],ivar[1],ivar[2]);
			}
			else{
				printf("ERROR: maximum number of variables is 3 (not supported, fix this).");
				exit(1);
			}
		}

		if(i==4){
			//tmin=(float *)malloc(sizeof(float));
			//tmax=(float *)malloc(sizeof(float));
			sscanf(line,"%f %f[ \t]+%s",tmin,tmax,dummyfl);
			//printf("parse: %f %f\n",*tmin,*tmax);
		}

		if(i==5){
			sscanf(line,"%s[ \t]+%s",outfl,dummyfl);
			//printf("parse: %s\n",outfl);
		}

		if(i==6){
			//nlag=(int *)malloc(sizeof(int));
			sscanf(line,"%d[ \t]+%s",nlag,dummyfl);
			//printf("parse: %d\n",*nlag);
		}

		if(i==7){
			//xlag=(float *)malloc(sizeof(float));
			sscanf(line,"%f[ \t]+%s",xlag,dummyfl);
			//printf("parse: %f\n",*xlag);
		}

		if(i==8){
			//xltol=(float *)malloc(sizeof(float));
			sscanf(line,"%f[ \t]+%s",xltol,dummyfl);
			//printf("parse: %f\n",*xltol);
		}

		if(i==9){
			//ndir=(int *)malloc(sizeof(int));
			sscanf(line,"%d[ \t]+%s",ndir,dummyfl);
			//printf("parse: %d\n",*ndir);
			//printf("line: %s",line);

			if(*ndir<1){
				printf("ERROR: ndir is too small. Check parameters.\n");
				exit(1);
			}

			flag=*ndir-1;

			azm=(float *)malloc((*ndir)*sizeof(float));
			*atol=malloc((*ndir)*sizeof(**atol));
			//bandwh=malloc((*ndir)*sizeof(float));
			*bandwh=malloc((*ndir)*sizeof(**bandwh));
			dip=(float *)malloc((*ndir)*sizeof(float));
			dtol=(float *)malloc((*ndir)*sizeof(float));
			//bandwd=malloc((*ndir)*sizeof(float));
			*bandwd=malloc((*ndir)*sizeof(**bandwd));
	
			j=0;
			while((fgets(line,BUFSIZE,fpparams))!=0 && j<*ndir){
				sscanf(line,"%f %f %f %f %f %f[ \t]+%s",&(azm[j]),&((*atol)[j]),&((*bandwh)[j]),&(dip[j]),&(dtol[j]),&((*bandwd)[j]),dummyfl);
				//printf("parse: %f %f %f %f %f %f\n",azm[j],(*atol)[j],(*bandwh)[j],dip[j],dtol[j],(*bandwd)[j]);
				if((*bandwh)[j]<0.0){
					printf("ERROR: Horizontal bandwidth is too small!.\n");
					exit(1);
				}
				if((*bandwd)[j]<0.0){
					printf("ERROR: Vertical bandwidth is too small!.\n");
					exit(1);
				}
				j++;
				i++;
				//printf("line: %s",line);
			}	
		}
		if(i==(10+flag)){
			//isill=(int *)malloc(sizeof(int));
			sscanf(line,"%d[ \t]+%s",isill,dummyfl);
			//printf("parse: %d\n",*isill);
		}

		if(i==(11+flag)){

			//nvarg=(int *)malloc(sizeof(int));
			sscanf(line,"%d[ \t]+%s",nvarg,dummyfl);
			//printf("parse: %d\n",*nvarg);
			//printf("line: %s",line);
			//mxdlv=(int *)malloc(sizeof(int));
			*mxdlv=(*ndir)*(*nlag+2)*(*nvarg);

			//dis=(double *)malloc(*mxdlv*(sizeof(double)));
			//gam=(double *)malloc(*mxdlv*(sizeof(double)));
			//hm= (double *)malloc(*mxdlv*(sizeof(double)));
			//tm= (double *)malloc(*mxdlv*(sizeof(double)));
			//hv= (double *)malloc(*mxdlv*(sizeof(double)));
			//tv= (double *)malloc(*mxdlv*(sizeof(double)));
			//np= (double *)malloc(*mxdlv*(sizeof(double)));
			//ivtail=(int *)malloc(*nvarg*(sizeof(int)));
			//ivhead=(int *)malloc(*nvarg*(sizeof(int)));
			//ivtype=(int *)malloc(*nvarg*(sizeof(int)));

			*dis=malloc(*mxdlv*(sizeof(**dis)));
			*gam=malloc(*mxdlv*(sizeof(**gam)));
			*hm =malloc(*mxdlv*(sizeof(**hm )));
			*tm =malloc(*mxdlv*(sizeof(**tm )));
			*hv =malloc(*mxdlv*(sizeof(**hv )));
			*tv =malloc(*mxdlv*(sizeof(**tv )));
			*np =malloc(*mxdlv*(sizeof(**np )));
			*ivtail=malloc(*nvarg*(sizeof(**ivtail)));
			*ivhead=malloc(*nvarg*(sizeof(**ivhead)));
			*ivtype=malloc(*nvarg*(sizeof(**ivtype)));


			j=0;
			while((fgets(line,BUFSIZE,fpparams))!=0 && j<*nvarg){
				sscanf(line,"%d %d %d[ \t]+%s",&((*ivtail)[j]),&((*ivhead)[j]),&((*ivtype)[j]),dummyfl);
				//printf("parse: %d %d %d\n",(*ivtail)[j],(*ivhead)[j],(*ivtype)[j]);
				if((*ivtype)[j]==9 || (*ivtype)[j]==10){
					printf("ERROR: Variogram types 9 and 10 not supported..\n");
					exit(1);
				}
				i++;
				j++;
				//printf("line: %s",line);
			}
			//MAXVAR=(int *)malloc(sizeof(int));
			*MAXVAR=nvar;
			flag3=1;
		}

		if(flag3==1) break;

		//printf("line: %s",line);
		if(i>0)i++;	
		if(strncmp(line,"START",5)==0)
			i++;
		

		/*sscanf(line,"%lf %d %lf",&vgval,&num,&lagval);
		expVariogram[i]=vgval;	
		npairs[i]=num;
		lagDistance[i]=lagval;	
		//printf("%f %d\n",vgval,num);
		*/
	}

	fclose(fpparams);


//
// Perform some quick error checking:
//
	if(*xltol<=0.0){
		printf("xltol is too small: resetting to xlag/2");
		*xltol = 0.5*(*xlag);
	}


	rewind(fpdata);
	//if( (fpdata = fopen(datafl,"r"))==NULL ){
	//	printf("ERROR: data file %s does not exist!\n",datafl);
	//	exit(1);
	//}			

//
// The data file exists so open the file and read in the header
// information. Initialize the storage that will be used to summarize
// the data found in the file:
//
	//printf("llega aca 1?\n");
	//while((fgets(line2,BUFSIZE,fpdata))!=0 && flag<=1){
	//vr=(float *)malloc(*maxdat * *MAXVAR * sizeof(float));
	*vr=malloc(*maxdat * *MAXVAR * sizeof(**vr));
	vrmin=(float *)malloc(*MAXVAR*sizeof(float));	
	vrmax=(float *)malloc(*MAXVAR*sizeof(float));	
	//sills=(double *)malloc(*MAXVAR*sizeof(double));	
	*sills=malloc(*MAXVAR*sizeof(**sills));	
	//*x=(float *)malloc(*maxdat*sizeof(float));
	//*y=(float *)malloc(*maxdat*sizeof(float));
	//*z=(float *)malloc(*maxdat*sizeof(float));
	*x=malloc(*maxdat * sizeof(**x));
	*y=malloc(*maxdat * sizeof(**y));
	*z=malloc(*maxdat * sizeof(**z));

	//printf("llega aca 2?\n");

	nvari=nvar;
	for(i=0;i<nvari;i++){
		num[i] = 0.0;
		avg[i] = 0.0;
		ssq[i] = 0.0;
	}

	char *ret;
	ret=fgets(line,BUFSIZE,fpdata);
	ret=fgets(line,BUFSIZE,fpdata);
	sscanf(line,"%d",&numcols);
	//printf("data: %d\n",numcols);
	flag=0;
	//nd=(int *)malloc(sizeof(int));
	//*maxdat=0;
	*nd=0;
	int testdat=0;
	while((fgets(line,BUFSIZE,fpdata))!=0){
		//if(flag==1){
		//	sscanf(line2,"%d",&numcols);
		//}
		if(flag>=numcols){
			if(nvari==1){
				if(ivar[0]==4){
					sscanf(line,"%f %f %f %f",&xtmp,&ytmp,&ztmp,&(var[0]));
				}
				else if(ivar[0]==5){
					sscanf(line,"%f %f %f %f %f",&xtmp,&ytmp,&ztmp,&dummy,&(var[0]));
				}
				//printf("data: %f %f %f %f\n",xtmp,ytmp,ztmp,var[0]);
				if(var[0]>=*tmin && var[0]<*tmax) testdat=1;
			}
			if(nvari==2){
				sscanf(line,"%f %f %f %f %f",&xtmp,&ytmp,&ztmp,&(var[0]),&(var[1]));
				//printf("data: %f %f %f %f %f\n",xtmp,ytmp,ztmp,var[0],var[1]);
				if(var[0]>=*tmin && var[0]<*tmax && var[1]>=*tmin && var[1]<*tmax) testdat=1;
			}
			if(nvari==3){
				sscanf(line,"%f %f %f %f %f %f",&xtmp,&ytmp,&ztmp,&(var[0]),&(var[1]),&(var[2]));
				//printf("data: %f %f %f %f %f %f\n",xtmp,ytmp,ztmp,var[0],var[1],var[2]);
				if(var[0]>=*tmin && var[0]<*tmax && var[1]>=*tmin && var[1]<*tmax && var[2]>=*tmin && var[2]<*tmax) testdat=1;
			}
			if(nvari>3){
				printf("ERROR: maximum number of variables is 3 (not supported, fix this).");
				exit(1);
			}

			if(testdat){
				*nd=*nd+1;

				for(iv=0;iv<nvar;iv++){
					j=iv;
					//j=ivar[iv]-4-nvar;
					//printf("%d\n",j);
					//printf("%f\t%f\n",(*vr)[*nd-1 + iv*(*maxdat)],var[j]);
					(*vr)[*nd-1 + iv*(*maxdat)] = var[j];
					//printf("%f\t%f\n",(*vr)[*nd-1 + iv*(*maxdat)],var[j]);
					if(var[j]>=*tmin && var[j]<*tmax){
						num[iv] = num[iv] + 1;
						avg[iv] = avg[iv] + (double)(var[j]);
						ssq[iv] = ssq[iv] + (double)(var[j]*var[j]);
						//printf("post: %f %d %f %f\n",var[j],num[iv],avg[iv],ssq[iv]);
					}
				}

				//printf("post2-pre: %f %f %f %d %d %d\n",x[*nd],y[*nd],z[*nd],ixl,iyl,izl);
				if(ixl-1<0)
					(*x)[*nd-1] = 0.0;
				else
					(*x)[*nd-1] = xtmp;
				//printf("llega aca 3?\n");

					//x[*nd] = var[ixl-1];
				if(iyl-1<0)
					(*y)[*nd-1] = 0.0;
				else
					(*y)[*nd-1] = ytmp;
				//printf("llega aca 4?\n");
					//y[*nd] = var[iyl-1];
				if(izl-1<0)
					(*z)[*nd-1] = 0.0;
				else
					(*z)[*nd-1] = ztmp;
				//printf("llega aca 5?\n");
					//z[*nd] = var[izl-1];
				
				//printf("post2-post: %f %f %f %d %d %d\n",x[*nd],y[*nd],z[*nd],ixl,iyl,izl);


			}

			// mover maxdat mas arriba, en una pasada previa por los datos
		}
		flag++;
		
	}
	fclose(fpdata);


//c
//c Compute the averages and variances as an error check for the user:
//c
//      do iv=1,nvar
//            sills(iv) = -999.
//            if(num(iv).gt.0) then
//                  avg(iv)   = avg(iv)/dble(num(iv))
//                  ssq(iv)   =(ssq(iv)/dble(num(iv)))-avg(iv)*avg(iv)
//                  sills(iv) = ssq(iv)
//                  write(*,*) 'Variable number ',iv
//                  write(*,*) '  Number   = ',num(iv)
//                  write(*,*) '  Average  = ',real(avg(iv))
//                  write(*,*) '  Variance = ',real(ssq(iv))
//            endif
//      end do
//c
//c Construct Indicator Variables if necessary:
//c
//      do ic=1,ncut
//            iv   = ivc(ic)
//            jv   = nvar + ic
//            ptot = 0.0
//            p1   = 0.0
//            do id=1,nd
//                  if(vr(id,iv).le.tmin.or.vr(id,iv).gt.tmax) then
//                        vr(id,jv) = tmin - EPSLON
//                  else
//                        if(indflag(ic).eq.1) then
//                              if(vr(id,iv).lt.cut(ic)) then
//                                    vr(id,jv) = 0.0
//                              else
//                                    vr(id,jv) = 1.0
//                              endif
//                              p1   = p1   + vr(id,jv)
//                              ptot = ptot + 1.0
//                        else
//                              vr(id,jv) = 0.0
//                              if(int(vr(id,iv)+0.5).eq.int(cut(ic)+0.5))
//     +                        vr(id,jv) = 1.0
//                              p1   = p1   + vr(id,jv)
//                              ptot = ptot + 1.0
//                        end if
//                  end if
//            end do
//            p1        = p1 / max(ptot,1.0)
//            sills(jv) = dble (p1*(1.0-p1))
//      end do
//c
//c Establish minimums and maximums:
//c
	for(i=0;i<*MAXVAR;i++){
		vrmin[i] =  1.0e21;
		vrmax[i] = -1.0e21;
	}
	ncut=0;
	int id=0;
	for(id=0;id<*nd;id++){
		for(iv=0;iv<nvar+ncut;iv++){
			if((*vr)[id + iv*(*maxdat)]>=*tmin && (*vr)[id + iv*(*maxdat)]<*tmax){
				if((*vr)[id + iv*(*maxdat)]<vrmin[iv]) vrmin[iv] = (*vr)[id + iv*(*maxdat)];
				if((*vr)[id + iv*(*maxdat)]>vrmax[iv]) vrmax[iv] = (*vr)[id + iv*(*maxdat)];
			}
		}
	}
	//printf("max=%f min=%f\n",vrmax[0],vrmin[0]);



//c
//c Loop over all the variograms to be computed:
//c
//      write(*,*)
//      do iv=1,nvarg
//c
//c Note the variogram type and the variables being used:
//c
//      it = abs(ivtype(iv))
//      if(it.eq. 1) title(1:24) = 'Semivariogram          :'
//      if(it.eq. 2) title(1:24) = 'Cross Semivariogram    :'
//      if(it.eq. 3) title(1:24) = 'Covariance             :'
//      if(it.eq. 4) title(1:24) = 'Correlogram            :'
//      if(it.eq. 5) title(1:24) = 'General Relative       :'
//      if(it.eq. 6) title(1:24) = 'Pairwise Relative      :'
//      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms:'
//      if(it.eq. 8) title(1:24) = 'Semimadogram           :'
//      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram:'
//      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram:'
//      write(title(25:64),100) names(ivtail(iv)),names(ivhead(iv))
// 100  format('  tail=',a12,' head=',a12)
//      write(*,101) iv,title(1:64)
// 101  format(' Variogram ',i2,1x,a64)
//c
//c Check for possible errors or inconsistencies:
//c
//      if(it.eq.2) then
//            if(ivtail(iv).eq.ivhead(iv)) write(*,201)
// 201        format('  WARNING: cross variogram with the same variable!')
//      else if(it.eq.5) then
//            if(ivtail(iv).ne.ivhead(iv)) write(*,501)
//            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
//     +            write(*,502)
//            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
//     +            write(*,502)
// 501        format('  WARNING: cross general relative variogram are',
//     +             ' difficult to interpret!')
// 502        format('  WARNING: there are both positive and negative',
//     +             ' values - lag mean could be zero!')
//      else if(it.eq.6) then
//            if(ivtail(iv).ne.ivhead(iv)) write(*,601)
//            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
//     +            write(*,602)
//            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
//     +            write(*,602)
// 601        format('  WARNING: cross pairwise relative variogram are',
//     +             ' difficult to interpret!')
// 602        format('  WARNING: there are both positive and negative',
//     +             ' values - pair means could be zero!')
//      else if(it.eq.7) then
//            if(ivtail(iv).ne.ivhead(iv)) write(*,701)
//            if(vrmin(ivtail(iv)).lt.0.0.or.vrmin(ivhead(iv)).lt.0.0)
//     +      write(*,702)
// 701        format('  WARNING: cross logarithmic variograms may be',
//     +             ' difficult to interpret!')
// 702        format('  WARNING: there are zero or negative',
//     +             ' values - logarithm undefined!')
//      else if(it.eq.8) then
//            if(ivtail(iv).ne.ivhead(iv)) write(*,801)
// 801        format('  WARNING: cross rodograms may be difficult to',
//     +             ' interpret!')
//      else if(it.eq.9) then
//            if(ivtail(iv).ne.ivhead(iv)) write(*,901)
// 901        format('  WARNING: cross madograms may be difficult to',
//     +             ' interpret!')
//      endif
//c
//c END Loop over all variograms:
//c
//      end do
//
//
//      goto 1001
//c
//c Error in an Input File Somewhere:
//c
// 98   stop 'ERROR in parameter file!'
// 99   stop 'ERROR in data file!'
//
//1001  print *,'Parameters ok.'
//
//cc
//cc Call gamv to compute the required variograms:
//cc
//c      call gamv
//




//
// Define the distance tolerance if it isn't already:
//
	if(*xltol<=0.0) *xltol = 0.5 * (*xlag);


//
// Define the angles and tolerance for each direction:
//

	//uvxazm=(float *)malloc(*ndir*sizeof(float));
	//uvyazm=(float *)malloc(*ndir*sizeof(float));
	//uvzdec=(float *)malloc(*ndir*sizeof(float));
	//uvhdec=(float *)malloc(*ndir*sizeof(float));
	//csatol=(float *)malloc(*ndir*sizeof(float));
	//csdtol=(float *)malloc(*ndir*sizeof(float));

	//*uvxazm=malloc(*ndir*sizeof(**uvxazm));
	//*uvyazm=malloc(*ndir*sizeof(**uvyazm));
	//*uvzdec=malloc(*ndir*sizeof(**uvzdec));
	//*uvhdec=malloc(*ndir*sizeof(**uvhdec));
	//*csatol=malloc(*ndir*sizeof(**csatol));
	//*csdtol=malloc(*ndir*sizeof(**csdtol));

	*uvxazm=malloc(100*sizeof(**uvxazm));
	*uvyazm=malloc(100*sizeof(**uvyazm));
	*uvzdec=malloc(100*sizeof(**uvzdec));
	*uvhdec=malloc(100*sizeof(**uvhdec));
	*csatol=malloc(100*sizeof(**csatol));
	*csdtol=malloc(100*sizeof(**csdtol));
	float azmuth,declin;

	for(id=0;id<*ndir;id++){
//
// The mathematical azimuth is measured counterclockwise from EW and
// not clockwise from NS as the conventional azimuth is:
//
            azmuth     = (90.0-azm[id])*PI/180.0;
            (*uvxazm)[id] = cosf(azmuth);
            (*uvyazm)[id] = sinf(azmuth);
            if((*atol)[id]<=0.0)
                  (*csatol)[id] = cosf(45.0*PI/180.0);
            else
                  (*csatol)[id] = cosf((*atol)[id]*PI/180.0);
            
//
// The declination is measured positive down from vertical (up) rather
// than negative down from horizontal:
//
            declin     = (90.0-dip[id])*PI/180.0;
            (*uvzdec)[id] = cosf(declin);
            (*uvhdec)[id] = sinf(declin);
            if(dtol[id]<=0.0)
                  (*csdtol)[id] = cosf(45.0*PI/180.0);
            else
                  (*csdtol)[id] = cosf(dtol[id]*PI/180.0);
            
	}

//
// Initialize the arrays for each direction, variogram, and lag:
//
	int nsiz = (*ndir)*(*nvarg)*(*nlag+2);
	//dismxs=(float *)malloc(sizeof(float));
	for(i=0;i<nsiz;i++){
		(*np )[i]  = 0.0;
		(*dis)[i] = 0.0;
		(*gam)[i] = 0.0;
		(*hm )[i]  = 0.0;
		(*tm )[i]  = 0.0;
		(*hv )[i]  = 0.0;
		(*tv )[i]  = 0.0;
	}

	*dismxs = (((float)(*nlag) + 0.5f - *EPSLON) * *xlag) * (((float)(*nlag) + 0.5f - *EPSLON) * *xlag) ;

	//irepo=(int *)malloc(sizeof(int));
	*irepo = MAX(1,MIN((*nd/10),1000));


//#ifdef _OPENMP	
//	if(*numThreads>1){
//		*reducedVariables=malloc(7 * *mxdlv * *numThreads * sizeof(**reducedVariables));
//	}
//#endif


	free(ivar);
	free(azm);
	free(dip);
	free(dtol);
	free(vrmin);
	free(vrmax);




	free(line);
	free(line2);

	//printf("nd=%d\n",*nd);

	//for(i=0;i<*maxdat;i++)
	//	printf("vr[%d]=%f\n",i,(*vr)[i]);

}


//int extractstatisticscwrapper_(
int gamv(
//      integer nd,irepo,maxdat,MAXVAR
	int *nd, int *irepo, int *maxdat, int *MAXVAR,
//      real x(maxdat),y(maxdat),z(maxdat)
	float **x, float **y, float **z,
//      real EPSLON
	float *EPSLON,
//      integer nlag
	int *nlag,
//      real xlag,xltol
	float *xlag, float *xltol,
//      integer mxdlv
	int *mxdlv,
//      real*8 np(mxdlv),dis(mxdlv),gam(mxdlv),hm(mxdlv),
//     + tm(mxdlv),hv(mxdlv),tv(mxdlv)
	double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv,
//      integer numThreads
	int *numThreads,
//      real*8 reducedVariables(7,mxdlv,numThreads)
	double **reducedVariables,
//      real dismxs,tmax,tmin
	float *dismxs, float *tmax, float *tmin,
//      integer ndir,nvarg
	int *ndir, int *nvarg,
//      real uvxazm(100),uvyazm(100),uvzdec(100),uvhdec(100)
	float **uvxazm, float **uvyazm, float **uvzdec, float **uvhdec,
//      real csatol(100),csdtol(100),bandwh(ndir),bandwd(ndir)
	float **csatol, float **csdtol, float **bandwh, float **bandwd,
//      real atol(ndir)
	float **atol,
//      integer ivtype(nvarg),ivtail(nvarg),ivhead(nvarg)
	int **ivtype, int **ivtail, int **ivhead,
//      real vr(maxdat,MAXVAR)
	float **vr	
	)
{

	int threadId,i,j,id,ii,il,it,iv,jj;
	float dx,dy,dz,dxs,dys,dzs,hs,h;
	int lagbeg,lagend,ilag;
	float band,dcazm,dcdec,dxy,gamma,vrh,vrhpr,vrt,vrtpr;
	int omni;


	//printf("Entering gamv calculation. nd=%d EPSLON=%f xlag=%f xltol=%f\n",*nd,*EPSLON,*xlag,*xltol);


	for(i=0;i<*nd;i++){
        	//if((int)(i/(*irepo)* (*irepo))==i) printf("   currently on %d of %d\n",i,*nd);
        	//if((int)(i/(*irepo)* (*irepo))==i) printf("   currently on %d of %d\n",i,*nd);
        	//printf("   currently on %d of %d\n",i,*nd);
		for(j=i;j<*nd;j++){

//
// Definition of the lag corresponding to the current pair:
//
			dx  = (*x)[j] - (*x)[i];
			dy  = (*y)[j] - (*y)[i];
			dz  = (*z)[j] - (*z)[i];
			dxs = dx*dx;
			dys = dy*dy;
			dzs = dz*dz;
			hs  = dxs + dys + dzs;
			//printf("i=%d\tj=%d\tx=(%f,%f)\ty=(%f,%f)\tz=(%f,%f)\tdx=%f\tdy=%f\tdz=%f\tdismxs=%f\n",i,j,(*x)[j],(*x)[i],(*y)[j],(*y)[i],(*z)[j],(*z)[i],dx,dy,dz,*dismxs);
			if(hs > *dismxs) continue;
			if(hs < 0.0) hs = 0.0;
			h   = sqrtf(hs);


//
// Determine which lag this is and skip if outside the defined distance
// tolerance:
//
			if(h<=*EPSLON){
				lagbeg = 1;
				lagend = 1;
			}
			else{
				lagbeg = -1;
				lagend = -1;
				for(ilag=2;ilag<=*nlag+2;ilag++){
					if(h>=(*xlag*(float)(ilag-2)-*xltol) && h<=(*xlag*(float)(ilag-2)+*xltol)){
						if(lagbeg<0) lagbeg = ilag; 
						lagend = ilag; 
					}
				}
				if(lagend<0) continue;
			}

			//printf("i=%d j=%d dx=%f dy=%f dz=%fh=%f lagbeg=%d lagend=%d\n",i,j,dx,dy,dz,h,lagbeg,lagend);




//
// Definition of the direction corresponding to the current pair. All
// directions are considered (overlapping of direction tolerance cones
// is allowed):
//
			for(id=0;id<*ndir;id++){
			//
			// Check for an acceptable azimuth angle:
			//
				dxy = sqrtf(MAX((dxs+dys),0.0));
				if(dxy<*EPSLON){
					dcazm = 1.0;
				}
				else{
					dcazm = (dx*(*uvxazm)[id]+dy*(*uvyazm)[id])/dxy;
				}
				if(fabs(dcazm)<(*csatol)[id]) continue;

//
// Check the horizontal bandwidth criteria (maximum deviation 
// perpendicular to the specified direction azimuth):
//
				band = (*uvxazm)[id]*dy - (*uvyazm)[id]*dx;
				if(fabs(band)>=(*bandwh)[id]) continue;

				//printf("dxy=%f\tdcazm=%f\tband=%f\tlagbed=%d\ta=%f\tb=%f\tc=%f\n",dxy,dcazm,band,lagbeg,(*uvhdec)[id],(*uvzdec)[id],(*csdtol)[id]);


//
// Check for an acceptable dip angle:
//
				if(dcazm<0.0) dxy = -dxy;
				//printf("dxy=%f\tdcazm=%f\tband=%f\n",dxy,dcazm,band);
				if(lagbeg==1)
					dcdec = 0.0;
				else{
					dcdec = (dxy*(*uvhdec)[id]+dz*(*uvzdec)[id])/h;
					if(fabs(dcdec)<(*csdtol)[id]) continue;
				}
				//printf("dxy=%f\tdcazm=%f\tband=%f\tdcdec=%f\n",dxy,dcazm,band,dcdec);
//
// Check the vertical bandwidth criteria (maximum deviation perpendicular
// to the specified dip direction):
//
				band = (*uvhdec)[id]*dz - (*uvzdec)[id]*dxy;
				if(fabs(band)>(*bandwd)[id]) continue;
//
// Check whether or not an omni-directional variogram is being computed:
//
				omni = 0;
				if((*atol)[id]>=90.0) omni = 1;
//
// This direction is acceptable - go ahead and compute all variograms:
//

				//printf("dxy=%f dcazm=%f uvxazm[0]=%f uvyazm[0]=%f band=%f dcdec=%f omni=%d csdtol[0]=%f\n",dxy,dcazm,(*uvxazm)[0],(*uvyazm)[0],band,dcdec,omni,(*csdtol)[0]);




//				fprintf(stdout,"dcazm=%f\tdcdec=%f\n",dcazm,dcdec);
				for(iv=0;iv<*nvarg;iv++){
//
// For this variogram, sort out which is the tail and the head value:
//
					it = (*ivtype)[iv];
					//printf("it=%d\n",it);
					//printf("i=%d,j=%d,id=%d,iv=%d,vr[%d]=%f\tvr[%d]=%f\n",i,j,id,iv,i,(*vr)[i],j,(*vr)[j]);
					if(dcazm>=0.0 && dcdec>=0.0){
						ii = (*ivtail)[iv]-1;
						vrh   = (*vr)[i+ii*(*maxdat)];
						ii = (*ivhead)[iv]-1;
						vrt   = (*vr)[j+ii*(*maxdat)];
						if(omni || it==2){
							ii    = (*ivhead)[iv]-1;
							vrtpr = (*vr)[i+ii*(*maxdat)];
							ii    = (*ivtail)[iv]-1;
							vrhpr = (*vr)[j+ii*(*maxdat)];
						}
					}
					else{
						ii = (*ivtail)[iv]-1;
						vrh   = (*vr)[j+ii*(*maxdat)];
						ii = (*ivhead)[iv]-1;
						vrt   = (*vr)[i+ii*(*maxdat)];
						if(omni || it==2){
							ii    = (*ivhead)[iv]-1;
							vrtpr = (*vr)[j+ii*(*maxdat)];
							ii    = (*ivtail)[iv]-1;
							vrhpr = (*vr)[i+ii*(*maxdat)];
						}
					}
					//printf("vrh=%f\tvrt=%f\n",vrh,vrt);



//
// Reject this pair on the basis of missing values:
//
					if(vrt<*tmin || vrh<*tmin || vrt>*tmax || vrh>*tmax) continue;
					if(it==2 && (vrtpr<*tmin || vrhpr<*tmin || vrtpr>*tmax || vrhpr>*tmax)) continue;
//
//             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
//
//
// The Semivariogram:
//


					if(it==1 || it==5 || it>=9){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;



							(*np)[ii]  = (*np)[ii]  + 1.0;
							(*dis)[ii] = (*dis)[ii] + (double)(h);
							(*tm)[ii]  = (*tm)[ii]  + (double)(vrt);
							(*hm)[ii]  = (*hm)[ii]  + (double)(vrh);
							(*gam)[ii] = (*gam)[ii] + (double)((vrh-vrt)*(vrh-vrt));
							
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									(*np)[ii]  = (*np)[ii]  + 1.0;
									(*dis)[ii] = (*dis)[ii] + (double)(h);
									(*tm)[ii]  = (*tm)[ii]  + (double)(vrtpr);
									(*hm)[ii]  = (*hm)[ii]  + (double)(vrhpr);
									(*gam)[ii] = (*gam)[ii] + (double)((vrhpr-vrtpr)*(vrhpr-vrtpr));
								}
							}

							//printf("ii=%d\tnp=%f\n",ii,(*np)[ii]);
							//printf("ii=%d\tdis=%f\n",ii,(*dis)[ii]);
							//printf("ii=%d\ttm=%f\n",ii,(*tm)[ii]);
							//printf("ii=%d\thm=%f\n",ii,(*hm)[ii]);
							//printf("ii=%d\tgam=%f\n",ii,(*gam)[ii]);
							//printf("ii=%d\ttv=%f\n",ii,(*tv)[ii]);
							//printf("ii=%d\thv=%f\n",ii,(*hv)[ii]);
							
						}
					}

//
// The Traditional Cross Semivariogram:
//
					else if(it==2){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							(*np)[ii]  = (*np)[ii]  + 1.0;
							(*dis)[ii] = (*dis)[ii] + (double)(h);
							(*tm)[ii]  = (*tm)[ii]  + (double)(0.5*(vrt+vrtpr));
							(*hm)[ii]  = (*hm)[ii]  + (double)(0.5*(vrh+vrhpr));
							(*gam)[ii] = (*gam)[ii] + (double)((vrhpr-vrh)*(vrt-vrtpr));
						}
					}
//
// The Covariance:
//
					else if(abs(it)==3){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							(*np)[ii]  = (*np)[ii]  + 1.0;
							(*dis)[ii] = (*dis)[ii] + (double)(h);
							(*tm)[ii]  = (*tm)[ii]  + (double)(vrt);
							(*hm)[ii]  = (*hm)[ii]  + (double)(vrh);
							(*gam)[ii] = (*gam)[ii] + (double)(vrh*vrt);
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									(*np)[ii]  = (*np)[ii]  + 1.0;
									(*dis)[ii] = (*dis)[ii] + (double)(h);
									(*tm)[ii]  = (*tm)[ii]  + (double)(vrtpr);
									(*hm)[ii]  = (*hm)[ii]  + (double)(vrhpr);
									(*gam)[ii] = (*gam)[ii] + (double)(vrhpr*vrtpr);
								}
							}
						}
					}
//
// The Correlogram:
//
					else if(abs(it)==4){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							(*np)[ii]  = (*np)[ii]  + 1.0;
							(*dis)[ii] = (*dis)[ii] + (double)(h);
							(*tm)[ii]  = (*tm)[ii]  + (double)(vrt);
							(*hm)[ii]  = (*hm)[ii]  + (double)(vrh);
							(*hv)[ii]  = (*hv)[ii]  + (double)(vrh*vrh);
							(*tv)[ii]  = (*tv)[ii]  + (double)(vrt*vrt);
							(*gam)[ii] = (*gam)[ii] + (double)(vrh*vrt);
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									(*np)[ii]  = (*np)[ii]  + 1.0;
									(*dis)[ii] = (*dis)[ii] + (double)(h);
									(*tm)[ii]  = (*tm)[ii]  + (double)(vrtpr);
									(*hm)[ii]  = (*hm)[ii]  + (double)(vrhpr);
									(*hv)[ii]  = (*hv)[ii]  + (double)(vrhpr*vrhpr);
									(*tv)[ii]  = (*tv)[ii]  + (double)(vrtpr*vrtpr);
									(*gam)[ii] = (*gam)[ii] + (double)(vrhpr*vrtpr);
								}
							}
						}
					}
//
// The Pairwise Relative:
//
					else if(it==6){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							if(abs(vrt+vrh)>*EPSLON){
								(*np)[ii]  = (*np)[ii]  + 1.0;
								(*dis)[ii] = (*dis)[ii] + (double)(h);
								(*tm)[ii]  = (*tm)[ii]  + (double)(vrt);
								(*hm)[ii]  = (*hm)[ii]  + (double)(vrh);
								gamma   = 2.0*(vrt-vrh)/(vrt+vrh);
								(*gam)[ii] = (*gam)[ii] + (double)(gamma*gamma);
							}
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									if(abs(vrtpr+vrhpr)>*EPSLON){
										(*np)[ii]  = (*np)[ii]  + 1.0;
										(*dis)[ii] = (*dis)[ii] + (double)(h);
										(*tm)[ii]  = (*tm)[ii]  + (double)(vrtpr);
										(*hm)[ii]  = (*hm)[ii]  + (double)(vrhpr);
										gamma   = 2.0*(vrt-vrh)/(vrt+vrh);
										(*gam)[ii] = (*gam)[ii] + (double)(gamma*gamma);
									}
								}
							}
						}
					}
//c
//c Variogram of Logarithms:
//c
//      else if(it.eq.7) then
//            do il=lagbeg,lagend
//               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
//               if(vrt.gt.EPSLON.and.vrh.gt.EPSLON) then
//                     np(ii)  = np(ii)  + 1.
//                     dis(ii) = dis(ii) + dble(h)
//                     tm(ii)  = tm(ii)  + dble(vrt)
//                     hm(ii)  = hm(ii)  + dble(vrh)
//                     gamma   = alog(vrt)-alog(vrh)
//                     gam(ii) = gam(ii) + dble(gamma*gamma)
//               endif
//               if(omni) then
//                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
//     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
//                     if(vrtpr.gt.EPSLON.and.vrhpr.gt.EPSLON) then
//                           np(ii)  = np(ii)  + 1.
//                           dis(ii) = dis(ii) + dble(h)
//                           tm(ii)  = tm(ii)  + dble(vrtpr)
//                           hm(ii)  = hm(ii)  + dble(vrhpr)
//                           gamma   = alog(vrt)-alog(vrh)
//                           gam(ii) = gam(ii) + dble(gamma*gamma)
//                     endif
//                     endif
//               endif
//            end do
//c
//c Madogram:
//c
//      else if(it.eq.8) then
//            do il=lagbeg,lagend
//               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
//               np(ii)  = np(ii)  + 1.
//               dis(ii) = dis(ii) + dble(h)
//               tm(ii)  = tm(ii)  + dble(vrt)
//               hm(ii)  = hm(ii)  + dble(vrh)
//               gam(ii) = gam(ii) + dble(abs(vrh-vrt))
//               if(omni) then
//                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
//     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
//                           np(ii)  = np(ii)  + 1.
//                           dis(ii) = dis(ii) + dble(h)
//                           tm(ii)  = tm(ii)  + dble(vrtpr)
//                           hm(ii)  = hm(ii)  + dble(vrhpr)
//                           gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
//                     endif
//               endif
//            end do
//      endif



				}			
			}		
		}
	}

	return 0;
//end routine
}


//int extractstatisticscwrapper_(
int gamvOMP(
//      integer nd,irepo,maxdat,MAXVAR
	int *nd, int *irepo, int *maxdat, int *MAXVAR,
//      real x(maxdat),y(maxdat),z(maxdat)
	float **x, float **y, float **z,
//      real EPSLON
	float *EPSLON,
//      integer nlag
	int *nlag,
//      real xlag,xltol
	float *xlag, float *xltol,
//      integer mxdlv
	int *mxdlv,
//      real*8 np(mxdlv),dis(mxdlv),gam(mxdlv),hm(mxdlv),
//     + tm(mxdlv),hv(mxdlv),tv(mxdlv)
	double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv,
//      integer numThreads
	int *numThreads,
//      real*8 reducedVariables(7,mxdlv,numThreads)
	double **reducedVariables,
//      real dismxs,tmax,tmin
	float *dismxs, float *tmax, float *tmin,
//      integer ndir,nvarg
	int *ndir, int *nvarg,
//      real uvxazm(100),uvyazm(100),uvzdec(100),uvhdec(100)
	float **uvxazm, float **uvyazm, float **uvzdec, float **uvhdec,
//      real csatol(100),csdtol(100),bandwh(ndir),bandwd(ndir)
	float **csatol, float **csdtol, float **bandwh, float **bandwd,
//      real atol(ndir)
	float **atol,
//      integer ivtype(nvarg),ivtail(nvarg),ivhead(nvarg)
	int **ivtype, int **ivtail, int **ivhead,
//      real vr(maxdat,MAXVAR)
	float **vr	
	)
{

	int threadId,i,j,id,ii,il,it,iv,jj;
	float dx,dy,dz,dxs,dys,dzs,hs,h,dxy;
	int lagbeg,lagend,ilag;
	float band,dcazm,dcdec,gamma,vrh,vrhpr,vrt,vrtpr;
	int omni;


	//printf("Entering gamv calculation. nd=%d EPSLON=%f xlag=%f xltol=%f\n",*nd,*EPSLON,*xlag,*xltol);

	double *nplocal,*dislocal,*gamlocal,*hmlocal,*tmlocal,*hvlocal,*tvlocal;
	double *reducedVariableslocal;
	
	reducedVariableslocal=(double *)calloc(7 * *mxdlv * *numThreads,sizeof(double));

#pragma omp parallel default(shared) \
	private(threadId,i,j,id,iv,il,it,ii,jj,\
		lagbeg,lagend,ilag,omni,\
		dx,dy,dz,dxs,dys,dzs,hs,h,dxy,\
		dcazm,band,dcdec,gamma,vrh,vrt,vrhpr,vrtpr,\
		nplocal,dislocal,gamlocal,hmlocal,tmlocal,hvlocal,tvlocal)
{
#ifdef _OPENMP
	threadId = omp_get_thread_num();
	//numThreads=omp_get_num_threads(); 
#else
	threadId = 0;
	//numThreads=0;
#endif

	nplocal=(double *)calloc(*mxdlv,sizeof(double));
	dislocal=(double *)calloc(*mxdlv,sizeof(double));
	gamlocal=(double *)calloc(*mxdlv,sizeof(double));
	hmlocal=(double *)calloc(*mxdlv,sizeof(double));
	tmlocal=(double *)calloc(*mxdlv,sizeof(double));
	hvlocal=(double *)calloc(*mxdlv,sizeof(double));
	tvlocal=(double *)calloc(*mxdlv,sizeof(double));

	//printf("threadId=%d, numThreads=%d\n",threadId,*numThreads);



#pragma omp for schedule(runtime) 
	for(i=0;i<*nd;i++){
        	//if((int)(i/(*irepo)* (*irepo))==i) printf("   currently on %d of %d\n",i,*nd);
        	//if((int)(i/(*irepo)* (*irepo))==i && threadId==0) printf("   currently on %d of %d\n",i,*nd);
        	//printf("   currently on %d of %d\n",i,*nd);
		for(j=i;j<*nd;j++){

//
// Definition of the lag corresponding to the current pair:
//
			dx  = (*x)[j] - (*x)[i];
			dy  = (*y)[j] - (*y)[i];
			dz  = (*z)[j] - (*z)[i];
			dxs = dx*dx;
			dys = dy*dy;
			dzs = dz*dz;
			hs  = dxs + dys + dzs;
			//printf("i=%d\tj=%d\tx=(%f,%f)\ty=(%f,%f)\tz=(%f,%f)\tdx=%f\tdy=%f\tdz=%f\tdismxs=%f\n",i,j,(*x)[j],(*x)[i],(*y)[j],(*y)[i],(*z)[j],(*z)[i],dx,dy,dz,*dismxs);
			if(hs > *dismxs) continue;
			if(hs < 0.0) hs = 0.0;
			h   = sqrtf(hs);


//
// Determine which lag this is and skip if outside the defined distance
// tolerance:
//
			if(h<=*EPSLON){
				lagbeg = 1;
				lagend = 1;
			}
			else{
				lagbeg = -1;
				lagend = -1;
				for(ilag=2;ilag<=*nlag+2;ilag++){
					if(h>=(*xlag*(float)(ilag-2)-*xltol) && h<=(*xlag*(float)(ilag-2)+*xltol)){
						if(lagbeg<0) lagbeg = ilag; 
						lagend = ilag; 
					}
				}
				if(lagend<0) continue;
			}

			//printf("i=%d j=%d dx=%f dy=%f dz=%fh=%f lagbeg=%d lagend=%d\n",i,j,dx,dy,dz,h,lagbeg,lagend);




//
// Definition of the direction corresponding to the current pair. All
// directions are considered (overlapping of direction tolerance cones
// is allowed):
//
			for(id=0;id<*ndir;id++){
			//
			// Check for an acceptable azimuth angle:
			//
				dxy = sqrtf(MAX((dxs+dys),0.0));
				if(dxy<*EPSLON){
					dcazm = 1.0;
				}
				else{
					dcazm = (dx*(*uvxazm)[id]+dy*(*uvyazm)[id])/dxy;
				}
				if(fabs(dcazm)<(*csatol)[id]) continue;

//
// Check the horizontal bandwidth criteria (maximum deviation 
// perpendicular to the specified direction azimuth):
//
				band = (*uvxazm)[id]*dy - (*uvyazm)[id]*dx;
				if(fabs(band)>=(*bandwh)[id]) continue;

				//printf("dxy=%f\tdcazm=%f\tband=%f\tlagbed=%d\ta=%f\tb=%f\tc=%f\n",dxy,dcazm,band,lagbeg,(*uvhdec)[id],(*uvzdec)[id],(*csdtol)[id]);


//
// Check for an acceptable dip angle:
//
				if(dcazm<0.0) dxy = -dxy;
				//printf("dxy=%f\tdcazm=%f\tband=%f\n",dxy,dcazm,band);
				if(lagbeg==1)
					dcdec = 0.0;
				else{
					dcdec = (dxy*(*uvhdec)[id]+dz*(*uvzdec)[id])/h;
					if(fabs(dcdec)<(*csdtol)[id]) continue;
				}
				//printf("dxy=%f\tdcazm=%f\tband=%f\tdcdec=%f\n",dxy,dcazm,band,dcdec);
//
// Check the vertical bandwidth criteria (maximum deviation perpendicular
// to the specified dip direction):
//
				band = (*uvhdec)[id]*dz - (*uvzdec)[id]*dxy;
				if(fabs(band)>(*bandwd)[id]) continue;
//
// Check whether or not an omni-directional variogram is being computed:
//
				omni = 0;
				if((*atol)[id]>=90.0) omni = 1;
//
// This direction is acceptable - go ahead and compute all variograms:
//

				//printf("dxy=%f dcazm=%f uvxazm[0]=%f uvyazm[0]=%f band=%f dcdec=%f omni=%d csdtol[0]=%f\n",dxy,dcazm,(*uvxazm)[0],(*uvyazm)[0],band,dcdec,omni,(*csdtol)[0]);




//				fprintf(stdout,"dcazm=%f\tdcdec=%f\n",dcazm,dcdec);
				for(iv=0;iv<*nvarg;iv++){
//
// For this variogram, sort out which is the tail and the head value:
//
					it = (*ivtype)[iv];
					//printf("it=%d\n",it);
					//printf("i=%d,j=%d,id=%d,iv=%d,vr[%d]=%f\tvr[%d]=%f\n",i,j,id,iv,i,(*vr)[i],j,(*vr)[j]);
					if(dcazm>=0.0 && dcdec>=0.0){
						ii = (*ivtail)[iv]-1;
						vrh   = (*vr)[i+ii*(*maxdat)];
						ii = (*ivhead)[iv]-1;
						vrt   = (*vr)[j+ii*(*maxdat)];
						if(omni || it==2){
							ii    = (*ivhead)[iv]-1;
							vrtpr = (*vr)[i+ii*(*maxdat)];
							ii    = (*ivtail)[iv]-1;
							vrhpr = (*vr)[j+ii*(*maxdat)];
						}
					}
					else{
						ii = (*ivtail)[iv]-1;
						vrh   = (*vr)[j+ii*(*maxdat)];
						ii = (*ivhead)[iv]-1;
						vrt   = (*vr)[i+ii*(*maxdat)];
						if(omni || it==2){
							ii    = (*ivhead)[iv]-1;
							vrtpr = (*vr)[j+ii*(*maxdat)];
							ii    = (*ivtail)[iv]-1;
							vrhpr = (*vr)[i+ii*(*maxdat)];
						}
					}
					//printf("vrh=%f\tvrt=%f\n",vrh,vrt);



//
// Reject this pair on the basis of missing values:
//
					if(vrt<*tmin || vrh<*tmin || vrt>*tmax || vrh>*tmax) continue;
					if(it==2 && (vrtpr<*tmin || vrhpr<*tmin || vrtpr>*tmax || vrhpr>*tmax)) continue;
//
//             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
//
//
// The Semivariogram:
//


					if(it==1 || it==5 || it>=9){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;


							nplocal[ii]  = nplocal[ii]  + 1.0;
							dislocal[ii] = dislocal[ii] + (double)(h);
							tmlocal[ii]  = tmlocal[ii]  + (double)(vrt);
							hmlocal[ii]  = hmlocal[ii]  + (double)(vrh);
							gamlocal[ii] = gamlocal[ii] + (double)((vrh-vrt)*(vrh-vrt));
							
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									nplocal[ii]  = nplocal[ii]  + 1.0;
									dislocal[ii] = dislocal[ii] + (double)(h);
									tmlocal[ii]  = tmlocal[ii]  + (double)(vrtpr);
									hmlocal[ii]  = hmlocal[ii]  + (double)(vrhpr);
									gamlocal[ii] = gamlocal[ii] + (double)((vrhpr-vrtpr)*(vrhpr-vrtpr));
								}
							}
							
							//printf("ii=%d\tnp=%f\n",ii,(*np)[ii]);
							//printf("ii=%d\tdis=%f\n",ii,(*dis)[ii]);
							//printf("ii=%d\ttm=%f\n",ii,(*tm)[ii]);
							//printf("ii=%d\thm=%f\n",ii,(*hm)[ii]);
							//printf("ii=%d\tgam=%f\n",ii,(*gam)[ii]);
							//printf("ii=%d\ttv=%f\n",ii,(*tv)[ii]);
							//printf("ii=%d\thv=%f\n",ii,(*hv)[ii]);
							
						}
					}
/*

//
// The Traditional Cross Semivariogram:
//
					else if(it==2){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							(*np)[ii]  = (*np)[ii]  + 1.0;
							(*dis)[ii] = (*dis)[ii] + (double)(h);
							(*tm)[ii]  = (*tm)[ii]  + (double)(0.5*(vrt+vrtpr));
							(*hm)[ii]  = (*hm)[ii]  + (double)(0.5*(vrh+vrhpr));
							(*gam)[ii] = (*gam)[ii] + (double)((vrhpr-vrh)*(vrt-vrtpr));
						}
					}
//
// The Covariance:
//
					else if(abs(it)==3){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							(*np)[ii]  = (*np)[ii]  + 1.0;
							(*dis)[ii] = (*dis)[ii] + (double)(h);
							(*tm)[ii]  = (*tm)[ii]  + (double)(vrt);
							(*hm)[ii]  = (*hm)[ii]  + (double)(vrh);
							(*gam)[ii] = (*gam)[ii] + (double)(vrh*vrt);
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									(*np)[ii]  = (*np)[ii]  + 1.0;
									(*dis)[ii] = (*dis)[ii] + (double)(h);
									(*tm)[ii]  = (*tm)[ii]  + (double)(vrtpr);
									(*hm)[ii]  = (*hm)[ii]  + (double)(vrhpr);
									(*gam)[ii] = (*gam)[ii] + (double)(vrhpr*vrtpr);
								}
							}
						}
					}
//
// The Correlogram:
//
					else if(abs(it)==4){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							(*np)[ii]  = (*np)[ii]  + 1.0;
							(*dis)[ii] = (*dis)[ii] + (double)(h);
							(*tm)[ii]  = (*tm)[ii]  + (double)(vrt);
							(*hm)[ii]  = (*hm)[ii]  + (double)(vrh);
							(*hv)[ii]  = (*hv)[ii]  + (double)(vrh*vrh);
							(*tv)[ii]  = (*tv)[ii]  + (double)(vrt*vrt);
							(*gam)[ii] = (*gam)[ii] + (double)(vrh*vrt);
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									(*np)[ii]  = (*np)[ii]  + 1.0;
									(*dis)[ii] = (*dis)[ii] + (double)(h);
									(*tm)[ii]  = (*tm)[ii]  + (double)(vrtpr);
									(*hm)[ii]  = (*hm)[ii]  + (double)(vrhpr);
									(*hv)[ii]  = (*hv)[ii]  + (double)(vrhpr*vrhpr);
									(*tv)[ii]  = (*tv)[ii]  + (double)(vrtpr*vrtpr);
									(*gam)[ii] = (*gam)[ii] + (double)(vrhpr*vrtpr);
								}
							}
						}
					}
//
// The Pairwise Relative:
//
					else if(it==6){
						for(il=lagbeg;il<=lagend;il++){
							ii = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;

							if(abs(vrt+vrh)>*EPSLON){
								(*np)[ii]  = (*np)[ii]  + 1.0;
								(*dis)[ii] = (*dis)[ii] + (double)(h);
								(*tm)[ii]  = (*tm)[ii]  + (double)(vrt);
								(*hm)[ii]  = (*hm)[ii]  + (double)(vrh);
								gamma   = 2.0*(vrt-vrh)/(vrt+vrh);
								(*gam)[ii] = (*gam)[ii] + (double)(gamma*gamma);
							}
							if(omni){
								if(vrtpr>=*tmin && vrhpr>=*tmin && vrtpr<*tmax && vrhpr<*tmax){
									if(abs(vrtpr+vrhpr)>*EPSLON){
										(*np)[ii]  = (*np)[ii]  + 1.0;
										(*dis)[ii] = (*dis)[ii] + (double)(h);
										(*tm)[ii]  = (*tm)[ii]  + (double)(vrtpr);
										(*hm)[ii]  = (*hm)[ii]  + (double)(vrhpr);
										gamma   = 2.0*(vrt-vrh)/(vrt+vrh);
										(*gam)[ii] = (*gam)[ii] + (double)(gamma*gamma);
									}
								}
							}
						}
					}
//c
//c Variogram of Logarithms:
//c
//      else if(it.eq.7) then
//            do il=lagbeg,lagend
//               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
//               if(vrt.gt.EPSLON.and.vrh.gt.EPSLON) then
//                     np(ii)  = np(ii)  + 1.
//                     dis(ii) = dis(ii) + dble(h)
//                     tm(ii)  = tm(ii)  + dble(vrt)
//                     hm(ii)  = hm(ii)  + dble(vrh)
//                     gamma   = alog(vrt)-alog(vrh)
//                     gam(ii) = gam(ii) + dble(gamma*gamma)
//               endif
//               if(omni) then
//                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
//     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
//                     if(vrtpr.gt.EPSLON.and.vrhpr.gt.EPSLON) then
//                           np(ii)  = np(ii)  + 1.
//                           dis(ii) = dis(ii) + dble(h)
//                           tm(ii)  = tm(ii)  + dble(vrtpr)
//                           hm(ii)  = hm(ii)  + dble(vrhpr)
//                           gamma   = alog(vrt)-alog(vrh)
//                           gam(ii) = gam(ii) + dble(gamma*gamma)
//                     endif
//                     endif
//               endif
//            end do
//c
//c Madogram:
//c
//      else if(it.eq.8) then
//            do il=lagbeg,lagend
//               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
//               np(ii)  = np(ii)  + 1.
//               dis(ii) = dis(ii) + dble(h)
//               tm(ii)  = tm(ii)  + dble(vrt)
//               hm(ii)  = hm(ii)  + dble(vrh)
//               gam(ii) = gam(ii) + dble(abs(vrh-vrt))
//               if(omni) then
//                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
//     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
//                           np(ii)  = np(ii)  + 1.
//                           dis(ii) = dis(ii) + dble(h)
//                           tm(ii)  = tm(ii)  + dble(vrtpr)
//                           hm(ii)  = hm(ii)  + dble(vrhpr)
//                           gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
//                     endif
//               endif
//            end do
//      endif





*/
				}			
			}		
		}
	}


//#ifdef _OPENMP
	//if(*numThreads>1){
		for(i=0;i<*mxdlv;i++){
			//printf("%d: hola %d\n",threadId,i);
			reducedVariableslocal[0+i*7+threadId*(*mxdlv * 7)]=dislocal[i];
			reducedVariableslocal[1+i*7+threadId*(*mxdlv * 7)]=gamlocal[i];
			reducedVariableslocal[2+i*7+threadId*(*mxdlv * 7)]=nplocal[i];
			reducedVariableslocal[3+i*7+threadId*(*mxdlv * 7)]=hmlocal[i];
			reducedVariableslocal[4+i*7+threadId*(*mxdlv * 7)]=tmlocal[i];
			reducedVariableslocal[5+i*7+threadId*(*mxdlv * 7)]=hvlocal[i];
			reducedVariableslocal[6+i*7+threadId*(*mxdlv * 7)]=tvlocal[i];
		}                
	//}        
//#endif

	free(nplocal);
	free(dislocal);
	free(gamlocal);
	free(hmlocal);
	free(tmlocal);
	free(hvlocal);
	free(tvlocal);

//end omp parallel
}

//#ifdef _OPENMP
	//if(*numThreads>1){
		for(i=0;i<*mxdlv;i++){
		      (*dis)[i]=0.0; 
		      (*gam)[i]=0.0;
		      (*np)[i]=0.0;  
		      (*hm)[i]=0.0;  
		      (*tm)[i]=0.0;  
		      (*hv)[i]=0.0;  
		      (*tv)[i]=0.0;
		}       
		for(ii=0;ii<*numThreads;ii++){
		   for(jj=0;jj<*mxdlv;jj++){
		      (*dis)[jj] = (*dis)[jj] + reducedVariableslocal[0+jj*7+ii*(*mxdlv * 7)];
		      (*gam)[jj] = (*gam)[jj] + reducedVariableslocal[1+jj*7+ii*(*mxdlv * 7)];
		      (*np)[jj]  = (*np)[jj]  + reducedVariableslocal[2+jj*7+ii*(*mxdlv * 7)];
		      (*hm)[jj]  = (*hm)[jj]  + reducedVariableslocal[3+jj*7+ii*(*mxdlv * 7)];
		      (*tm)[jj]  = (*tm)[jj]  + reducedVariableslocal[4+jj*7+ii*(*mxdlv * 7)];
		      (*hv)[jj]  = (*hv)[jj]  + reducedVariableslocal[5+jj*7+ii*(*mxdlv * 7)];
		      (*tv)[jj]  = (*tv)[jj]  + reducedVariableslocal[6+jj*7+ii*(*mxdlv * 7)];
		   }
		}
		
	//}
//#endif


	free(reducedVariableslocal);

	return 0;
//end routine
}





int gamvAverages(int *ndir, int *nvarg, int *nlag, int *isill, 
			int **ivtype, int **ivtail, int **ivhead, double **sills,
			float *EPSLON,
			double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv){

	int id,iv,il,i,it,iii;
	double rnum,htave,variance;

	for(id=0;id<*ndir;id++){
	for(iv=0;iv<*nvarg;iv++){
	for(il=1;il<=*nlag+2;il++){
		i = (id)*(*nvarg)*((*nlag)+2)+(iv)*((*nlag)+2)+il -1;
		if((*np)[i]<=0.0) continue;
		rnum   = (*np)[i];
		(*dis)[i] = (*dis)[i] / (double)(rnum);
		(*gam)[i] = (*gam)[i] / (double)(rnum);
		(*hm)[i]  = (*hm)[i]  / (double)(rnum);
		(*tm)[i]  = (*tm)[i]  / (double)(rnum);
		(*hv)[i]  = (*hv)[i]  / (double)(rnum);
		(*tv)[i]  = (*tv)[i]  / (double)(rnum);
		it = (*ivtype)[iv];
//
// Attempt to standardize:
//
		if(*isill==1){
			if((*ivtail)[iv]==(*ivhead)[iv]){
				iii = (*ivtail)[iv]-1;
				if((it==1 || it>=9) && (*sills)[iii]>0.0)
					(*gam)[i] = (*gam)[i] / (*sills)[iii];
			}
		}
//
//  1. report the semivariogram rather than variogram
//  2. report the cross-semivariogram rather than variogram
//  3. the covariance requires "centering"
//  4. the correlogram requires centering and normalizing
//  5. general relative requires division by lag mean
//  6. report the semi(pairwise relative variogram)
//  7. report the semi(log variogram)
//  8. report the semi(madogram)
//
		if(it==1 || it==2)
			(*gam)[i] = 0.5 * (*gam)[i];
		else if(abs(it)==3){
			(*gam)[i] = (*gam)[i] - (*hm)[i]*(*tm)[i];
			if(it<0){
				if((*sills)[(*ivtail)[iv]]<0.0 || (*sills)[(*ivhead)[iv]]<0.0)
					(*gam)[i] = -999.0;
				else{
					variance =  sqrt((*sills)[(*ivtail)[iv]])*sqrt((*sills)[(*ivhead)[iv]]) ;
					(*gam)[i] = variance - (*gam)[i];
				}
			}
		}
		else if(abs(it)==4){
			(*hv)[i]  = (*hv)[i]-(*hm)[i]*(*hm)[i];
			if((*hv)[i]<0.0) (*hv)[i] = 0.0;
			(*hv)[i]  = sqrt((*hv)[i]);
			(*tv)[i]  = (*tv)[i]-(*tm)[i]*(*tm)[i];
			if((*tv)[i]<0.0) (*tv)[i] = 0.0;
			(*tv)[i]  = sqrt((*tv)[i]);
			if(((*hv)[i]*(*tv)[i])<*EPSLON)
				(*gam)[i] = 0.0;
			else
				(*gam)[i] =((*gam)[i]-(*hm)[i]*(*tm)[i])/((*hv)[i]*(*tv)[i]);
			
			if(it<0) (*gam)[i] = 1.0 - (*gam)[i];
//
// Square "hv" and "tv" so that we return the variance:
//
			(*hv)[i]  = (*hv)[i]*(*hv)[i];
			(*tv)[i]  = (*tv)[i]*(*tv)[i];
		}
		else if(it==5){
			htave  = 0.5*((*hm)[i]+(*tm)[i]);
			htave  = htave   *   htave;
			if(htave<*EPSLON)
				(*gam)[i] = 0.0;
			else
				(*gam)[i] = (*gam)[i]/htave;
		}
		else if(it>=6)
			(*gam)[i] = 0.5 * (*gam)[i];
		
	}
	}
	}

}


int gamvFreeMemory(
		float **x, float **y, float **z,
		double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv,
		float **uvxazm, float **uvyazm, float **uvzdec, float **uvhdec,
		float **csatol, float **csdtol, float **bandwh, float **bandwd,
		float **atol,
		int **ivtype, int **ivtail, int **ivhead,
		float **vr,
		double **sills,
		double **reducedVariables,
		int *numThreads
		){

	free(*x);free(*y);free(*z);
	free(*np);free(*dis);free(*gam);free(*hm);free(*tm);free(*hv);free(*tv);
	free(*uvxazm);free(*uvyazm);free(*uvzdec);free(*uvhdec);
	free(*csatol);free(*csdtol);free(*bandwh);free(*bandwd);
	free(*atol);
	free(*ivtype);free(*ivtail);free(*ivhead);
	free(*vr);
	free(*sills);
#ifdef _OPENMP
	if(*numThreads>1) free(*reducedVariables);
#endif
}



