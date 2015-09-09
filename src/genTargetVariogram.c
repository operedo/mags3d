#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mags3d.h>

int genTargetVariogram(int nx, int ny, int nz, TYPE xlo, TYPE ylo, TYPE zlo, TYPE h, TYPE a, TYPE *imageZero,int useNscore){
	//int nx = atoi(argv[1]);
	//int ny = atoi(argv[2]);
	//int nz = atoi(argv[3]);
	//TYPE xlo = (TYPE)atof(argv[4]);
	//TYPE ylo = (TYPE)atof(argv[5]);
	//TYPE zlo = (TYPE)atof(argv[6]);
	//TYPE h = (TYPE)atof(argv[7]);
	//TYPE a = (TYPE)atof(argv[8]);

//	float *x,*y,*z,*vr;

	int bufferNodes=ceil(a/h);
	int neighRadius=bufferNodes;
	int neighSide=2*neighRadius+1;

	int nxExt=nx+2*bufferNodes;
	int nyExt=ny+2*bufferNodes;
	int nzExt=nz+2*bufferNodes;


	FILE *fpoctave,*fpgslib,*fpdistanceweight;

/*
	x=(float *)malloc(nx*ny*nz*sizeof(float));
	y=(float *)malloc(nx*ny*nz*sizeof(float));
	z=(float *)malloc(nx*ny*nz*sizeof(float));
	vr=(float *)malloc(nx*ny*nz*sizeof(float));
*/
	fpoctave = fopen("targetimage.m","w");
	fpgslib  = fopen("targetimage.dat","w");
	//fpgslib  = fopen("gslibout.nscore","w");
	fpdistanceweight = fopen("targetdistanceweights.dat","w");
	//fpgslib  = fopen("gslibouttmp.dat","w");

//#ifdef OCTAVE
	fprintf(fpoctave,"image=zeros(%d,%d,%d);\n",nx,ny,nz);
//#else
//#ifdef GSLIB
	fprintf(fpgslib,"Test\n4\nX\nY\nZ\ndata\n");
//#endif
//#endif

	latticeParams params;
	params.nx=nx;
	params.ny=ny;
	params.nz=nz;
	params.xlo=xlo;
	params.ylo=ylo;
	params.zlo=zlo;
	params.h=h;

	//TYPE *image;
	//image = genRandomImage(nxExt,nyExt,nzExt);
	//int initialAddress = 	bufferNodes*nyExt*nxExt + 
	//			bufferNodes*nxExt + 
	//			bufferNodes; 
	//TYPE *imageZero = image+initialAddress;

	int i,j,k,ni,nj,nk;
	TYPE sum=0.0;
	TYPE xcoord1,ycoord1,zcoord1;
	TYPE xcoord2,ycoord2,zcoord2;
	TYPE weight;
	TYPE r, rsquared, asquared, asquaredinv, cons;

	asquared = a*a;
	asquaredinv=1.0/asquared;
	//cons = pow((4.0/((TYPE)(neighRadius*neighRadius)*PI)),0.75);
	cons = pow((4.0/(asquared*PI)),0.75);

	TYPE *kernelWeight = (TYPE *)malloc(sizeof(TYPE)*neighSide*neighSide*neighSide);
	//printf("%d\n",neighSide);
	getCoordinates(0,0,0,params,&xcoord1,&ycoord1,&zcoord1);
	//printf("weight=zeros(%d,%d,%d);\n",neighSide,neighSide,neighSide);

	for(nk=-neighRadius;nk<=neighRadius;nk++){
	   for(nj=-neighRadius;nj<=neighRadius;nj++){
	      for(ni=-neighRadius;ni<=neighRadius;ni++){
	         getCoordinates(ni,nj,nk,params,&xcoord2,&ycoord2,&zcoord2);
	         rsquared = 	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
	  			(ycoord1-ycoord2)*(ycoord1-ycoord2)+
	  			(zcoord1-zcoord2)*(zcoord1-zcoord2);
		 r = sqrt(rsquared);
	         *(kernelWeight+
				(ni+neighRadius)+
				(nj+neighRadius)*neighSide+
				(nk+neighRadius)*neighSide*neighSide) 
				= cons * exp(-2.0*rsquared*asquaredinv); 	// gaussian 3D
				//= r<=a*0.5? cons : 0.0 ; 			// circular 2D
				//= (1.0/(r+2.0)) * (1.0/sqrt(2.0*PI*a))*exp(-(r+2.0)/a);	// exponetial 3D
				//= (cons /(r+1.0)) * exp(-r/a);	// exponetial 3D normalized

                 //printf("weight(%d,%d,%d)=%f;\n",(ni+neighRadius)+1,
		 //		(nj+neighRadius)+1,
		 //		(nk+neighRadius)+1,
                 //              cons * exp(-2.0*rsquared*asquaredinv) 
                 //);
                 //printf("weight[%d]=%f;\n",
		 //       	(ni+neighRadius)+
		 //       	(nj+neighRadius)*neighSide+
		 //       	(nk+neighRadius)*neighSide*neighSide,
                 //              cons * exp(-2.0*rsquared*asquaredinv) 
                 //);
                 if(nk<=0 && nj<=0 && ni<=0)
                    fprintf(fpdistanceweight,"%f %f\n",r,
			*(kernelWeight+
			(ni+neighRadius)+
			(nj+neighRadius)*neighSide+
			(nk+neighRadius)*neighSide*neighSide));
	      }
	   }
	}

	for(k=0;k<nz;k++){
	  for(j=0;j<ny;j++){
	    for(i=0;i<nx;i++){
              sum=0.0;
              getCoordinates(i,j,k,params,&xcoord1,&ycoord1,&zcoord1);
	      for(nk=-neighRadius;nk<=neighRadius;nk++){
	        for(nj=-neighRadius;nj<=neighRadius;nj++){
	          for(ni=-neighRadius;ni<=neighRadius;ni++){
                    //getCoordinates(i+ni,j+nj,k+nk,params,&xcoord2,&ycoord2,&zcoord2);
                    //rsquared = 	(xcoord1-xcoord2)*(xcoord1-xcoord2)+
		    //		(ycoord1-ycoord2)*(ycoord1-ycoord2)+
		    //		(zcoord1-zcoord2)*(zcoord1-zcoord2);
                    ////weight = cons * exp(-2.0*rsquared/((TYPE)(neighRadius*neighRadius)));
                    //weight = cons * exp(-2.0*rsquared*asquaredinv);
                    weight= *(kernelWeight+
				(ni+neighRadius)+
				(nj+neighRadius)*neighSide+
				(nk+neighRadius)*neighSide*neighSide);
 
	            sum=sum + weight * 
                    		*(imageZero+(
					(i+ni)+
					(j+nj)*(nxExt)+
					(k+nk)*(nxExt)*(nyExt)
					)
				);
	          }
	        }
	      }
//	            sum= *(imageZero+(
//					(i)+
//					(j)*(nxExt)+
//					(k)*(nxExt)*(nyExt)
//					)
//				);
//#ifdef OCTAVE
              fprintf(fpoctave,"image(%d,%d,%d)= %f;\n",i+1,j+1,k+1,sum);
//#else
//#ifdef GSLIB
/*
		x[i+j*nx+k+nx*ny]=(float)xcoord1;
		y[i+j*nx+k+nx*ny]=(float)ycoord1;
		z[i+j*nx+k+nx*ny]=(float)zcoord1;
		vr[i+j*nx+k+nx*ny]=(float)sum;
*/
              fprintf(fpgslib,"%f %f %f %f\n",xcoord1,ycoord1,zcoord1,sum);
//#endif
//#endif

	    }
	  }
	}

	free(kernelWeight);
	//freeRandomImage(image);

/*
	int nd=nx*ny*nz;
	int maxdat=nd;
	int MAXVAR=1;	
	nscore(&nd,&maxdat,&MAXVAR,x,y,z,vr);
	for(k=0;k<nz;k++)
		for(j=0;j<ny;j++)
			for(i=0;i<nx;i++)
				fprintf(fpgslib,"%f %f %f %f\n",x[i+j*nx+k+nx*ny],y[i+j*nx+k+nx*ny],z[i+j*nx+k+nx*ny],vr[i+j*nx+k+nx*ny]);

*/
	fclose(fpoctave);
	fclose(fpgslib);
	fclose(fpdistanceweight);

/*
	free(x);
	free(y);
	free(z);
	free(vr);
*/

	int sysret;
	//sysret=system("../gslib90/gamv gamv.par > gslib.log 2>&1");
	//sysret=system("../gslib90/gam-io gam.par > gslib.log 2>&1");

	if(useNscore){
		sysret=system("../../gslib90/nscore nscore.par > /dev/null 2>&1");
		//sysret=system("../gslib90/gam-io gam64x64x1.par > gslib.log 2>&1");
		//sysret=system("../../gslib90/gam-io gam128x128x1_nscore.par > gslib.log 2>&1");
		//sysret=system("../../gslib90/gamv gamv_128x128x1_nscore.par > gslib.log 2>&1");
		sysret=system("../../gslib90/gamv gamv_nscore.par > gslib.log 2>&1");
	}
	else{
		//sysret=system("../../gslib90/gam-io gam128x128x1.par > gslib.log 2>&1");
		//sysret=system("../../gslib90/gamv gamv_128x128x1.par > gslib.log 2>&1");
		sysret=system("../../gslib90/gamv gamv.par > gslib.log 2>&1");
	}
	//sysret=system("../gslib90/gam-io gam256x256x1.par > gslib.log 2>&1");
	//sysret=system("../gslib90/gam-io gam128x128x32.par > gslib.log 2>&1");
	sysret=system("tail -n +2 targetimage.out | wc -l > tmplinenum");
	sysret=system("tail -n +2 targetimage.out | awk '{ print $3, $4, $2 }' > tmpdata");
	sysret=system("cat tmplinenum tmpdata > targetvariogram.dat; rm tmplinenum; rm tmpdata");

	return 0;
}
