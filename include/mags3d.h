#define TYPE double
#define PI 3.14159265358979323846264338327
#define TWO_PI 6.2831853071795864769252866
#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a <= _b ? _a : _b; })
#define MAXDISTVALS 1000

typedef struct latticeParams_{
	int nx;
	int ny;
	int nz;
	TYPE xlo;
	TYPE ylo;
	TYPE zlo;
	TYPE h;
} latticeParams;


int nweight;
TYPE *weight;

void quickSort(TYPE *a, int lb, int ub); 
int existValue(TYPE *array, int sizearray, TYPE value);
double generateGaussianNoise();
TYPE * genRandomImage(int nx, int ny, int nz);
int freeRandomImage(TYPE *image);
//int getCoordinates(int index, latticeParams params, TYPE *xcoord, TYPE *ycoord, TYPE *zcoord);
int getCoordinates(int i, int j, int k, latticeParams params, TYPE *xcoord, TYPE *ycoord, TYPE *zcoord);
TYPE costFunction(int nlags, TYPE *expVariogram, int *npairs,TYPE *image,int nx,int ny,int nz,int bufferNodes);
TYPE costFunctionGSLIB(int nlags,TYPE *expVariogram, int *npairs,int iter,int useNscore);
TYPE costFunctionFast(int nsiz, double *gam, TYPE *expVariogram, int iter);
int genMovingAverageImage(	
				TYPE *imageZero, 
				TYPE* kernelWeight, 
				int nx, 
				int ny, 
				int nz, 
				int nxExt, 
				int nyExt, 
				int nzExt,
				int neighRadius, 
				int neighSide, 
				latticeParams params,
				int iter,
				int report
			);

int genMovingAverageImageFast(	
				TYPE *imageZero, 
				TYPE* kernelWeight, 
				int nx, 
				int ny, 
				int nz, 
				int nxExt, 
				int nyExt, 
				int nzExt,
				int neighRadius, 
				int neighSide, 
				latticeParams params,
				int iter,
				int report,
				int *nd, 
				int *MAXVAR,
				int *maxdat,
				float **x, 
				float **y, 
				float **z,
				float **vr
			);

void nscore(
				int *nd, 
				int *MAXVAR,
				int *maxdat,
				float **x, 
				float **y, 
				float **z,
				float **vr
			);

int genTargetVariogram(	
				int nx, 
				int ny, 
				int nz, 
				TYPE xlo,
				TYPE ylo,
				TYPE zlo,
				TYPE h,
				TYPE a,
				TYPE *imageZero,
				int useNscore
			);

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
			);

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
			);

int finiteDifferenceGSLIB(TYPE* weights, int nweights, int nlags, TYPE *expVariogram, int *npairs,int iter);

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
	);

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
	);


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
	);




int gamvAverages(int *ndir, int *nvarg, int *nlag, int *isill, 
			int **ivtype, int **ivtail, int **ivhead, double **sills,
			float *EPSLON,
			double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv);

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
		);


int gamvCUDA(
	int *nd, int *irepo, int *maxdat, int *MAXVAR,
	float **x, float **y, float **z,
	float *EPSLON,
	int *nlag,
	float *xlag, float *xltol,
	int *mxdlv,
	double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv,
	int *numThreads,
	double **reducedVariables,
	float *dismxs, float *tmax, float *tmin,
	int *ndir, int *nvarg,
	float **uvxazm, float **uvyazm, float **uvzdec, float **uvhdec,
	float **csatol, float **csdtol, float **bandwh, float **bandwd,
	float **atol,
	int **ivtype, int **ivtail, int **ivhead,
	float **vr
	);


int gamvCUDAoptimized(
	int *nd, int *irepo, int *maxdat, int *MAXVAR,
	float **x, float **y, float **z,
	float *EPSLON,
	int *nlag,
	float *xlag, float *xltol,
	int *mxdlv,
	double **np, double **dis, double **gam, double **hm, double **tm, double **hv, double **tv,
	int *numThreads,
	double **reducedVariables,
	float *dismxs, float *tmax, float *tmin,
	int *ndir, int *nvarg,
	float **uvxazm, float **uvyazm, float **uvzdec, float **uvhdec,
	float **csatol, float **csdtol, float **bandwh, float **bandwd,
	float **atol,
	int **ivtype, int **ivtail, int **ivhead,
	float **vr
	);

int printGSLIB(	
				int nx, 
				int ny, 
				int nz, 
				latticeParams params,
				int iter,
				int report,
				float **vr,
				char *filename
			);

