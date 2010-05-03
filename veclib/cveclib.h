#ifndef VECLIB_H
#define VECLIB_H
/*****************************************************************************
 *                           V E C L I B . H
 *
 * $Id$
 *****************************************************************************/

#include <stdarg.h>

#include <cfemdef.h>

#ifndef BLAS			/* Assume vendor-supplied FORTRAN BLAS       */
#define BLAS 3			/* library will be available.  If not,       */
#endif                          /* define BLAS to be 0 before this point.    */

#ifndef LAPACK			/* Assume vendor-supplied FORTRAN LAPACK     */
#define LAPACK 1	        /* library will be available.  If not,       */
#endif                          /* define LAPACK to be 0 before this point.  */

#define NVREG  8		/* Registers reserved for FORTRAN linkage.   */

extern integer _vecIreg[NVREG];
extern char    _vecCreg[NVREG];
extern float   _vecSreg[NVREG];
extern double  _vecDreg[NVREG];

#define SQR(a)     ((a) * (a))
#define SGN(a)     ((a) < 0.0 ?    -1.0 : 1.0)
#define MIN(a, b)  ((a) < (b) ?     (a) : (b))
#define MAX(a, b)  ((a) > (b) ?     (a) : (b))
#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))

#ifndef STR_MAX
#define STR_MAX 2048
#endif

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif

#ifndef TWOPI
#define TWOPI  6.28318530717958647692
#endif

#ifndef PI_180
#define PI_180 0.01745329251994329576
#endif

#define EPSm3   1.0e-3
#define EPSm4   1.0e-4
#define EPSm5   1.0e-5
#define EPSm6   1.0e-6
#define EPSSP   6.0e-7
#define EPSm7   1.0e-7
#define EPSm8   1.0e-8
#define EPSm12  1.0e-12
#define EPSDP   6.0e-14
#define EPSm14  1.0e-14
#define EPSm20  1.0e-20
#define EPSm30  1.0e-30
 
typedef struct  {float  Re, Im;}         complex;
typedef struct  {double Re, Im;}         zomplex;
enum    err_lev {WARNING, ERROR, REMARK};


/* ------------------------------------------------------------------------- *
 * UTILITIES:
 * ------------------------------------------------------------------------- */

char buf[STR_MAX];

void    message (const char *routine, const char *txt, int level);
FILE*   efopen  (const char *file,    const char *mode);
double  dclock  (void);
float   sclock  (void);

void    printDvector (FILE *fp,  integer width, integer prec,
		      integer ntot,  integer inc,   integer nfield, ...);
void    printIvector (FILE *fp,  integer width,           
		      integer ntot,  integer inc,   integer nfield, ...);
void    printSvector (FILE *fp,  integer width, integer prec,
		       integer ntot,  integer inc,   integer nfield, ...);

/* ------------------------------------------------------------------------- *
 * MEMORY MANAGEMENT:
 * ------------------------------------------------------------------------- */

#define tempVector(v, n) double *v = dvector (0, n-1)
#define freeVector(v)    free (v)

complex   *cvector  (integer nl, integer nh);
complex  **cmatrix  (integer rl, integer rh, integer cl, integer ch);
complex ***c3matrix (integer rl, integer rh, integer cl, integer ch,
		     integer dl, integer dh);

double    *dvector  (integer nl, integer nh);
double   **dmatrix  (integer rl, integer rh, integer cl, integer ch);
double  ***d3matrix (integer rl, integer rh, integer cl, integer ch,
		     integer dl, integer dh);

float     *svector  (integer nl, integer nh);
float    **smatrix  (integer rl, integer rh, integer cl, integer ch);
float   ***s3matrix (integer rl, integer rh, integer cl, integer ch,
		     integer dl, integer dh);

integer   *ivector  (integer nl, integer nh);
integer  **imatrix  (integer rl, integer rh, integer cl, integer ch);
integer ***i3matrix (integer rl, integer rh, integer cl, integer ch,
		     integer dl, integer dh);

zomplex   *zvector  (integer nl, integer nh);
zomplex  **zmatrix  (integer rl, integer rh, integer cl, integer ch);
zomplex ***z3matrix (integer rl, integer rh, integer cl, integer ch,
		     integer dl, integer dh);

void freeCvector  (complex   *v, integer nl);
void freeCmatrix  (complex  **m, integer nrl, integer ncl);
void freeC3matrix (complex ***t, integer nrl, integer ncl, integer ndl);

void freeDvector  (double    *v, integer nl);
void freeDmatrix  (double   **m, integer nrl, integer ncl);
void freeD3matrix (double  ***t, integer nrl, integer ncl, integer ndl);

void freeSvector  (float     *v, integer nl);
void freeSmatrix  (float    **m, integer nrl, integer ncl);
void freeS3matrix (float   ***t, integer nrl, integer ncl, integer ndl);

void freeIvector  (integer   *v, integer nl);
void freeImatrix  (integer  **m, integer nrl, integer ncl);
void freeI3matrix (integer ***t, integer nrl, integer ncl, integer ndl);

void freeZvector  (zomplex   *v, integer nl);
void freeZmatrix  (zomplex  **m, integer nrl, integer ncl);
void freeZ3matrix (zomplex ***t, integer nrl, integer ncl, integer ndl);


/* ------------------------------------------------------------------------- *
 * MATHEMATICAL PRIMITIVES:
 * ------------------------------------------------------------------------- */

void dcopy (integer n,
	    const double*  x, integer incx, double*  y, integer incy);
void icopy (integer n, 
	    const integer* x, integer incx, integer* y, integer incy);
void scopy (integer n,
	    const float*   x, integer incx, float*   y, integer incy);

void dfill (integer n, double  alpha, double*  x, integer incx);
void ifill (integer n, integer alpha, integer* x, integer incx);
void sfill (integer n, float   alpha, float*   x, integer incx);

void dneg (integer n, double*  x, integer incx);
void ineg (integer n, integer* x, integer incx);
void sneg (integer n, float*   x, integer incx);

void dvneg (integer n,
	    const double*  x, integer incx, double*  y, integer incy);
void ivneg (integer n,
	    const integer* x, integer incx, integer* y, integer incy);
void svneg (integer n,
	    const float*   x, integer incx, float*   y, integer incy);

void dvsgn (integer n,
	    const double*  x, integer incx, double*  y, integer incy);
void ivsgn (integer n,
	    const integer* x, integer incx, integer* y, integer incy);
void svsgn (integer n,
	    const float*   x, integer incx, float*   y, integer incy);

void dsadd (integer n, double  alpha, const double*  x, integer incx, 
	                                    double*  y, integer incy);
void isadd (integer n, integer alpha, const integer* x, integer incx,
	                                    integer* y, integer incy);
void ssadd (integer n, float   alpha, const float*   x, integer incx,
	                                    float*   y, integer incy);

void dspow (const integer n, const double alpha,
	    const double* x, integer incx, double* y, integer incy);
void sspow (const integer n, const float  alpha,
	    const float*  x, integer incx, float*  y, integer incy);

void dvadd (integer n, const double*  x, integer incx,
	    const double*  y, integer incy, double*  z, integer incz);
void ivadd (integer n, const integer* x, integer incx,
	    const integer* y, integer incy, integer* z, integer incz);
void svadd (integer n, const float*   x, integer incx,
	    const float*   y, integer incy, float*   z, integer incz);

void dssub (integer n, double  alpha, const double*  x, integer incx, 
	                                    double*  y, integer incy);
void issub (integer n, integer alpha, const integer* x, integer incx,
	                                    integer* y, integer incy);
void sssub (integer n, float   alpha, const float*   x, integer incx,
	                                    float*   y, integer incy);

void dvsub (integer n, const double*  x, integer incx,
	    const double*  y, integer incy, double*  z, integer incz);
void ivsub (integer n, const integer* x, integer incx,
	    const integer* y, integer incy, integer* z, integer incz);
void svsub (integer n, const float*   x, integer incx,
	    const float*   y, integer incy, float*   z, integer incz);

void dsmul (integer n, double  alpha, const double*  x, integer incx,
	                                    double*  y, integer incy);
void ismul (integer n, integer alpha, const integer* x, integer incx,
	                                    integer* y, integer incy);
void ssmul (integer n, float   alpha, const float*   x, integer incx,
	                                    float*   y, integer incy);

void dvmul (integer n, const double*  x, integer incx,
	    const double*  y, integer incy, double*  z, integer incz);
void ivmul (integer n, const integer* x, integer incx,
	    const integer* y, integer incy, integer* z, integer incz);
void svmul (integer n, const float*   x, integer incx,
	    const float*   y, integer incy, float*   z, integer incz);

void dsdiv (integer n, double  alpha, const double*  x, integer incx,
	                                    double*  y, integer incy);
void isdiv (integer n, integer alpha, const integer* x, integer incx,
	                                    integer* y, integer incy);
void ssdiv (integer n, float   alpha, const float*   x, integer incx,
	                                    float*   y, integer incy);

void dvrecp (integer n,
	     const double* x, integer incx, double* y, integer incy);
void svrecp (integer n,
	     const float*  x, integer incx, float*  y, integer incy);

void dvdiv (integer n, const double*  x, integer incx,
	    const double*  y, integer incy, double*  z, integer incz);
void svdiv (integer n, const float*   x, integer incx,
	    const float*   y, integer incy, float*   z, integer incz);

void dzero (integer n, double*  x, integer incx);
void izero (integer n, integer* x, integer incx);
void szero (integer n, float*   x, integer incx);


/* ------------------------------------------------------------------------- *
 * OTHER MATHEMATICAL FUNCTIONS:
 * ------------------------------------------------------------------------- */

void dvabs (integer n, const double*  x, integer incx,
                             double*  y, integer incy);
void ivabs (integer n, const integer* x, integer incx,
	                     integer* y, integer incy);
void svabs (integer n, const float*   x, integer incx,
	                     float*   y, integer incy);

void dvamax (integer n, const double*  x, integer incx,
	     const double*  y, integer incy, double*  z, integer incz);
void ivamax (integer n, const integer* x, integer incx,
	     const integer* y, integer incy, integer* z, integer incz);
void svamax (integer n, const float*   x, integer incx,
	     const float*   y, integer incy, float*   z, integer incz);

void dvexp (integer n, const double* x, integer incx, double* y, integer incy);
void svexp (integer n, const float*  x, integer incx, float*  y, integer incy);

void dvlg10 (integer n, const double* x, integer incx,
	                      double* y, integer incy);
void svlg10 (integer n, const float*  x, integer incx,
	                      float*  y, integer incy);

void dvlog (integer n, const double* x, integer incx, double* y, integer incy);
void svlog (integer n, const float*  x, integer incx, float*  y, integer incy);

void dvatan (integer n, const double* x, integer incx,
	                      double* y, integer incy);
void svatan (integer n, const float*  x, integer incx,
	                      float*  y, integer incy);

void dvatn2 (integer n, const double* x, integer incx,
	     const double* y, integer incy, double* z, integer incz);
void svatn2 (integer n, const float*  x, integer incx,
	     const float*  y, integer incy, float*  z, integer incz);

void dvcos (integer n, const double* x, integer incx, double* y, integer incy);
void svcos (integer n, const float*  x, integer incx, float*  y, integer incy);

void dvsin (integer n, const double* x, integer incx, double* y, integer incy);
void svsin (integer n, const float*  x, integer incx, float*  y, integer incy);

void dvsqrt (integer n, const double* x, integer incx,
	                      double* y, integer incy);
void svsqrt (integer n, const float*  x, integer incx,
	                      float*  y, integer incy);

void   raninit  (integer flag);
double dranu    (void);
float  sranu    (void);
double drang    (void);
float  srang    (void);
double dnormal  (double mean, double sdev);
float  snormal  (float  mean, float  sdev);
void   dvrandom (integer n, double* x, integer incx);
void   svrandom (integer n, float*  x, integer incx);
void   dvgauss  (integer n, double* x, integer incx);
void   dsgauss  (integer n, float*  x, integer incx);
void   dvnormal (integer n, double mean, double sdev, double* x, integer incx);
void   svnormal (integer n, float  mean, float  sdev, float*  x, integer incx);

integer ispow2   (integer k);
integer irpow2   (integer k);
void zpreft (integer k,                          zomplex *wtab, integer sign);
void cpreft (integer k,                          complex *wtab, integer sign);
void zfft   (integer n, zomplex *x,
	     integer l, const zomplex *wtab, integer dir);
void cfft   (integer n, complex *x,
	     integer l, const complex *wtab, integer dir);
void dzfft  (integer n, zomplex *x,
	     integer l, const zomplex *wtab, integer dir);
void scfft  (integer n, complex *x,
	     integer l, const complex *wtab, integer dir);
void zpfft  (integer n, zomplex *x,
	     integer l, const zomplex *wtab, integer dir);
void spfft  (integer n, complex *x,
	     integer l, const complex *wtab, integer dir);

void dvhypot (integer n, const double* x, integer incx,
	      const double* y, integer incy, double* z, integer incz);
void svhypot (integer n, const float*  x, integer incx,
	      const float*  y, integer incy, float*  z, integer incz);
void dvmag   (integer n, const double* w, integer incw, 
	      const double* x, integer incx, const double* y,
	      integer incy, double* z, integer incz);
void svmag   (integer n, const float* w, integer incw, 
	      const float*  x, integer incx, const float*  y,
	      integer incy, float*  z, integer incz);

void dvpow (integer n, const double* x, integer incx,
	    const double* y, integer incy, double* z, integer incz);
void svpow (integer n, const float*  x, integer incx,
	    const float*  y, integer incy, float*  z, integer incz);


/* ------------------------------------------------------------------------- *
 * TRIAD OPERATIONS:
 * ------------------------------------------------------------------------- */

void dsvmvt (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvmvt (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvpvt (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvpvt (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvtsp (integer n, double alpha, double beta,
	     const double* x, integer incx, double* y, integer incy);
void ssvtsp (integer n, float  alpha, float  beta,
	     const float*  x, integer incx, float*  y, integer incy);

void dsvtvm (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvtvm (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvtvp (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvtvp (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvvmt (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvvmt (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvvpt (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvvpt (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvvtm (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvvtm (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvvtp (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvvtp (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dsvvtt (integer n, double alpha, const double* x, integer incx,
                                      const double* y, integer incy,
	                                    double* z, integer incz);
void ssvvtt (integer n, float  alpha, const float*  x, integer incx,
                                      const float*  y, integer incy,
	                                    float*  z, integer incz);

void dvvmvt (integer n, const double* w, integer incw, 
	                const double* x, integer incx,
	                const double* y, integer incy,
	                      double* z, integer incz);
void svvmvt (integer n, const float*  w, integer incw,
                        const float*  x, integer incx,
                        const float*  y, integer incy,
	                      float*  z, integer incz);

void dvvpvt (integer n, const double* w, integer incw,
	                const double* x, integer incx,
                        const double* y, integer incy,
	                      double* z, integer incz);
void svvpvt (integer n, const float*  w, integer incw,
	                const float*  x, integer incx,
                        const float*  y, integer incy,
	                      float*  z, integer incz);

void dvvpvt (integer n, const double* w, integer incw,
	                const double* x, integer incx,
                        const double* y, integer incy,
	                      double* z, integer incz);
void svvpvt (integer n, const float*  w, integer incw,
	                const float*  x, integer incx,
                        const float*  y, integer incy,
	                      float*  z, integer incz);

void dvvtvp (integer n, const double* w, integer incw,
	                const double* x, integer incx,
                        const double* y, integer incy,
	                      double* z, integer incz);
void svvtvp (integer n, const float*  w, integer incw,
	                const float*  x, integer incx,
                        const float*  y, integer incy,
	                      float*  z, integer incz);

void dvvtvm (integer n, const double* w, integer incw,
	                const double* x, integer incx,
                        const double* y, integer incy,
	                      double* z, integer incz);
void svvtvm (integer n, const float*  w, integer incw,
	                const float*  x, integer incx,
                        const float*  y, integer incy,
	                      float*  z, integer incz);

void dvvvtm (integer n, const double* w, integer incw,
	                const double* x, integer incx,
                        const double* y, integer incy,
	                      double* z, integer incz);
void svvvtm (integer n, const float*  w, integer incw,
	                const float*  x, integer incx,
                        const float*  y, integer incy,
	                      float*  z, integer incz);

void dvvvtt (integer n, const double* w, integer incw,
	                const double* x, integer incx,
                        const double* y, integer incy,
	                      double* z, integer incz);
void svvvtt (integer n, const float*  w, integer incw,
	                const float*  x, integer incx,
                        const float*  y, integer incy,
	                      float*  z, integer incz);


/* ------------------------------------------------------------------------- *
 * RELATIONAL PRIMITIVE OPERATIONS:
 * ------------------------------------------------------------------------- */

void iseq (integer n, integer alpha,
	   const integer* x, integer incx, integer *y, integer incy);
void dsge (integer n, double  alpha,
	   const double*  x, integer incx, integer *y, integer incy);
void isge (integer n, integer alpha,
	   const integer* x, integer incx, integer *y, integer incy);
void ssge (integer n, float   alpha,
	   const float*   x, integer incx, integer *y, integer incy);
void dsle (integer n, double  alpha,
	   const double*  x, integer incx, integer *y, integer incy);
void isle (integer n, integer alpha,
	   const integer* x, integer incx, integer *y, integer incy);
void ssle (integer n, float   alpha,
	   const float*   x, integer incx, integer *y, integer incy);
void dslt (integer n, double  alpha,
	   const double*  x, integer incx, integer *y, integer incy);
void islt (integer n, integer alpha,
	   const integer* x, integer incx, integer *y, integer incy);
void sslt (integer n, float   alpha,
	   const float*   x, integer incx, integer *y, integer incy);
void dsne (integer n, double  alpha,
	   const double*  x, integer incx, integer *y, integer incy);
void isne (integer n, integer alpha,
	   const integer* x, integer incx, integer *y, integer incy);
void ssne (integer n, float   alpha,
	   const float*   x, integer incx, integer *y, integer incy);


/* ------------------------------------------------------------------------- *
 * REDUCTION FUNCTIONS:
 * ------------------------------------------------------------------------- */

double  dsum   (integer n, const double*  x, integer incx);
integer isum   (integer n, const integer* x, integer incx);
float   ssum   (integer n, const float*   x, integer incx);
integer idmax  (integer n, const double*  x, integer incx);
integer iimax  (integer n, const integer* x, integer incx);
integer ismax  (integer n, const float*   x, integer incx);
integer idmin  (integer n, const double*  x, integer incx);
integer iimin  (integer n, const integer* x, integer incx);
integer ismin  (integer n, const float*   x, integer incx);
integer icount (integer n, const integer* x, integer incx);
integer ifirst (integer n, const integer* x, integer incx);
integer lany   (integer n, const integer* x, integer incx);
integer lisame (integer n, const integer* x, integer incx,
		           const integer* y, integer incy);
integer ldsame (integer n, const double*  x, integer incx,
		           const double*  y, integer incy);
integer lssame (integer n, const float*   x, integer incx,
		           const float*   y, integer incy);

/* ------------------------------------------------------------------------- *
 * CONVERSION PRIMITIVES:
 * ------------------------------------------------------------------------- */

void vdble  (integer n, const float*   x, integer incx,
	                      double*  y, integer incy);
void vsngl  (integer n, const double*  x, integer incx,
	                      float*   y, integer incy);

void dvfloa (integer n, const integer* x, integer incx,
	                      double*  y, integer incy);
void svfloa (integer n, const integer* x, integer incx,
	                      float*   y, integer incy);

integer iformat(void);
void format (char*);
void dbrev  (integer n, const double*  x, integer incx,
	                      double*  y, integer incy);
void ibrev  (integer n, const integer* x, integer incx,
	                      integer* y, integer incy);
void sbrev  (integer n, const float*   x, integer incx,
	                      float*   y, integer incy);


/* ------------------------------------------------------------------------- *
 * MISCELLANEOUS FUNCTIONS:
 *
 * NB: xmxm & xmxv operations are replaced by macros that call BLAS xgemm,
 * xgemv routines (generally faster).
 * ------------------------------------------------------------------------- */

void dscatr (integer n, double*  x, integer *y, double*  z);
void iscatr (integer n, integer* x, integer *y, integer* z);
void sscatr (integer n, float*   x, integer *y, float*   z);

void dgathr (integer n, double*  x, integer *y, double*  z);
void igathr (integer n, integer* x, integer *y, integer* z);
void sgathr (integer n, float*   x, integer *y, float*   z);

void dramp (integer n, double  alpha, double  beta, double*  x, integer incx);
void iramp (integer n, integer alpha, integer beta, integer* x, integer incx);
void sramp (integer n, float   alpha, float   beta, float*   x, integer incx);

void dclip (integer n, const  double alpha,  const double   beta,
	    const double* x,  integer incx,  double*  y, integer incy);
void iclip (integer n, const  integer alpha, const integer  beta,
	    const integer* x, integer incx,  integer* y, integer incy);
void sclip (integer n, const  float alpha,   const float    beta,
	    const float* x,  integer incx,   float*   y, integer incy);

void dclipup (integer n, const  double alpha,
	      const double* x,  integer incx,  double*  y, integer incy);
void iclipup (integer n, const  integer alpha,
	      const integer* x, integer incx,  integer* y, integer incy);
void sclipup (integer n, const  float alpha,
	      const float* x,  integer incx,   float*   y, integer incy);

void dclipdn (integer n, const  double alpha,
	      const double* x,  integer incx,  double*  y, integer incy);
void iclipdn (integer n, const  integer alpha,
	      const integer* x, integer incx,  integer* y, integer incy);
void sclipdn (integer n, const  float alpha,
	      const float* x,  integer incx,   float*   y, integer incy);

void diclip (integer n, const  double alpha,  const double   beta,
	     const double* x,  integer incx,  double*  y, integer incy);
void iiclip (integer n, const  integer alpha, const integer  beta,
	     const integer* x, integer incx,  integer* y, integer incy);
void siclip (integer n, const  float alpha,   const float    beta,
	     const float* x,  integer incx,   float*   y, integer incy);

void dcndst (integer n, const  double*  x, integer incx,
	     const integer* y, integer incy, double*  z, integer incz);
void icndst (integer n, const  integer* x, integer incx,
	     const integer* y, integer incy, integer* z, integer incz);
void scndst (integer n, const  float*   x, integer incx,
	     const integer* y, integer incy, float*   z, integer incz);

void dmask (integer n, const double*  w, integer incw,
	               const double*  x, integer incx,
	               const integer* y, integer incy,
	                     double*  z, integer incz);
void imask (integer n, const integer* w, integer incw,
	               const integer* x, integer incx,
	               const integer* y, integer incy,
	                     integer* z, integer incz);
void smask (integer n, const float*   w, integer incw,
	               const float*   x, integer incx,
	               const integer* y, integer incy,
	                     float*   z, integer incz);

void dvpoly (integer n, const double* x, integer incx, integer m,
	                const double* c, integer incc, 
		              double* y, integer incy);
void svpoly (integer n, const float*  x, integer incx, integer m,
	                const float*  c, integer incc, 
		              float*  y, integer incy);

double dpoly   (integer n, double x, const double* xp, const double* yp);
float  spoly   (integer n, float  x, const float*  xp, const float*  yp);
void   dpolint (const double* xa, const double* ya, integer n,
	              double x,        double* y,  double* dy);
void   spolint (const float*  xa, const float*  ya, integer n,
	              float   x,        float*  y, float*  dy);

void   dspline  (integer n, double yp1, double ypn,
	         const double* x, const double* y,  double* y2);
void   sspline  (integer n, float yp1, float ypn,
	         const float*  x, const float*  y,  float*  y2);
double dsplint  (integer n, double x,
	         const  double* xa, const double* ya, const double* y2a);
float  ssplint  (integer n, float  x,
	         const  float*  xa, const float*  ya, const float*  y2a);
double dsplquad (const double* xa, const double* ya,   const double* y2a,
		 const integer n,      const double xmin,  const double xmax);
float  ssplquad (const float*  xa, const float*  ya,   const float*  y2a,
		 const integer n,      const float  xmin,  const float  xmax);

void dmxm (double* A, integer nra, double* B, integer nca,
	   double* C, integer ncb);
void smxm (float*  A, integer nra, float*  B, integer nca,
	   float*  C, integer ncb);

void dmxv (double* A, integer nra, double* B, integer nca, double* C);
void smxv (float*  A, integer nra, float*  B, integer nca, float*  C);

void dmxva (double* A, integer iac, integer iar, double* B, integer ib,
	    double* C, integer ic,  integer nra, integer nca);
void smxva (float*  A, integer iac, integer iar, float*  B, integer ib,
	    float*  C, integer ic,  integer nra, integer nca);

void dvvtvvtp (integer n, const double* v, integer incv,
          	          const double* w, integer incw,
	                  const double* x, integer incx,
	                  const double* y, integer incy,
	                        double* z, integer incz);
void svvtvvtp (integer n, const float*  v, integer incv,
	                  const float*  w, integer incw,
	                  const float*  x, integer incx,
	                  const float*  y, integer incy,
	                        float*  z, integer incz);
void dvvtvvtm (integer n, const double* v, integer incv,
          	          const double* w, integer incw,
	                  const double* x, integer incx,
	                  const double* y, integer incy,
	                        double* z, integer incz);
void svvtvvtm (integer n, const float*  v, integer incv,
	                  const float*  w, integer incw,
	                  const float*  x, integer incx,
	                  const float*  y, integer incy,
	                        float*  z, integer incz);

void dsvvttvp (integer n, const double  alpha,
		          const double* w, integer incw,
		          const double* x, integer incx,
		          const double* y, integer incy,
		                double* z, integer incz);
void ssvvttvp (integer n, const float   alpha,
	                  const float*  w, integer incw,
		          const float*  x, integer incx,
		          const float*  y, integer incy,
		                float*  z, integer incz);  

/*****************************************************************************
 * Prototypes and macros for vendor-supplied FORTRAN libraries follow.
 *****************************************************************************/


/* ------------------------------------------------------------------------- *
 * BLAS level 1 selection.
 *
 * Note that the regular (FORTRAN) BLAS routines idamax & isamax have their
 * return values reduced by 1 to conform with usual C practice.
 * ------------------------------------------------------------------------- */

#if BLAS >= 1

double ddot_  (integer* n, double* x, integer* incx, double* y, integer* incy);
double dasum_ (integer* n, double* x, integer* incx);
double dnrm2_ (integer* n, double* x, integer* incx);

void drotg_  (double* a, double* b, double* c, double* s);
void drot_   (integer* n, double* x, integer* incx, double* y, integer* incy,
	      double* c, double* s);
void dswap_  (integer* n, double* x, integer* incx, double* y, integer* incy);
void dscal_  (integer* n, double* alpha, double* x, integer* incx);
void daxpy_  (integer* n, double* alpha, double* x, integer* incx, 
		                         double* y, integer* incy);
integer idamax_ (integer* n, double* x, integer* incx);

float sdot_  (integer* n, float*  x, integer* incx, float*  y, integer* incy);
float sasum_ (integer* n, float*  x, integer* incx);
float snrm2_ (integer* n, float*  x, integer* incx);

void srotg_  (float*  a, float*  b, float*  c, float*  s);
void srot_   (integer* n, float *x, integer* incx, float *y, integer* incy,
	      float *c, float *s);
void sswap_  (integer* n, float *x, integer* incx, float *y, integer* incy);
void sscal_  (integer* n, float *alpha, float *x, integer* incx);
void saxpy_  (integer* n, float *alpha, float *x, integer* incx, 
		                        float *y, integer* incy);
integer isamax_ (integer* n, float *x, integer* incx);


#define drotg(a, b, c, s) ( drotg_(a, b, c, s) )
#define dswap(n, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   dswap_(_vecIreg, x ,_vecIreg+1, y, _vecIreg+2) )
#define dscal(n, alpha, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecDreg[0]=alpha, \
   dscal_(_vecIreg, _vecDreg, x, _vecIreg+1) )
#define daxpy(n, alpha, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, _vecDreg[0]=alpha, \
   daxpy_(_vecIreg, _vecDreg, x, _vecIreg+1, y, _vecIreg+2) )
#define ddot(n, x, incx, y, incy) \
  (_vecIreg[0]=n,_vecIreg[1]=incx,_vecIreg[2]=incy,\
   ddot_ (_vecIreg, x, _vecIreg+1, y, _vecIreg+2) )
#define dasum(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, dasum_(_vecIreg, x, _vecIreg+1) )
#define dnrm2(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, dnrm2_(_vecIreg, x, _vecIreg+1) )
#define drot(n, x, incx, y, incy, c, s) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   _vecDreg[0]=c, _vecDreg[1]=s, \
   drot_(_vecIreg, x, _vecIreg+1, y, _vecIreg+2, _vecDreg, _vecDreg+1) )
#define idamax(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, idamax_(_vecIreg, x, _vecIreg+1)-1 )

#define srotg(a, b, c, s) ( srotg_(a, b, c, s) )
#define sswap(n, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   sswap_(_vecIreg, x ,_vecIreg+1, y, _vecIreg+2) )
#define sscal(n, alpha, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecSreg[0]=alpha, \
   sscal_(_vecIreg, _vecSreg, x, _vecIreg+1) )
#define saxpy(n, alpha, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, _vecSreg[0]=alpha, \
   saxpy_(_vecIreg, _vecSreg, x, _vecIreg+1, y, _vecIreg+2) )
#define sdot(n, x, incx, y, incy) \
  (_vecIreg[0]=n,_vecIreg[1]=incx,_vecIreg[2]=incy,\
   sdot_ (_vecIreg, x, _vecIreg+1, y, _vecIreg+2) )
#define sasum(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, sasum_(_vecIreg, x, _vecIreg+1) )
#define snrm2(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, snrm2_(_vecIreg, x, _vecIreg+1) )
#define srot(n, x, incx, y, incy, c, s) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   _vecSreg[0]=c, _vecSreg[1]=s,                      \
   srot_(_vecIreg, x, _vecIreg+1, y, _vecIreg+2, _vecSreg, _vecSreg+1) )
#define isamax(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, isamax_(_vecIreg, x, _vecIreg+1)-1 )

#endif

/* ------------------------------------------------------------------------- *
 * BLAS level 2 selection.
 * ------------------------------------------------------------------------- */

#if BLAS >= 2

void dgemv_ (char *trans, integer* m, integer* n, double* alpha,
	     double* a, integer* lda,
	     double* x, integer* incx, double* beta, double* y, integer* incy);
void sgemv_ (char *trans, integer* m, integer* n, float*  alpha,
	     float*  a, integer* lda,
	     float*  x, integer* incx, float*  beta, float*  y, integer* incy);

void dger_ (integer* m, integer* n, double* alpha, double* x, integer* incx,
	    double* y, integer* incy, double* a, integer* lda);
void sger_ (integer* m, integer* n, float*  alpha, float*  x, integer* incx,
	    float*  y, integer* incy, float*  a, integer* lda);

void dspmv_(char *uplo, integer* n, double* alpha, double* ap, double* x,
	    integer* incx, double* beta, double* y, integer* incy);
void sspmv_(char *uplo, integer* n, float*  alpha, float*  ap, float*  x,
	    integer* incx, float*  beta, float*  y, integer* incy);

#define dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)           \
  (_vecCreg[0]=trans,_vecIreg[0]=m,_vecIreg[1]=n,_vecIreg[2]=lda, \
   _vecIreg[3]=incx,_vecIreg[4]=incy,_vecDreg[0]=alpha,           \
   _vecDreg[1]=beta,                                              \
   dgemv_(_vecCreg,_vecIreg,_vecIreg+1,_vecDreg,a,_vecIreg+2,x,   \
	  _vecIreg+3,_vecDreg+1,y,_vecIreg+4))
#define sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)           \
  (_vecCreg[0]=trans,_vecIreg[0]=m,_vecIreg[1]=n,_vecIreg[2]=lda, \
   _vecIreg[3]=incx,_vecIreg[4]=incy,_vecSreg[0]=alpha,           \
   _vecSreg[1]=beta,                                              \
   sgemv_(_vecCreg,_vecIreg,_vecIreg+1,_vecSreg,a,_vecIreg+2,x,   \
	  _vecIreg+3,_vecSreg+1,y,_vecIreg+4))

#define dspmv(uplo,n,alpha,a,x,incx,beta,y,incy)        \
  (_vecCreg[0]=uplo,_vecIreg[0]=n,_vecIreg[1]=incx,     \
   _vecIreg[2]=incy,_vecDreg[0]=alpha,_vecDreg[1]=beta, \
   dspmv_(_vecCreg,_vecIreg,_vecDreg,a,x,               \
	  _vecIreg+1,_vecDreg+1,y,_vecIreg+2))
#define sspmv(uplo,n,alpha,a,x,incx,beta,y,incy)        \
  (_vecCreg[0]=uplo,_vecIreg[0]=n,_vecIreg[1]=incx,     \
   _vecIreg[2]=incy,_vecSreg[0]=alpha,_vecSreg[1]=beta, \
   sspmv_(_vecCreg,_vecIreg,_vecSreg,a,x,               \
	  _vecIreg+1,_vecSreg+1,y,_vecIreg+2))

#define dger(m,n,alpha,x,incx,y,incy,a,lda)           \
  (_vecIreg[0]=m,_vecIreg[1]=n,_vecDreg[0]=alpha,     \
   _vecIreg[2]=incx,_vecIreg[3]=incy,_vecIreg[4]=lda, \
   dger_(_vecIreg,_vecIreg+1,_vecDreg,x,_vecIreg+2,y,_vecIreg+3,a,_vecIreg+4))
#define sger(m,n,alpha,x,incx,y,incy,a,lda)           \
  (_vecIreg[0]=m,_vecIreg[1]=n,_vecSreg[0]=alpha,     \
   _vecIreg[2]=incx,_vecIreg[3]=incy,_vecIreg[4]=lda, \
   dger_(_vecIreg,_vecIreg+1,_vecSreg,x,_vecIreg+2,y,_vecIreg+3,a,_vecIreg+4))

#define dmxv(A,nra,x,nca,y)                                       \
  (_vecCreg[0]='T',_vecIreg[0]=nra,_vecIreg[1]=nca,_vecIreg[2]=1, \
   _vecDreg[0]=1.0,_vecDreg[1]=0.0,                               \
   dgemv_(_vecCreg,_vecIreg+1,_vecIreg,_vecDreg,A,_vecIreg+1,     \
	  x,_vecIreg+2,_vecDreg+1,y,_vecIreg+2))
#define smxv(A,nra,x,nca,y)\
  (_vecCreg[0]='T',_vecIreg[0]=nra,_vecIreg[1]=nca,_vecIreg[2]=1, \
   _vecSreg[0]=1.0,_vecSreg[1]=0.0,                               \
   sgemv_(_vecCreg,_vecIreg+1,_vecIreg,_vecSreg,A,_vecIreg+1,     \
	  x,_vecIreg+2,_vecSreg+1,y,_vecIreg+2))

#endif

/* ------------------------------------------------------------------------- *
 * BLAS level 3 selection.
 * ------------------------------------------------------------------------- */

#if BLAS == 3

void dgemm_ (char *ta, char *tb, integer* m, integer* n, integer* k,
	     double* alpha,
	     double* a, integer* lda,
	     double* b, integer* ldb, double* beta,
	     double* c, integer* ldc);
void sgemm_ (char *ta, char *tb, integer* m, integer* n, integer* k,
	     float *alpha,
	     float *a, integer* lda, float *b,
	     integer* ldb, float *beta,
	     float *c, integer* ldc);
 
#define dgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)           \
  (_vecCreg[0]=ta,_vecCreg[1]=tb,_vecIreg[0]=m,_vecIreg[1]=n,     \
   _vecIreg[2]=k,_vecIreg[3]=lda,_vecIreg[4]=ldb,_vecIreg[5]=ldc, \
   _vecDreg[0]=alpha,_vecDreg[1]=beta,                            \
   dgemm_(_vecCreg,_vecCreg+1,_vecIreg,_vecIreg+1,_vecIreg+2,     \
	  _vecDreg,a,_vecIreg+3,b,_vecIreg+4,_vecDreg+1,c,        \
	  _vecIreg+5))
#define sgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)           \
  (_vecCreg[0]=ta,_vecCreg[1]=tb,_vecIreg[0]=m,_vecIreg[1]=n,     \
   _vecIreg[2]=k,_vecIreg[3]=lda,_vecIreg[4]=ldb,_vecIreg[5]=ldc, \
   _vecSreg[0]=alpha,_vecSreg[1]=beta,                            \
   sgemm_(_vecCreg,_vecCreg+1,_vecIreg,_vecIreg+1,_vecIreg+2,     \
	  _vecSreg,a,_vecIreg+3,b,_vecIreg+4,_vecSreg+1,c,        \
	  _vecIreg+5))

#define dmxm(A,nra,B,nca,C,ncb)                                     \
  (_vecCreg[0]='N',_vecCreg[1]='N',_vecIreg[0]=nra,_vecIreg[1]=nca, \
   _vecIreg[2]=ncb,_vecDreg[0]=1.0,_vecDreg[1]=0.0,                 \
   dgemm_(_vecCreg,_vecCreg+1,_vecIreg+2,_vecIreg,_vecIreg+1,       \
	  _vecDreg,B,_vecIreg+2,A,_vecIreg+1,_vecDreg+1,C,_vecIreg+2))
#define smxm(A,nra,B,nca,C,ncb)                                     \
  (_vecCreg[0]='N',_vecCreg[1]='N',_vecIreg[0]=nra,_vecIreg[1]=nca, \
   _vecIreg[2]=ncb,_vecSreg[0]=1.0,_vecSreg[1]=0.0,                 \
   dgemm_(_vecCreg,_vecCreg+1,_vecIreg+2,_vecIreg,_vecIreg+1,       \
	  _vecSreg,B,_vecIreg+2,A,_vecIreg+1,_vecSreg+1,C,_vecIreg+2))


#endif

/* ------------------------------------------------------------------------- *
 * LAPACK selection.
 * ------------------------------------------------------------------------- */

#if LAPACK == 1

/* Factor a general matrix. */

void dgetrf_ (integer* m, integer* n, double* a, integer* lda,
	      integer* ipiv, integer* info);
void sgetrf_ (integer* m, integer* n, float*  a, integer* lda,
	      integer* ipiv, integer* info);

#define dgetrf(m,n,a,lda,ipiv,info)               \
  (_vecIreg[0]=m, _vecIreg[1]=n, _vecIreg[2]=lda, \
   dgetrf_(_vecIreg, _vecIreg+1, a, _vecIreg+2, ipiv, info))
#define sgetrf(m,n,a,lda,ipiv,info)               \
  (_vecIreg[0]=m, _vecIreg[1]=n, _vecIreg[2]=lda, \
   sgetrf_(_vecIreg, _vecIreg+1, a, _vecIreg+2, ipiv, info))

/* Invert a general matrix from factors. */

void dgetri_ (integer* n,    double* a,     integer* lda,   integer* ipiv,  
	      double* work, integer* lwork, integer* info);
void sgetri_ (integer* n,    float*  a,     integer* lda,   integer* ipiv,
	      float*  work, integer* lwork, integer* info);

#define dgetri(n,a,lda,ipiv,work,lwork,info)          \
  (_vecIreg[0]=n, _vecIreg[1]=lda, _vecIreg[2]=lwork, \
   dgetri_(_vecIreg, a, _vecIreg+1, ipiv, work, _vecIreg+2, info))
#define sgetri(n,a,lda,ipiv,work,lwork,info)          \
  (_vecIreg[0]=n, _vecIreg[1]=lda, _vecIreg[2]=lwork, \
   sgetri_(_vecIreg, a, _vecIreg+1, ipiv, work, _vecIreg+2, info))

/* Factor and solve band-symmetric matrix problem. */

void dpbtrf_ (char *uplo, integer* n, integer* kd,
	      double* ab, integer* ldab, integer* info);
void spbtrf_ (char *uplo, integer* n, integer* kd,
	      float*  ab, integer* ldab, integer* info);

#define dpbtrf(uplo,n,kd,ab,ldab,info)                                \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=ldab, \
   dpbtrf_(_vecCreg, _vecIreg, _vecIreg+1, ab, _vecIreg+2, info))
#define spbtrf(uplo,n,kd,ab,ldab,info)                                \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=ldab, \
   spbtrf_(_vecCreg, _vecIreg, _vecIreg+1, ab, _vecIreg+2, info))

void dpbtrs_ (char *uplo, integer* n, integer* kd, integer* nrhs, double* ab,
	      integer* ldab, double* b, integer* ldb, integer* info);
void spbtrs_ (char *uplo, integer* n, integer* kd, integer* nrhs, float*  ab,
	      integer* ldab, float*  b, integer* ldb, integer* info);

#define dpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)                      \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=nrhs,  \
   _vecIreg[3]=ldab, _vecIreg[4]=ldb,                                  \
   dpbtrs_(_vecCreg, _vecIreg, _vecIreg+1, _vecIreg+2, ab, _vecIreg+3, \
	   b, _vecIreg+4, info))
#define spbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)                      \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=nrhs,  \
   _vecIreg[3]=ldab, _vecIreg[4]=ldb,                                  \
   spbtrs_(_vecCreg, _vecIreg, _vecIreg+1, _vecIreg+2, ab, _vecIreg+3, \
	   b, _vecIreg+4, info))

/* Factor and solve packed-symmetric matrix problem. */

void dpptrf_(char *uplo, integer* n, double* ap, integer* info);
void spptrf_(char *uplo, integer* n, double* ap, integer* info);

#define dpptrf(uplo,n,ap,info) \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, dpptrf_(_vecCreg, _vecIreg, ap, info))
#define spptrf(uplo,n,ap,info) \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, spptrf_(_vecCreg, _vecIreg, ap, info))

void dpptrs_(char *uplo, integer* n, integer* nrhs, double* ap,
	     double* b, integer* ldb, integer* info);
void spptrs_(char *uplo, integer* n, integer* nrhs, float*  ap,
	     float*  b, integer* ldb, integer* info);

#define dpptrs(uplo,n,nrhs,ap,b,ldb,info)                              \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=nrhs, _vecIreg[2]=ldb, \
   dpptrs_(_vecCreg, _vecIreg, _vecIreg+1, ap, b, _vecIreg+2, info))
#define spptrs(uplo,n,nrhs,ap,b,ldb,info)                              \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=nrhs, _vecIreg[2]=ldb, \
   spptrs_(_vecCreg, _vecIreg, _vecIreg+1, ap, b, _vecIreg+2, info))

#endif

#endif

