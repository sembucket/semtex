#ifndef VECLIB_H
#define VECLIB_H

/* VECLIB
 *
 * $Revision$
 *
 * Author:  R. D. Henderson
 *
 * This is the number-crunching library for Prism.  This library also 
 * contains code segments from the BLAS and LAPACK libraries that are needed
 * for manipulating linear systems.
 * ------------------------------------------------------------------------- */

#define NVREG  8   /* Reserved registers for FORTRAN linkage */

extern char   _vlib_creg [NVREG];
extern int    _vlib_ireg [NVREG];
extern float  _vlib_sreg [NVREG]; 
extern double _vlib_dreg [NVREG]; 

typedef struct zcmplx {  /* Type definitions for complex numbers */
  double r;
  double i;
} zcomplex;

#ifndef __OSF1__
typedef struct cmplx {
  float r;
  float i;
} complex;
#endif

#include "veclib/blas.h"
#include "veclib/lapack.h"
 
/* ----------------- SECTION 1: Mathematical Primatives -------------------- */

double drand  (void);

void   ifill  (int n, int    a, int     *x, int incx);
void   dfill  (int n, double a, double  *x, int incx);

void   dsadd  (int n, double a, double  *x, int incx, double  *y, int incy);
void   dsmul  (int n, double a, double  *x, int incx, double  *y, int incy);
void   dsdiv  (int n, double a, double  *x, int incx, double  *y, int incy);

void   izero  (int n, int     *x, int incx);
void   dzero  (int n, double  *x, int incx);
void   sneg   (int n, float   *x, int incx);
void   dneg   (int n, double  *x, int incx);
void   dvrand (int n, double  *x, int incx);

void   dvneg  (int n, double  *x, int incx, double  *y, int incy);
void   dvrecp (int n, double  *x, int incx, double  *y, int incy);
void   dvabs  (int n, double  *x, int incx, double  *y, int incy);
void   dvlg10 (int n, double  *x, int incx, double  *y, int incy);
void   dvcos  (int n, double  *x, int incx, double  *y, int incy);
void   dvsin  (int n, double  *x, int incx, double  *y, int incy);
void   dvsqrt (int n, double  *x, int incx, double  *y, int incy);
void   dvexp  (int n, double  *x, int incx, double  *y, int incy);

void   dvamax (int n, double  *x, int incx, double  *y, int incy,
	              double  *z, int incz);
void   dvdiv  (int n, double  *x, int incx, double  *y, int incy,
	              double  *z, int incz);
void   dvadd  (int n, double  *x, int incx, double  *y, int incy,
	              double  *z, int incz);
void   dvsub  (int n, double  *x, int incx, double  *y, int incy,
	              double  *z, int incz);
void   dvmul  (int n, double  *x, int incx, double  *y, int incy,
	              double  *z, int incz);


/* -------------------- SECTION 2: Triad Operations ------------------------ */

void   dsvtvm (int n, double a, double *x, int incx, double *y, int incy,
	                        double *z, int incz);
void   dsvtvp (int n, double a, double *x, int incx, double *y, int incy,
	                        double *z, int incz);
void   dsvvmt (int n, double a, double *x, int incx, double *y, int incy,
	                        double *z, int incz);
void   dsvvpt (int n, double a, double *x, int incx, double *y, int incy,
	                        double *z, int incz);
void   dsvvtm (int n, double a, double *x, int incx, double *y, int incy,
	                        double *z, int incz);
void   dsvvtp (int n, double a, double *x, int incx, double *y, int incy,
	                        double *z, int incz);
void   dsvtsp (int n, double a, double b, double *x, int incx, 
	                                  double *y, int incy);

void   dvvtvp (int n, double *w, int incw, double *x, int incx, 
	              double *y, int incy, double *z, int incz);
void   dvvtvm (int n, double *w, int incw, double *x, int incx, 
	              double *y, int incy, double *z, int incz);
void   dvvpvt (int n, double *w, int incw, double *x, int incx, 
	              double *y, int incy, double *z, int incz);
void   dvvmvt (int n, double *w, int incw, double *x, int incx, 
	              double *y, int incy, double *z, int incz);
void   dvvvtm (int n, double *w, int incw, double *x, int incx, 
	              double *y, int incy, double *z, int incz);

/* ------------------ SECTION 3: Relational Primitives --------------------- */

void   iseq   (int n, int    *a, int    *x, int incx, int *y, int incy);
void   dsle   (int n, double *a, double *x, int incx, int *y, int incy);
void   dslt   (int n, double *a, double *x, int incx, int *y, int incy);
void   dsne   (int n, double *a, double *x, int incx, int *y, int incy);

/* ------------------- SECTION 4: Reduction Functions ---------------------- */

int    isum   (int n, int    *x, int incx);
double dsum   (int n, double *x, int incx);
double dvpoly (int n, double *x, int incx, int m, double *c, int incc,
	              double *y, int incy);
int    idmax  (int n, double *x, int incx);
int    ismax  (int n, float  *x, int incx);
int    iimax  (int n, int    *x, int incx);
int    idmin  (int n, double *x, int incx);
int    ismin  (int n, float  *x, int incx);
int    icount (int n, int    *x, int incx);
int    ifirst (int n, int    *x, int incx);

/* ------------------ SECTION 5: Conversion Primitives --------------------- */

void   dvfloa (int n, int      *x, int incx, double *y, int incy);
void   svfloa (int n, int      *x, int incx, float  *y, int incy);
void   vsngl  (int n, double   *x, int incx, float  *y, int incy);
void   zvreal (int n, zcomplex *x, int incx, double *y, int incy);
void   zvimag (int n, zcomplex *x, int incx, double *y, int incy);
void   zvcmplx(int n, double   *x, int incx, double *y, int incy,
	              zcomplex *z, int incz);

void   dbrev  (int n, double *x, int incx, double *y, int incy);
void   sbrev  (int n, float  *x, int incx, float  *y, int incy);
void   ibrev  (int n, int    *x, int incx, int    *y, int incy);

/* ----------------- SECTION 6: Miscellaneous Functions -------------------- */

void   dscatr (int n, double  *x, int *y, double  *z);
void   sscatr (int n, float   *x, int *y, float   *z);
void   dgathr (int n, double  *x, int *y, double  *z);
void   sgathr (int n, float   *x, int *y, float   *z);
void   dramp  (int n, double  *a, double *b, double *x, int incx);
void   iramp  (int n, int     *a, int    *b, int    *x, int incx);
void   dcndst (int n, double  *x, int incx, int *y, int incy, 
	              double  *z, int incz);
void   dmask  (int n, double  *w, int incw, double *x, int incx, 
	              int     *y, int incy, double *z, int incz);
double dpoly  (int n, double xp, double *x, double *y);
void   polint (double *x, double *y, int n, double xp, double *v, double *err);
void   spline (int n, double yp1, double ypn, double *x, double *y,double *y2);

double splint (int n, double x, double *xp, double *yp, double *ypp);
double dclock (void);

void   realft (int n, double *f, int dir);
#if !defined(HPUX)
void   fftdf  (double *f, int n, int nskip, int m, int mskip, int isign,
	       double *trig, int irev, int *bitrev, int ierr, int ireal);
#endif

void   mxm    (double *a, int nra, double *b, int nca, double *c, int ncb);
void   mxv    (double *a, int nra, double *b, int nca, double *ab);
void   mxva   (double *a, int, int, double *b, int, 
	       double *c, int, int, int);

/* --------------- SECTION 7: Memory Management Functions ------------------ */

int      *ivector     (int rmin, int rmax);  
float    *svector     (int rmin, int rmax);  
double   *dvector     (int rmin, int rmax);          

int     **imatrix     (int rmin, int rmax, int cmin, int cmax); 
float   **smatrix     (int rmin, int rmax, int cmin, int cmax); 
double  **dmatrix     (int rmin, int rmax, int cmin, int cmax); 

void     free_dvector (double *v, int rmin);
void     free_svector (float  *v, int rmin);
void     free_ivector (int    *v, int rmin);

void     free_dmatrix (double **mat, int rmin, int cmin);
void     free_smatrix (float  **mat, int rmin, int cmin);
void     free_imatrix (int    **mat, int rmin, int cmin);

void     vecerr       (char *sec, int msgid);

/* Macro Translations */

#define fftdf(f,n,nskip,m,mskip,isign,trig,irev,bitrev,ierr,ireal) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=nskip,_vlib_ireg[2]=m,_vlib_ireg[3]=mskip,\
   _vlib_ireg[4]=isign,_vlib_ireg[5]=irev,_vlib_ireg[6]=ireal,\
   fftdf_(f,_vlib_ireg,_vlib_ireg+1,_vlib_ireg+2,_vlib_ireg+3,_vlib_ireg+4,\
	  trig,_vlib_ireg+5,bitrev,&ierr,_vlib_ireg+6))

#if defined(CRAY)
#
  double  SECOND  (void); 
# define  dclock  SECOND
# define  fftdf_  FFTDF
#
# define tempVector(v,n)      auto double v[n]
# define freeVector(v)        /* automatic */
#
#endif

#if defined(HPUX)
# define fftdf_ fftdf
#endif

#ifndef  tempVector
#define  tempVector(v,n)      double *v = dvector(0,n-1)
#define  freeVector(v)        free (v)
#endif

#endif
