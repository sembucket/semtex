#ifndef VECLIB_H
#define VECLIB_H
///////////////////////////////////////////////////////////////////////////////
// Veclib.h:  C++ wrappers for veclib subroutine calls.
//
// Veclib is described in the iPSC/2 Programmer's Reference Manual.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cfemdef.h>

extern "C" {

  /* -- MATHEMATICAL PRIMITIVES */

  void  dcopy  (integer n,
		const double*  x, integer incx, double*  y, integer incy);
  void  icopy  (integer n,
		const integer* x, integer incx, integer* y, integer incy);
  void  scopy  (integer n,
		const float*   x, integer incx, float*   y, integer incy);
  
  void  dfill  (integer n, double  alpha, double*  x, integer incx);
  void  ifill  (integer n, integer alpha, integer* x, integer incx);
  void  sfill  (integer n, float   alpha, float*   x, integer incx);
  
  void  dneg   (integer n, double*  x, integer incx);
  void  ineg   (integer n, integer* x, integer incx);
  void  sneg   (integer n, float*   x, integer incx);
  
  void  dvneg  (integer n,
		const double*  x, integer incx, double*  y, integer incy);
  void  ivneg  (integer n,
		const integer* x, integer incx, integer* y, integer incy);
  void  svneg  (integer n,
		const float*   x, integer incx, float*   y, integer incy);
  
  void  dvsgn  (integer n,
		const double*  x, integer incx, double*  y, integer incy);
  void  ivsgn  (integer n,
		const integer* x, integer incx, integer* y, integer incy);
  void  svsgn  (integer n,
		const float*   x, integer incx, float*   y, integer incy);
  
  void  dsadd  (integer n, double  alpha,
		const double*  x, integer incx, double*  y, integer incy);
  void  isadd  (integer n, integer alpha,
		const integer* x, integer incx, integer* y, integer incy);
  void  ssadd  (integer n, float  alpha,
		const float*   x, integer incx, float*   y, integer incy);
  
  void  dspow  (const integer n, const double alpha,
		const double* x, integer incx, double* y, integer incy);
  void  sspow  (const integer n, const float  alpha,
		const float*  x, integer incx, float*  y, integer incy);
  
  void  dvadd  (integer n, const double*  x, integer incx,
		const double*  y, integer incy, double*  z, integer incz);
  void  ivadd  (integer n, const integer* x, integer incx, 
		const integer* y, integer incy, integer* z, integer incz);
  void  svadd  (integer n, const float*  x, integer incx,
		const float*   y, integer incy, float*   z, integer incz);
  
  void  dssub  (integer n, double  alpha,
		const double*  x, integer incx, double*  y, integer incy);
  void  issub  (integer n, integer alpha,
		const integer* x, integer incx, integer* y, integer incy);
  void  sssub  (integer n, float   alpha,
		const float*   x, integer incx, float*   y, integer incy);
  
  void  dvsub  (integer n, const double*  x, integer incx,
		const double*  y, integer incy, double*  z, integer incz);
  void  ivsub  (integer n, const integer* x, integer incx,
		const integer* y, integer incy, integer* z, integer incz);
  void  svsub  (integer n, const float*   x, integer incx,
		const float*   y, integer incy, float*   z, integer incz);
  
  void  dsmul  (integer n, double alpha,
		const double*  x, integer incx, double*  y, integer incy);
  void  ismul  (integer n, integer alpha,
		const integer* x, integer incx, integer* y, integer incy);
  void  ssmul  (integer n, float  alpha,
		const float*   x, integer incx, float*   y, integer incy);
  
  void  dvmul  (integer n, const double* x, integer incx,
		const double*  y, integer incy, double*  z, integer incz);
  void  ivmul  (integer n, const integer* x, integer incx,
		const integer* y, integer incy, integer* z, integer incz);
  void  svmul  (integer n, const float*   x, integer incx,
		const float*   y, integer incy, float*   z, integer incz);
  
  void  dsdiv  (integer n, double  alpha,
		const double*  x, integer incx, double*  y, integer incy);
  void  isdiv  (integer n, integer alpha,
		const integer* x, integer incx, integer* y, integer incy);
  void  ssdiv  (integer n, float   alpha,
		const float*   x, integer incx, float*   y, integer incy);
  
  void  dvrecp (integer n,
		const double* x, integer incx, double* y, integer incy);
  void  svrecp (integer n,
		const float*  x, integer incx, float*  y, integer incy);
  
  void  dvdiv  (integer n, const double* x, integer incx,
		const double*  y, integer incy, double*  z, integer incz);
  void  svdiv  (integer n, const float*  x, integer incx,
		const float*   y, integer incy, float*   z, integer incz);
  
  void  dzero  (integer n, double*  x, integer incx);
  void  izero  (integer n, integer* x, integer incx);
  void  szero  (integer n, float*   x, integer incx);
  
  /* -- OTHER MATHEMATICAL FUNCTIONS */
  
  void    dvabs    (integer n,
		    const double*  x, integer incx, double*  y, integer incy);
  void    ivabs    (integer n,
		    const integer* x, integer incx, integer* y, integer incy);
  void    svabs    (integer n,
		    const float*   x, integer incx, float*   y, integer incy);
  
  void    dvamax   (integer n, const double*  x, integer incx,
		    const double*  y, integer incy, double*  z, integer incz);
  void    ivamax   (integer n, const integer* x, integer incx,
		    const integer* y, integer incy, integer* z, integer incz);
  void    svamax   (integer n, const float*   x, integer incx,
		    const float*   y, integer incy, float*   z, integer incz);
  
  void    dvexp    (integer n,
		    const double* x, integer incx, double* y, integer incy);
  void    svexp    (integer n,
		    const float*  x, integer incx, float*  y, integer incy);
  
  void    dvlg10   (integer n,
		    const double* x, integer incx, double* y, integer incy);
  void    svlg10   (integer n,
		    const float*  x, integer incx, float*  y, integer incy);
  
  void    dvlog    (integer n,
		    const double* x, integer incx, double* y, integer incy);
  void    svlog    (integer n,
		    const float*  x, integer incx, float*  y, integer incy);
  
  void    dvatan   (integer n,
		    const double* x, integer incx, double* y, integer incy);
  void    svatan   (integer n,
		    const float*  x, integer incx, float*  y, integer incy);
  
  void    dvatn2   (integer n, const double* x, integer incx,
		    const double* y, integer incy, double* z, integer incz);
  void    svatn2   (integer n, const float*  x, integer incx,
		    const float*  y, integer incy, float*  z, integer incz);
  
  void    dvcos    (integer n,
		    const double* x, integer incx, double* y, integer incy);
  void    svcos    (integer n,
		    const float*  x, integer incx, float*  y, integer incy);
  
  void    dvsin    (integer n,
		    const double* x, integer incx, double* y, integer incy);
  void    svsin    (integer n,
		    const float*  x, integer incx, float*  y, integer incy);
  
  void    dvsqrt   (integer n,
		    const double* x, integer incx, double* y, integer incy);
  void    svsqrt   (integer n,
		    const float*  x, integer incx, float*  y, integer incy);
  
  void    raninit  (integer flag);
  double  dranu    (void);
  float   sranu    (void);
  double  dnormal  (double mean, double sdev);
  float   snormal  (float  mean, float  sdev);
  void    dvrandom (integer n, double* x, integer incx);
  void    svrandom (integer n, float*  x, integer incx);
  void    dvnormal (integer n, double mean, double sdev,
		    double* x, integer incx);
  void    svnormal (integer n, float  mean, float  sdev, 
		    float*  x, integer incx);

  void    dvhypot  (integer n, const double* x, integer incx,
		    const double* y, integer incy, double* z, integer incz);
  void    svhypot  (integer n, const float*  x, integer incx,
		    const float*  y, integer incy, float*  z, integer incz);
  void    dvmag    (integer n, const double* w, integer incw, 
		    const double* x, integer incx, const double* y,
		    integer incy, double* z, integer incz);
  void    svmag    (integer n, const float* w, integer incw, 
		    const float*  x, integer incx, const float*  y,
		    integer incy, float*  z, integer incz);

  void    dvpow    (integer n, const double* x, integer incx,
		    const double* y, integer incy, double* z, integer incz);
  void    svpow    (integer n, const float*  x, integer incx,
		    const float*  y, integer incy, float*  z, integer incz);

  
  /* -- TRIAD OPERATIONS */

  void   dsvmvt (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvmvt (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dsvpvt (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvpvt (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dsvtsp (integer n, double alpha, double beta,
		 const double* x, integer incx, double* y, integer incy);
  void   ssvtsp (integer n, float  alpha, float  beta,
		 const float*  x, integer incx, float*  y, integer incy);
  
  void   dsvtvm (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvtvm (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dsvtvp (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvtvp (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dsvvmt (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvvmt (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dsvvpt (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvvpt (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dsvvtm (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvvtm (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dsvvtp (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvvtp (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);

  void   dsvvtt (integer n, double alpha, const double* x, integer incx,
		 const double* y, integer incy, double* z, integer incz);
  void   ssvvtt (integer n, float  alpha, const float*  x, integer incx,
		 const float*  y, integer incy, float*  z, integer incz);
  
  void   dvvmvt (integer n, const double* w, integer incw, const double* x,
		 integer incx, const double* y, integer incy,
		 double* z, integer incz);
  void   svvmvt (integer n, const float*  w, integer incw, const float*  x,
		 integer incx, const float*  y, integer incy,
		 float*  z, integer incz);
  
  void   dvvpvt (integer n, const double* w, integer incw, const double* x,
		 integer incx, const double* y, integer incy,
		 double* z, integer incz);
  void   svvpvt (integer n, const float*  w, integer incw, const float*  x,
		 integer incx, const float*  y, integer incy,
		 float*  z, integer incz);

  void   dvvtvm (integer n, const double* w, integer incw, const double* x,
		 integer incx, const double* y, integer incy,
		 double* z, integer incz);
  void   svvtvm (integer n, const float*  w, integer incw, const float*  x,
		 integer incx, const float*  y, integer incy,
		 float*  z, integer incz);
  
  void   dvvtvp (integer n, const double* w, integer incw, const double* x,
		 integer incx, const double* y, integer incy,
		 double* z, integer incz);
  void   svvtvp (integer n, const float*  w, integer incw, const float*  x,
		 integer incx, const float*  y, integer incy,
		 float*  z, integer incz);
  
  void   dvvvtm (integer n, const double* w, integer incw, const double* x,
		 integer incx, const double* y, integer incy,
		 double* z, integer incz);
  void   svvvtm (integer n, const float*  w, integer incw, const float*  x,
		 integer incx, const float*  y, integer incy,
		 float*  z, integer incz);
  
  /* -- RELATIONAL PRIMITIVE OPERATIONS */

  void   iseq (integer n, integer    alpha,
	       const  integer* x, integer incx, integer *y, integer incy);
  void   dsge (integer n, double alpha,
	       const  double* x, integer incx, integer *y, integer incy);
  void   isge (integer n, integer    alpha,
	       const  integer* x, integer incx, integer *y, integer incy);
  void   ssge (integer n, float  alpha,
	       const  float*  x, integer incx, integer *y, integer incy);
  void   dsle (integer n, double alpha,
	       const  double* x, integer incx, integer *y, integer incy);
  void   isle (integer n, integer    alpha,
	       const  integer* x, integer incx, integer *y, integer incy);
  void   ssle (integer n, float  alpha,
	       const  float*  x, integer incx, integer *y, integer incy);
  void   dslt (integer n, double alpha,
	       const  double* x, integer incx, integer *y, integer incy);
  void   islt (integer n, integer    alpha,
	       const  integer* x, integer incx, integer *y, integer incy);
  void   sslt (integer n, float  alpha,
	       const  float*  x, integer incx, integer *y, integer incy);
  void   dsne (integer n, double alpha,
	       const  double* x, integer incx, integer *y, integer incy);
  void   isne (integer n, integer    alpha,
	       const  integer* x, integer incx, integer *y, integer incy);
  void   ssne (integer n, float  alpha,
	       const  float*   x, integer incx, integer *y, integer incy);
  
  /* -- REDUCTION FUNCTIONS */

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

  /* -- CONVERSION PRIMITIVES */
  
  void   vdble   (integer n,
		  const float*   x, integer incx, double* y, integer incy);
  void   vsngl   (integer n,
		  const double*  x, integer incx, float*  y, integer incy);
  void   dvfloa  (integer n,
		  const integer* x, integer incx, double* y, integer incy);
  void   svfloa  (integer n,
		  const integer* x, integer incx, float*  y, integer incy);


  integer iformat (void);
  void    format  (char*);
  void    dbrev   (integer n,
		   const double*  x, integer incx, double*  y, integer incy);
  void    ibrev   (integer n,
		   const integer* x, integer incx, integer* y, integer incy);
  void    sbrev   (integer n,
		   const float*   x, integer incx, float*   y, integer incy);
  
  /* -- MISCELLANEOUS FUNCTIONS */

  double dclock  (void);
  float  sclock  (void);
  
  void   dscatr  (integer n, const double*  x, const integer* y, double*  z);
  void   iscatr  (integer n, const integer* x, const integer* y, integer* z);
  void   sscatr  (integer n, const float*   x, const integer* y, float*   z);
  
  void   dgathr  (integer n, const double*  x, const integer* y, double*  z);
  void   igathr  (integer n, const integer* x, const integer* y, integer* z);
  void   sgathr  (integer n, const float*   x, const integer* y, float*   z);

  void   dgathr_scatr (integer n, const double*  w, 
		       const integer* x, const integer* y, double*  z);
  void   igathr_scatr (integer n, const integer* w,
		       const integer* x, const integer* y, integer* z);
  void   sgathr_scatr (integer n, const float*   w,
		       const integer* x, const integer* y,  float*  z);
  
  void   dscatr_sum (integer n,
		     const double*  x, const integer *y, double*  z);
  void   iscatr_sum (integer n,
		     const integer* x, const integer* y, integer* z);
  void   sscatr_sum (integer n,
		     const float*   x, const integer *y, float*   z);
  
  void   dgathr_sum (integer n,
		     const double*  x, const integer* y, double* z);
  void   igathr_sum (integer n,
		     const integer* x, const integer* y, integer* z);
  void   sgathr_sum (integer n,
		     const float*   x, const integer* y, float*   z);

  void   dgathr_scatr_sum (integer n, const double*  w, 
			   const integer* x, const integer* y, double*  z);
  void   igathr_scatr_sum (integer n, const integer* w,
			   const integer* x, const integer* y, integer* z);
  void   sgathr_scatr_sum (integer n, const float*   w,
			   const integer* x, const integer* y, float*   z);
  
  void   dramp   (integer n, double  alpha, double  beta,
		  double*  x, integer incx);
  void   iramp   (integer n, integer alpha, integer beta,
		  integer* x, integer incx);
  void   sramp   (integer n, float   alpha, float   beta,
		  float*   x, integer incx);
  
  void   dcndst  (integer n, const  double*  x, integer incx,
		  const integer* y, integer incy, double*  z, integer incz);
  void   icndst  (integer n, const  integer* x, integer incx,
		  const integer* y, integer incy, integer* z, integer incz);
  void   scndst  (integer n, const  float*   x, integer incx,
		  const integer* y, integer incy, float*   z, integer incz);
  
  void   dmask   (integer n,
		  const double*  w, integer incw,
		  const double*  x, integer incx,
		  const integer* y, integer incy,
		        double*  z, integer incz);
  void   imask   (integer n,
		  const integer* w, integer incw,
		  const integer* x, integer incx,
		  const integer* y, integer incy,
		        integer* z, integer incz);
  void   smask   (integer n,
		  const float*   w, integer incw,
		  const float*   x, integer incx,
		  const integer* y, integer incy,
		        float*   z, integer incz);

  void dclip (integer n, const  double alpha,  const double   beta,
	      const double* x,  integer incx,  const double*  y, integer incy);
  void iclip (integer n, const  integer alpha, const integer  beta,
	      const integer* x, integer incx,  const integer* y, integer incy);
  void sclip (integer n, const  float alpha,   const float    beta,
	      const float* x,  integer incx,   const float*   y, integer incy);

  void dclipup (integer n, const double alpha, const double* x,
		integer incx,  const double* y, integer incy);
  void iclipup (integer n, const integer alpha, const integer* x,
		integer incx, const integer* y, integer incy);
  void sclipup (integer n, const float alpha, const float* x,
		integer incx, const float* y, integer incy);

  void dclipdn (integer n, const double alpha, const double* x,
		integer incx,  const double* y, integer incy);
  void iclipdn (integer n, const integer alpha, const integer* x,
		integer incx, const integer* y, integer incy);
  void sclipdn (integer n, const float alpha, const float* x,
		integer incx, const float* y, integer incy);

  void diclip (integer n, const  double alpha,  const double   beta,
	       const double* x,  integer incx, const double*  y, integer incy);
  void iiclip (integer n, const  integer alpha, const integer  beta,
	       const integer* x, integer incx, const integer* y, integer incy);
  void siclip (integer n, const  float alpha,   const float    beta,
	       const float* x,  integer incx,  const float*   y, integer incy);
  
  void   dvpoly  (integer n, const double* x, integer incx, integer m,
		             const double* c, integer incc, 
		                   double* y, integer incy);
  void   svpoly  (integer n, const float*  x, integer incx, integer m,
		             const float*  c, integer incc, 
		                   float*  y, integer incy);
  
  double dpoly   (integer n, double x, const double* xp, const double* yp);
  float  spoly   (integer n, float  x, const float*  xp, const float*  yp);
  void   dpolint (const double* xa, const double* ya, integer n,
		  double x,        double* y,  double* dy);
  void   spolint (const float*  xa, const float*  ya, integer n,
		  float  x,        float*  y, float*  dy);
  
  void   dspline  (integer n, double yp1, double ypn,
		   const double* x, const double* y,  double* y2);
  void   sspline  (integer n, float yp1, float ypn,
		   const float*  x, const float*  y,  float*  y2);
  double dsplint  (integer n, double x,
		   const  double* xa, const double* ya, const double* y2a);
  float  ssplint  (integer n, float  x,
		   const  float*  xa, const float*  ya, const float*  y2a);
  double dsplquad (const double* xa, const double* ya,   const double* y2a,
		   const integer n,  const double xmin,  const double xmax);
  float  ssplquad (const float*  xa, const float*  ya,   const float*  y2a,
		   const integer n,  const float  xmin,  const float  xmax);

  void   dvvtvvtp (integer n, const double* v, integer incv,
          	              const double* w, integer incw,
	                      const double* x, integer incx,
	                      const double* y, integer incy,
	                            double* z, integer incz);
  void   svvtvvtp (integer n, const float*  v, integer incv,
	                      const float*  w, integer incw,
	                      const float*  x, integer incx,
	                      const float*  y, integer incy,
	                            float*  z, integer incz); 

  void   dvvtvvtm (integer n, const double* v, integer incv,
          	              const double* w, integer incw,
	                      const double* x, integer incx,
	                      const double* y, integer incy,
	                            double* z, integer incz);
  void   svvtvvtm (integer n, const float*  v, integer incv,
	                      const float*  w, integer incw,
	                      const float*  x, integer incx,
	                      const float*  y, integer incy,
	                            float*  z, integer incz); 

  void   dsvvttvp (integer n, const double  alpha,
          	              const double* w, integer incw,
	                      const double* x, integer incx,
	                      const double* y, integer incy,
	                            double* z, integer incz);
  void   ssvvttvp (integer n, const float   alpha,
	                      const float*  w, integer incw,
	                      const float*  x, integer incx,
	                      const float*  y, integer incy,
	                            float*  z, integer incz);  
}


class Veclib {
 public:

  // -- ZERO-OFFSET 2D MATRIX ADDRESSES:

  static inline integer row_major (integer i, integer j, integer n)
  { return j + i * n; }
  static inline integer col_major (integer i, integer j, integer n)
  { return i + j * n; }


  // -- MATHEMATICAL PRIMITIVES:

  static void copy (integer n, const double*  x, integer incx,
                                     double*  y, integer incy)
  { dcopy (n, x, incx, y, incy); }
  static void copy (integer n, const integer* x, integer incx,
                                     integer* y, integer incy)
  { icopy (n, x, incx, y, incy); }
  static void copy (integer n, const float*   x, integer incx,
                                     float*   y, integer incy)
  { scopy (n, x, incx, y, incy); } 

  
  static void fill (integer n, double  alpha, double*  x, integer incx)
  { dfill (n, alpha, x, incx); }
  static void fill (integer n, integer alpha, integer* x, integer incx)
  { ifill (n, alpha, x, incx); }
  static void fill (integer n, float   alpha, float*   x, integer incx)
  { sfill (n, alpha, x, incx); }


  static void neg (integer n, double*  x, integer incx)
  { dneg (n, x, incx); }
  static void neg (integer n, integer* x, integer incx)
  { ineg (n, x, incx); }
  static void neg (integer n, float*   x, integer incx)
  { sneg (n, x, incx); }


  static void vneg (integer n,
		    const double*  x, integer incx, double* y, integer incy)
  { dvneg (n, x, incx, y, incy); }
  static void vneg (integer n,
		    const integer* x, integer incx, integer* y, integer incy)
  { ivneg (n, x, incx, y, incy); }
  static void vneg (integer n,
		    const float*   x, integer incx, float*   y, integer incy)
  { svneg (n, x, incx, y, incy); }


  static void vsgn (integer n,
		    const double*  x, integer incx, double* y, integer incy)
  { dvsgn (n, x, incx, y, incy); }
  static void vsgn (integer n,
		    const integer* x, integer incx, integer* y, integer incy)
  { ivsgn (n, x, incx, y, incy); }
  static void vsgn (integer n,
		    const float*   x, integer incx, float*   y, integer incy)
  { svsgn (n, x, incx, y, incy); }

  
  static void sadd (integer n, double  alpha, const double*  x, integer incx, 
		                                   double*  y, integer incy)
  { dsadd (n, alpha, x, incx, y, incy); }
  static void sadd (integer n, integer alpha, const integer* x, integer incx,
		                                    integer* y, integer incy)
  { isadd (n, alpha, x, incx, y, incy); }
  static void sadd (integer n, float  alpha,  const float*   x, integer incx,
		                                    float*   y, integer incy)
  { ssadd (n, alpha, x, incx, y, incy); }


  static void spow (const integer n, const double alpha,
		    const double* x, integer incx, double* y, integer incy)
  { dspow (n, alpha, x, incx, y, incy); }
  static void spow (const integer n, const float alpha,
		    const float *x, integer incx, float *y, integer incy)
  { sspow (n, alpha, x, incx, y, incy); }


  static void vadd (integer n, const double*  x, integer incx,
		               const double*  y, integer incy,
                                     double*  z, integer incz)
  { dvadd (n, x, incx, y, incy, z, incz); }
  static void vadd (integer n, const integer* x, integer incx,
		               const integer* y, integer incy,
                                     integer* z, integer incz)
  { ivadd (n, x, incx, y, incy, z, incz); }
  static void vadd (integer n, const float*   x, integer incx,
		               const float*   y, integer incy,
                                     float*   z, integer incz)
  { svadd (n, x, incx, y, incy, z, incz); }

 
  static void ssub (integer n, double  alpha, const double*  x, integer incx, 
		                                    double*  y, integer incy)
  { dssub (n, alpha, x, incx, y, incy); }
  static void ssub (integer n, integer alpha, const integer* x, integer incx,
		                                    integer* y, integer incy)
  { issub (n, alpha, x, incx, y, incy); }
  static void ssub (integer n, float   alpha, const float*   x, integer incx,
		                                    float*   y, integer incy)
  { sssub (n, alpha, x, incx, y, incy); }

  
  static void vsub (integer n, const double*  x, integer incx,
		               const double*  y, integer incy,
                                     double*  z, integer incz)
  { dvsub (n, x, incx, y, incy, z, incz); }
  static void vsub (integer n, const integer* x, integer incx,
		               const integer* y, integer incy,
                                     integer* z, integer incz)
  { ivsub (n, x, incx, y, incy, z, incz); }
  static void vsub (integer n, const float*   x, integer incx,
		               const float*   y, integer incy,
                                     float*   z, integer incz)
  { svsub (n, x, incx, y, incy, z, incz); }

  
  static void  smul  (integer n, double  alpha, const double*  x, integer incx,
		                                      double*  y, integer incy)
  { dsmul (n, alpha, x, incx, y, incy); }
  static void  smul  (integer n, integer alpha, const integer* x, integer incx,
		                                      integer* y, integer incy)
  { ismul (n, alpha, x, incx, y, incy); }
  static void  smul  (integer n, float  alpha,  const float*   x, integer incx,
		                                      float*   y, integer incy)
  { ssmul (n, alpha, x, incx, y, incy); }

  
  static void vmul (integer n, const double*  x, integer incx,
		               const double*  y, integer incy,
                                     double*  z, integer incz)
  { dvmul (n, x, incx, y, incy, z, incz); } 
  static void vmul (integer n, const integer* x, integer incx,
		               const integer* y, integer incy,
                                     integer* z, integer incz)
  { ivmul (n, x, incx, y, incy, z, incz); }
  static void vmul (integer n, const float*   x, integer incx,
		               const float*   y, integer incy,
                                     float*   z, integer incz)
  { svmul (n, x, incx, y, incy, z, incz); }


  static void sdiv (integer n, double  alpha, const double*  x, integer incx,
		                                    double*  y, integer incy)
  { dsdiv (n, alpha, x, incx, y, incy); }
  static void sdiv (integer n, integer alpha, const integer* x, integer incx,
		                                    integer* y, integer incy)
  { isdiv (n, alpha, x, incx, y, incy); }
  static void  sdiv (integer n, float  alpha, const float*   x, integer incx,
		                                    float*   y, integer incy)
  { ssdiv (n, alpha, x, incx, y, incy); }

  
  static void vrecp (integer n, const double* x, integer incx,
                                      double* y, integer incy)
  { dvrecp (n, x, incx, y, incy); }
  static void vrecp (integer n, const float*  x, integer incx,
                                      float*  y, integer incy)
  { svrecp (n, x, incx, y, incy); }
  

  static void vdiv (integer n, const double*  x, integer incx,
		               const double*  y, integer incy,
		                     double*  z, integer incz)
  { dvdiv (n, x, incx, y, incy, z, incz); }
  static void vdiv (integer n, const float*   x, integer incx,
		               const float*   y, integer incy,
		                     float*   z, integer incz)
  { svdiv (n, x, incx, y, incy, z, incz); }

  
  static void zero (integer n, double*  x, integer incx)
  { dzero (n, x, incx); }
  static void zero (integer n, integer* x, integer incx)
  { izero (n, x, incx); }
  static void zero (integer n, float*   x, integer incx)
  { szero (n, x, incx); }


  // -- OTHER MATHEMATICAL FUNCTIONS:
  
  static void vabs (integer n, const double*  x, integer incx,
                                     double*  y, integer incy)
  { dvabs (n, x, incx, y, incy); }
  static void vabs (integer n, const integer* x, integer incx,
                                     integer* y, integer incy)
  { ivabs (n, x, incx, y, incy); }
  static void vabs (integer n, const float*   x, integer incx,
                                     float*   y, integer incy)
  { svabs (n, x, incx, y, incy); }


  static void vamax (integer n, const double*  x, integer incx,
		                const double*  y, integer incy,
		                      double*  z, integer incz)
  { dvamax (n, x, incx, y, incy, z, incz); }
  static void vamax (integer n, const integer* x, integer incx,
		                const integer* y, integer incy,
                                      integer* z, integer incz)
  { ivamax (n, x, incx, y, incy, z, incz); }
  static void vamax (integer n, const float*   x, integer incx,
		                const float*   y, integer incy,
                                      float*   z, integer incz)
  { svamax (n, x, incx, y, incy, z, incz); }

  
  static void vexp (integer n, const double* x, integer incx,
                                     double* y, integer incy)
  { dvexp (n, x, incx, y, incy); }
  static void vexp (integer n, const float*  x, integer incx,
		                     float*  y, integer incy)
  { svexp (n, x, incx, y, incy); }

  
  static void vlg10 (integer n, const double* x, integer incx,
		                      double* y, integer incy)
  { dvlg10 (n, x, incx, y, incy); }
  static void vlg10 (integer n, const float*  x, integer incx,
		                      float*  y, integer incy)
  { svlg10 (n, x, incx, y, incy); }


  static void vlog (integer n, const double* x, integer incx,
                                     double* y, integer incy)
  { dvlog (n, x, incx, y, incy); }
  static void vlog (integer n, const float*  x, integer incx,
                                     float*  y, integer incy)
  { svlog (n, x, incx, y, incy); }


  static void vatan (integer n, const double* x, integer incx,
	                              double* y, integer incy)
  { dvatan (n, x, incx, y, incy); }
  static void vatan (integer n, const float*  x, integer incx,
                                      float*  y, integer incy)
  { svatan (n, x, incx, y, incy); } 


  static void vatn2 (integer n, const double* x, integer incx,
		                const double* y, integer incy,
                                      double* z, integer incz)
  { dvatn2 (n, x, incx, y, incy, z, incz); }
  static void vatn2 (integer n, const float*  x, integer incx,
		                const float*  y, integer incy,
                                      float*  z, integer incz)
  { svatn2 (n, x, incx, y, incy, z, incz); }


  static void vcos (integer n, const double* x, integer incx,
		                     double* y, integer incy)
  { dvcos (n, x, incx, y, incy); }
  static void vcos (integer n, const float*  x, integer incx,
		                     float*  y, integer incy)
  { svcos (n, x, incx, y, incy); }


  static void vsin (integer n, const double* x, integer incx,
		                     double* y, integer incy)
  { dvsin (n, x, incx, y, incy); }
  static void vsin (integer n, const float*  x, integer incx,
		                     float*  y, integer incy)
  { svsin (n, x, incx, y, incy); }


  static void vsqrt (integer n, const double* x, integer incx,
		                      double* y, integer incy)
  { dvsqrt (n, x, incx, y, incy); }
  static void vsqrt (integer n, const float*  x, integer incx,
		                      float*  y, integer incy)
  { svsqrt (n, x, incx, y, incy); }


  static void ranInit (integer flag)
  { raninit (flag); }
  static double dranu ()         
  { return dranu (); }
  static float  sranu ()         
  { return sranu (); }
  static double normal (double mean, double sdev)
  { return dnormal (mean, sdev); }
  static float  normal (float  mean, float  sdev)
  { return snormal (mean, sdev); }
  static void  vrandom (integer n, double* x, integer incx)
  { dvrandom (n, x, incx); }
  static void  vrandom (integer n, float*  x, integer incx)
  { svrandom (n, x, incx); }
  static void  vnormal (integer n, double mean, double sdev,
			double* x, integer incx)
  { dvnormal (n, mean, sdev, x, incx); }
  static void  vnormal (integer n, float  mean, float  sdev,
			float*  x, integer incx)
  { svnormal (n, mean, sdev, x, incx); }


  static void vhypot (integer n, const double* x, integer incx,
		                 const double* y, integer incy,
                                       double* z, integer incz)
  { dvhypot (n, x, incx, y, incy, z, incz); }
  static void vhypot (integer n, const float*  x, integer incx,
		                 const float*  y, integer incy,
                                       float*  z, integer incz)
  { svhypot (n, x, incx, y, incy, z, incz); }
  static void vmag (integer n, const double* w, integer incw,
		    const double* x, integer incx,
		    const double* y, integer incy,
		          double* z, integer incz)
  { dvmag (n, w, incw, x, incx, y, incy, z, incz); }
  static void vmag (integer n, const float*  w, integer incw,
		    const float*  x, integer incx,
		    const float*  y, integer incy,
		          float*  z, integer incz)
  { svmag (n, w, incw, x, incx, y, incy, z, incz); }

  static void vpow (integer n, const double* x, integer incx,
		    const double* y, integer incy, double* z, integer incz)
  { dvpow (n, x, incx, y, incy, z, incz); } 
  static void vpow (integer n, const float*  x, integer incx,
		    const float*  y, integer incy, float*  z, integer incz)
  { svpow (n, x, incx, y, incy, z, incz); } 


  // -- TRIAD OPERATIONS:
  
  static void svmvt (integer n, double alpha, const double*  x, integer incx,
		                              const double*  y, integer incy,
		                                    double*  z, integer incz)
  { dsvmvt (n, alpha, x, incx, y, incy, z, incz); }
  static void svmvt (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
                                                    float*  z, integer incz)
  { ssvmvt (n, alpha, x, incx, y, incy, z, incz); }

 
  static void svpvt (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
                                                double* z, integer incz)
  { dsvpvt (n, alpha, x, incx, y, incy, z, incz); }
  static void svpvt (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
                                                    float*  z, integer incz)
  { ssvpvt (n, alpha, x, incx, y, incy, z, incz); }


  static void svtsp (integer n, double alpha, double beta,
		     const double* x, integer incx,
		           double* y, integer incy)
  { dsvtsp (n, alpha, beta, x, incx, y, incy); }
  static void svtsp (integer n, float  alpha, float  beta,
		     const float*  x, integer incx,
		           float*  y, integer incy)
  { ssvtsp (n, alpha, beta, x, incx, y, incy); }
 
 
  static void svtvm (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
		                                    double* z, integer incz)
  { dsvtvm (n, alpha, x, incx, y, incy, z, incz); }
  static void svtvm (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
                                                    float*  z, integer incz)
  { ssvtvm (n, alpha, x, incx, y, incy, z, incz); }


  static void svtvp (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
                                                    double* z, integer incz)
  { dsvtvp (n, alpha, x, incx, y, incy, z, incz); }
  static void svtvp (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
                                                    float*  z, integer incz)
  { ssvtvp (n, alpha, x, incx, y, incy, z, incz); }

  
  static void svvmt (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
                                                    double* z, integer incz)
  { dsvvmt (n, alpha, x, incx, y, incy, z, incz); }
  static void svvmt (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
		                                    float*  z, integer incz)
  { ssvvmt (n, alpha, x, incx, y, incy, z, incz); }

 
  static void svvpt (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
		                                    double* z, integer incz)
  { dsvvpt (n, alpha, x, incx, y, incy, z, incz); }
  static void svvpt (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
                                                    float*  z, integer incz)
  { ssvvpt (n, alpha, x, incx, y, incy, z, incz); }


  static void svvtm (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
		                                    double* z, integer incz)
  { dsvvtm (n, alpha, x, incx, y, incy, z, incz); }
  static void svvtm (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
		                                    float*  z, integer incz)
  { ssvvtm (n, alpha, x, incx, y, incy, z, incz); }

 
  static void svvtp (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
		                                    double* z, integer incz)
  { dsvvtp (n, alpha, x, incx, y, incy, z, incz); }
  static void svvtp (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
		                                    float*  z, integer incz)
  { ssvvtp (n, alpha, x, incx, y, incy, z, incz); }


  static void svvtt (integer n, double alpha, const double* x, integer incx,
		                              const double* y, integer incy,
		                                    double* z, integer incz)
  { dsvvtt (n, alpha, x, incx, y, incy, z, incz); }
  static void svvtt (integer n, float  alpha, const float*  x, integer incx,
		                              const float*  y, integer incy,
		                                    float*  z, integer incz)
  { ssvvtt (n, alpha, x, incx, y, incy, z, incz); }


  static void vvmvt (integer n, const double* w, integer incw,
		                const double* x, integer incx,
		                const double* y, integer incy,
		                      double* z, integer incz)
  { dvvmvt (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvmvt (integer n, const float*  w, integer incw,
		                const float*  x, integer incx,
		                const float*  y, integer incy,
		                      float*  z, integer incz)
  { svvmvt (n, w, incw, x, incx, y, incy, z, incz); }

  
  static void vvpvt (integer n, const double* w, integer incw,
		                const double* x, integer incx,
		                const double* y, integer incy,
                                      double* z, integer incz)
  { dvvpvt (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvpvt (integer n, const float*  w, integer incw,
		                const float*  x, integer incx,
		                const float*  y, integer incy,
		                      float*  z, integer incz)
  { svvpvt (n, w, incw, x, incx, y, incy, z, incz); }

 
  static void vvtvp (integer n, const double* w, integer incw,
		                const double* x, integer incx,
		                const double* y, integer incy,
		                      double* z, integer incz)
  { dvvtvp (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvp (integer n, const float*  w, integer incw,
		                const float*  x, integer incx,
		                const float*  y, integer incy,
		                      float*  z, integer incz)
  { svvtvp (n, w, incw, x, incx, y, incy, z, incz); }

 
  static void vvtvm (integer n, const double* w, integer incw,
		                const double* x, integer incx,
		                const double* y, integer incy,
		                      double* z, integer incz)
  { dvvtvm (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvm (integer n, const float*  w, integer incw,
		                const float*  x, integer incx,
		                const float*  y, integer incy,
		                      float*  z, integer incz)
  { svvtvm (n, w, incw, x, incx, y, incy, z, incz); }

  
  static void vvvtm (integer n, const double* w, integer incw,
		                const double* x, integer incx,
		                const double* y, integer incy,
		                      double* z, integer incz)
  { dvvvtm (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvvtm (integer n, const float*  w, integer incw,
		                const float*  x, integer incx,
		                const float*  y, integer incy,
		                      float*  z, integer incz)
  { svvvtm (n, w, incw, x, incx, y, incy, z, incz); }

 
  // -- RELATIONAL PRIMITIVE OPERATIONS:
  
  static void seq (integer n, integer alpha, const integer* x, integer incx,
                                                   integer* y, integer incy)
  { iseq (n, alpha, x, incx, y, incy); }


  static void sge (integer n, double  alpha, const double* x, integer incx,
		                                  integer* y, integer incy)
  { dsge (n, alpha, x, incx, y, incy); }
  static void sge (integer n, integer alpha, const integer* x, integer incx,
                                                   integer* y, integer incy)
  { isge (n, alpha, x, incx, y, incy); }
  static void sge (integer n, float   alpha, const float*   x, integer incx,
		                                   integer* y, integer incy)
  { ssge (n, alpha, x, incx, y, incy); }


  static void sle (integer n, double  alpha,  const double*  x, integer incx,
		                                    integer* y, integer incy)
  { dsle (n, alpha, x, incx, y, incy); }
  static void sle (integer n, integer alpha, const integer* x, integer incx,
                                                   integer* y, integer incy)
  { isle (n, alpha, x, incx, y, incy); }
  static void sle (integer n, float   alpha, const float*   x, integer incx,
		                                   integer* y, integer incy)
  { ssle (n, alpha, x, incx, y, incy); }


  static void slt (integer n, double  alpha, const double*  x, integer incx,
		                                   integer* y, integer incy)
  { dslt (n, alpha, x, incx, y, incy); }
  static void slt (integer n, integer alpha, const integer* x, integer incx,
		                                   integer* y, integer incy)
  { islt (n, alpha, x, incx, y, incy); }
  static void slt (integer n, float   alpha, const float*   x, integer incx,
		                                   integer* y, integer incy)
  { sslt (n, alpha, x, incx, y, incy); }


  static void sne (integer n, double  alpha, const double*  x, integer incx,
		                                   integer* y, integer incy)
  { dsne (n, alpha, x, incx, y, incy); }
  static void sne (integer n, integer alpha, const integer* x, integer incx,
		                                   integer* y, integer incy)
  { isne (n, alpha, x, incx, y, incy); }
  static void sne (integer n, float   alpha, const float*   x, integer incx,
		                                   integer* y, integer incy)
  { ssne (n, alpha, x, incx, y, incy); } 
 

  // -- REDUCTION FUNCTIONS:
  
  static double  sum (integer n, const double*  x, integer incx)
  { return dsum (n, x, incx); }
  static integer sum (integer n, const integer* x, integer incx)
  { return isum (n, x, incx); }
  static float   sum (integer n, const float*   x, integer incx)
  { return ssum (n, x, incx); }


  static integer imax (integer n, const double*  x, integer incx)
  { return idmax (n, x, incx); }
  static integer imax (integer n, const integer* x, integer incx)
  { return iimax (n, x, incx); }
  static integer imax (integer n, const float*   x, integer incx)
  { return ismax (n, x, incx); }


  static integer imin (integer n, const double*  x, integer incx)
  { return idmin (n, x, incx); }
  static integer imin (integer n, const integer* x, integer incx)
  { return iimin (n, x, incx); }
  static integer imin (integer n, const float*   x, integer incx)
  { return ismin (n, x, incx); }


  static integer count (integer n, const integer *x, integer incx)
  { return icount (n, x, incx); }
  static integer first (integer n, const integer *x, integer incx)
  { return ifirst (n, x, incx); }
  static integer any   (integer n, const integer *x, integer incx)
  { return lany  (n, x, incx); }


  static integer same (integer n, const integer* x, integer incx,
                                  const integer* y, integer incy)
  { return lisame (n, x, incx, y, incy); }
  static integer same (integer n, const double*  x, integer incx,
                                  const double*  y, integer incy)
  { return ldsame (n, x, incx, y, incy); }
  static integer same (integer n, const float*   x, integer incx,
                                  const float*   y, integer incy)
  { return lssame (n, x, incx, y, incy); }

    
  // -- CONVERSION PRIMITIVES:
  
  static void vdble (integer n, const float*  x, integer incx,
		                      double* y, integer incy)
  { vdble (n, x, incx, y, incy); }
  static void vsngl (integer n, const double* x, integer incx,
		                      float*  y, integer incy)
  { vsngl (n, x, incx, y, incy); }


  static void vfloa (integer n, const integer* x, integer incx,
                                      double*  y, integer incy)
  { dvfloa (n, x, incx, y, incy); }
  static void vfloa (integer n, const integer* x, integer incx,
                                      float*   y, integer incy)
  { svfloa (n, x, incx, y, incy); }


  static integer testFormat ()
  { return iformat (); }
  static void describeFormat (char* s)
  { format (s); }
  static void brev (integer n, const double*  x, integer incx,
		                     double*  y, integer incy)
  { dbrev (n, x, incx, y, incy); }
  static void brev (integer n, const integer* x, integer incx,
                                     integer* y, integer incy)
  { ibrev (n, x, incx, y, incy); }
  static void brev (integer n, const float*   x, integer incx,
                                     float*   y, integer incy)
  { sbrev (n, x, incx, y, incy); }


   // -- MISCELLANEOUS FUNCTIONS:

  static double clock ()
  { return dclock () ; }
  

  static void scatr (integer n, const double*  x, const integer *y, double*  z)
  { dscatr (n, x, y, z); }
  static void scatr (integer n, const integer* x, const integer *y, integer* z)
  { iscatr (n, x, y, z); }
  static void scatr (integer n, const float*   x, const integer *y, float*   z)
  { sscatr (n, x, y, z); }


  static void gathr (integer n, const double*  x, const integer *y, double*  z)
  { dgathr (n, x, y, z); }
  static void gathr (integer n, const integer* x, const integer *y, integer* z)
  { igathr (n, x, y, z); }
  static void gathr (integer n, const float*   x, const integer *y, float*   z)
  { sgathr (n, x, y, z); }


  static void gathr_scatr (integer n, const double*  w, const integer* x,
			              const integer* y,       double  *z)
  { dgathr_scatr (n, w, x, y, z); }
  static void gathr_scatr (integer n, const integer* w, const integer* x,
			              const integer* y,       integer* z)
  { igathr_scatr (n, w, x, y, z); }
  static void gathr_scatr (integer n, const float*   w, const integer* x,
			              const integer* y,       float   *z)
  { sgathr_scatr (n, w, x, y, z); }

  
  static void scatr_sum (integer n,
			 const double*  x, const integer *y, double*  z)
  { dscatr_sum (n, x, y, z); }
  static void scatr_sum (integer n,
			 const integer* x, const integer *y, integer* z)
  { iscatr_sum (n, x, y, z); }
  static void scatr_sum (integer n,
			 const float*   x, const integer *y, float*   z)
  { sscatr_sum (n, x, y, z); }


  static void gathr_sum (integer n,
			 const double*  x, const integer *y, double*  z)
  { dgathr_sum (n, x, y, z); }
  static void gathr_sum (integer n,
			 const integer* x, const integer *y, integer* z)
  { igathr_sum (n, x, y, z); }
  static void gathr_sum (integer n,
			 const float*   x, const integer *y, float*   z)
  { sgathr_sum (n, x, y, z); }


  static void gathr_scatr_sum (integer n, const double*  w, const integer* x,
			                  const integer* y,       double*  z)
  { dgathr_scatr_sum (n, w, x, y, z); }
  static void gathr_scatr_sum (integer n, const integer* w, const integer* x,
			                  const integer* y,       integer* z)
  { igathr_scatr_sum (n, w, x, y, z); }
  static void gathr_scatr_sum (integer n, const float*   w, const integer* x,
			                  const integer* y,       float   *z)
  { sgathr_scatr_sum (n, w, x, y, z); }


  static void ramp (integer n, double  alpha, double  beta,
		    double*  x, integer incx)
  { dramp (n, alpha, beta, x, incx); }
  static void ramp (integer n, integer alpha, integer beta,
		    integer* x, integer incx)
  { iramp (n, alpha, beta, x, incx); }
  static void ramp (integer n, float   alpha, float   beta,
		    float*   x, integer incx)
  { sramp (n, alpha, beta, x, incx); }

  static void clip (integer n, const double alpha,  const double    beta,
		    const double* x, integer incx,
 		          double* y, integer incy)
  { dclip (n, alpha, beta, x, incx, y, incy); }
  static void clip (integer n, const  integer alpha, const integer  beta,
		    const integer* x, integer incx,
		          integer* y, integer incy)
  { iclip (n, alpha, beta, x, incx, y, incy); }
  static void clip (integer n, const  float alpha,   const float    beta,
		    const float* x,  integer incx,
		          float* y, integer incy)
  { sclip (n, alpha, beta, x, incx, y, incy); }

  static void clipup (integer n, const double alpha,
		      const double* x, integer incx,
		            double* y, integer incy)
  { dclipup (n, alpha, x, incx, y, incy); }
  static void clipup (integer n, const  integer alpha,
		      const integer* x, integer incx,
		            integer* y, integer incy)
  { iclipup (n, alpha, x, incx, y, incy); }
  static void clipup (integer n, const  float alpha,
		      const float* x, integer incx,
		            float* y, integer incy)
  { sclipup (n, alpha, x, incx, y, incy); }

  static void clipdn (integer n, const double alpha,
		    const double* x, integer incx,
		          double* y, integer incy)
  { dclipdn (n, alpha, x, incx, y, incy); }
  static void clipdn (integer n, const  integer alpha,
		    const integer* x, integer incx,
		          integer* y, integer incy)
  { iclipdn (n, alpha, x, incx, y, incy); }
  static void clipdn (integer n, const  float alpha,
		    const float* x,  integer incx,
		          float* y, integer incy)
  { sclipdn (n, alpha, x, incx, y, incy); }

  static void iclip (integer n, const double alpha,  const double    beta,
		      const double* x, integer incx,
		            double* y, integer incy)
  { diclip (n, alpha, beta, x, incx, y, incy); }
  static void iclip (integer n, const  integer alpha, const integer  beta,
		     const integer* x, integer incx,
		           integer* y, integer incy)
  { iiclip (n, alpha, beta, x, incx, y, incy); }
  static void iclip (integer n, const float alpha,   const float     beta,
		     const float* x,  integer incx,
		           float* y,  integer incy)
  { siclip (n, alpha, beta, x, incx, y, incy); }

  static void cndst (integer n, const double*  x, integer incx,
		                const integer* y, integer incy,
		                      double*  z, integer incz)
  { dcndst (n, x, incx, y, incy, z, incz); }
  static void cndst (integer n, const integer* x, integer incx,
		                const integer* y, integer incy,
		                      integer* z, integer incz)
  { icndst (n, x, incx, y, incy, z, incz); }
  static void cndst (integer n, const float*   x, integer incx,
		                const integer* y, integer incy,
		                      float*   z, integer incz)
  { scndst (n, x, incx, y, incy, z, incz); }
 

  static void mask (integer n, const double*  w, integer incw,
		               const double*  x, integer incx,
		               const integer* y, integer incy,
		                     double*  z, integer incz)
  { dmask (n, w, incw, x, incx, y, incy, z, incz); }
  static void mask (integer n, const integer* w, integer incw,
		               const integer* x, integer incx,
		               const integer* y, integer incy,
		                     integer *z, integer incz)
  { imask (n, w, incw, x, incx, y, incy, z, incz); }
  static void mask (integer n, const float*   w, integer incw,
		           const float*   x, integer incx,
		           const integer* y, integer incy, 
		                 float*   z, integer incz)
  { smask (n, w, incw, x, incx, y, incy, z, incz); }

  
  static void vpoly (integer n, const double* x, integer incx, integer m,
		                const double* c, integer incc, 
		                      double* y, integer incy)
  { dvpoly (n, x, incx, m, c, incc, y, incy); }
  static void vpoly (integer n, const float*  x, integer incx, integer m,
		                const float*  c, integer incc, 
		                      float*  y, integer incy)
  { svpoly (n, x, incx, m, c, incc, y, incy); }  


  static double poly (integer n, double x, const double* xp, const double* yp)
  { return dpoly (n, x, xp, yp); }
  static float  poly (integer n, float  x, const float*  xp, const float*  yp)
  { return spoly (n, x, xp, yp); }


  static void polint (const double* xa, const double* ya, integer n,
		            double  x,        double* y,  double* dy)
  { dpolint (xa, ya, n, x, y, dy); }
  static void polint (const float*  xa, const float*  ya, integer n,
		            float   x,        float*  y, float*  dy)
  { spolint (xa, ya, n, x, y, dy); }
  

  static void spline (integer n, double yp1, double ypn,
		      const double* x,  const double* y, double* y2)
  { dspline (n, yp1, ypn, x,  y, y2); }
  static void spline (integer n, float yp1, float ypn,
		      const float*  x, const float*  y,  float*  y2)
  { sspline (n, yp1, ypn, x,  y, y2); }


  static double splint (integer n, double x,
			const  double* xa, const double* ya, const double* y2a)
  { return dsplint (n, x, xa, ya, y2a); }
  static float  splint (integer n, float  x,
			const  float*  xa, const float*  ya, const float*  y2a)
  { return ssplint (n, x, xa, ya, y2a); }


  static double splquad (const double* xa, const double* ya, const double* y2a,
			 const integer n,  const double xmn, const double  xmx)
  { return dsplquad (xa, ya, y2a, n, xmn, xmx); }
  static float  splquad (const float*  xa, const float*  ya, const float*  y2a,
			 const integer n,  const float  xmn, const float   xmx)
  { return ssplquad (xa, ya, y2a, n, xmn, xmx); }


  static void vvtvvtp (integer n, const double* v, integer incv,
          	                  const double* w, integer incw,
	                          const double* x, integer incx,
	                          const double* y, integer incy,
	                                double* z, integer incz)
  { dvvtvvtp (n, v, incv, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvvtp (integer n, const float*  v, integer incv,
	                          const float*  w, integer incw,
	                          const float*  x, integer incx,
	                          const float*  y, integer incy,
		                        float*  z, integer incz)
  { svvtvvtp (n, v, incv, w, incw, x, incx, y, incy, z, incz); }

  static void vvtvvtm (integer n, const double* v, integer incv,
          	                  const double* w, integer incw,
	                          const double* x, integer incx,
	                          const double* y, integer incy,
	                                double* z, integer incz)
  { dvvtvvtm (n, v, incv, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvvtm (integer n, const float*  v, integer incv,
	                          const float*  w, integer incw,
	                          const float*  x, integer incx,
	                          const float*  y, integer incy,
		                        float*  z, integer incz)
  { svvtvvtm (n, v, incv, w, incw, x, incx, y, incy, z, incz); }


  static void svvttvp (integer n, const double  alpha,
          	                  const double* w, integer incw,
	                          const double* x, integer incx,
	                          const double* y, integer incy,
	                                double* z, integer incz)
  { dsvvttvp (n, alpha, w, incw, x, incx, y, incy, z, incz); }
  static void svvttvp (integer n, const float   alpha,
	                          const float*  w, integer incw,
	                          const float*  x, integer incx,
	                          const float*  y, integer incy,
		                  float*  z, integer incz)
  { ssvvttvp (n, alpha, w, incw, x, incx, y, incy, z, incz); }

};

#endif
