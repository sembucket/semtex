#ifndef BLAS_H
#define BLAS_H

/* BLAS: Basic Linear Algebra Subroutines
 *
 * $Revision$
 *
 * Author: R. D. Henderson
 *
 * This file defines a set of macros to translate call to the BLAS from C to
 * FORTRAN linkage.  It uses a set of registers provided by veclib for the
 * translation from call-by-value to call-by-reference.
 * ------------------------------------------------------------------------- */

/* Level 1 */
 
void   icopy  (int n, int      *x, int incx, int      *y, int incy);
void   dcopy  (int n, double   *x, int incx, double   *y, int incy);
void   zcopy  (int n, zcomplex *x, int incx, zcomplex *y, int incy);
void   dscal  (int n, double alpha, double *x, int incx);
void   daxpy  (int n, double alpha, double *x, int incx, double *y, int incy);
void   dswap  (int n, double *x, int incx, double *y, int incy);

float  sdot   (int n, float  *x, int incx, float  *y, int incy);
double ddot   (int n, double *x, int incx, double *y, int incy);
double dasum  (int n, double *x, int incx);
double dnrm2  (int n, double *x, int incx);
int    idamax (int n, double *x, int incx);  /* range = [0,n-1] */

/* Level 2 */

void   dgemv  (char  t, int m, int n, double alpha, double *a, int lda, 
	          double *x, int incx, double beta, double *y, int incy); 

/* Level 3 */ 

void   dgemm  (char ta, char tb, int m, int n, int k, 
	       double alpha, double *a, int lda, double *b, int ldb,
	       double beta , double *c, int ldc);

/* Translations */

double ddot_  (int *n, double *x, int *incx, double *y, int *incy);
double dasum_ (int *n, double *x, int *incx);
double dnrm2_ (int *n, double *x, int *incx);

void dcopy_ (int *n, double *x, int *incx, double *y, int *incy);
void dswap_ (int *n, double *x, int *incx, double *y, int *incy);
void dscal_ (int *n, double *d, double *x, int *incx);
void daxpy_ (int *n, double *d, double *x, int *incx, double *y, int *incy);
int idamax_ (int *n, double *x, int *incx);

#define ddot(n,x,incx,y,incy) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,_vlib_ireg[2]=incy,\
   ddot_(_vlib_ireg,x,_vlib_ireg+1,y,_vlib_ireg+2))
#define dasum(n,x,incx)\
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,dasum_(_vlib_ireg,x,_vlib_ireg+1))
#define dnrm2(n,x,incx) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,dnrm2_(_vlib_ireg,x,_vlib_ireg+1))
#define dcopy(n,x,incx,y,incy) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,_vlib_ireg[2]=incy,\
   dcopy_(_vlib_ireg,x,_vlib_ireg+1,y,_vlib_ireg+2))
#define dswap(n,x,incx,y,incy) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,_vlib_ireg[2]=incy,\
   dswap_(_vlib_ireg,x,_vlib_ireg+1,y,_vlib_ireg+2))
#define dscal(n,alpha,x,incx) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,_vlib_dreg[0]=alpha,\
   dscal_(_vlib_ireg,_vlib_dreg,x,_vlib_ireg+1))
#define daxpy(n,alpha,x,incx,y,incy) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,_vlib_ireg[2]=incy,_vlib_dreg[0]=alpha,\
   daxpy_(_vlib_ireg,_vlib_dreg,x,_vlib_ireg+1,y,_vlib_ireg+2))
#define idamax(n,x,incx) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=incx,idamax_(_vlib_ireg,x,_vlib_ireg+1)-1)

void    dgemv_(char *t, int *m, int *n, double *alpha, double *a, int *lda,
	       double *x, int *incx, double *beta, double *y, int *incy);
void    sgemv_(char *t, int *m, int *n, float  *alpha, float  *a, int *lda,
	       float  *x, int *incx, float  *beta, float  *y, int *incy);
#define dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) \
  (_vlib_creg[0]=trans,_vlib_ireg[0]=m,_vlib_ireg[1]=n,_vlib_ireg[2]=lda,\
   _vlib_ireg[3]=incx,_vlib_ireg[4]=incy,_vlib_dreg[0]=alpha,\
   _vlib_dreg[1]=beta,\
   dgemv_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,_vlib_dreg,a,_vlib_ireg+2,x,\
	  _vlib_ireg+3,_vlib_dreg+1,y,_vlib_ireg+4))
#define sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) \
  (_vlib_creg[0]=trans,_vlib_ireg[0]=m,_vlib_ireg[1]=n,_vlib_ireg[2]=lda,\
   _vlib_ireg[3]=incx,_vlib_ireg[4]=incy,_vlib_sreg[0]=alpha,\
   _vlib_sreg[1]=beta,\
   sgemv_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,_vlib_sreg,a,_vlib_ireg+2,x,\
	  _vlib_ireg+3,_vlib_sreg+1,y,_vlib_ireg+4))

void    dgemm_(char *ta, char *tb, int *m, int *n, int *k, double *alpha,
	       double *a, int *lda, double *b, int *ldb, double *beta,
	       double *c, int *ldc);
#define dgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)\
  (_vlib_creg[0]=ta,_vlib_creg[1]=tb,_vlib_ireg[0]=m,_vlib_ireg[1]=n,\
   _vlib_ireg[2]=k,_vlib_ireg[3]=lda,_vlib_ireg[4]=ldb,_vlib_ireg[5]=ldc,\
   _vlib_dreg[0]=alpha,_vlib_dreg[1]=beta,\
   dgemm_(_vlib_creg,_vlib_creg+1,_vlib_ireg,_vlib_ireg+1,_vlib_ireg+2,\
	  _vlib_dreg,a,_vlib_ireg+3,b,_vlib_ireg+4,_vlib_dreg+1,c,\
	  _vlib_ireg+5))


#if defined(CRAY)
#
# include <fortran.h>
#
double SDOT  (int *n, double *x, int *incx, double *y, int *incy);
double SASUM (int *n, double *x, int *incx);
double SNRM2 (int *n, double *x, int *incx);

void SCOPY (int *n, double *x, int *incx, double *y, int *incy);
void SSWAP (int *n, double *x, int *incx, double *y, int *incy);
void SSCAL (int *n, double *d, double *x, int *incx);
void SAXPY (int *n, double *d, double *x, int *incx, double *y, int *incy);
int ISAMAX (int *n, double *x, int *incx);
#
# define  ddot_   SDOT
# define  dasum_  SASUM
# define  dnrm2_  SNRM2
# define  dswap_  SSWAP
# define  dscal_  SSCAL
# define  daxpy_  SAXPY
# define  idamax_ ISAMAX
# define  dcopy_  SCOPY
#
# define dgemv_(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)\
  SGEMV (_cptofcd(trans,1),m,n,alpha,a,lda,x,incx,beta,y,incy)
#
# define dgemm_(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)\
  SGEMM (_cptofcd(ta,1),_cptofcd(tb,1),m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#endif

#endif
