#ifndef  LAPACK_H
#define  LAPACK_H

/* LAPACK: translation macros
 *
 * $Revision$
 *
 * Author: R. D. Henderson
 *
 * The following are just wrappers for calls to LAPACK routines written in
 * Fortran.  You have to mess around with these to get the linkage write,
 * depending on what system you compile for.  Registers for storing call-
 * by-reference values are provided by veclib.
 *
 * IMPORTANT NOTE: Arrays in C are stored in row-major order, while arrays
 * in Fortran are stored in column-major order.  You are responsible for
 * passing a sensible (i.e. Fortran-formatted) array into these routines.
 * There is no automatic transpose of array elements.
 * ------------------------------------------------------------------------ */

/* Macros */

#define dpttrf(n,d,e,info)\
  (_vlib_ireg[0]=n,dpttrf_(_vlib_ireg,d,e,&info))
#define dpttrs(n,nrhs,d,e,b,ldb,info)\
  (_vlib_ireg[0]=n,_vlib_ireg[1]=nrhs,_vlib_ireg[2]=ldb,\
   dpttrs_(_vlib_ireg,_vlib_ireg+1,d,e,b,_vlib_ireg+2,&info))

#define dpbtrf(uplo,n,kd,ab,ldab,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=kd,_vlib_ireg[2]=ldab,\
   dpbtrf_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ab,_vlib_ireg+2,&info))
#define spbtrf(uplo,n,kd,ab,ldab,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=kd,_vlib_ireg[2]=ldab,\
   spbtrf_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ab,_vlib_ireg+2,&info))

#define dpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=kd,_vlib_ireg[2]=nrhs,\
   _vlib_ireg[3]=ldab,_vlib_ireg[4]=ldb,\
   dpbtrs_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,_vlib_ireg+2,ab,_vlib_ireg+3,\
	   b,_vlib_ireg+4,&info))
#define spbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=kd,_vlib_ireg[2]=nrhs,\
   _vlib_ireg[3]=ldab,_vlib_ireg[4]=ldb,\
   spbtrs_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,_vlib_ireg+2,ab,_vlib_ireg+3,\
	   b,_vlib_ireg+4,&info))

#define dpbcon(uplo,n,kd,ab,ldab,anorm,rcond,work,iwork,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=kd,_vlib_ireg[2]=ldab,\
   _vlib_dreg[0]=anorm,\
   dpbcon_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ab,_vlib_ireg+2,_vlib_dreg,\
	   &rcond,work,iwork,&info))
#define spbcon(uplo,n,kd,ab,ldab,anorm,rcond,work,iwork,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=kd,_vlib_ireg[2]=ldab,\
   _vlib_sreg[0]=anorm,\
   spbcon_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ab,_vlib_ireg+2,_vlib_sreg,\
	   &rcond,work,iwork,&info))

#define dpptrf(uplo,n,ap,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,\
   dpptrf_(_vlib_creg,_vlib_ireg,ap,&info))
#define spptrf(uplo,n,ap,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,\
   spptrf_(_vlib_creg,_vlib_ireg,ap,&info))

#define dpptrs(uplo,n,nrhs,ap,b,ldb,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=nrhs,_vlib_ireg[2]=ldb,\
   dpptrs_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ap,b,_vlib_ireg+2,&info))
#define spptrs(uplo,n,nrhs,ap,b,ldb,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=nrhs,_vlib_ireg[2]=ldb,\
   spptrs_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ap,b,_vlib_ireg+2,&info))

#define dppcon(uplo,n,ap,anorm,rcond,work,iwork,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_dreg[0]=anorm,\
   dppcon_(_vlib_creg,_vlib_ireg,ap,_vlib_dreg,&rcond,work,iwork,&info))
#define sppcon(uplo,n,ap,anorm,rcond,work,iwork,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_sreg[0]=anorm,\
   sppcon_(_vlib_creg,_vlib_ireg,ap,_vlib_sreg,&rcond,work,iwork,&info))

#define dsptrf(uplo,n,ap,ipiv,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,\
   dsptrf_(_vlib_creg,_vlib_ireg,ap,ipiv,&info))
#define ssptrf(uplo,n,ap,ipiv,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,\
   ssptrf_(_vlib_creg,_vlib_ireg,ap,ipiv,&info))

#define dsptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=nrhs,_vlib_ireg[2]=ldb,\
   dsptrs_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ap,ipiv,b,_vlib_ireg+2,&info))
#define ssptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=nrhs,_vlib_ireg[2]=ldb,\
   ssptrs_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,ap,ipiv,b,_vlib_ireg+2,&info))

#define dspcon(uplo,n,ap,ipiv,anorm,rcond,work,iwork,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_dreg[0]=anorm,\
   dspcon_(_vlib_creg,_vlib_ireg,ap,ipiv,_vlib_dreg,&rcond,work,iwork,&info))
#define sspcon(uplo,n,ap,ipiv,anorm,rcond,work,iwork,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_sreg[0]=anorm,\
   sspcon_(_vlib_creg,_vlib_ireg,ap,ipiv,_vlib_sreg,&rcond,work,iwork,&info))

#define dgbtrf(m,n,kl,ku,ab,ldab,ipiv,info)\
  (_vlib_ireg[0]=m,_vlib_ireg[1]=n,_vlib_ireg[2]=kl,_vlib_ireg[3]=ku,\
   _vlib_ireg[4]=ldab,\
   dgbtrf_(_vlib_ireg,_vlib_ireg+1,_vlib_ireg+2,_vlib_ireg+3,ab,_vlib_ireg+4,\
	   ipiv,&info))
#define dgbtrs(uplo,n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)\
  (_vlib_creg[0]=uplo,_vlib_ireg[0]=n,_vlib_ireg[1]=kl,_vlib_ireg[2]=ku,\
   _vlib_ireg[3]=nrhs,_vlib_ireg[4]=ldab,_vlib_ireg[5]=ldb,\
   dgbtrs_(_vlib_creg,_vlib_ireg,_vlib_ireg+1,_vlib_ireg+2,_vlib_ireg+3,ab,\
	   _vlib_ireg+4,ipiv,b,_vlib_ireg+5,&info))

/* Simple Drivers for Linear Equations */

#define dgesv(n,nrhs,a,lda,ipiv,b,ldb,info) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=nrhs,_vlib_ireg[2]=lda,_vlib_ireg[3]=ldb,\
   dgesv_(_vlib_ireg,_vlib_ireg+1,a,_vlib_ireg+2,ipiv,b,_vlib_ireg+3,&info))
#define sgesv(n,nrhs,a,lda,ipiv,b,ldb,info) \
  (_vlib_ireg[0]=n,_vlib_ireg[1]=nrhs,_vlib_ireg[2]=lda,_vlib_ireg[3]=ldb,\
   sgesv_(_vlib_ireg,_vlib_ireg+1,a,_vlib_ireg+2,ipiv,b,_vlib_ireg+3,&info))

/* Expert Drivers for Linear Equations */

#define dpbsvx(fact,uplo,n,kd,nrhs,ab,ldab,afb,ldafb,equed,s,b,ldb,x,ldx,\
	       rcond,ferr,berr,work,iwork,info)\
  (_vlib_creg[0]=fact,_vlib_creg[1]=uplo ,_vlib_creg[2]=equed,\
   _vlib_ireg[0]=n   ,_vlib_ireg[1]=kd   ,_vlib_ireg[2]=nrhs,\
   _vlib_ireg[3]=ldab,_vlib_ireg[4]=ldafb,_vlib_ireg[5]=ldb,\
   _vlib_ireg[6]=ldx,\
   dpbsvx_(_vlib_creg,_vlib_creg+1,_vlib_ireg,_vlib_ireg+1,_vlib_ireg+2,ab,\
	   _vlib_ireg+3,afb,_vlib_ireg+4,_vlib_creg+2,s,b,_vlib_ireg+5,x,\
	   _vlib_ireg+6,&rcond,ferr,berr,work,iwork,&info))

#define dppsvx(fact,uplo,n,nrhs,ap,afp,equed,s,b,ldb,x,ldx,\
	       rcond,ferr,berr,work,iwork,info)\
  (_vlib_creg[0]=fact,_vlib_creg[1]=uplo ,_vlib_creg[2]=equed,\
   _vlib_ireg[0]=n   ,_vlib_ireg[1]=nrhs ,_vlib_ireg[2]=ldb,\
   _vlib_ireg[3]=ldx,\
   dppsvx_(_vlib_creg,_vlib_creg+1,_vlib_ireg,_vlib_ireg+1,ap,afp,\
	   _vlib_creg+2,s,b,_vlib_ireg+2,x,_vlib_ireg+3,\
	   &rcond,ferr,berr,work,iwork,&info))

/* Simple Drivers for Standard Eigenvalue problems */

#define dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info) \
  ( _vlib_creg[0]=jobvl, _vlib_creg[1]=jobvr, _vlib_ireg[0]=n, \
    _vlib_ireg[1]=lda,   _vlib_ireg[2]=ldvl,  _vlib_ireg[3]=ldvr, \
    _vlib_ireg[4]=lwork, \
      dgeev_(_vlib_creg, _vlib_creg+1, _vlib_ireg, a, _vlib_ireg+1, \
	     wr, wi, vl, _vlib_ireg+2, vr, _vlib_ireg+3, work, \
	     _vlib_ireg+4, &info))

#define dgegv(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,\
              vl,ldvl,vr,ldvr,work,lwork,info) \
  ( _vlib_creg[0]=jobvl, _vlib_creg[1]=jobvr, _vlib_ireg[0]=n, \
    _vlib_ireg[1]=lda,   _vlib_ireg[2]=ldb,   _vlib_ireg[3]=ldvl,\
    _vlib_ireg[4]=ldvr,  _vlib_ireg[5]=lwork, \
      dgegv_(_vlib_creg, _vlib_creg+1, _vlib_ireg, a, _vlib_ireg+1, b,\
             _vlib_ireg+2, alphar, alphai, beta, vl, _vlib_ireg+3, vr,\
             _vlib_ireg+4, work, _vlib_ireg+5, &info))
		 
#define dsbev(jobz,uplo,n,m,AB,ldab,w,z,ldz,work,info) \
  ( _vlib_creg[0]=jobz,  _vlib_creg[1]=uplo, _vlib_ireg[0]=n,\
    _vlib_ireg[1]=m   ,  _vlib_ireg[2]=ldab, _vlib_ireg[3]=ldz,\
     dsbev_(_vlib_creg,_vlib_creg+1,_vlib_ireg,_vlib_ireg+1,AB,_vlib_ireg+2,\
	   w, z, _vlib_ireg+3, work, &info) )

#define dspev(jobz,uplo,n,AB,w,z,ldz,work,info) \
  ( _vlib_creg[0]=jobz,  _vlib_creg[1]=uplo, _vlib_ireg[0]=n,\
    _vlib_ireg[1]=ldz,\
     dspev_(_vlib_creg,_vlib_creg+1,_vlib_ireg, AB, w, z,_vlib_ireg+1,\
	   work, &info) )

#if defined(CRAY)
#
# include<fortran.h>
#
# define  spbtrs_ SPBTRS
# define  spptrs_ SPPTRS
# define  dgesv_  SGESV
#
# define  dpbtrf_(uplo,n,kd,ab,ldab,info)\
  SPBTRF(_cptofcd(uplo,1),n,kd,ab,ldab,info)
# define  dpbtrs_(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)\
  SPBTRS(_cptofcd(uplo,1),n,kd,nrhs,ab,ldab,b,ldb,info)
# define  dpbcon_(uplo,n,kd,ab,ldab,anorm,rcond,work,iwork,info)\
  SPBCON(_cptofcd(uplo,1),n,kd,ab,ldab,anorm,rcond,work,iwork,info)
# define  dpptrf_(uplo,n,ap,info)\
  SPPTRF(_cptofcd(uplo,1),n,ap,info)
# define  dpptrs_(uplo,n,nrhs,ap,b,ldb,info)\
  SPPTRS(_cptofcd(uplo,1),n,nrhs,ap,b,ldb,info)
# define  dppcon_(uplo,n,ap,anorm,rcond,work,iwork,info)\
  SPPCON(_cptofcd(uplo,1),n,ap,anorm,rcond,work,iwork,info)
#endif

/* End of LAPACK routines */

#endif

