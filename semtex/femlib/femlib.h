#ifndef FEMLIB_H
#define FEMLIB_H
/*****************************************************************************
 *         FUNCTION PROTOTYPES FOR ROUTINES IN LIBRARY LIBFEM.A
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

/* -- Routines from initial.c: */

void    yy_initialize (void);
double  yy_interpret  (const char*);

void    yy_vec_init   (const char*, const char*);
void    yy_vec_interp (const integer, ...);

void    yy_help       (void);
void    yy_show       (void);
integer yy_dump       (char*, const integer);

/* -- Routines from polyops.c: */

void   dermat_g (const integer, const double*, const integer,
		 const double*, double**, double**);
void   intmat_g (const integer, const double*, const integer,
		 const double*, double**, double**);
void   dermat_k (const integer, const double*, double**, double**);

void   jacg     (const integer, const double, const double, double*);
void   jacgr    (const integer, const double, const double, double*);
void   jacgl    (const integer, const double, const double, double*);

void   zwgl     (double*, double*, const integer);
void   zwgrl    (double*, double*, const integer);
void   zwgll    (double*, double*, const integer);

double pnleg    (const double,  const integer);
double pndleg   (const double,  const integer);
double pnd2leg  (const double,  const integer);
void   legtr2d  (const integer, const double*, double*);
void   legtr2i  (const integer, const double*, double*);

void   dgll     (const integer, const double*, double**, double**);

void   uniknot  (const integer, double*);

integer    quadComplete (const integer, const integer);

/* -- Routines from operators.c: */

void  dQuadOps(const integer   rule, /* input: quadrature rule: GL or LL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   nq  , /* input: number of quadrature points   */
	       const double**  kp  , /* pointer to knot point storage        */
	       const double**  qp  , /* pointer to quadrature point storage  */
	       const double**  qw  , /* pointer to quadrature weight storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  );/* pointer to transposed deriv matrix   */

void  dMeshOps(const integer   old , /* input: element basis: STD or GLL     */
	       const integer   new , /* input: desired basis: STD or GLL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   ni  , /* input: number of interpolant points  */
	       const double**  mesh, /* pointer to interpolated mesh storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  );/* pointer to transposed deriv matrix   */

void  sQuadOps(const integer   rule, /* input: quadrature rule: GL or LL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   nq  , /* input: number of quadrature points   */
	       const float**   kp  , /* pointer to knot point storage        */
	       const float**   qp  , /* pointer to quadrature point storage  */
	       const float**   qw  , /* pointer to quadrature weight storage */
	       const float***  in  , /* pointer to interpolation matrix      */
	       const float***  it  , /* pointer to transposed interp matrix  */
	       const float***  dr  , /* pointer to derivative matrix         */
	       const float***  dt  );/* pointer to transposed deriv matrix   */

void  sMeshOps(const integer   old , /* input: element basis: STD or GLL     */
	       const integer   new , /* input: desired basis: STD or GLL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   ni  , /* input: number of interpolant points  */
	       const float**   mesh, /* pointer to interpolated mesh storage */
	       const float***  in  , /* pointer to interpolation matrix      */
	       const float***  it  , /* pointer to transposed interp matrix  */
	       const float***  dr  , /* pointer to derivative matrix         */
	       const float***  dt  );/* pointer to transposed deriv matrix   */

void dIntpOps (const integer basis,  /* element basis: STD or GLL            */
	       const integer np   ,  /* number of knot points                */
	       const double  r    ,  /* location of r in [-1, 1]             */
	       const double  s    ,  /* location of s in [-1, 1]             */
	       double*       inr  ,  /* 1D shape function at r               */
	       double*       ins  ,  /* 1D shape function at s               */
	       double*       dvr  ,  /* 1D shape function derivative at r    */
	       double*       dvs  ); /* 1D shape function derivative at s    */

void sIntpOps (const integer basis,
	       const integer np   ,
	       const float   r    ,
	       const float   s    ,
	       float*        inr  ,
	       float*        ins  ,
	       float*        dvr  ,
	       float*        dvs  );

void LegCoefGL (const integer  np,  /* input: number of points for Leg polys */
		const double** cd,  /* output: pointer to table of coeffs.   */
		const float**  cs); /* output: single-precision coeffs.      */

/* -- Routines from mapping.c: */

void edgemaps (const integer np, const integer dim,
	       integer** emap, integer** pmap);

/* -- Routines from family.c: */

void iadopt   (const integer, integer**);
void dadopt   (const integer, double** );
void sadopt   (const integer, float**  );

void iabandon (integer**);
void dabandon (double** );
void sabandon (float**  );

integer  FamilySize (integer*, integer*, integer*);

/* -- Routines from RCM.f: */

void genrcm_ (integer*, integer*, integer*, integer*, integer*, integer*);
#define genrcm(neqns, xadj, adjncy, perm, mask, xls) \
(_alpIreg[0] = neqns, genrcm_(_alpIreg, xadj, adjncy, perm, mask, xls))

void fnroot_ (integer*, integer*, integer*, integer*,
	      integer*, integer*, integer*);
#define fnroot(root, xadj, adncy, mask, nlvl, xls, ls) \
(_alpIreg[0] = root, _alpIreg[1] = nlvl,                   \
 fnroot_(_alpIreg, xadj, adjncy, mask, _alpIreg + 1, xls, ls))

void rcm_    (integer*, integer*, integer*, integer*,
	      integer*, integer*, integer*);
#define rcm(root, xadj, adjncy, mask, perm, ccsize, deg)  \
(_alpIreg[0] = root, rcm_(_alpIreg, xadj, adjncy, mask, perm, ccsize, deg))

/* -- Routines from fftpack.f (NETLIB/FFTPACK): */

void srffti_ (integer*, float*);
void srfftf_ (integer*, float*, float*);
void srfftb_ (integer*, float*, float*);

#define srffti(n,wsave)   (_alpIreg[0]=n, srffti_ (_alpIreg, wsave))
#define srfftf(n,r,wsave) (_alpIreg[0]=n, srfftf_ (_alpIreg, r, wsave))
#define srfftb(n,r,wsave) (_alpIreg[0]=n, srfftb_ (_alpIreg, r, wsave))

void drffti_ (integer*, double*);
void drfftf_ (integer*, double*, double*);
void drfftb_ (integer*, double*, double*);

#define drffti(n,wsave)   (_alpIreg[0]=n, drffti_ (_alpIreg, wsave))
#define drfftf(n,r,wsave) (_alpIreg[0]=n, drfftf_ (_alpIreg, r, wsave))
#define drfftb(n,r,wsave) (_alpIreg[0]=n, drfftb_ (_alpIreg, r, wsave))

/* -- Routines from canfft.f (Canuto FFT routines): */

void factor_ (integer*, integer*, integer*);

#define factor(n, nfac, ifac) (_alpIreg[0]=n, factor_ (_alpIreg, nfac, ifac))

void spreft_ (integer*, integer*, integer*, float*);
void smrcft_ (float*, integer*, integer*, float*,
	      integer*, integer*, float*, integer*);

#define spreft(n,nfac,ifac,trig)  \
(_alpIreg[0]=n, spreft_ (_alpIreg, nfac, ifac, trig))
#define smrcft(v, np, nz, w, nfac, ifac, trig, sign)  \
(_alpIreg[0]=np, _alpIreg[1]=nz, _alpIreg[2]=nfac, _alpIreg[3]=sign,  \
smrcft_(v, _alpIreg, _alpIreg+1, w, _alpIreg+2, ifac, trig, _alpIreg+3))

void dpreft_ (integer*, integer*, integer*, double*);
void dmrcft_ (double*, integer*, integer*, double*,
	      integer*, integer*, double*, integer*);

#define dpreft(n,nfac,ifac,trig)  \
(_alpIreg[0]=n, dpreft_ (_alpIreg, nfac, ifac, trig))
#define dmrcft(v, np, nz, w, nfac, ifac, trig, sign)  \
(_alpIreg[0]=np, _alpIreg[1]=nz, _alpIreg[2]=nfac, _alpIreg[3]=sign,  \
dmrcft_(v, _alpIreg, _alpIreg+1, w, _alpIreg+2, ifac, trig, _alpIreg+3))

/* -- Routines from temfftx.f (Temperton FFT routines): */

void prf235_ (integer*, integer*, integer*, integer*, integer*);
#define prf235(n,ip,iq,ir,ipqr2) (prf235_(n,ip,iq,ir,ipqr2))

void dsetpf_ (double*, integer*, integer*, integer*, integer*);
void dmpfft_ (double*, double*,  integer*, integer*, integer*,
	      integer*, integer*, double*, integer*);

#define dsetpf(trig,n,ip,iq,ir)                               \
(_alpIreg[0]=n,_alpIreg[1]=ip,_alpIreg[2]=iq, _alpIreg[3]=ir, \
dsetpf_(trig,_alpIreg,_alpIreg+1,_alpIreg+2,_alpIreg+3))
#define dmpfft(v,w,np,nz,ip,iq,ir,trig,isign)                 \
(_alpIreg[0]=np,_alpIreg[1]=nz,_alpIreg[2]=ip,                \
_alpIreg[3]=iq, _alpIreg[4]=ir,_alpIreg[5]=isign,             \
dmpfft_(v,w,_alpIreg,_alpIreg+1,_alpIreg+2,                   \
_alpIreg+3,_alpIreg+4,trig,_alpIreg+5))

void ssetpf_ (float*, integer*, integer*, integer*, integer*);
void smpfft_ (float*, float*, integer*, integer*, integer*,
	      integer*, integer*, float*, integer*);

#define ssetpf(trig,n,ip,iq,ir)                               \
(_alpIreg[0]=n,_alpIreg[1]=ip,_alpIreg[2]=iq, _alpIreg[3]=ir, \
ssetpf_(trig,_alpIreg,_alpIreg+1,_alpIreg+2,_alpIreg+3))
#define smpfft(v,w,np,nz,ip,iq,ir,trig,isign)                 \
(_alpIreg[0]=np,_alpIreg[1]=nz,_alpIreg[2]=ip,                \
_alpIreg[3]=iq, _alpIreg[4]=ir,_alpIreg[5]=isign,             \
smpfft_(v,w,_alpIreg,_alpIreg+1,_alpIreg+2,                   \
_alpIreg+3,_alpIreg+4,trig,_alpIreg+5))

/* -- Routines from matops.F: */

void dgrad2_ (double*,double*,double*,double*,double*,double*,
	      integer*,integer*);
void sgrad2_ (float*, float*, float*, float*, float*, float*,
	      integer*, integer*);

#define dgrad2(u,v,ur,vs,dv,dt,np,nel)                          \
(_alpIreg[0]=np,_alpIreg[1]=nel,dgrad2_(u,v,ur,vs,dv,dt,_alpIreg,_alpIreg+1))
#define sgrad2(u,v,ur,vs,dv,dt,np,nel)                          \
(_alpIreg[0]=np,_alpIreg[1]=nel,sgrad2_(u,v,ur,vs,dv,dt,_alpIreg,_alpIreg+1))

#if defined(SX)

/* -- Routines from NEC FFT library: floating precision depends on library. */

void rftfax_ (integer*,integer*,real*);
void rfft_   (real*,real*,real*,integer*,integer*,integer*,real*);

#define rftfax(n,ifax,trigs)                        \
(_alpIreg[0]=n,rftfax_(_alpIreg,ifax,trigs))
#define rfft(r,w,trigs,ifax,n,l,xnorm)              \
(_alpIreg[0]=n,_alpIreg[1]=l,_alpDreg[0]=xnorm,     \
rfft_(r,w,trigs,ifax,_alpIreg,_alpIreg+1,_alpDreg))

#endif

/* -- Routines from fourier.c */

void sDFTr (float*,  const integer, const integer, const integer);
void dDFTr (double*, const integer, const integer, const integer);

/* -- Routines from filter.c */

void bvdFilter (const integer,const integer,const integer, const real, real*);

/* -- Routines from message.c: */

void message_init      (int*, char***);
void message_stop      ();
void message_sync      ();
void message_dsend     (double*,  const integer, const integer);
void message_drecv     (double*,  const integer, const integer);
void message_ssend     (float*,   const integer, const integer);
void message_srecv     (float*,   const integer, const integer);
void message_isend     (integer*, const integer, const integer);
void message_irecv     (integer*, const integer, const integer);
void message_dexchange (double*,  const integer, const integer,const integer);
void message_sexchange (float*,   const integer, const integer,const integer);
void message_iexchange (integer*, const integer, const integer,const integer);

#endif
