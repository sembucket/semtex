#ifndef C_FEMLIB_H
#define C_FEMLIB_H
/*****************************************************************************
 *         FUNCTION PROTOTYPES FOR ROUTINES IN LIBRARY LIBFEM.A
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>

/* -- Routines from initial.c: */

void    yy_initialize (void);
double  yy_interpret  (const char*);

void    yy_vec_init   (const char*, const char*);
void    yy_vec_interp (const integer, ...);

void    yy_help       (void);
void    yy_show       (void);

/* -- Routines from polyops.c: */

void   dermat_g (const integer, const double*, const integer,
		 const double*, double**, double**);
void   intmat_g (const integer, const double*, const integer,
		 const double*, double**, double**);
void   dermat_k (const integer, const double*, double**, double**);


void   JACG     (const integer, const double, const double, double*);
void   JACGR    (const integer, const double, const double, double*);
void   JACGL    (const integer, const double, const double, double*);

void   ZWGL     (double*, double*, const integer);
void   ZWGRL    (double*, double*, const integer);
void   ZWGLL    (double*, double*, const integer);
void   ZWGLJ    (double*, double*, const double, const double, const integer);

void   DGLL     (const integer, const double*, double**, double**);

double PNLEG    (const double,  const integer);
double PNDLEG   (const double,  const integer);
double PND2LEG  (const double,  const integer);
double PNMOD    (const double, const integer);

void   uniknot  (const integer, double*);

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

void dIntpOps (const integer basis,  /* element basis: STD or GLL            */
	       const integer np   ,  /* number of knot points                */
	       const double  r    ,  /* location of r in [-1, 1]             */
	       const double  s    ,  /* location of s in [-1, 1]             */
	       double*       inr  ,  /* 1D shape function at r               */
	       double*       ins  ,  /* 1D shape function at s               */
	       double*       dvr  ,  /* 1D shape function derivative at r    */
	       double*       dvs  ); /* 1D shape function derivative at s    */
#if 0
void dglldpc (const integer  np,    /* input: number of points for Leg polys */
	      const double** cd);   /* output: pointer to table of coeffs.   */

void dglldpt (const integer  np,    /* input:  number of points for DLT      */
	      const double** fw,    /* output: 1D forward transform matrix   */
	      const double** ft,    /* output: transpose of fw               */
	      const double** bw,    /* output: 1D inverse transform matrix   */
	      const double** bt,    /* output: transpose of bw               */
	      const double** fu,    /* output: 2D forward transform matrix   */
	      const double** bu);   /* output: 2D inverse transform matrix   */
#endif

/* -- Routines from mapping.c: */

void edgemaps (const integer np, const integer dim,
	       integer** emap, integer** pmap);

/* -- Routines from family.c: */

void iadopt   (const integer, integer**);
void dadopt   (const integer, double** );
void sadopt   (const integer, float**  );

void iabandon (integer**);
void dabandon (double**);
void sabandon (float**);

integer FamilySize (integer*, integer*, integer*);

/* -- Routines from RCM.f: */

void genrcm_ (integer*, integer*, integer*, integer*, integer*, integer*);
#define genrcm(neqns, xadj, adjncy, perm, mask, xls) \
(_vecIreg[0] = neqns, genrcm_(_vecIreg, xadj, adjncy, perm, mask, xls))

void fnroot_ (integer*, integer*, integer*, integer*,
	      integer*, integer*, integer*);
#define fnroot(root, xadj, adncy, mask, nlvl, xls, ls) \
(_vecIreg[0] = root, _vecIreg[1] = nlvl,                   \
 fnroot_(_vecIreg, xadj, adjncy, mask, _vecIreg + 1, xls, ls))

void rcm_    (integer*, integer*, integer*, integer*,
	      integer*, integer*, integer*);
#define rcm(root, xadj, adjncy, mask, perm, ccsize, deg)  \
(_vecIreg[0] = root, rcm_(_vecIreg, xadj, adjncy, mask, perm, ccsize, deg))

/* -- Routines from fftpack.f (NETLIB/FFTPACK): */

#if 1
void drffti_ (integer*, double*, integer*);
void drfftf_ (integer*, double*, double*, double*, integer*);
void drfftb_ (integer*, double*, double*, double*, integer*);

#define drffti(n,wa,ifac) \
  (_vecIreg[0]=n, drffti_ (_vecIreg, wa, ifac))
#define drfftf(n,c,ch,wa,ifac) \
  (_vecIreg[0]=n, drfftf_ (_vecIreg, c, ch, wa, ifac))
#define drfftb(n,c,ch,wa,ifac) \
  (_vecIreg[0]=n, drfftb_ (_vecIreg, c, ch, wa, ifac))
#else
void drffti_ (integer*, double*);
void drfftf_ (integer*, double*, double*);
void drfftb_ (integer*, double*, double*);

#define drffti(n,wsave)   (_vecIreg[0]=n, drffti_ (_vecIreg, wsave))
#define drfftf(n,r,wsave) (_vecIreg[0]=n, drfftf_ (_vecIreg, r, wsave))
#define drfftb(n,r,wsave) (_vecIreg[0]=n, drfftb_ (_vecIreg, r, wsave))
#endif
/* -- Routines from canfft.f (Canuto FFT routines): */

void factor_ (integer*, integer*, integer*);

#define factor(n, nfac, ifac) (_vecIreg[0]=n, factor_ (_vecIreg, nfac, ifac))

#define dpreft(n,nfac,ifac,trig)  \
(_vecIreg[0]=n, dpreft_ (_vecIreg, nfac, ifac, trig))
#define dmrcft(v, np, nz, w, nfac, ifac, trig, sign)  \
(_vecIreg[0]=np, _vecIreg[1]=nz, _vecIreg[2]=nfac, _vecIreg[3]=sign,  \
dmrcft_(v, _vecIreg, _vecIreg+1, w, _vecIreg+2, ifac, trig, _vecIreg+3))

/* -- Routines from temfftd.f (Temperton FFT routines): */

void prf235_ (integer*, integer*, integer*, integer*, integer*);
#define prf235(n,ip,iq,ir,ipqr2) (prf235_(n,ip,iq,ir,ipqr2))

void dsetpf_ (double*, integer*, integer*, integer*, integer*);
void dmpfft_ (double*, double*,  integer*, integer*, integer*,
	      integer*, integer*, double*, integer*);

#define dsetpf(trig,n,ip,iq,ir)                               \
(_vecIreg[0]=n,_vecIreg[1]=ip,_vecIreg[2]=iq, _vecIreg[3]=ir, \
dsetpf_(trig,_vecIreg,_vecIreg+1,_vecIreg+2,_vecIreg+3))
#define dmpfft(v,w,np,nz,ip,iq,ir,trig,isign)                 \
(_vecIreg[0]=np,_vecIreg[1]=nz,_vecIreg[2]=ip,                \
_vecIreg[3]=iq, _vecIreg[4]=ir,_vecIreg[5]=isign,             \
dmpfft_(v,w,_vecIreg,_vecIreg+1,_vecIreg+2,                   \
_vecIreg+3,_vecIreg+4,trig,_vecIreg+5))

/* -- Routines from matops.F: */

void dgrad2_ (double*,double*,double*,double*,double*,double*,
	      integer*,integer*);

#define dgrad2(u,v,ur,vs,dv,dt,np,nel)                          \
(_vecIreg[0]=np,_vecIreg[1]=nel,dgrad2_(u,v,ur,vs,dv,dt,_vecIreg,_vecIreg+1))

#if defined(SX)

/* -- Routines from NEC FFT library: floating precision depends on library. */

void rftfax_ (integer*,integer*,real*);
void rfft_   (real*,real*,real*,integer*,integer*,integer*,real*);

#define rftfax(n,ifax,trigs)                        \
(_vecIreg[0]=n,rftfax_(_vecIreg,ifax,trigs))
#define rfft(r,w,trigs,ifax,n,l,xnorm)              \
(_vecIreg[0]=n,_vecIreg[1]=l,_vecDreg[0]=xnorm,     \
rfft_(r,w,trigs,ifax,_vecIreg,_vecIreg+1,_vecDreg))

#endif

/* -- Routines from fourier.c */

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
