#ifndef FEMLIB_H
#define FEMLIB_H
///////////////////////////////////////////////////////////////////////////////
// C++ function prototypes for routines in library libfem.a
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cfemdef.h>
#include <polylib.h>

extern "C" {

extern char buf[STR_MAX];

// -- Routines from initial.y:

void  yy_initialize (void);
real  yy_interpret  (const char*);

void  yy_vec_init   (const char*, const char*);
void  yy_vec_interp (const integer, ...);

void  yy_help       (void);
void  yy_show       (void);

// -- Routines from polyops.c:

void  dermat_g (const integer, const real*, const integer,
		const real*, real**, real**);
void  intmat_g (const integer, const real*, const integer,
		const real*, real**, real**);
void  uniknot  (const integer, double*);

#if 0
real  pnleg    (const real, const integer);
real  pnmod    (const real, const integer);	
#endif

// -- Routines from operators.c:

void zquad (const real**  point ,  // Quadrature points.
	    const real**  weight,  // Quadrature weights.
	    const real**  dv    ,  // Derivative operator at points.
	    const real**  dt    ,  // (Transpose) derivative operator.
	    const integer np    ,  // Input: Number of quadrature points.
	    const char    rule  ,  // Input: 'G'auss, 'R'adau, or 'L'obatto.
	    const real    alpha ,  // Input: Jacobi polynomial constant.
	    const real    beta  ); // Input: Jacobi polynomial constant.

void proj (const real**  IN    ,  // Interpolant operator matrix
	   const real**  IT    ,  // Transposed interpolant operator matrix
	   const integer nfrom ,  // Input: Number of "from" points
	   const char    rulefr,  // Input: 'G', 'R', 'L', or 'U', "from"
	   const real    alphfr,  // Input: Jacobi constant, "from"
	   const real    betafr,  // Input: Jacobi constant, "from"
	   const integer nto   ,  // Input: Number of "to" points
	   const char    ruleto,  // Input: 'G', 'R', 'L', or 'U', "to"
	   const real    alphto,  // Input: Jacobi constant, "to"
	   const real    betato); // Input: Jacobi constant, "to"

void intp (real*         inr   ,  // 1D shape function at r
	   real*         ins   ,  // 1D shape function at s
	   real*         dvr   ,  // 1D shape function derivative at r
	   real*         dvs   ,  // 1D shape function derivative at s
	   const integer nr    ,  // Input: number of points, r-direction
	   const char    basisr,  // Input: 'G', 'R', 'L', r-direction
	   const real    alphar,  // Input: Jacobi constant, r-direction
	   const real    betar ,  // Input: Jacobi constant, r-direction
	   const integer ns    ,  // Input: number of points, s-direction
	   const char    basiss,  // Input: 'G', 'R', 'L', s-direction
	   const real    alphas,  // Input: Jacobi constant, s-direction
	   const real    betas ,  // Input: Jacobi constant, s-direction
	   const real    r     ,  // Input: location of r in [-1, 1]
	   const real    s     ); // Input: location of s in [-1, 1]


void dQuadOps (const integer rule,  // input: quadrature rule: GL or LL
	       const integer np  ,  // input: number of knot points
	       const integer nq  ,  // input: number of quadrature points
	       const real**  kp  ,  // pointer to knot point storage
	       const real**  qp  ,  // pointer to quadrature point storage
	       const real**  qw  ,  // pointer to quadrature weight storage
	       const real*** in  ,  // pointer to interpolation matrix
	       const real*** it  ,  // pointer to transposed interp matrix
	       const real*** dr  ,  // pointer to collocation deriv matrix
	       const real*** dt  ); // pointer to transposed  deriv matrix

void dMeshOps (const integer oldb,  // input: old element basis: STD or GLL
	       const integer newb,  // input: new/desired basis: STD or GLL
	       const integer np  ,  // input: number of knot points
	       const integer ni  ,  // input: number of interpolant points
	       const real**  mesh,  // pointer to interpolated mesh storage
	       const real*** in  ,  // pointer to interpolation matrix
	       const real*** it  ,  // pointer to transposed interp matrix
	       const real*** dr  ,  // pointer to interpolation deriv matrix
	       const real*** dt  ); // pointer to transposed    deriv matrix

void dIntpOps (const integer basis,  // element basis: STD or GLL
	       const integer np   ,  // number of knot points
	       const real    r    ,  // location of r in [-1, 1]
	       const real    s    ,  // location of s in [-1, 1]
	       real*         inr  ,  // 1D shape function at r
	       real*         ins  ,  // 1D shape function at s
	       real*         dvr  ,  // 1D shape function derivative at r
	       real*         dvs  ); // 1D shape function derivative at s

void dglldpc (const integer  np,     // input:  number of points for Leg polys
	      const real**   cd);    // output: pointer to table of coeffs

void dglldpt (const integer np,     // input:  number of points for DLT
	      const real**  fw,     // output: 1D forward transform matrix
	      const real**  ft,     // output: transpose of fw
	      const real**  bw,     // output: 1D inverse transform matrix
	      const real**  bt,     // output: transpose of bw
	      const real**  fu,     // output: 2D forward transform matrix
	      const real**  bu);    // output: 2D inverse transform matrix

void dglmdpc (const integer np,     // input:  no. of points for expansions.
	      const real**  cd);    // output: pointer to table of coeffs

void dglmdpt (const integer np,     // input:  number of points for DPT
	      const real**  fw,     // output: 1D forward transform matrix
	      const real**  ft,     // output: transpose of fw
	      const real**  bw,     // output: 1D inverse transform matrix
	      const real**  bt,     // output: transpose of bw
	      const real**  fu,     // output: 2D forward transform matrix
	      const real**  bu);    // output: 2D inverse transform matrix

// -- Routines from mapping.c:

void edgemaps (const integer np, const integer dim, 
	       integer** emap, integer** pmap);

// -- Routines from family.c:

void iadopt   (const integer, integer**);
void dadopt   (const integer, real**);
void sadopt   (const integer, float**);

void iabandon (integer**);
void dabandon (real**);
void sabandon (float**);

integer  FamilySize (integer*, integer*, integer*);

// -- Routines from RCM.f (NETLIB/SPARSPAK):

void F77NAME(genrcm) (const integer&, integer*, integer*, integer*,
		      integer*, integer*);
void F77NAME(fnroot) (integer&, integer*, integer*, integer*,
		      integer&, integer*, integer*);
void F77NAME(rcm)    (const integer&, integer*, integer*, integer*,
		      integer*, integer&, integer*);

// -- Routines from fftpack.f (NETLIB/FFTPACK):

#if 1
void F77NAME(drffti) (const integer&,real*,integer*);
void F77NAME(drfftf) (const integer&,real*,real*,const real*,const integer*);
void F77NAME(drfftb) (const integer&,real*,real*,const real*,const integer*);
#else
void F77NAME(drffti) (const integer&, real*);
void F77NAME(drfftf) (const integer&, real*, const real*);
void F77NAME(drfftb) (const integer&, real*, const real*);
#endif

// -- Routines from canfft.f:

void F77NAME(factor) (const integer&, integer&, integer*);

void F77NAME(dpreft) (const integer&, integer&, integer*, real*);
void F77NAME(dmrcft) (real*, const integer&, const integer&, real*,
		      const integer&, integer*, real*, const integer&);
void F77NAME(dfft1)  (real*, real*, const integer&, const integer&,
		      const integer*, const integer&, const real*,
		      const integer&);

// -- Routines from temfft.f:

void F77NAME(prf235) (integer&, integer&, integer&, integer&, integer&);

void F77NAME(dsetpf) (const real*, const integer&, const integer&,
		      const integer&, const integer&);
void F77NAME(dmpfft) (real*, real*, const integer&, const integer&,
		      const integer&, const integer&, const integer&, 
		      const real*, const integer&);
void F77NAME(dgpfa)  (real*, real*, const real*, const integer&,
		      const integer&, const integer&, const integer&,
		      const integer&, const integer&, const integer&,
		      const integer&);

// -- Routines from fourier.c:

void dDFTr (real*, const integer, const integer, const integer);

// -- Routines from filter.c

void bvdFilter (const integer,const real,const real, const real, real*);

// -- Routines from message.c:

void message_init      (int*, char***);
void message_stop      ();
void message_sync      ();
void message_dsend     (real*,    const integer, const integer);
void message_drecv     (real*,    const integer, const integer);
void message_ssend     (float*,   const integer, const integer);
void message_srecv     (float*,   const integer, const integer);
void message_isend     (integer*, const integer, const integer);
void message_irecv     (integer*, const integer, const integer);
void message_dexchange (real*,    const integer, const integer,const integer);
void message_sexchange (float*,   const integer, const integer,const integer);
void message_iexchange (integer*, const integer, const integer,const integer);

// -- Routines from matops.F:

void F77NAME(dgrad2) (const real*, const real*, real*, real*, const real*,
		      const real*, const integer&, const integer&,
		      const integer&);
void F77NAME(dtpr2d) (const real*, real*, real*, const real*, const real*,
		      const integer&, const integer&, const integer&);
}

class Femlib {
public:

  static void initialize (int* argc, char*** argv)
    { message_init (argc, argv); }
  static void prepVec (const char* v, const char* f)
    { yy_vec_init (v, f); }
  static void finalize () 
    { message_stop (); }
  static void synchronize ()
    { message_sync(); }

#if 0
  // -- Don't know why this won't compile OK but:
  static void parseVec (const integer n ... )
    { yy_vec_interp (n ... ); }
#else
  #define Femlib__parseVec yy_vec_interp
#endif

  static void value (const char* s, const real p)
    { sprintf (buf, "%s = %.17g", s, p); yy_interpret (buf); }
  static void ivalue (const char* s, const integer p)
    { sprintf (buf, "%s = %1d", s, p); yy_interpret (buf); }
  static real value (const char* s)
    { return yy_interpret (s); }
  static integer ivalue (const char* s)
    { return static_cast<integer>(yy_interpret (s)); }
  
  static void equispacedMesh (const integer np, real* z)
    { uniknot (np, z); }

#if 0
  static void     LagrangeInt (const integer N,  const real* zero,
                               const integer I,  const real* x   ,
                               real**        IN, real**      IT  )
    { intmat_g (N, zero, I, x, IN, IT); }
  static void     LagrangeDer (const integer N,  const real* zero,
                               const integer I,  const real* x   ,
                               real**      IN, real**      IT  )
    { dermat_g (N, zero, I, x, IN, IT); }
  static void     GLLzw     (const integer N, real* z, real* w)
    { zwgll (z, w, N); }
  static real   LegendreVal (const integer N, const real x)
    { return pnleg (x, N); }
  static real   ModalVal    (const integer N, const real x)
    { return pnmod (x, N); }
  
#endif
#if 0
  static void quad (const integer rule, const integer np  ,
		    const integer nq  , const real**  kp  ,
		    const real**  qp  , const real**  qw  ,
		    const real*** in  , const real*** it  ,
		    const real*** dr  , const real*** dt  ) 
    { dQuadOps (rule, np, nq, kp, qp, qw, in, it, dr, dt); }

  static void mesh (const integer oldb, const integer newb,
		    const integer np  , const integer ni  ,
		    const real**  mesh,
		    const real*** in  , const real*** it  ,
		    const real*** dr  , const real*** dt  )
    { dMeshOps (oldb, newb, np, ni, mesh, in, it, dr, dt); }

  static void interp (const integer basis, const integer np,
		      const real r       , const real    s ,
		      real* inr, real* ins, real* dvr, real* dvs)
    { dIntpOps (basis, np, r, s, inr, ins, dvr, dvs); }
#endif
  static void quadrature (const real**  point ,
			  const real**  weight,
			  const real**  dv    ,
			  const real**  dt    ,
			  const integer np    ,
			  const char    rule  ,
			  const real    alpha ,
			  const real    beta  )
  { zquad (point, weight, dv, dt, np, rule, alpha, beta); }

  static void projection (const real**  IN    ,
			  const real**  IT    ,
			  const integer nfrom ,
			  const char    rulefr,
			  const real    alphfr,
			  const real    betafr,
			  const integer nto   ,
			  const char    ruleto,
			  const real    alphto,
			  const real    betato)
  { proj (IN, IT, nfrom,rulefr,alphfr,betafr, nto,ruleto,alphto,betato); }

  static void interpolation (real*         inr   ,
			     real*         ins   ,
			     real*         dvr   ,
			     real*         dvs   ,
			     const integer nr    ,
			     const char    basisr,
			     const real    alphar,
			     const real    betar ,
			     const integer ns    ,
			     const char    basiss,
			     const real    alphas,
			     const real    betas ,
			     const real    r     ,
			     const real    s     )
  { intp (inr,ins,dvr,dvs,nr,basisr,alphar,betar,ns,basiss,alphas,betas,r,s); }


  static void legCoef (const integer n,
		       const real**  c)
    { dglldpc (n, c); }

  static void legTran (const integer n,
		       const real**  f1, const real** ft,
		       const real**  i1, const real** it,
		       const real**  f2, const real** i2)
    { dglldpt (n, f1, ft, i1, it, f2, i2); }

  static void modCoef (const integer n,
		       const real**  c)
    { dglmdpc (n, c); }

  static void modTran (const integer n,
		       const real**  f1, const real** ft,
		       const real**  i1, const real** it,
		       const real**  f2, const real** i2)
    { dglmdpt (n, f1, ft, i1, it, f2, i2); }

  static void buildMaps (const integer np, const integer dim,
			 integer** e, integer** i)
    { edgemaps (np, dim, e, i); }

  static void adopt (const integer np, integer** v)
    { iadopt (np, v); }
  static void adopt (const integer np, real** v)
    { dadopt (np, v); }
  static void adopt (const integer np, float** v)
    { sadopt (np, v); }
  static void abandon (integer** v)
    { iabandon (v); }
  static void abandon (real** v)
    { dabandon (v); }
  static void abandon (float** v)
    { sabandon (v); }
  static integer  fwords (integer* ni, integer* nd, integer* ns)
    { return FamilySize (ni, nd, ns); }

  static void genrcm (const integer& n, integer* x, integer* a,
		      integer* p, integer* m, integer* l) 
    { F77NAME(genrcm) (n, x, a, p, m, l); }
  static void fnroot (integer& r, integer* x, integer* a, integer* m,
		      integer& n, integer* l, integer* p) 
    { F77NAME(fnroot) (r, x, a, m, n, l, p); }
  static void rcm (const integer& n, integer* x, integer* a, integer* m,
		   integer* p, integer& c, integer* l)
    { F77NAME(rcm) (n, x, a, m, p, c, l); }

  static void rffti (const integer& n , real* w, integer* i)
    { F77NAME(drffti) (n, w, i); }
  static void rfftf (const integer& n, real* c, real* ch,
	const real* wa, const integer* ifac)
    { F77NAME(drfftf) (n, c, ch, wa, ifac); }
  static void rfftb (const integer& n, real* c, real* ch,
  	const real* wa, const integer* ifac) 
    { F77NAME(drfftb) (n, c, ch, wa, ifac); }

  static void DFTr (real*  data, const integer nz, const integer np,
		    const integer sign)
    { dDFTr (data, nz, np, sign); }

  static void preft  (const integer& n, integer& nfax, integer* ifax,
		      real* trig) 
    { F77NAME(dpreft) (n, nfax, ifax, trig); }
  static void mrcft  (real* v, const integer& np, const integer& nz,
		      real* w, const integer& nfx, integer* ifx,
		      real* trig, const integer& sign)
    { F77NAME(dmrcft) (v, np, nz, w, nfx, ifx, trig, sign); }
  static void primes23 (const integer& n, integer& np, integer* fact)
    { F77NAME(factor) (n, np, fact); }
  static void fft1 (real* a, real* c, const integer& n,
		    const integer& nfax,
                    const integer* ifax, const integer& isign,
		    const real* trig, const integer& len)
    { F77NAME(dfft1) (a,c,n,nfax,ifax,isign,trig,len); }  

  static void primes235 (integer& n,
			 integer&ip, integer&iq, integer& ir, integer& ipqr2)
    { F77NAME(prf235) (n, ip, iq, ir, ipqr2); }
  static void setpf (const real* t, const integer& n, const integer& ip,
		     const integer& iq, const integer& ir)
    { F77NAME(dsetpf) (t, n, ip, iq, ir); }
  static void mpfft (real* v, real* w, const integer& np,
		     const integer& nz,
		     const integer& ip, const integer& iq, const integer& ir,
		     const real* trig, const integer& isign)
    { F77NAME(dmpfft) (v, w, np, nz, ip, iq, ir, trig, isign); }
  static void gpfa (real* a, real* b, const real* trig,
		    const integer& inc, const integer& jump, const integer& n,
		    const integer& ip, const integer& iq, const integer& ir,
		    const integer& lot, const integer& isign)
    { F77NAME(dgpfa) (a, b, trig, inc, jump, n, ip, iq, ir, lot, isign); }

  static void send  (real*    data, const integer N, const integer tgt)
    { message_dsend (data, N, tgt); }
  static void recv  (real*    data, const integer N, const integer src)
    { message_drecv (data, N, src); }
  static void send  (float*   data, const integer N, const integer tgt)
    { message_ssend (data, N, tgt); }
  static void recv  (float*   data, const integer N, const integer src)
    { message_srecv (data, N, src); }
  static void send  (integer* data, const integer N, const integer tgt)
    { message_isend (data, N, tgt); }
  static void recv  (integer* data, const integer N, const integer src)
    { message_irecv (data, N, src); }
  static void exchange  (real* data, const integer nZ,
			 const integer nP, const integer sign)
    { message_dexchange (data, nZ, nP, sign); }
  static void exchange  (float* data, const integer nZ,
			 const integer nP, const integer sign)
    { message_sexchange (data, nZ, nP, sign); }
  static void exchange  (integer* data,  const integer nZ,
			 const integer nP, const integer sign)
    { message_iexchange (data, nZ, nP, sign); }

  static void grad2 (const real* x, const real* y, real* xr, real* ys,
		     const real* dv, const real* dt,
		     const integer& nr, const integer& ns, const integer& nel)
    { F77NAME(dgrad2) (x, y, xr, ys, dv, dt, nr, ns, nel); }

  static void tpr2d (const real* x, real* y, real* t,
		     const real* dv, const real* dt,
		     const integer& nr, const integer& ns, const integer& nel)
    { F77NAME(dtpr2d) (x, y, t, dv, dt, nr, ns, nel); }

  static void erfcFilter (const integer N, const real p, const real s,
			  const real a, real* f)
    { bvdFilter (N, p, s, a, f); }

};

#endif
