/*****************************************************************************
 * operators.c:  operators for mesh points, derivatives, quadratures.
 *
 * Copyright (c) Hugh Blackburn 1994--2000.
 *
 * All 2D matrices have row-major ordering.
 * All routines are double precision.
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <femdef.h>
#include <alplib.h>
#include <femlib.h>

typedef struct dquadopr {	/* ---- quadrature operator information  --- */
  integer          rule    ;	/* quadrature rule: GL or LL                 */
  integer          np      ;	/* number of interpolant knot points         */
  integer          nq      ;	/* number of quadrature points               */
  double*          knot    ;	/* Lagrange knots    on [-1, 1]              */
  double*          quad    ;	/* Quadrature points on [-1, 1]              */
  double*          weight  ;	/* Quadrature weights                        */
  double**         interp  ;	/* Interpolant operator: knots-->quadrature  */
  double**         interpT ;	/* Transpose of the above                    */
  double**         deriv   ;	/* Derivative operator:  knots-->quadrature  */
  double**         derivT  ;	/* Transpose of the above                    */
  struct dquadopr* next    ;	/* link to next one                          */
} dQuadOpr;			/* ----------------------------------------- */

static dQuadOpr* dQhead = 0;


void dQuadOps (const integer   rule, /* input: quadrature rule: GL or LL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   nq  , /* input: number of quadrature points   */
	       const double**  kp  , /* pointer to knot point storage        */
	       const double**  qp  , /* pointer to quadrature point storage  */
	       const double**  qw  , /* pointer to quadrature weight storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  ) /* pointer to transposed deriv matrix   */
/* ------------------------------------------------------------------------- *
 * Maintain/return QUADRATURE operators for finite elements with GLL BASIS
 * FUNCTIONS.  Quadrature rule may be at nodes (LL, nq = np enforced) or at
 * Gauss(-Legendre) points (GL).
 *
 * Operators are defined on the interval [-1, 1].
 *
 * If the required operators are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 * ------------------------------------------------------------------------- */
{
  char               routine[] = "dQuadOps";
  register integer   found = 0;
  register dQuadOpr* p;

  for (p = dQhead; p; p=p->next) {
    found = p->rule == rule && p->np == np && p->nq == nq;
    if (found) break;
  }

  if (!found) {

    if (rule != LL && rule != GL)
      message (routine, "unrecognized quadrature rule", ERROR);

    p = (dQuadOpr *) calloc (1, sizeof (dQuadOpr));
    if (dQhead) p -> next = dQhead;
    dQhead = p;

    p -> rule = rule;
    p -> np   = np;
    p -> nq   = (rule == GL) ? nq : np;

    if (rule == LL && np != nq)
      message (routine, "np != nq in LL rule...enforcing equality", WARNING);

    if (rule == LL) {

      p -> knot    = dvector (0, np-1);
      p -> quad    = p -> knot;
      p -> weight  = dvector (0, np-1);
      p -> interp  = (double**) 0;	/* No interpolation needed. */
      p -> interpT = (double**) 0;
      p -> deriv   = dmatrix (0, np-1, 0, np-1);
      p -> derivT  = dmatrix (0, np-1, 0, np-1);

      zwgll (p->knot, p->weight, np);
      dgll  (np, p->knot, p->deriv, p->derivT);

    } else {

      p -> knot    = dvector (0, np-1);
      p -> quad    = dvector (0, nq-1);
      p -> weight  = dvector (0, nq-1);
      p -> interp  = dmatrix (0, nq-1, 0, np-1);
      p -> interpT = dmatrix (0, np-1, 0, nq-1);
      p -> deriv   = dmatrix (0, nq-1, 0, np-1);
      p -> derivT  = dmatrix (0, np-1, 0, nq-1);

      jacgl    (np-1, 0.0, 0.0, p->knot);     /* Knots at Lobatto points.    */
      zwgl     (p->quad, p->weight, nq);      /* Quadrature at Gauss points. */
      intmat_g (np, p->knot, nq, p->quad, p->interp, p->interpT);
      dermat_g (np, p->knot, nq, p->quad, p->deriv,  p->derivT );

    }
  }

  /* -- p now points to valid storage: return requested operators. */

  if (kp) *kp = (const double* ) p -> knot;
  if (qp) *qp = (const double* ) p -> quad;
  if (qw) *qw = (const double* ) p -> weight;
  if (in) *in = (const double**) p -> interp;
  if (it) *it = (const double**) p -> interpT;
  if (dr) *dr = (const double**) p -> deriv;
  if (dt) *dt = (const double**) p -> derivT;
}


typedef struct dmeshopr {	/* ------- mesh operator information ------- */
  integer          oldbasis;	/* STD or GLL                                */
  integer          newbasis;	/* STD or GLL                                */
  integer          np      ;	/* Number of points on original mesh         */
  integer          ni      ;	/* Number of points in interpolated mesh     */
  double*          mesh    ;	/* Location of points on interpolated mesh   */
  double**         interp  ;	/* Interpolant operator: original-->new      */
  double**         interpT ;	/* Transpose of the above                    */
  double**         deriv   ;	/* Derivative operator: original-->new       */
  double**         derivT  ;	/* Transpose of the above                    */
  struct dmeshopr* next    ;	/* link to next one                          */
} dMeshOpr;			/* ----------------------------------------- */

static dMeshOpr* dMhead = 0;


void dMeshOps (const integer   old , /* input: element basis: STD or GLL     */
	       const integer   new , /* input: desired basis: STD or GLL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   ni  , /* input: number of interpolant points  */
	       const double**  mesh, /* pointer to interpolated mesh storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  ) /* pointer to transposed deriv matrix   */
/* ------------------------------------------------------------------------- *
 * Maintain/return INTERPOLATION operators for STD and GLL meshes.
 * Operators are defined on the interval [-1, 1].
 *
 * It is legal to ask for interpolation onto the same mesh as input, in
 * which case the mesh storage and interpolation matrices are NULL.
 *
 * If the required operators are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 * ------------------------------------------------------------------------- */
{
  char               routine[] = "dMeshOps";
  register integer   found = 0;
  register dMeshOpr* p;
  double*            oldmesh;

  for (p = dMhead; p; p = p->next) {
    found = p->oldbasis == old
      &&    p->newbasis == new 
      &&    p->np       == np
      &&    p->ni       == ni;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */

    p = (dMeshOpr *) calloc (1, sizeof (dMeshOpr));
    if (dMhead) p -> next = dMhead;
    dMhead = p;

    p -> oldbasis = old;
    p -> newbasis = new;
    p -> np       = np;
    p -> ni       = ni;

    oldmesh       = dvector (0, np-1);
    p -> mesh     = dvector (0, ni-1);
    p -> deriv    = dmatrix (0, ni-1, 0, np-1);
    p -> derivT   = dmatrix (0, np-1, 0, ni-1);
    
    if (old == new && np == ni) {	/* Meshes are the same; */
      p -> interp  = (double**) 0;      /* no interpolation.    */
      p -> interpT = (double**) 0;
    } else {
      p -> interp  = dmatrix (0, ni-1, 0, np-1);
      p -> interpT = dmatrix (0, np-1, 0, ni-1);
    }   

    if (old == GLL && new == GLL) {

      if (np == ni) {  /* Meshes are the same. */
	jacgl (np-1, 0.0, 0.0, p->mesh);
	dgll  (np, p->mesh, p->deriv, p->derivT);
      } else {
	jacgl    (np-1, 0.0, 0.0, oldmesh);
	jacgl    (ni-1, 0.0, 0.0, p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
      }

    } else if (old == GLL && new == STD) {

      jacgl    (np-1, 0.0, 0.0,  oldmesh);
      uniknot  (ni,              p->mesh);
      intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
      dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );

    } else if (old == STD && new == STD) {

      if (np == ni) {  /* Meshes are the same. */
	uniknot  (np, p->mesh);
	dermat_k (np, p->mesh, p->deriv, p->derivT);
      } else {
	uniknot  (np, oldmesh);
	uniknot  (ni, p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
      }

    } else if (old == STD && new == GLL) {

      uniknot  (np,              oldmesh);
      jacgl    (ni-1, 0.0, 0.0,  p->mesh);
      intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
      dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );

    } else

      message (routine, "basis function unrecognized as STD or GLL", ERROR);

    freeDvector (oldmesh, 0);
  }

  /* p now points to valid storage: return requested operators. */

  if (mesh) *mesh = (const double* ) p -> mesh;
  if (in)   *in   = (const double**) p -> interp;
  if (it)   *it   = (const double**) p -> interpT;
  if (dr)   *dr   = (const double**) p -> deriv;
  if (dt)   *dt   = (const double**) p -> derivT;
}


void dIntpOps (const integer basis,  /* element basis: STD or GLL            */
	       const integer np   ,  /* number of knot points                */
	       const double  r    ,  /* location of r in [-1, 1]             */
	       const double  s    ,  /* location of s in [-1, 1]             */
	       double*       inr  ,  /* 1D shape function at r               */
	       double*       ins  ,  /* 1D shape function at s               */
	       double*       dvr  ,  /* 1D shape function derivative at r    */
	       double*       dvs  )  /* 1D shape function derivative at s    */
/* ------------------------------------------------------------------------- *
 * Return the interpolation/derivative vectors for tensor-product shape
 * functions evaluated at canonical coordinates (r, s), each in [-1. 1].
 *
 * Take no action for empty vectors.
 * ------------------------------------------------------------------------- */
{
  const double* kp;
  double        x[1], **dv, **dt;

  dv = dmatrix (0, 0, 0, np - 1);
  dt = dmatrix (0, np - 1, 0, 0);

  dQuadOps (basis, np, np, &kp, 0, 0, 0, 0, 0, 0);

  x[0] = r;
  if (inr) {
    intmat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, inr, 1);
  }
  if (dvr) {
    dermat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, dvr, 1);
  }

  x[0] = s;
  if (ins) {
    intmat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, ins, 1);
  }
  if (dvs) {
    dermat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, dvs, 1);
  }

  freeDmatrix (dv, 0, 0);
  freeDmatrix (dt, 0, 0);
}


typedef struct legcoef {	/* ---- Table for GLL Legendre transform --- */
  integer         np  ;		/* Number of mesh points                     */
  double*         dtab;		/* (np+1)*np table of polynomials & weights  */
  struct legcoef* next;		/* link to next one                          */
} legCoef;			/* ----------------------------------------- */

static legCoef* lChead = 0;


void dglldpc (const integer  np,
	      const double** cd)
/* ------------------------------------------------------------------------- *
 * Return pointers to look-up tables of Legendre polynomials and
 * weights, based on the Gauss--Lobatto nodes.  The tables can be used
 * in calculating discrete Legendre transforms.
 *
 * Tables contain values of Legendre polynomials evaluated at GLL
 * quadrature points, and weights.  Storage is row-major: the first np
 * points contain values of the zeroth Legendre polynomial at the GLL
 * points (these are identically 1.0), the second np points the values
 * of the first polynomial, etc.  So rows index polynomial, while
 * columns index spatial location. 
 *
 * NB: this is the opposite ordering used for dglmdpc, below.
 *
 * The final np points (the np-th row of the table, for 0-based
 * indexing) contains the weights 1/gamma_k (see Canuto et al.,
 * 2.3.13).
 * -------------------------------------------------------------------------
 * */
{
  register integer  found = 0;
  register legCoef* p;

  for (p = lChead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register integer i, j, k;
    const integer    nm = np - 1;
    const double*    z;

    p = (legCoef*) calloc (1, sizeof (legCoef));
    if (lChead) p -> next = lChead;
    lChead = p;

    p -> np   = np;
    p -> dtab = dvector (0, np * (np + 1) - 1);

    dQuadOps (LL, np, np, &z, 0, 0, 0, 0, 0, 0);

    for (i = 0; i < np; i++) {
      k = i + np * np;
      p -> dtab[k] = (i < nm) ?  0.5*(i+i+1) : 0.5*nm;
      for (j = 0; j < np; j++) {
	k = j + i * np;
	p -> dtab[k] = pnleg (z[j], i);
      }
    }
  }

  /* p now points to valid storage: return requested operators. */

  if (cd) *cd = (const double*) p -> dtab;
}


typedef struct legtran {	/* ---- GLL Legendre transform matrices  --- */
  integer         np  ;		/* Number of mesh points                     */
  double*         FW  ;		/* np*np forward Legendre transform matrix.  */
  double*         FT  ;		/* Transpose of FW.                          */
  double*         BW  ;		/* np*np inverse Legendre transform matrix.  */
  double*         BT  ;		/* Transpose of BW.                          */
  double*         UF  ;		/* np^2*np^2 forward transform matrix.       */
  double*         UB  ;		/* np^2*np^2 inverse transform matrix.       */
  struct legtran* next;		/* link to next one                          */
} legTran;			/* ----------------------------------------- */

static legTran* lThead = 0;


void dglldpt (const integer  np,
	      const double** fw,
	      const double** ft,
	      const double** bw,
	      const double** bt,
	      const double** fu,
	      const double** bu)
/* ------------------------------------------------------------------------- *
 * Return pointers to matrices for 1D and 2D Legendre polynomial
 * transforms of np data, based on the Gauss--Lobatto nodes.  See the
 * definitions of forward and inverse polynomial transforms in Canuto
 * et al., Sections 2.2.3, 2.2.13.
 *
 * The 1D discrete transform is defined as (summation implied)
 *
 *   Ak = Ck Wj Uk Lkj,
 *
 * with inverse
 *
 *   Uj = Ak Lkj,
 *
 * where Lij is the ith Legendre polynomial evaluated at the jth
 * quadrature point, Ck correspond to 1/\gamma_k in Canuto et al., and
 * Wj are the quadrature weights for the corresponding points.
 *
 * The equivalent 2D tensor-product (or sum-factorisation) relationships are
 *                                            t
 *   Aij = Ci Cj Wp Lip Wq Upq Ljq = Rip Upq Rqj
 *                                    t
 *   Uij = Lpi Apq Lqj             = Lip Apq Lqj
 *
 * where Rip = Ci Wp Lip.
 *
 * Instead of these tensor-product forms, we can also "unroll" the
 * matrix multiplies above to be a single matrix-vector product in
 * each case, so that
 *
 *   Aij = Fijpq Upq,  Uij = Bijpq Apq.
 *
 * This routine supplies the transform matrices.  If any pointer is
 * passed in as NULL, we don't return the corresponding value.
 * Matrices are supplied 1D, with row-major ordering.
 * ------------------------------------------------------------------------- */
{
  register integer  found = 0;
  register legTran* p;

  for (p = lThead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register integer i, j, k, l, r, s;
    const integer    np2 = np * np;
    const double     *tab, *w;
    double           ci;

    p = (legTran*) calloc (1, sizeof (legTran));
    if (lThead) p -> next = lThead;
    lThead = p;

    p -> np = np;
    p -> FW = dvector (0, 4 * np2 + 2 * np2*np2 - 1);
    p -> FT = p -> FW + np2;
    p -> BW = p -> FT + np2;
    p -> BT = p -> BW + np2;
    p -> UF = p -> BT + np2;
    p -> UB = p -> UF + np2*np2;

    dglldpc  (np, &tab);
    dQuadOps (LL, np, np, 0, 0, &w, 0, 0, 0, 0);

    /* -- Create forward & inverse 1D transform matrices. */

    for (i = 0; i < np; i++) {
      ci = tab[np2 + i];
      for (j = 0; j < np; j++) {
	p->FW[i*np + j] = ci * w[j] * tab[i*np +j];
	p->BW[i*np + j] = tab[j*np+ i];
      }
    }

    /* -- And their transposes. */

    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++) {
	p->FT[i*np + j] = p->FW[j*np + i];
	p->BT[i*np + j] = p->BW[j*np + i];
      }

    /* -- Manufacture 2D forward, inverse DLT matrices. */

    for (k = 0, i = 0; i < np; i++)
      for (j = 0; j < np; j++, k++)
	for (l = 0, r = 0; r < np; r++)
	  for (s = 0; s < np; s++, l++) {
	    p->UF[k*np2 + l] = p->FW[i*np + r] * p->FT[s*np + j];
	    p->UB[k*np2 + l] = p->BW[i*np + r] * p->BT[s*np + j];
	  }
  }

  /* p now points to valid storage: return requested operators. */

  if (fw) *fw = (const double*) p -> FW;
  if (ft) *ft = (const double*) p -> FT;
  if (bw) *bw = (const double*) p -> BW;
  if (bt) *bt = (const double*) p -> BT;
  if (fu) *fu = (const double*) p -> UF;
  if (bu) *bu = (const double*) p -> UB;
}


typedef struct modcoef {	/* -- Table for modal expansion transform -- */
  integer         np  ;		/* Number of mesh points                     */
  double*         dtab;		/* np*np table of basis function values      */
  struct modcoef* next;		/* link to next one                          */
} modCoef;			/* ----------------------------------------- */

static modCoef* mChead = 0;


void dglmdpc (const integer  np,
	      const double** cd)
/* ------------------------------------------------------------------------- *
 * Return pointers to look-up tables of modal expansion functions and
 * weights, based on the Gauss--Lobatto nodes.  The tables can be used
 * in calculating discrete Legendre transforms.
 *
 * Tables contain values of modal expansion functions evaluated at GLL
 * quadrature points, and the associated quadrature weights.  Storage
 * is column-major: the first np points contain values of the all the
 * modal basis function at the first GLL point, the second np points
 * the values of the all the basis functions at the second point, etc.
 * So rows index quadrature point, while columns index polynomial.
 *
 * NB: this ordering is THE OPPOSITE of that used for dglldpc, but it
 * conforms to that used by Karniadakis & Sherwin.
 *
 * The modal expansion functions are a set of hierarchical functions
 * based on the Jacobi polynomials.
 *
 * p_0(z) = 0.5  (1 + z)
 *
 * p_1(z) = 0.5  (1 - z)
 *                                1,1
 * p_n(z) = 0.25 (1 + z) (1 - z) J   (z)   n >= 2.
 *                                n-2
 *
 * where J is a Jacobi polynomial.
 * ------------------------------------------------------------------------- */
{
  register integer  found = 0;
  register modCoef* p;

  for (p = mChead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register integer i, j, k;
    const double     *z;

    p = (modCoef*) calloc (1, sizeof (modCoef));
    if (mChead) p -> next = mChead;
    mChead = p;

    p -> np   = np;
    p -> dtab = dvector (0, np * np - 1);

    dQuadOps (LL, np, np, &z, 0, 0, 0, 0, 0, 0);

    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++) {
	k = i + j * np;
	p -> dtab[k] = pnmod (z[j], i);
      }
  }

  /* -- p now points to valid storage: return requested operators. */

  if (cd) *cd = (const double*) p -> dtab;
}


typedef struct modtran {	/* - GL modal expansion transform matrices - */
  integer         np  ;		/* Number of mesh points                     */
  double*         FW  ;		/* np*np forward modal transform matrix.     */
  double*         FT  ;		/* Transpose of FW.                          */
  double*         BW  ;		/* np*np inverse modal transform matrix.     */
  double*         BT  ;		/* Transpose of BW.                          */
  double*         UF  ;		/* np^2*np^2 forward transform matrix.       */
  double*         UB  ;		/* np^2*np^2 inverse transform matrix.       */
  struct modtran* next;		/* link to next one                          */
} modTran;			/* ----------------------------------------- */

static modTran* mThead = 0;


void dglmdpt (const integer  np,
	      const double** fw,
	      const double** ft,
	      const double** bw,
	      const double** bt,
	      const double** fu,
	      const double** bu)
/* ------------------------------------------------------------------------- *

 * Return pointers to matrices for 1D and 2D modal expansion
 * transforms of np data, based on the Gauss--Lobatto nodes.  The
 * "modal" expansion functions are a set of hierarchical basis
 * functions closely associated with the GLL Lagrange basis. The
 * transform equations are just the standard "Normal form" projections
 * of one basis onto another.  (For a Legendre tranform -- see dglldpt
 * -- the equations are simplified by the discrete orthogonality of
 * all the basis functions.)  Here, not all (but most) of the basis
 * functions are orthogonal.  The treatment and intended use of the
 * transform matrices is similar to that for dglldpt.
 *
 * The 1D discrete forward transform is defined as
 *
 *   ^     t     -1 t                                  ^
 *   u = [B W B ]  B  W u  = F u,  with inverse  u = B u
 *   ~                  ~      ~                 ~     ~
 *
 * where u_i is a vector of function values evaluated at the ith
 * quadrature point, B_ij is the jth basis function evaluated at the
 * ith quadrature point and W_i are the associated quadrature weights.
 * The vector of coefficients of the new basis functions is u^.
 *
 * The equivalent 2D tensor-product (sum-factorisation) relationships are
 *
 *   ^              t                 ^    t
 *   uij = Fip upq Fqj,     uij = Bip upq Bqj
 *
 * Instead of these tensor-product forms, we can also "unroll" the
 * matrix multiplies above to be a single matrix-vector product in
 * each case, so that
 *   ^                                  ^
 *   uij = Fijpq upq,       uij = Bijpq upq.
 *
 * This form, although involving more operations, can actually be
 * faster on vector architectures (depending on np).
 *
 * This routine supplies the transform matrices.  If any pointer is
 * passed in as NULL, we don't return the corresponding value.
 * Matrices are supplied 1D, with row-major ordering.
 *
 * Refs:
 *
 * eq.(2.16), R.D. Henderson, "Adaptive Spectral Element Methods", in
 * "High-Order Methods for Computational Physics", eds T.J. Barth &
 * H. Deconinck, Springer, 1999.
 *
 * eq.(2.40), G.E. Karniadakis & S.J. Sherwin, "Spectral/hp Element
 * Methods for CFD", Oxford, 1999.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "dglmdpt";

  register integer  found = 0;
  register modTran* p;

  for (p = mThead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register integer i, j, k, l, r, s;
    const integer    np2 = np * np;
    const double     *B, *W;
    double           *work, *BtW, *BtWB, *rwrk;
    integer          *iwrk, info;

    p = (modTran*) calloc (1, sizeof (modTran));
    if (mThead) p -> next = mThead;
    mThead = p;

    p -> np = np;
    p -> FW = dvector (0, 4 * np2 + 2 * np2*np2 - 1);
    p -> FT = p -> FW + np2;
    p -> BW = p -> FT + np2;
    p -> BT = p -> BW + np2;
    p -> UF = p -> BT + np2;
    p -> UB = p -> UF + np2*np2;

    dglmdpc  (np, &B);
    dQuadOps (LL, np, np, 0, 0, &W, 0, 0, 0, 0);

    iwrk = ivector (0, np - 1);
    work = dvector (0, 4 * np2 - 1);
    BtW  = work + np2;
    BtWB = BtW  + np2;
    rwrk = BtWB + np2;

    /* -- Create matrices BtW & (symmetric) BtWB. */
    
    for (i = 0; i < np; i++)
      dsmul (np, W[i], B + i * np, 1, BtW + i, np);

    dmxm (BtW, np, (double*) B, np, BtWB, np);

    /* -- Invert BtWB. */

    dgetrf (np, np, BtWB, np, iwrk, &info);
    if (info) message (routine, "matrix BtWB has singular factor", ERROR);

    dgetri (np, BtWB, np, iwrk, rwrk, 2*np2, &info);
    if (info) message (routine, "matrix BtWB is singular", ERROR);

    /* -- Create forward (FW) & inverse (BW) 1D transform matrices. */

    dmxm  (BtWB, np, BtW, np, p->FW, np);
    dcopy (np2, B, 1,         p->BW,  1);

    /* -- And their transposes. */

    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++) {
	p->FT[i*np + j] = p->FW[j*np + i];
	p->BT[i*np + j] = p->BW[j*np + i];
      }

    /* -- Manufacture 2D forward, inverse DPT matrices. */

    for (k = 0, i = 0; i < np; i++)
      for (j = 0; j < np; j++, k++)
	for (l = 0, r = 0; r < np; r++)
	  for (s = 0; s < np; s++, l++) {
	    p->UF[k*np2 + l] = p->FW[i*np + r] * p->FT[s*np + j];
	    p->UB[k*np2 + l] = p->BW[i*np + r] * p->BT[s*np + j];
	  }

    freeIvector (iwrk, 0);
    freeDvector (work, 0);
  }

  /* p now points to valid storage: return requested operators. */

  if (fw) *fw = (const double*) p -> FW;
  if (ft) *ft = (const double*) p -> FT;
  if (bw) *bw = (const double*) p -> BW;
  if (bt) *bt = (const double*) p -> BT;
  if (fu) *fu = (const double*) p -> UF;
  if (bu) *bu = (const double*) p -> UB;
}
