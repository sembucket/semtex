///////////////////////////////////////////////////////////////////////////////
// pseudo.C: (2D) Pseudospectral multiplication with dealiasing.
//
// Copyright (c) 2003 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "dns.h"


void daprod (const int   nz,
	     const real* u ,
	     const real* v ,
	     real*       w )
// ---------------------------------------------------------------------------
// Compute w = u * v, but with dealiasing on in-plane products. Inputs
// u & v are assumed to have usual AuxField structure.
//
// Given a vectors {u}_{N x 1}, {v}_{N x 1}, we first project these
// onto the longer vectors {u^+}_{3N/2 x 1}, {v^+}_{3N/2 x 1} by
// interpolation of the Lagrange interpolants, via matrix
// multiplication.
//
//   {u^+}_{3N/2 x 1} = I^+_{3N/2 x N} {u}_{N x 1}
//   {v^+}_{3N/2 x 1} = I^+_{3N/2 x N} {v}_{N x 1}
//
// Then we form the product
//
//   {w^+} = {u^+} * {w^+} where * represents vector multiply.
//
// This product is not polluted by aliasing; the task is to project
// back to the smaller space only those polynomial coefficients which
// are appropriate to the space of {N-1}-degree polynomials. This is
// done by projecting to a space of orthogonal polynomials, filtering
// and projecting back. This is done via the discrete polynomial
// transform, with a (diagonal) filter matrix L, which has 1 as the
// diagonal entry for the first N rows, 0 in all other locations.
//
//   {w^+_f} = B L B^{-1} {w^+},
//
// where all matrices are 3N/2 x 3N/2, and the vectors are length
// 3N/2.  The Bs are matrices of basis functions (Legendre
// polynomials, given the knots and weights we are using).
//
// Finally we project back to the space of GLL Lagrange polynomials:
//
//   {w}_{N x 1} = I^-_{N x 3N/2} {w^+_f}_{3N/2 x 1}
//
// The matrices I^- B L B^{-1} may be premultiplied.
//
// The above description is for 1D products. Here we extend the
// treatment to 2D, using tensor-product forms.
// ---------------------------------------------------------------------------
{
  const int         N   = Geometry::nP();
  const int         N2  = sqr (N);
  const int         Np  = (3*N + 1) / 2;
  const int         Np2 = sqr(Np);
  const int         nP  = Geometry::planeSize();
  const int         nel = Geometry::nElmt();
  int               i, j;
  static real       *Fn, *Ft;
  static const real **In, **It;
  
  if (!Fn) {			// -- Set up matrices that do the work.

    Fn = new real [Np * N];
    Ft = new real [N * Np];

    const real   **Im, *B, *Bm;
    vector<real> L(Np), F1 (Np2), F2 (Np2);

    // -- Retrieve the interpolation and basis function matrices.

    Femlib::mesh    (GLL, GLL, N, Np, 0, &In, &It, 0, 0);
    Femlib::mesh    (GLL, GLL, Np, N, 0, &Im, 0,   0, 0);
    Femlib::legTran (Np, &B, 0, &Bm, 0, 0, 0);

    // -- Create the filter vector.
    
    Veclib::fill (N, 1.0, &L[0], 1);
    Veclib::zero (Np-N,   &L[N], 1);

    // -- Create the (non-transpose) matrix to do the filtering
    //    and back-projection, Fn.

    for (i = 0; i < Np; i++)
      Veclib::smul (Np, L[i], B + i*Np, 1, &F1[0] + i*Np, 1);
    
    Blas::mxm (Bm,   Np, &F1[0], Np, &F2[0], Np);
    Blas::mxm (Im[0], N, &F2[0], Np, &Fn[0], Np);

    // -- And its transpose, Ft.

    for (i = 0; i < N; i++)
      for (j = 0; j < Np; j++)
	Ft[Veclib::row_major(j, i, N)] = Fn[Veclib::row_major (i, j, Np)];
  }

#if 1
  // -- Here the matrices are applied on an element-by-element basis.

  const real        *uplane, *vplane, *uel, *vel;
  real              *up, *vp, *wp, *tp, *wplane, *wel;
  static real       *work = new real [4 * Np2];

  up = &work[0];
  vp = up + Np2;
  wp = vp + Np2;
  tp = wp + Np2;

  for (i = 0; i < nz; i++) {

    uplane = u + i*nP;
    vplane = v + i*nP;
    wplane = w + i*nP;

    for (j = 0; j < nel; j++) {

      uel = uplane + j*N2;
      vel = vplane + j*N2;
      wel = wplane + j*N2;
      
      // -- Interpolate u & v.

      Blas::mxm (In[0], Np, uel,   N, tp, N );
      Blas::mxm (tp,    Np, It[0], N, up, Np);
      Blas::mxm (In[0], Np, vel,   N, tp, N );
      Blas::mxm (tp,    Np, It[0], N, vp, Np);

      // -- Form product wp.

      Veclib::vmul (Np2, up, 1, vp, 1, wp, 1);

      // -- Filter/interpolate back.

      Blas::mxm (Fn, N, wp, Np, tp, Np);
      Blas::mxm (tp, N, Ft, Np, wel, N);
    }
  }
#else
  // -- Here on a plane-at-once basis.

  const int   npad  = Np2*nel;
  static real *work = new real [4*npad + nP];
  const real  *uplane, *vplane;
  real        *uplus, *vplus, *wplus, *wplane, *temp;

  uplus = work;
  vplus = uplus + npad;
  wplus = vplus + npad;
  temp  = wplus + npad;
  

  for (i = 0; i < nz; i++) {
    uplane = u + i*nP;
    vplane = v + i*nP;
    wplane = w + i*nP;

    Femlib::tpr2d (uplane, uplus, temp, In[0], It[0], Np, N, nel);
    Femlib::tpr2d (vplane, vplus, temp, In[0], It[0], Np, N, nel);
    Veclib::vmul  (Np2*nel, uplus, 1, vplus, 1, wplus, 1);
    Femlib::tpr2d (wplus, wplane, temp, Fn, Ft, N, Np, nel);
  }

#endif
}














