///////////////////////////////////////////////////////////////////////////////
// nonlinear.C: compute nonlinear velocity terms and rates of strain in
// physical space, then transform, so returning the transform of
//
//   -N(u) - div(SGSS) + body force
//
// Copyright (c) 2000 Hugh Blackburn
//
// NB: the past pressure field is destroyed here.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "les.h"

static void transform    (const TransformKind, const ExchangeKind, real*);
static void realGradient (const AuxField*, real*, const int);


void nonLinear (Domain*       D ,
		matrix<real>& Ut,
		vector<real>& ff)
// ---------------------------------------------------------------------------
// Compute nonlinear and SGSS terms, N(u) and \tau_ij.
//
// On entry, D contains the Fourier-transforms of the old velocity
// (and pressure) fields, the lowest levels of Us & Uf are free.
// On exit, the lowest levels of Us contain the Fourier transforms
// of the old velocity fields and lowest levels of Uf contain the 
// Fourier transforms of the corresponding nonlinear terms plus
// body force terms, minus divergence of SGSS.
//
// Nonlinear terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 ~ ~
//           N  = -0.5 ( u . grad u + div uu )
//           ~           ~        ~       ~~
//
// i.e. in Cartesian component form
//
//           N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ),
//            i           j    i      j      i j      j
//
// while for cylindrical coordinates
//
//           Nx = -0.5 {ud(u)/dx + vd(u)/dy +  d(uu)/dx +
//                 1/y [wd(u)/dz + d(uw)/dz + d(yvu)/dy      ]}
//           Ny = -0.5 {ud(v)/dx + vd(v)/dy +  d(uv)/dx +
//                 1/y [wd(v)/dz + d(vw)/dz + d(yvv)/dy - 2ww]}
//           Nz = -0.5 {ud(w)/dx + vd(w)/dy +  d(uw)/dx +
//                 1/y [wd(w)/dz + d(ww)/dz + d(yvw)/dy + 2wv]}.
//
// Data are transformed to physical space for most of the operations,
// but for gradients in the Fourier direction the data must be
// transferred back to Fourier space.
//
// NB: no dealiasing.
// ---------------------------------------------------------------------------
{
  const integer nZ   = Geometry::nZ();
  const integer nZP  = Geometry::nZProc();
  const integer nP   = Geometry::planeSize();
  const integer nPP  = Geometry::nBlock();
  const integer nPR  = Geometry::nProc();
  const integer nTot = Geometry::nTotProc();
  const real    refV = Femlib::value ("KINVIS");
  const real    EPS  = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;
  Field*        meta = D -> u[0]; // -- A handle to use for Field operations.
  integer       i, j, ij;

  if (Geometry::system() == Geometry::Cylindrical)
    message ("nonLinear", "no cylindrical coordinate version yet", ERROR);

  // -- Create names for local computations.

  vector<real*> u (26);
  real**        Sr  = u();
  real**        St  = Sr + 6;
  real**        Sm  = St + 6;
  real**        Ua  = Sm + 2;
  real**        Us  = Ua + 3;
  real**        Ud  = Us + 3;
  real**        Nl  = Ud + 3;
  real*         tmp = D -> udat[3];

  // -- Set pointers into supplied workspace, Ut.

  for (i = 0; i < 6; i++) Sr [i] = Ut (     i);
  for (i = 0; i < 6; i++) St [i] = Ut ( 6 + i);
  for (i = 0; i < 3; i++) Ua [i] = Ut (12 + i);
  for (i = 0; i < 2; i++) Sm [i] = Ut (15 + i);
  for (i = 0; i < 3; i++) Us [i] = Ut (17 + i);
  for (i = 0; i < 3; i++) Nl [i] = Ut (20 + i);

  for (i = 0; i < 3; i++) Ud [i] = D->udat(i);

  // -- Transform and filter velocities.

  for (i = 0; i < 3; i++) {
    Veclib::zero (nTot, Nl[i], 1);
    Veclib::copy (nTot, Ud[i], 1, Us[i], 1);
    Veclib::copy (nTot, Ud[i], 1, Ua[i], 1);
    lowpass      (Ud[i]);
    transform    (INVERSE, FULL, Ud[i]);
    transform    (INVERSE, FULL, Ua[i]);
  }

  // -- Form off-diagonal terms for RoS tensors, nonconservative Nl terms.
  
  for (i = 0; i < 3; i++)
    for (j = i + 1; j < 3; j++) {
      ij = i + j - 1;
      Veclib::copy   (nTot, Us[i], 1, Sr[ij], 1);
      meta->gradient (nZP, nP, Sr[ij], j);
      Veclib::copy   (nTot, Sr[ij], 1, St[ij], 1);
      lowpass        (St[ij]);
      transform      (INVERSE, FULL, Sr[ij]);
      transform      (INVERSE, FULL, St[ij]);
      Veclib::vvtvp  (nTot, Ua[j], 1, Sr[ij], 1, Nl[i], 1, Nl[i], 1);
    }
  
  for (i = 0; i < 3; i++)
    for (j = i + 1; j < 3; j++) {
      Veclib::copy   (nTot, Us[j], 1, Sr[i], 1);
      meta->gradient (nZP, nP, Sr[i], i);
      Veclib::copy   (nTot, Sr[i], 1, St[i], 1);
      lowpass        (St[i]);
      transform      (INVERSE, FULL, St[i]);
      transform      (INVERSE, FULL, Sr[i]);
      Veclib::vvtvp  (nTot, Ua[i], 1, Sr[i], 1, Nl[j], 1, Nl[j], 1);
    }
  
  for (i = 0; i < 3; i++) {
    Veclib::svvtt (nTot, 0.5, Sr[i], 1, Sr[i+3], 1, Sr[i+3], 1);
    Veclib::svvtt (nTot, 0.5, St[i], 1, St[i+3], 1, St[i+3], 1);
  }

  // -- Form diagonal terms for RoS tensors, nonconservative Nl terms.

  for (i = 0; i < 3; i++) {
    Veclib::copy   (nTot, Us[i], 1, Sr[i], 1);
    meta->gradient (nZP, nP, Sr[i], i);
    Veclib::copy   (nTot, Sr[i], 1, St[i], 1);
    lowpass        (St[i]);
    transform      (INVERSE, FULL, St[i]);
    transform      (INVERSE, FULL, Sr[i]);
    Veclib::vvtvp  (nTot, Ua[i], 1, Sr[i], 1, Nl[i], 1, Nl[i], 1);
  }

  // -- Form strain rate magnitude estimates |S|, |S~|.

  Veclib::zero (nTot, Sm[0], 1);
  Veclib::zero (nTot, Sm[1], 1);

  for (i = 0; i < 6; i++) {
    Veclib::vvtvp (nTot, Sr[i], 1, Sr[i], 1, Sm[0], 1, Sm[0], 1);
    Veclib::vvtvp (nTot, St[i], 1, St[i], 1, Sm[1], 1, Sm[1], 1);
  }

  Veclib::vsqrt (nTot, Sm[0], 1, Sm[0], 1);
  Veclib::vsqrt (nTot, Sm[1], 1, Sm[1], 1);

  // -- Delta^2 |S|.

  Veclib::zero (nP, tmp, 1);
  meta -> lengthScale (tmp);
  Veclib::vmul (nP, tmp, 1, tmp, 1, tmp, 1);

  for (i = 0; i < nZP; i++) {
    Veclib::vmul  (nP,      tmp, 1, Sm[0]+i*nP, 1, Sm[0]+i*nP, 1);
    Veclib::svvtt (nP, 4.0, tmp, 1, Sm[1]+i*nP, 1, Sm[1]+i*nP, 1);
  }

  // -- Accumulate LijMij & MijMij terms, conservative Nl terms.

  for (i = 0; i < 6; i++) {
    Veclib::vmul (nTot, Sm[0], 1, Sr[i], 1, Sr[i], 1);
    Veclib::vmul (nTot, Sm[1], 1, St[i], 1, St[i], 1);
  }

  real* L   = Sm[0]; real* M   = Sm[1];
  real* Num = Us[0]; real* Den = Us[1];
  Veclib::zero (nTot, Num, 1);
  Veclib::zero (nTot, Den, 1);

  // -- Diagonal terms.

  for (i = 0; i < 3; i++) {
    Veclib::vmul  (nTot, Ua[i], 1, Ua[i], 1, tmp, 1);
    Veclib::copy  (nTot, tmp, 1, L, 1);
    realGradient  (meta, tmp, i);
    Veclib::vadd  (nTot, tmp, 1, Nl[i], 1, Nl[i], 1);
    transform     (FORWARD, FULL, L);
    lowpass       (L);
    transform     (INVERSE, FULL, L);
    Veclib::vvvtm (nTot, L, 1, Ud[i], 1, Ud[i], 1, L, 1);
    Veclib::copy  (nTot, Sr[i], 1, M, 1);
    transform     (FORWARD, FULL, M);
    lowpass       (M);
    transform     (INVERSE, FULL, M);
    Veclib::vsub  (nTot, St[i], 1, M, 1, M, 1);
    Veclib::vvtvp (nTot, L, 1, M, 1, Num, 1, Num, 1);
    Veclib::vvtvp (nTot, M, 1, M, 1, Den, 1, Den, 1);
  }

  // -- Off-diagonal terms.

  for (i = 0; i < 3; i++)
    for (j = i + 1; j < 3; j++) {
      ij = i + j - 1;
      Veclib::vmul    (nTot, Ua[i], 1, Ua[j], 1, tmp, 1);
      Veclib::copy    (nTot, tmp, 1, L, 1);
      Veclib::copy    (nTot, tmp, 1, M, 1);
      realGradient    (meta, M,   i);
      realGradient    (meta, tmp, j);
      Veclib::vadd    (nTot, M,   1, Nl[j], 1, Nl[j], 1);
      Veclib::vadd    (nTot, tmp, 1, Nl[i], 1, Nl[i], 1);
      transform       (FORWARD, FULL, L);
      lowpass         (L);
      transform       (INVERSE, FULL, L);
      Veclib::vvvtm   (nTot, L, 1, Ud[i], 1, Ud[j], 1, L, 1);
      Veclib::copy    (nTot, Sr[j], 1, M, 1);
      transform       (FORWARD, FULL, M);
      lowpass         (M);
      transform       (INVERSE, FULL, M);
      Veclib::vsub    (nTot, St[j], 1, M, 1, M, 1);
      Veclib::svvttvp (nTot, 2.0, L, 1, M, 1, Num, 1, Num, 1);
      Veclib::svvttvp (nTot, 2.0, M, 1, M, 1, Den, 1, Den, 1);
    }

  // -- Here is the dynamic estimate, L = -2 Cs^2 = LijMij/MijMij.
  //    NB: WE WILL NEED SOME AVERAGING OF Cs HERE.
  
  Veclib::vdiv (nTot, Num, 1, Den, 1, L, 1);

  // -- Subtract off our spatially-constant reference viscosity
  //    (disguised as KINIVS), note factor of 2.

  Veclib::sadd (nTot, 2.0*refV, L, 1, L, 1);

  // -- Create SGSS \tau_ij = -2 (Cs^2 Delta^2 |S| - refVisc) Sij.

  for (i = 0; i < 6; i++)
    Veclib::vmul (nTot, L, 1, Sr[i], 1, Sr[i], 1);

  // -- Normalise skewsymmetric nonlinear terms.

  for (i = 0; i < 3; i++)
    Blas::scal (nTot, -0.5, Nl[i], 1);

#if !defined (NOMODEL)
  // -- Subtract divergence of SGSS from nonlinear terms.

  for (i = 0; i < 3; i++) {	// -- Diagonal terms.
    gradient     (meta, Sr[i], i);
    Veclib::vsub (nTot, Nl[i], 1, Sr[i], 1, Nl[i], 1);
  }

  for (i = 0; i < 3; i++)	// -- Off-diagonal terms.
    for (j = i + 1; j < 3; j++) {
      ij = i + j - 1;
      Veclib::copy (nTot, Sr[j], 1, tmp, 1);
      realGradient (meta, Sr[j], j);
      realGradient (meta, tmp,   i);
      Veclib::vsub (nTot, Nl[i], 1, Sr[j], 1, Nl[i], 1);
      Veclib::vsub (nTot, Nl[j], 1, tmp,   1, Nl[j], 1);
    }
#endif

  // -- Direct stiffness summation.

  for (i = 0; i < 3; i++)
    meta -> smooth (nZP, Nl[i]);

  // -- Fourier transform velocities and nonlinear terms.

  for (i = 0; i < 3; i++) {
    Veclib::copy (nTot, Ua[i], 1, Us[i], 1);
    transform    (FORWARD, FULL,  Us[i]);
    transform    (FORWARD, FULL,  Nl[i]);
  }

  // -- Add on body force terms.

  ROOTONLY
    for (i = 0; i < 3; i++)
      if (fabs (ff[i]) > EPS)
	Veclib::sadd (Geometry::nPlane(), ff[i], Nl[i], 1, Nl[i], 1);
}


static void transform (const TransformKind SIGN,
		       const ExchangeKind  EXCH,
		       real*               data)
// ---------------------------------------------------------------------------
// Fourier transform a block of data in place with half or full
// exchange, according to flags.  No dealiasing.
// ---------------------------------------------------------------------------
{
  const integer nZ  = Geometry::nZ();
  const integer nZP = Geometry::nZProc();
  const integer nP  = Geometry::planeSize();
  const integer nPP = Geometry::nBlock();
  const integer nPR = Geometry::nProc();

  if (SIGN == FORWARD) {
    Femlib::exchange   (data, nZP, nP,  FORWARD);
    Femlib::DFTr       (data, nZ,  nPP, FORWARD);
    if (EXCH == FULL) 
      Femlib::exchange (data, nZP, nP,  INVERSE);
  } else {
    if (EXCH == FULL) 
      Femlib::exchange (data, nZP, nP,  FORWARD);
    Femlib::DFTr       (data, nZ,  nPP, INVERSE);
    Femlib::exchange   (data, nZP, nP,  INVERSE);
  }
}


static void realGradient (const AuxField* meta, 
			  real*           ui  ,
			  const int       xj  )
// ---------------------------------------------------------------------------
// Carry out gradient operation on data, special case for homogeneous
// direction.
// ---------------------------------------------------------------------------
{
  const integer nZ  = Geometry::nZ();
  const integer nZP = Geometry::nZProc();
  const integer nP  = Geometry::planeSize();
  const integer nPP = Geometry::nBlock();
      
  if (xj == 2) {
    transform        (FORWARD, HALF, ui);
    meta -> gradient (nZ, nPP, ui, xj);
    transform        (INVERSE, HALF, ui);
  } else
    meta -> gradient (nZP, nP, ui, xj);
}
