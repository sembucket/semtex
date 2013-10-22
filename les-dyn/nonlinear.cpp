///////////////////////////////////////////////////////////////////////////////
// nonlinear.C: compute nonlinear velocity terms and rates of strain in
// physical space, then Fourier transform, so returning the transform of
//
//   -N(u) - div(SGSS) + body force
//
// Copyright (c) 2000--2001 Hugh Blackburn
//
// NB: the past pressure field is destroyed here.
//
// In order to return the Smagorinsky length scale instead of the
// dynamic estimate (i.e. so that the Smagorinsky non-dynamic SGSS
// model is used), define SMAG during compilation. The dynamic estimate
// is still computed, but overridden.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "les.h"

static void transform       (const TransformKind, const ExchangeKind, real_t*);
static void realGradient    (const AuxField*, real_t*, const int_t);
static void complexGradient (const AuxField*, real_t*, const int_t);


void nonLinear (Domain*        D ,
		SumIntegrator* S ,
		vector<real_t*>& Ut,
		vector<real_t>&  ff)
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
// Data are transformed to physical space for most of the operations,
// but for gradients in the Fourier direction the data must be
// transferred back to Fourier space.
//
// NB: no dealiasing.
// ---------------------------------------------------------------------------
{
  const int_t nZ   = Geometry::nZ();
  const int_t nZP  = Geometry::nZProc();
  const int_t nP   = Geometry::planeSize();
  const int_t nPP  = Geometry::nBlock();
  const int_t nPR  = Geometry::nProc();
  const int_t nTot = Geometry::nTotProc();
  const int_t CYL  = Geometry::system() == Geometry::Cylindrical;

  Field*      meta = D -> u[0]; // -- A handle to use for Field operations.
  int_t       i, j, ij;

  // -- Create names for local computations.

  vector<real_t*> u (26);
  real_t**        Sr  = &u[0];	// -- Strain rate tensor, unfiltered velocity.
  real_t**        St  = Sr + 6;	// -- Strain rate tensor, filtered velocity.
  real_t**        Ua  = St + 6;	// -- Physical space velocities.
  real_t**        Sm  = Ua + 3;	// -- Strain rate magnitude.
  real_t**        Us  = Sm + 2;	// -- Top level of velocity storage.
  real_t**        Nl  = Us + 3;	// -- Top level of nonlinear/forcing storage.
  real_t**        Ud  = Nl + 3;	// -- Domain velocity field, filtered.
  real_t*         tmp = D -> udat[3]; // -- NB: pressure field gets destroyed.

  // -- Set pointers into supplied workspace, Ut.

  for (i = 0; i < 6; i++) Sr [i] = Ut[     i];
  for (i = 0; i < 6; i++) St [i] = Ut[ 6 + i];
  for (i = 0; i < 3; i++) Ua [i] = Ut[12 + i];
  for (i = 0; i < 2; i++) Sm [i] = Ut[15 + i];
  for (i = 0; i < 3; i++) Us [i] = Ut[17 + i];
  for (i = 0; i < 3; i++) Nl [i] = Ut[20 + i];

  for (i = 0; i < 3; i++) Ud [i] = D->udat [i];

  // -- Zero workspace areas.

  Veclib::zero (6*nTot, Sr[0], 1);
  Veclib::zero (6*nTot, St[0], 1);

  // -- Compute the dynamic mixing length estimate L_mix^2 in Ut[15],
  //    the strain rate magnitude |S| in Ut[16], along with the strain
  //    rate tensor in Ut[0-5].  Also compute the nonlinear terms.

  dynamic (D, Ut);

  real_t* Lmix2 = Sm[0];
  real_t* RoS   = Sm[1];

  const real_t refvis = Femlib::value ("KINVIS");
  const real_t molvis = Femlib::value ("REFVIS");

#if !defined(SMAG)
  // -- Temporal filtering of L_mix^2.

  S -> update (Lmix2);
#endif

  // -- Compose the turbulent eddy viscosity in Sm[0]: nut = L_mix^2 |S|.

  real_t* nut = Sm[0];

  Veclib::vmul(nTot, Lmix2, 1, RoS, 1, nut, 1);

#if 1				// -- Diagnostic printout.
  ROOTONLY
  if (!(D->step % (int) Femlib::value ("IO_CFL")))
    cout << "-- Eddy/molecular viscosity"
	 << ", min: "
	 << nut[Veclib::imin (nTot, nut, 1)]/molvis
	 << ", max: "
	 << nut[Veclib::imax (nTot, nut, 1)]/molvis
	 << endl;
#endif

  // -- Add molecular viscosity to get total viscosity, then clip and smooth.

  Veclib::sadd   (nTot, molvis, nut, 1, nut, 1);
  Veclib::clipup (nTot, 0.0,    nut, 1, nut, 1);
  meta -> smooth (nZP, nut);

  // -- Subtract spatially-constant reference viscosity.

  Veclib::sadd (nTot, -refvis, nut, 1, nut, 1);

  // -- Create SGSS \tau_ij = -2 (L_mix^2 |S| + molvisc - refVisc) Sij.

  for (i = 0; i < 6; i++)
    Veclib::svvtt (nTot, -2.0, nut, 1, Sr[i], 1, Sr[i], 1);

#if !defined (NOMODEL)
  // -- Define NOMODEL, set REFVIS = KINVIS to run as *aliased* dns.

  // -- Subtract divergence of SGSS from nonlinear terms.

  if (CYL) {			// -- Cylindrical geometry.

    // -- Diagonal stress-divergence terms.
      
    realGradient (meta, Sr[0], 0);                      // -- Txx.
    Veclib::vsub (nTot, Nl[0], 1, Sr[0], 1, Nl[0], 1);

    Veclib::copy (nTot, Sr[1], 1, tmp, 1);              // -- Tyy.
    meta -> divY (nZP,  tmp);
    realGradient (meta, Sr[1], 1);
    Veclib::vsub (nTot, Nl[1], 1, tmp,   1, Nl[1], 1);
    Veclib::vsub (nTot, Nl[1], 1, Sr[1], 1, Nl[1], 1);

    Veclib::copy (nTot, Sr[2], 1, tmp, 1);              // -- Tzz.
    meta -> divY (nZP,  tmp);
    Veclib::vadd (nTot, Nl[1], 1, tmp,   1, Nl[1], 1);
    realGradient (meta, Sr[2], 2);
    meta -> divY (nZP,  Sr[2]);
    Veclib::vsub (nTot, Nl[2], 1, Sr[2], 1, Nl[2], 1);

    // -- Off-diagonal terms.

    Veclib::copy (nTot, Sr[3], 1, Sr[0], 1);            // -- Txy, Tyx.
    Veclib::copy (nTot, Sr[3], 1, Sr[1], 1);
    realGradient (meta, Sr[0], 0);
    realGradient (meta, Sr[1], 1);
    meta -> divY (nZP,  Sr[3]);
    Veclib::vsub (nTot, Nl[1], 1, Sr[0], 1, Nl[1], 1);
    Veclib::vsub (nTot, Nl[0], 1, Sr[1], 1, Nl[0], 1);
    Veclib::vsub (nTot, Nl[0], 1, Sr[3], 1, Nl[0], 1);

    Veclib::copy (nTot, Sr[4], 1, tmp, 1);              // -- Txz, Tzx.
    realGradient (meta, Sr[4], 2);
    realGradient (meta, tmp,   0);
    meta -> divY (nZP,  Sr[4]);
    Veclib::vsub (nTot, Nl[0], 1, Sr[4], 1, Nl[0], 1);
    Veclib::vsub (nTot, Nl[2], 1, tmp,   1, Nl[2], 1);

    Veclib::copy (nTot, Sr[5], 1, Sr[0], 1);            // -- Tyz, Tzy.
    Veclib::copy (nTot, Sr[5], 1, Sr[1], 1);
    realGradient (meta, Sr[0], 2);
    meta -> divY (nZP,  Sr[0]);
    realGradient (meta, Sr[1], 1);
    meta -> divY (nZP,  Sr[5]);
    Veclib::vsub (nTot, Nl[1], 1, Sr[0], 1, Nl[1], 1);
    Veclib::vsub (nTot, Nl[2], 1, Sr[1], 1, Nl[2], 1);
    Blas::axpy   (nTot, -2.0,     Sr[5], 1, Nl[2], 1);

  } else {			// -- Cartesian geometry.

    for (i = 0; i < 3; i++) {	// -- Diagonal terms.
      realGradient (meta, Sr[i], i);
      Veclib::vsub (nTot, Nl[i], 1, Sr[i], 1, Nl[i], 1);
    }

    for (i = 0; i < 3; i++)	// -- Off-diagonal terms.
      for (j = i + 1; j < 3; j++) {
	ij = 3 + i + j - 1;
	Veclib::copy (nTot, Sr[ij], 1, tmp, 1);
	realGradient (meta, Sr[ij], j);
	realGradient (meta, tmp,    i);
	Veclib::vsub (nTot, Nl[i], 1, Sr[ij], 1, Nl[i], 1);
	Veclib::vsub (nTot, Nl[j], 1, tmp,    1, Nl[j], 1);
      }
  }

#endif

  // -- Fourier transform velocities and nonlinear terms.

  for (i = 0; i < 3; i++) {
    Veclib::copy (nTot, Ua[i], 1, Us[i], 1);
    transform    (FORWARD, FULL,  Us[i]);
    transform    (FORWARD, FULL,  Nl[i]);
  }

  // -- Add on body force terms.

  ROOTONLY
    for (i = 0; i < 3; i++)
      if (fabs (ff[i]) > EPSDP)
	Veclib::sadd (Geometry::nPlane(), ff[i], Nl[i], 1, Nl[i], 1);

  // -- Direct stiffness summation.

  for (i = 0; i < 3; i++) meta -> smooth (nZP, Nl[i]);
}


void dynamic (Domain*          D ,
	      vector<real_t*>& Ut,
	      const bool       NL)
// ---------------------------------------------------------------------------
// Compute nonlinear and SGSS terms: N(u) and eddy viscosity.
//
// On entry, D contains the Fourier-transforms of the old velocity
// (and pressure) fields, the lowest levels of Us & Uf are free.  On
// exit, the lowest levels of Us the old velocity fields and lowest
// levels of Uf contain the nonlinear terms, both in physical space,
// while other areas hold the mixing length, the rate of strain
// magnitude and components (again, all in physical space).
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
// in cylindrical coordinates
//
//           Nx = -0.5 {ud(u)/dx + vd(u)/dy +  d(uu)/dx + d(vu)/dy +
//                 1/y [wd(u)/dz + d(uw)/dz + vu      ]}
//           Ny = -0.5 {ud(v)/dx + vd(v)/dy +  d(uv)/dx + d(vv)/dy +
//                 1/y [wd(v)/dz + d(vw)/dz + vv - 2ww]}
//           Nz = -0.5 {ud(w)/dx + vd(w)/dy +  d(uw)/dx + d(vw)/dy +
//                 1/y [wd(w)/dz + d(ww)/dz + 3wv     ]}
//
// Data are kept in physical space for most of the operations, but for
// gradients in the Fourier direction the data must be transferred
// back to Fourier space.
//
// If flag NL is set, then we compute the nonlinear terms and return
// also the scaled (dynamic estimate) mixing length
//
//          Cs^2 \Delta^2 = L_mix^2 = -0.5 <LijMij>/<MijMji>
//
// in Ut[15] = Sm[0].  If NL isn't set, do not compute the nonlinear
// terms and return the dynamic estimate of eddy viscosity
//
//          \nu_t = Cs^2 \Delta^2 |S| = L_mix^2 |S|
//
// This coding approach is ugly but means we can use the same code
// both for normal running and for computing the eddy viscosity
// estimate as a stand-alone diagnostic.
// 
// In either case, the unfiltered strain rate tensor terms are also
// returned, in the Ut[0-5], together with the unfiltered strain rate
// magnitude |S| in Ut[19] = Us[2].
//
// NB: no dealiasing occurs, either for serial or parallel
// calculations: they should return the same values.
// ---------------------------------------------------------------------------
{
  const int_t    nZ   = Geometry::nZ();
  const int_t    nZP  = Geometry::nZProc();
  const int_t    nP   = Geometry::planeSize();
  const int_t    nPP  = Geometry::nBlock();
  const int_t    nPR  = Geometry::nProc();
  const int_t    nTot = Geometry::nTotProc();
  const int_t    CYL  = Geometry::system() == Geometry::Cylindrical;

  Field*         meta = D -> u[0]; // -- A handle for Field operations.
  register int_t i, j, ij;

  // -- Create names for local computations.

  vector<real_t*> u (26);
  real_t**        Sr  = &u[0];	// -- Strain rate tensor, unfiltered velocity.
  real_t**        St  = Sr + 6;	// -- Strain rate tensor, filtered velocity.
  real_t**        Ua  = St + 6;	// -- Physical space velocities.
  real_t**        Sm  = Ua + 3;	// -- Strain rate magnitude.
  real_t**        Us  = Sm + 2;	// -- Top level of velocity storage.
  real_t**        Nl  = Us + 3;	// -- Top level of nonlinear/forcing storage.
  real_t**        Ud  = Nl + 3;	// -- Domain velocity field, filtered.
  real_t*         tmp = D -> udat[3];

  // -- Set pointers into supplied workspace, Ut.

  for (i = 0; i < 6; i++) Sr [i] = Ut [     i];
  for (i = 0; i < 6; i++) St [i] = Ut [ 6 + i];
  for (i = 0; i < 3; i++) Ua [i] = Ut [12 + i];
  for (i = 0; i < 2; i++) Sm [i] = Ut [15 + i];
  for (i = 0; i < 3; i++) Us [i] = Ut [17 + i];
  for (i = 0; i < 3; i++) Nl [i] = Ut [20 + i];

  for (i = 0; i < 3; i++) Ud [i] = D->udat [i];

  // -- Transform and filter velocities.

  for (i = 0; i < 3; i++) {
    if (NL) Veclib::zero (nTot, Nl[i], 1);
    Veclib::copy (nTot, Ud[i], 1, Us[i], 1);
    Veclib::copy (nTot, Ud[i], 1, Ua[i], 1);
    lowpass      (Ud[i]);
    transform    (INVERSE, FULL, Ud[i]);
    transform    (INVERSE, FULL, Ua[i]);
  }

  // -- Form off-diagonal terms for RoS tensors, nonconservative Nl terms.
  
  for (i = 0; i < 3; i++)
    for (j = i + 1; j < 3; j++) { // -- Superdiagonal terms.
      ij = 3 + i + j - 1;
      Veclib::copy    (nTot, Us[i], 1, Sr[ij], 1);
      complexGradient (meta, Sr[ij], j);
      Veclib::copy    (nTot, Sr[ij], 1, St[ij], 1);
      lowpass         (St[ij]);
      transform       (INVERSE, FULL, Sr[ij]);
      transform       (INVERSE, FULL, St[ij]);
      if (NL) {
	if (CYL) {
	  switch (ij) { 
	  case 3:
	    Veclib::vmul    (nTot, Ua[1], 1, Sr[3], 1, Nl[0], 1);
	    break;
	  case 4:
	    Veclib::vmul    (nTot, Ua[2], 1, Sr[4], 1, tmp, 1);
	    Veclib::vvtvp   (nTot, Ua[0], 1, Ua[1], 1, tmp, 1, tmp, 1);
	    meta -> divY    (nZP, tmp);
	    Veclib::vadd    (nTot, Nl[0], 1, tmp, 1, Nl[0], 1);
	    break;
	  case 5:
	    Veclib::vmul    (nTot, Ua[2], 1, Sr[5], 1,     Nl[1], 1);
	    Veclib::vvtvp   (nTot, Ua[1], 1, Ua[1], 1,     Nl[1], 1, Nl[1], 1);
	    Veclib::svvttvp (nTot, -2.0, Ua[2],1, Ua[2],1, Nl[1], 1, Nl[1], 1);
	    meta -> divY    (nZP, Nl[1]);
	    break;
	  }
	} else {
	  Veclib::vvtvp (nTot, Ua[j], 1, Sr[ij], 1, Nl[i], 1, Nl[i], 1);
	}
      }
      if (CYL) {
	if (ij == 5) {
	  Veclib::vsub (nTot, Sr[ij], 1, Ua[2], 1, Sr[ij], 1);
	  Veclib::vsub (nTot, St[ij], 1, Ud[2], 1, St[ij], 1);
	}  
	if (ij > 3) {
	  meta -> divY (nZP, Sr[ij]);
	  meta -> divY (nZP, St[ij]);
	}
      }
    }  
  
  for (i = 0; i < 3; i++)
    for (j = i + 1; j < 3; j++) { // -- Subdiagonal terms.
      ij = i + j - 1;
      Veclib::copy    (nTot, Us[j],  1, Sr[ij], 1);
      complexGradient (meta, Sr[ij], i);
      Veclib::copy    (nTot, Sr[ij], 1, St[ij], 1);
      lowpass         (St[ij]);
      transform       (INVERSE, FULL, St[ij]);
      transform       (INVERSE, FULL, Sr[ij]);
      if (NL) Veclib::vvtvp (nTot, Ua[i], 1, Sr[ij], 1, Nl[j], 1, Nl[j], 1);
    }
  
  for (i = 0; i < 3; i++) {
    Veclib::svvpt (nTot, 0.5, Sr[i], 1, Sr[i+3], 1, Sr[i+3], 1);
    Veclib::svvpt (nTot, 0.5, St[i], 1, St[i+3], 1, St[i+3], 1);
  }

  // -- Form diagonal terms for RoS tensors, nonconservative Nl terms.

  for (i = 0; i < 3; i++) {
    Veclib::copy    (nTot, Us[i], 1, Sr[i], 1);
    complexGradient (meta, Sr[i], i);
    Veclib::copy    (nTot, Sr[i], 1, St[i], 1);
    lowpass         (St[i]);
    transform       (INVERSE, FULL, St[i]);
    transform       (INVERSE, FULL, Sr[i]);
    if (NL) {
      if (CYL && i == 2) {
	Veclib::vmul    (nTot,      Ua[2], 1, Sr[2], 1, tmp, 1);
	Veclib::svvttvp (nTot, 3.0, Ua[1], 1, Ua[2], 1, tmp, 1, tmp, 1);
	meta -> divY    (nZP, tmp);
	Veclib::vadd    (nTot,      tmp, 1, Nl[2], 1, Nl[2], 1);
      } else
	Veclib::vvtvp   (nTot, Ua[i], 1, Sr[i], 1, Nl[i], 1, Nl[i], 1);
    }
    if (CYL && i == 2) {
      Veclib::vadd (nTot, Ua[1], 1, Sr[2], 1, Sr[2], 1);
      Veclib::vadd (nTot, Ud[1], 1, St[2], 1, St[2], 1);
      meta -> divY (nZP, Sr[2]);
      meta -> divY (nZP, St[2]);
    }
  }

  // -- Form strain rate magnitude estimates |S|, |S~|: sqrt (2 Sij Sji).

  Veclib::zero  (nTot,           Sm[0], 1);
  Veclib::zero  (nTot,           Sm[1], 1);
  for (i = 3; i < 6; i++) {	// -- Off-diagonal.
    Veclib::vvtvp (nTot, Sr[i], 1, Sr[i], 1, Sm[0], 1, Sm[0], 1);
    Veclib::vvtvp (nTot, St[i], 1, St[i], 1, Sm[1], 1, Sm[1], 1);
  }
  Blas::scal    (nTot, 2.0,      Sm[0], 1);
  Blas::scal    (nTot, 2.0,      Sm[1], 1);
  for (i = 0; i < 3; i++) {	// -- Diagonal.
    Veclib::vvtvp (nTot, Sr[i], 1, Sr[i], 1, Sm[0], 1, Sm[0], 1);
    Veclib::vvtvp (nTot, St[i], 1, St[i], 1, Sm[1], 1, Sm[1], 1);
  }
  Blas::scal    (nTot, 2.0,      Sm[0], 1);
  Blas::scal    (nTot, 2.0,      Sm[1], 1);
  Veclib::vsqrt (nTot, Sm[0], 1, Sm[0], 1);
  Veclib::vsqrt (nTot, Sm[1], 1, Sm[1], 1);

  // -- Accumulate LijMij & MijMij terms, conservative Nl terms.
  //
  //    First, multiply LAMBDA_M^2 |S~| through S~, where the token
  //    LAMBDA_M is the assumed ratio of length scales on the coarse
  //    mesh to those on the fine mesh (e.g. 2).
  //  
  //    And save a copy of |S| (called RoS for Rate of Strain) in Us[2].

  for (i = 0; i < 6; i++)
    Veclib::svvtt
      (nTot, Femlib::value ("LAMBDA_M^2"), Sm[1], 1, St[i], 1, St[i], 1);

  real_t* RoS = Us[2];
  Veclib::copy (nTot, Sm[0], 1, RoS, 1);

  real_t* L   = Sm[0]; real_t* M   = Sm[1];
  real_t* Num = Us[0]; real_t* Den = Us[1];

  Veclib::zero (nTot, Num, 1);
  Veclib::zero (nTot, Den, 1);

  // -- Diagonal terms.

  for (i = 0; i < 3; i++) {
    Veclib::vmul   (nTot, Ua[i], 1, Ua[i], 1, tmp, 1);
    Veclib::copy   (nTot, tmp, 1, L, 1);
    if (NL) {
      realGradient (meta, tmp, i);
      if (CYL && i == 2) meta -> divY (nZP, tmp);
      Veclib::vadd (nTot, tmp, 1, Nl[i], 1, Nl[i], 1);
    }
    transform      (FORWARD, FULL, L);
    lowpass        (L);
    transform      (INVERSE, FULL, L);
    Veclib::vvvtm  (nTot, L, 1, Ud[i], 1, Ud[i], 1, L, 1);

    Veclib::vmul   (nTot, Sr[i], 1, RoS, 1, M, 1);
    transform      (FORWARD, FULL, M);
    lowpass        (M);
    transform      (INVERSE, FULL, M);
    Veclib::vsub   (nTot, St[i], 1, M, 1, M, 1);

    Veclib::vvtvp  (nTot, L, 1, M, 1, Num, 1, Num, 1);
    Veclib::vvtvp  (nTot, M, 1, M, 1, Den, 1, Den, 1);
  }

  // -- Off-diagonal terms.

  for (i = 0; i < 3; i++)
    for (j = i + 1; j < 3; j++) {
      ij = 3 + i + j - 1;
      Veclib::vmul    (nTot, Ua[i], 1, Ua[j], 1, tmp, 1);
      Veclib::copy    (nTot, tmp, 1, L, 1);
      Veclib::copy    (nTot, tmp, 1, M, 1);
      if (NL) {
	realGradient  (meta, M,   i);
	realGradient  (meta, tmp, j);
	if (CYL && j == 2) meta -> divY (nZP, tmp);
	Veclib::vadd  (nTot, M,   1, Nl[j], 1, Nl[j], 1);
	Veclib::vadd  (nTot, tmp, 1, Nl[i], 1, Nl[i], 1);
      }
      transform       (FORWARD, FULL, L);
      lowpass         (L);
      transform       (INVERSE, FULL, L);
      Veclib::vvvtm   (nTot, L, 1, Ud[i], 1, Ud[j], 1, L, 1);

      Veclib::vmul    (nTot, Sr[ij], 1, RoS, 1, M, 1);
      transform       (FORWARD, FULL, M);
      lowpass         (M);
      transform       (INVERSE, FULL, M);
      Veclib::vsub    (nTot, St[ij], 1, M, 1, M, 1);

      Veclib::svvttvp (nTot, 2.0, L, 1, M, 1, Num, 1, Num, 1);
      Veclib::svvttvp (nTot, 2.0, M, 1, M, 1, Den, 1, Den, 1);
    }

  // -- Averaging of Num = <LijMij> and Den = <MijMij>, for stability.

#if 1 // -- Homogeneous average.

  Femlib::exchange (Num, nZP, nP, FORWARD);
  Femlib::exchange (Den, nZP, nP, FORWARD);
  for (i = 1; i < nZ; i++) {
    Veclib::vadd (nPP, Num + i * nPP, 1, Num, 1, Num, 1);
    Veclib::vadd (nPP, Den + i * nPP, 1, Den, 1, Den, 1);
  }
  Blas::scal (nPP, 1.0/nZ, Num, 1);
  Blas::scal (nPP, 1.0/nZ, Den, 1);
  for (i = 1; i < nZ; i++) {
    Veclib::copy (nPP, Num, 1, Num + i * nPP, 1);
    Veclib::copy (nPP, Den, 1, Den + i * nPP, 1);
  }
  Femlib::exchange (Num, nZP, nP, INVERSE);
  Femlib::exchange (Den, nZP, nP, INVERSE);
#endif

#if 0 // -- Elemental average.

  const int_t nEl = Geometry::nElmt();
  const int_t nP2 = Geometry::nTotElmt();

  for (i = 0; i < nEl; i++) {
    Num[ij] = Veclib::sum (nP2, Num + i * nP2, 1) / nP2;
    Den[ij] = Veclib::sum (nP2, Den + i * nP2, 1) / nP2;
    Veclib::fill (nP2, Num[ij], Num + i * nP2, 1);
    Veclib::fill (nP2, Den[ij], Den + i * nP2, 1);
  }
  for (i = 1; i < nZP; i++) {
    Veclib::copy (nP, Num, 1, Num + i * nP, 1);
    Veclib::copy (nP, Den, 1, Den + i * nP, 1);
  }
#endif

  // -- Nudge denominator in case it's very small.

  Veclib::sadd (nTot, EPSDP, Den, 1, Den, 1);

  // -- Here is the dynamic estimate, L = L_mix^2 = -0.5 <LijMij>/<MijMij>.

  Veclib::vdiv (nTot, Num, 1, Den, 1, L, 1);
  Blas::scal   (nTot, -0.5, L, 1);

#if defined (SMAG)		// -- For debugging.

  // -- Substitute a length scale based on the Smagorinsky constant.
  //    L = L_mix^2 = C_s^2 \Delta^2.

  meta -> lengthScale (tmp);
  Blas::scal   (nP, Femlib::value ("C_SMAG"), tmp, 1);
  Veclib::vmul (nP, tmp, 1, tmp, 1, tmp, 1);
  
  for (i = 0; i < nZP; i++)
    Veclib::copy (nP, tmp, 1, L + i * nP, 1);

#endif

  if (NL)
    // -- Normalise skewsymmetric nonlinear terms.

    for (i = 0; i < 3; i++) Blas::scal (nTot, -0.5, Nl[i], 1);

  else {
    // -- We are calculating eddy viscosity from L, make sure also we
    //    put transformed velocities back in D.  L is in physical space.
    //    \nu_t = L_mix^2 |S|.

    Veclib::vmul (nTot, RoS, 1, L, 1, L, 1);

    for (i = 0; i < 3; i++) {
      meta -> smooth (nZP, L);
      Veclib::copy (nTot, Ua[i], 1, Us[i], 1);
      transform    (FORWARD, FULL,  Us[i]);
    }
  }

  // -- Copy |S| to Sm[1], so we have L_mix^2 in Sm[0], |S| in Sm[1].

  Veclib::copy (nTot, RoS, 1, Sm[1], 1);
}


static void transform (const TransformKind SIGN,
		       const ExchangeKind  EXCH,
		       real_t*             data)
// ---------------------------------------------------------------------------
// Fourier transform a block of data in place with half or full
// exchange, according to flags.  No dealiasing.
// ---------------------------------------------------------------------------
{
  const int_t nZ  = Geometry::nZ();
  const int_t nZP = Geometry::nZProc();
  const int_t nP  = Geometry::planeSize();
  const int_t nPP = Geometry::nBlock();

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
			  real_t*         ui  ,
			  const int_t     xj  )
// ---------------------------------------------------------------------------
// Carry out gradient operation on data, special case for homogeneous
// direction.
// ---------------------------------------------------------------------------
{
  const int_t nZ  = Geometry::nZ();
  const int_t nZP = Geometry::nZProc();
  const int_t nP  = Geometry::planeSize();
  const int_t nPP = Geometry::nBlock();
      
  if (xj == 2) {
    transform        (FORWARD, HALF, ui);
    meta -> gradient (nZ, nPP, ui, xj);
    transform        (INVERSE, HALF, ui);
  } else
    meta -> gradient (nZP, nP, ui, xj);
}


static void complexGradient (const AuxField* meta, 
			     real_t*         ui  ,
			     const int_t     xj  )
// ---------------------------------------------------------------------------
// Carry out gradient operation on data, special case for homogeneous
// direction.
// ---------------------------------------------------------------------------
{
  const int_t nZ  = Geometry::nZ();
  const int_t nZP = Geometry::nZProc();
  const int_t nP  = Geometry::planeSize();
  const int_t nPP = Geometry::nBlock();
      
  if (xj == 2) {
    Femlib::exchange (ui, nZP, nP, FORWARD);
    meta -> gradient (nZ, nPP, ui, xj);
    Femlib::exchange (ui, nZP, nP, INVERSE);
  } else
    meta -> gradient (nZP, nP, ui, xj);
}
