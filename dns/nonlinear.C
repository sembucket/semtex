///////////////////////////////////////////////////////////////////////////////
// nonlinear.C
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// Compute nonlinear terms in Navier--Stokes equations.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "dns.h"


void nonLinear (Domain*       D ,
		AuxField**    Us,
		AuxField**    Uf,
		vector<real>& ff)
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u) + ff.
//
// Here N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS and ff is a vector of body force per unit mass.
//
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms N(u) in skew-symmetric form are
//                 ~ ~
//           N  = -0.5 ( u . grad u + div uu )
//           ~           ~        ~       ~~
//
// i.e., in Cartesian component form
//
//           N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//            i           j    i      j      i j      j
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
// If STOKES is defined for compilation, the nonlinear terms are set to zero.
//
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// Define ALIAS to force Fourier-aliased nonlinear terms (serial only).
// Define SKEW to get skew-symmetric form (default is now non-conservative).
// ---------------------------------------------------------------------------
{
  const int NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int NCOM = D -> nField() - 1;	// -- Number of velocity components.
  int       i, j;

#if defined(STOKES)

  for (i = 0; i < NCOM; i++) {
    *Uf[i] = 0.0;
    AuxField::swapData (D -> u[i], Us[i]);
    ROOTONLY if (fabs (ff[i]) > EPSDP) {
      Veclib::fill (nP, -ff[i], tmp, 1);
      if (i < 2) master -> mulY (1, tmp);
      N[i] -> addToPlane (0, tmp);
    }
  }

#else

#if defined (ALIAS)
  const int         nZ32   = Geometry::nZProc();
#else
  const int         nZ32   = Geometry::nZ32();
#endif 
  const int         nZ     = Geometry::nZ();
  const int         nZP    = Geometry::nZProc();
  const int         nP     = Geometry::planeSize();
  const int         nPP    = Geometry::nBlock();
  const int         nPR    = Geometry::nProc();
  const int         nTot   = Geometry::nTotProc();
  const int         nTot32 = nZ32 * nP;
  vector<real>      work ((2 * NCOM + 1) * nTot32);
  vector<real*>     u32 (NCOM);
  vector<real*>     n32 (NCOM);
  vector<AuxField*> U   (NCOM);
  vector<AuxField*> N   (NCOM);
  Field*            master = D -> u[0];
  real*             tmp    = &work[0] + 2 * NCOM * nTot32;

  Veclib::zero ((2 * NCOM + 1) * nTot32, &work[0], 1); // -- Catch-all cleanup.

  for (i = 0; i < NCOM; i++) {
    u32[i] = &work[0] +  i         * nTot32;
    n32[i] = &work[0] + (i + NCOM) * nTot32;

    AuxField::swapData (D -> u[i], Us[i]);

    U[i] = Us[i];
    U[i] -> transform32 (INVERSE, u32[i]);
    N[i] = Uf[i];
  }

  if (Geometry::cylindrical()) {

    for (i = 0; i < NCOM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

#if defined (SKEW)
      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);
#endif

      if (NCOM == 3) {

#if defined (SKEW)
	if (i == 1)
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);
#else
	if (i == 1)
	  Veclib::svvttvp (nTot32, -1.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  1.0, u32[2], 1, u32[1], 1,      n32[2], 1);
#endif

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);

#if defined (SKEW)
	  Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
#endif
	}
      }

      if (i == 2) master -> divY (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	if (i < 2) master -> mulY (nZ32, tmp);
	Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

#if defined (SKEW)
      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	if (i < 2) master -> mulY (nZ32, tmp);
	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }
#endif

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (FORWARD, n32[i]);
#if 0
      master -> smooth (N[i]);
#endif
#if defined (SKEW)
      ROOTONLY if (fabs (ff[i]) > EPSDP) {
	Veclib::fill (nP, -2.0*ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	N[i] -> addToPlane (0, tmp);
      }
      *N[i] *= -0.5;
#else
      ROOTONLY if (fabs (ff[i]) > EPSDP) {
	Veclib::fill (nP, -ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	N[i] -> addToPlane (0, tmp);
      }
      *N[i] *= -1.0;
#endif

    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NCOM; i++) {
      for (j = 0; j < NDIM; j++) {
      
	// -- Perform n_i += u_j d(u_i) / dx_j.

	Veclib::copy (nTot32, u32[i], 1, tmp,  1);
	if (j == 2) {
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

#if defined (SKEW)
	// -- Perform n_i += d(u_i u_j) / dx_j.

	Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
	if (j == 2) {
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
#endif	
      }

      // -- Transform to Fourier space, smooth, add forcing.
      
      N[i]   -> transform32 (FORWARD, n32[i]);
#if 0
      master -> smooth (N[i]);
#endif
#if defined (SKEW)
      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
#else
      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -ff[i]);
      *N[i] *= -1.0;
#endif
    }
  }

#endif
}
