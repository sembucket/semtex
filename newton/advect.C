///////////////////////////////////////////////////////////////////////////////
// advect.C: different forms of advection terms for Navier--Stokes problems
// in primitive variables.
//
// Copyright (c) 2002 <--> $Date$, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "newt.h"

void nonlinear (Domain*         D ,
		AuxField**      Us,
		AuxField**      Uf,
		vector<real_t>& ff)
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u),
// where N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS.
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
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// Define ALIAS to force Fourier-aliased nonlinear terms (serial only).
// Define CONV  to get convective as opposed to default skew-symmetric form.
// ---------------------------------------------------------------------------
{
  const int NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int NCOM = D -> nField() - 1;	// -- Number of velocity components.
  int       i, j;
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

#if !defined (CONV)
      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);
#endif

      if (NCOM == 3) {

#if !defined (CONV)
	if (i == 1)
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);
#else
	if (i == 1)
	  Veclib::svvttvp (nTot32, -1.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::vmul    (nTot32, u32[2], 1, u32[1], 1, n32[2], 1);
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

#if !defined (CONV)
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
	if (i <  2) master -> mulY (nZ32, tmp);
	Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

#if !defined (CONV)
      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	if (i <  2) master -> mulY (nZ32, tmp);
	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }
#endif

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (FORWARD, n32[i]);

#if !defined (CONV)
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

#if !defined (CONV)
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

#if !defined (CONV)
      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
#else
      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -ff[i]);
      *N[i] *= -1.0;
#endif
    }
  }
}


void linear (Domain*         D ,
	     AuxField**      Us,
	     AuxField**      Uf,
	     vector<real_t>& ff)
// ---------------------------------------------------------------------------
// Compute linearised (forcing) terms in Navier--Stokes equations: L(u).
//
// Velocity field data areas of D and first level of Us are swapped,
// then the next stage of linear forcing terms D(u) are computed from
// velocity fields and left in the first level of Uf.
//
// Linearised terms L(u) are computed in convective form
//                 
//   L  = - ( U.grad u + u.grad U )
//   ~        ~      ~   ~      ~
//
// where U (the "base flow") and u (the "perturbation") are separate
// storage in D.  They both have the same spatial structure.
//
// In Cartesian coordinates:
//
//   L  = - ( U  d(u )/dx  +  u   d(U )/dx)
//    i       j    i    j      i    j   j
//
// in cylindrical coordinates
//
//  -Lx = Udu/dx + udU/dx + Vdu/dy + vdU/dy + 1/y [Wdu/dz + wdU/dz]
//  -Ly = Udv/dx + udV/dx + Vdv/dy + vdV/dy + 1/y [Wdv/dz + wdV/dz - 2wW]
//  -Lz = Udw/dx + udW/dx + Vdw/dy + vdW/dy + 1/y [Wdw/dz + wdW/dz + vW + Vw]
//
// ---------------------------------------------------------------------------
{ 
  const int_t       NCOM   = D -> nField() - 1;
  const int_t       NDIM   = Geometry::nDim();
  const int_t       nZ     = Geometry::nZ();
  const int_t       nZP    = Geometry::nZProc();
  const int_t       nP     = Geometry::planeSize();
  const int_t       nPP    = Geometry::nBlock();
  const int_t       nPR    = Geometry::nProc();
  const int_t       nTot   = Geometry::nTotProc();
#if defined (ALIAS)
  const int_t       nZ32   = Geometry::nZProc();
#else
  const int_t       nZ32   = Geometry::nZ32();
#endif 
  const int_t       nTot32 = nZ32 * nP;
  int_t             i, j;
  vector<real_t>    work((3 * NCOM + 1) * nTot32);
  vector<real_t*>   u32(NCOM), U32(NCOM), L32(NCOM);
  vector<AuxField*> u  (NCOM), U  (NCOM), L  (NCOM);
  Field*            master = D -> u[0];
  real_t*           tmp    = &work[0] + 3 * NCOM * nTot32;

  Veclib::zero ((3 * NCOM + 1) * nTot32, &work[0], 1);

  for (i = 0; i < NCOM; i++) {
    u32[i] = &work[0] +  i                * nTot32;
    U32[i] = &work[0] + (i + NCOM)        * nTot32;
    L32[i] = &work[0] + (i + NCOM + NCOM) * nTot32;

    AuxField::swapData (D -> u[i], Us[i]);

    u[i] = Us[i];
    u[i] -> transform32 (INVERSE, u32[i]);
    U[i] = D -> U[i];
    U[i] -> transform32 (INVERSE, U32[i]);
    L[i] = Uf[i];
  }

  if (Geometry::cylindrical()) {

    for (i = 0; i < NCOM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (NCOM == 3) {
	if (i == 1) {
	  Veclib::svvttvp (nTot32, 2., u32[2],1, U32[2],1, L32[1],1, L32[1],1);
	} else if (i == 2) {
	  Veclib::vvvtm   (nTot32, L32[2],1, u32[1],1, U32[2],1, L32[2],1);
	  Veclib::vvvtm   (nTot32, L32[2],1, U32[1],1, u32[2],1, L32[2],1);
	}

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vvvtm      (nTot32, L32[i], 1, U32[2], 1, tmp, 1, L32[i], 1);

	  Veclib::copy       (nTot32, U32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vvvtm      (nTot32, L32[i], 1, u32[2], 1, tmp, 1, L32[i], 1);
	}
      }

      if (i == 2) master -> divY (nZ32, L32[i]);

      // -- 2D convective derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, U32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	if      (i  < 2 && j  < 2) master -> mulY (nZ32, tmp);
	else if (i == 2 && j == 2) master -> divY (nZ32, tmp);
	Veclib::vvvtm      (nTot32, L32[i], 1, u32[j], 1, tmp, 1, L32[i], 1);

	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	if ( i < 2) master -> mulY (nZ32, tmp);
	Veclib::vvvtm      (nTot32, L32[i], 1, U32[j], 1, tmp, 1, L32[i], 1);
      }

    }
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NCOM; i++)
      for (j = 0; j < NDIM; j++) {

	// -- Perform L_i -= U_j d(u_i) / dx_j.

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
	Veclib::vvvtm (nTot32, L32[i], 1, U32[j], 1, tmp, 1, L32[i], 1);

	// -- Perform L_i -= u_j d(U_i) / dx_j.

	Veclib::copy (nTot32, U32[i], 1, tmp,  1);
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
	Veclib::vvvtm (nTot32, L32[i], 1, u32[j], 1, tmp, 1, L32[i], 1);
      }
  }
  
  for (i = 0; i < NCOM; i++) {
    L[i] -> transform32 (FORWARD, L32[i]);
    if (Geometry::cylindrical()) {
      ROOTONLY if (fabs (ff[i]) > EPSDP) {
	Veclib::fill (nP, ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	L[i] -> addToPlane (0, tmp);
      }
    } else
      ROOTONLY if (fabs (ff[i]) > EPSDP) L[i] -> addToPlane (0, ff[i]);
  }
}
