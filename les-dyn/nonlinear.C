///////////////////////////////////////////////////////////////////////////////
// nonlinear.C: compute nonlinear type terms in physical space.
//
// Copyright (c) 2000 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "les.h"


void nonLinear (Domain*       D ,
		AuxField***   Us,
		AuxField***   Uf,
		matrix<real>& Ut)
// ---------------------------------------------------------------------------
// Compute nonlinear + forcing terms in Navier--Stokes equations
//                       N(u) + div(2*EV*S) + ff.
//
// Here N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS, div(2*EV*S) is the divergence of the non-constant
// component of the eddy-viscosity SGS terms and ff is a vector of body
// force per unit mass.
//
// On entry, D contains the Fourier-transforms of the old velocity
// (and pressure) fields, the lowest levels of Us & Uf are free.  On
// exit, the velocity field storage areas of D are free, the zeroth
// level of Us contains the old velocities and the zeroth level of Uf
// contains the most recent explicit forcing terms.  Velocity field
// data areas of D and first level of Us are swapped, then the next
// stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
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
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// NB: dealiasing is not used for multiprocessor operation.
// ---------------------------------------------------------------------------
{
  integer           i, j;
  const real        EPS    = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;
  const integer     nZ     = Geometry::nZ();
  const integer     nZP    = Geometry::nZProc();
  const integer     nP     = Geometry::planeSize();
  const integer     nPP    = Geometry::nBlock();
  const integer     nPR    = Geometry::nProc();
  const integer     nTot   = nZP * nP;
  const integer     nZ32   = (nPR > 1) ? nZP : (3 * nZ) >> 1;
  const integer     nTot32 = nZ32 * nP;
  vector<real>      work ((2 * DIM + 1) * nTot32);
  vector<real*>     u32 (DIM);
  vector<real*>     n32 (DIM);
  vector<AuxField*> U   (DIM);
  vector<AuxField*> N   (DIM);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * DIM * nTot32;

  vector<real*>     Sr (6), St(6), Ua(3), Sm(2);

  // -- Set pointers into supplied workspace.

  for (i = 0; i < 6; i++) { Sr (i) = Ut (i);      St (i) = Ut (6+i); }
  for (i = 0; i < 3; i++)   Ua (i) = Ut (12 + i);
  for (i = 0; i < 2; i++)   Sm (i) = Ut (15 + i);

  // -- Transform and filter velocities.

  for (i = 0; i < 3; i++) {
    *Us[i][0] = *D -> u(i);
    D -> u(i) -> transform32 (-1, Ua(i));
    lowpass (D -> udat(i));
    D -> u(i) -> transform   (-1);
  }

  for (i = 0; i < DIM; i++) {
    u32[i] = work() +  i        * nTot32;
    n32[i] = work() + (i + DIM) * nTot32;
  }

  // -- Start with contribution from divergence of SGS.

  EV -> transform32 (-1, tmp);
  Blas::scal (nTot32, -4.0, tmp, 1); // -- 2 = -0.5 * -4 ... see below.

  if (CYL) {			// -- Cylindrical coordinates.

    // -- 2D stress-divergence terms.

    Us[0][0] -> transform32 (-1, n32[0]);
    Veclib::vmul            (nTot32,  tmp, 1, n32[0], 1, n32[0], 1);
    master   -> gradient    (nZ32, nP, n32[0], 0);

    Us[1][0] -> transform32 (-1, n32[1]);
    Veclib::vmul            (nTot32,  tmp, 1, n32[1], 1, n32[1], 1);
    master   -> mulR        (nZ32, n32[1]);
    master   -> gradient    (nZ32, nP, n32[1], 1);
    master   -> divR        (nZ32, n32[1]);

    Uf[0][0] -> transform32 (-1, u32[0]);
    Veclib::vmul            (nTot32,  tmp, 1, u32[0], 1, u32[0], 1);
    Veclib::copy            (nTot32,          u32[0], 1, u32[1], 1);
    master   -> mulR        (nZ32, u32[0]);
    master   -> gradient    (nZ32, nP, u32[0], 1);
    master   -> divR        (nZ32, u32[0]);
    master   -> gradient    (nZ32, nP, u32[1], 0);
    Veclib::vadd            (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
    Veclib::vadd            (nTot32, u32[1], 1, n32[1], 1, n32[1], 1);

    // -- 3D stress-divergence terms.

    if (C3D) {		

      Us[2][0] -> transform32 (-1, n32[2]);
      Veclib::vmul            (nTot32, tmp, 1, n32[2], 1, n32[2], 1);
      Veclib::copy            (nTot32,         n32[2], 1, u32[0], 1);
      master -> divR          (nZ32, u32[0]);
      Veclib::vsub            (nTot32, n32[1], 1, u32[0], 1, n32[1], 1);
      Femlib::exchange        (n32[2], nZ32,        nP, +1);
      Femlib::DFTr            (n32[2], nZ32 * nPR, nPP, +1);
      Veclib::zero            (nTot32 - nTot, n32[2] + nTot, 1);
      master -> gradient      (nZ, nPP, n32[2], 2);
      Femlib::DFTr            (n32[2], nZ32 * nPR, nPP, -1);
      Femlib::exchange        (n32[2], nZ32,        nP, -1);
      master -> divR          (nZ32, n32[2]);

      Uf[1][0] -> transform32 (-1, u32[0]);
      Veclib::vmul            (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy            (nTot32,         u32[0], 1, u32[1], 1);
      Femlib::exchange        (u32[0], nZ32,        nP, +1);
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, +1);
      Veclib::zero            (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient      (nZ, nPP, u32[0], 2);
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, -1);
      Femlib::exchange        (u32[0], nZ32,        nP, -1);
      master -> divR          (nZ32, u32[0]);
      master -> gradient      (nZ32, nP, u32[1], 0);
      Veclib::vadd            (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
      Veclib::vadd            (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);

      Uf[2][0] -> transform32 (-1, u32[0]);
      Veclib::vmul            (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy            (nTot32,         u32[0], 1, u32[1], 1);
      Veclib::copy            (nTot32,         u32[0], 1, u32[2], 1);
      Femlib::exchange        (u32[0], nZ32,        nP, +1);      
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, +1);
      Veclib::zero            (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient      (nZ, nPP, u32[0], 2);
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, -1);
      Femlib::exchange        (u32[0], nZ32,        nP, -1);      
      master -> divR          (nZ32, u32[0]);
      master -> gradient      (nZ32, nP, u32[1], 1);
      master -> divR          (nZ32, u32[2]);
      Veclib::vadd            (nTot32, u32[0], 1, n32[1], 1, n32[1], 1);
      Veclib::vadd            (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);
      Blas::axpy              (nTot32, 2.0, u32[2], 1,       n32[2], 1);
    }

  } else {			// -- Cartesian coordinates.

    // -- Diagonal stress-divergence terms.

    for (i = 0; i < DIM; i++) {

      Us[i][0] -> transform32 (-1, n32[i]);
      Veclib::vmul (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

      if (i == 2) {
	Femlib::exchange   (n32[2], nZ32,        nP, +1);      
	Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, +1);
	Veclib::zero       (nTot32 - nTot, n32[2] + nTot, 1);
	master -> gradient (nZ, nPP, n32[2], 2);
	Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, -1);
	Femlib::exchange   (n32[2], nZ32,        nP, -1);      
      } else
	master -> gradient (nZ32, nP, n32[i], i);
    }

    // -- Off-diagonal stress-divergence terms.

    for (i = 0; i < DIM; i++)
      for (j = i + 1; j < DIM; j++) {

	Uf[i + j - 1][0] -> transform32 (-1, u32[0]);
	Veclib::vmul                    (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
	Veclib::copy                    (nTot32,         u32[0], 1, u32[1], 1);

	// -- Super-diagonal.

	if (j == 2) {
	  Femlib::exchange   (u32[0], nZ32,        nP, +1);
	  Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, u32[0] + nTot, 1);
	  master -> gradient (nZ, nPP, u32[0], 2);
	  Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, -1);
	  Femlib::exchange   (u32[0], nZ32,        nP, -1);
	} else
	  master -> gradient (nZ32, nP, u32[0], j);
	Veclib::vadd (nTot32, u32[0], 1, n32[i], 1, n32[i], 1);

	// -- Sub-diagonal.

	master -> gradient (nZ32, nP, u32[1], i);
	Veclib::vadd       (nTot32, u32[1], 1, n32[j], 1, n32[j], 1);
      }
  }

  // -- Nonlinear terms.

  for (i = 0; i < DIM; i++) {
    AuxField::swapData (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    N[i] = Uf[i][0];
    U[i] -> transform32 (-1, u32[i]);
    if (CYL) N[i] -> mulR (nZ32, n32[i]);
  }
  
  if (CYL) {			// -- Cylindrical coordinates.

    for (i = 0; i < DIM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (DIM == 3) {
	if      (i == 1)	// -- Centripetal.
	  Veclib::svvttvp (nTot32, -2.0,u32[2],1,u32[2],1,n32[1],1,n32[1],1);
	else if (i == 2)	// -- Coriolis.
	  Veclib::svvttvp (nTot32,  2.0,u32[2],1,u32[1],1,n32[2],1,n32[2],1);

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::exchange   (tmp, nZ32,        nP, -1);
	  Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);
	
	  Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::exchange   (tmp, nZ32,        nP, -1);
	  Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
	}
      }

      // -- Terms made with product of radius.
      
      Veclib::vmul       (nTot32, u32[1], 1, u32[i], 1, tmp,  1);
      master -> mulR     (nZ32, tmp);
      master -> gradient (nZ32, nP, tmp, 1);
      Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      master -> divR     (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vvtvp      (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Remaining conservative derivative.
      
      Veclib::vmul       (nTot32, u32[0], 1, u32[i], 1, tmp, 1);
      master -> gradient (nZ32, nP, tmp, 0);
      Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (+1, n32[i]);
      master -> smooth (N[i]);
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
      
	// -- Nonconservative contribution,  n_i += u_j d(u_i) / dx_j.

	Veclib::copy (nTot32, u32[i], 1, tmp,  1);
	if (j == 2) {
	  Femlib::exchange   (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::exchange   (tmp, nZ32,        nP, -1);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

	// -- Conservative contribution,  n_i += d(u_i u_j) / dx_j.

	Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
	if (j == 2) {
	  Femlib::exchange   (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::exchange   (tmp, nZ32,        nP, -1);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Transform to Fourier space, smooth, add forcing.
      
      N[i]   -> transform32 (+1, n32[i]);
      master -> smooth (N[i]);
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
    }
  }
}
