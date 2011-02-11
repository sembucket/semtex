///////////////////////////////////////////////////////////////////////////////
// nonlinear.C
//
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u) + ff.
//
// Here N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS and ff is a vector of body force per unit mass.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <dns.h>


void skewSymmetric (Domain*         D ,
		    AuxField**      Us,
		    AuxField**      Uf,
		    vector<real_t>& ff)
// ---------------------------------------------------------------------------
// Velocity field data areas of D (which on entry contain velocity
// data from the previous timestep) and first level of Us are swapped
// (so that subsequently Us stores the old velocity data), then the
// next stage of nonlinear forcing terms N(u) are computed from
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
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in reference[3].
//
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// Define ALIAS to force Fourier-aliased nonlinear terms (serial only
// -- these are never dealiased in parallel execution).
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int_t NCOM = D -> nField() - 1;	// -- Number of velocity components.
  const int_t nP   = Geometry::planeSize();

#if defined (ALIAS)
  const int_t       nZ32   = Geometry::nZProc();
#else
  const int_t       nZ32   = Geometry::nZ32();
#endif 
  const int_t       nZ     = Geometry::nZ();
  const int_t       nZP    = Geometry::nZProc();

  const int_t       nPP    = Geometry::nBlock();
  const int_t       nPR    = Geometry::nProc();
  const int_t       nTot   = Geometry::nTotProc();
  const int_t       nTot32 = nZ32 * nP;

  vector<real_t*>   u32 (NCOM), n32 (NCOM);
  vector<AuxField*> U   (NCOM), N   (NCOM);
  Field*            master = D -> u[0];

  vector<real_t> work ((2 * NCOM + 1) * nTot32);
  real_t*        tmp  = &work[0] + 2 * NCOM * nTot32;
  int_t          i, j;

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

      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);

      if (NCOM == 3) {

	if (i == 1)
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);

	  Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
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

      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);

	if (i <  2) master -> mulY (nZ32, tmp);

	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Transform to Fourier space, smooth, add forcing.

      N[i] -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) {
	Veclib::fill (nP, -2.0*ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	N[i] -> addToPlane (0, tmp);
      }

      *N[i] *= -0.5;
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
      }

      // -- Transform to Fourier space, smooth, add forcing.
      
      N[i] -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
    }
  }
}


void altSkewSymmetric (Domain*         D ,
		       AuxField**      Us,
		       AuxField**      Uf,
		       vector<real_t>& ff)
// ---------------------------------------------------------------------------
// This routine is very similar to skewSymmetric(), but here either
// the convective (u.grad(u)) or conservative (div(uu)) parts of the
// Skew symmetric form are computed on alternating timesteps, thus
// approximating skew symmetry on average. This form seems nearly
// identical in stability properties to the full skew symmetric form,
// but has a lower operation count.
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int_t NCOM = D -> nField() - 1;	// -- Number of velocity components.
  const int_t nP   = Geometry::planeSize();

#if defined (ALIAS)
  const int_t       nZ32   = Geometry::nZProc();
#else
  const int_t       nZ32   = Geometry::nZ32();
#endif 
  const int_t       nZ     = Geometry::nZ();
  const int_t       nZP    = Geometry::nZProc();

  const int_t       nPP    = Geometry::nBlock();
  const int_t       nPR    = Geometry::nProc();
  const int_t       nTot   = Geometry::nTotProc();
  const int_t       nTot32 = nZ32 * nP;

  vector<real_t*>   u32 (NCOM), n32 (NCOM);
  vector<AuxField*> U   (NCOM), N   (NCOM);
  Field*            master = D -> u[0];

  vector<real_t> work ((2 * NCOM + 1) * nTot32);
  real_t*        tmp  = &work[0] + 2 * NCOM * nTot32;
  int_t          i, j;

  static int     toggle = 1; 	// -- Switch between u.grad(u) and div(uu).

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

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i < NCOM; i++) {

	// -- Terms involving azimuthal derivatives and frame components.

	if (NCOM == 3) {
	  if (i == 1)
	    Veclib::svvttvp (nTot32, -1.0, u32[2],1,u32[2],1,n32[1],1,n32[1],1);
	  if (i == 2)
	    Veclib::vmul    (nTot32, u32[2], 1, u32[1], 1, n32[2], 1);

	  if (nZ > 2) {
	    Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	    Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	    Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	    Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	    master -> gradient (nZ, nPP, tmp, 2);
	    Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	    Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	    Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i],1);
	  }
	}

	if (i == 2) master -> divY (nZ32, n32[i]);

	// -- 2D convective derivatives.

	for (j = 0; j < 2; j++) {
	  Veclib::copy (nTot32, u32[i], 1, tmp, 1);
	  master -> gradient (nZ32, nP, tmp, j);

	  if (i <  2) master -> mulY (nZ32, tmp);

	  Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
	}

	// -- Transform to Fourier space, smooth, add forcing.

	N[i] -> transform32 (FORWARD, n32[i]);
	master -> smooth (N[i]);

	ROOTONLY if (fabs (ff[i]) > EPSDP) {
	  Veclib::fill (nP, -ff[i], tmp, 1);
	  if (i < 2) master -> mulY (1, tmp);
	  N[i] -> addToPlane (0, tmp);
	}
	*N[i] *= -1.0;
      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NCOM; i++) {

	// -- Terms involving azimuthal derivatives and frame components.

	if (i == 0)
	  Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
	if (i == 1)
	  Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);

	if (NCOM == 3) {

	  if (i == 1)
	    Veclib::svvttvp (nTot32, -1., u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	  if (i == 2)
	    Veclib::svvtt   (nTot32,  2., u32[2], 1, u32[1], 1,      n32[2], 1);

	  if (nZ > 2) {
	    Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	    Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	    Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	    Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	    master -> gradient (nZ, nPP, tmp, 2);
	    Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	    Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	    Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
	  }
	}

	if (i == 2) master -> divY (nZ32, n32[i]);

	// -- 2D conservative derivatives.
     
	for (j = 0; j < 2; j++) {
	  Veclib::vmul (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	  master -> gradient (nZ32, nP, tmp, j);

	  if (i <  2) master -> mulY (nZ32, tmp);

	  Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
	}

	// -- Transform to Fourier space, smooth, add forcing.

	N[i] -> transform32 (FORWARD, n32[i]);
	master -> smooth (N[i]);

	ROOTONLY if (fabs (ff[i]) > EPSDP) {
	  Veclib::fill (nP, -ff[i], tmp, 1);
	  if (i < 2) master -> mulY (1, tmp);
	  N[i] -> addToPlane (0, tmp);
	}

	*N[i] *= -1.0;
      }
    }

  } else {			// -- Cartesian coordinates.

    if (toggle) { // -- Convective component u.grad(u).

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
	    
	}

	// -- Transform to Fourier space, smooth, add forcing.
	
	N[i] -> transform32 (FORWARD, n32[i]);
	master -> smooth (N[i]);
	
	ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -ff[i]);

	*N[i] *= -1.0;
      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NCOM; i++) {
	for (j = 0; j < NDIM; j++) {

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
	}

	// -- Transform to Fourier space, smooth, add forcing.
      
	N[i] -> transform32 (FORWARD, n32[i]);
	master -> smooth (N[i]);

	ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -ff[i]);
	*N[i] *= -1.0;
      }
    }
  }
 
  toggle = 1 - toggle;
}


void convective (Domain*         D ,
		AuxField**      Us,
		AuxField**      Uf,
		vector<real_t>& ff)
// ---------------------------------------------------------------------------
// Nonlinear terms N(u) in convective form are
//                 ~ ~
//           N  = - u . grad u
//           ~      ~        ~
//
// This has a similar operation count to altSkewSymmetric() but seems
// less stable.
//
// i.e., in Cartesian component form
//
//           N  = -  u  d(u ) / dx 
//            i       j    i      j
//
// in cylindrical coordinates
//
//           Nx = -{ud(u)/dx + vd(u)/dy + 1/y [wd(u)/dz]}
//           Ny = -{ud(v)/dx + vd(v)/dy + 1/y [wd(v)/dz - ww]}
//           Nz = -{ud(w)/dx + vd(w)/dy + 1/y [wd(w)/dz + wv]}
//
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in reference[3].
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int_t NCOM = D -> nField() - 1;	// -- Number of velocity components.
  const int_t nP   = Geometry::planeSize();

#if defined (ALIAS)
  const int_t       nZ32   = Geometry::nZProc();
#else
  const int_t       nZ32   = Geometry::nZ32();
#endif 
  const int_t       nZ     = Geometry::nZ();
  const int_t       nZP    = Geometry::nZProc();

  const int_t       nPP    = Geometry::nBlock();
  const int_t       nPR    = Geometry::nProc();
  const int_t       nTot   = Geometry::nTotProc();
  const int_t       nTot32 = nZ32 * nP;

  vector<real_t*>   u32 (NCOM), n32 (NCOM);
  vector<AuxField*> U   (NCOM), N   (NCOM);
  Field*            master = D -> u[0];

  vector<real_t> work ((2 * NCOM + 1) * nTot32);
  real_t*        tmp  = &work[0] + 2 * NCOM * nTot32;
  int_t          i, j;

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

      if (NCOM == 3) {
	if (i == 1)
	  Veclib::svvttvp (nTot32, -1.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::vmul    (nTot32, u32[2], 1, u32[1], 1, n32[2], 1);

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);
	}
      }

      if (i == 2) master -> divY (nZ32, n32[i]);

      // -- 2D convective derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);

	if (i <  2) master -> mulY (nZ32, tmp);

	Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Transform to Fourier space, smooth, add forcing.

      N[i] -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) {
	Veclib::fill (nP, -ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	N[i] -> addToPlane (0, tmp);
      }
      *N[i] *= -1.0;

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

      }

      // -- Transform to Fourier space, smooth, add forcing.
      
      N[i] -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -ff[i]);

      *N[i] *= -1.0;
    }
  }
}


void Stokes (Domain*         D ,
	     AuxField**      Us,
	     AuxField**      Uf,
	     vector<real_t>& ff)
// ---------------------------------------------------------------------------
// Stokes flow by definition does not have nonlinear terms but we
// still need to zero storage areas and add in body force ff.
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int_t NCOM = D -> nField() - 1;	// -- Number of velocity components.
  const int_t nP   = Geometry::planeSize();

#if defined (ALIAS)
  const int_t       nZ32   = Geometry::nZProc();
#else
  const int_t       nZ32   = Geometry::nZ32();
#endif 
  const int_t       nZ     = Geometry::nZ();
  const int_t       nZP    = Geometry::nZProc();

  const int_t       nPP    = Geometry::nBlock();
  const int_t       nPR    = Geometry::nProc();
  const int_t       nTot   = Geometry::nTotProc();
  const int_t       nTot32 = nZ32 * nP;

  vector<real_t*>   u32 (NCOM), n32 (NCOM);
  vector<AuxField*> U   (NCOM), N   (NCOM);
  Field*            master = D -> u[0];

  vector<real_t> work ((2 * NCOM + 1) * nTot32);
  real_t*        tmp  = &work[0] + 2 * NCOM * nTot32;
  int_t          i, j;

  for (i = 0; i < NCOM; i++) {
    *Uf[i] = 0.0;
    AuxField::swapData (D -> u[i], Us[i]);
    ROOTONLY if (fabs (ff[i]) > EPSDP) {
      Veclib::fill (nP, ff[i], tmp, 1);
      if (i < 2 && Geometry::cylindrical()) master -> mulY (1, tmp);
      Uf[i] -> addToPlane (0, tmp);
    }
  }
}
