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


void nonlinear (Domain*         D ,
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
// Nonlinear terms N(u) in robust skew-symmetric form are (the average
// of the so-called non-conservative and conservative ways of
// computing the nonlinear terms, theoretically equivalent if the
// velocity field is divergence-free):
//
//           N  = -0.5 ( u . grad u + div uu )
//           ~           ~        ~       ~~
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
// The formulation of nonlinear terms used here is so-called
// "alternating skew symmetric" method (first documented by Bob Kerr)
// which uses the non-conservative and conservative forms of the
// nonlinear terms on alternating timesteps. This has shown in testing
// to be as robust as the full skew symmetric method but costs half as
// much. For simplicity it is now (as of early 2013) the only form
// implemented.
//
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// As of January 2013, dealiasing in Z of serial 3D computations has
// been removed in the interests of simplicity and maintainability.
// (The idea is that more Z planes could be used, with parallel
// operation if needs be.)
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int_t NCOM = D -> nField() - 1;	// -- Number of velocity components.
  const int_t nP   = Geometry::planeSize();
  const int_t nZ   = Geometry::nZ();
  const int_t nZP  = Geometry::nZProc();
  const int_t nTot = Geometry::nTotProc();

  vector<AuxField*> U (NCOM), N (NCOM); 
  Field*            master = D -> u[0];
  int_t             i, j;

  static int               toggle = 1; // -- Switch u.grad(u) or div(uu).
  static real_t*           work;
  static vector<AuxField*> Uphys (NCOM);
  static AuxField*         tmp;

  if (!work) { 			// -- First time through.
    work = new real_t [static_cast<size_t>((NCOM + 1) * nTot)];
    for (i = 0; i < NCOM; i++)
      Uphys[i] = new AuxField (work + i * nTot, nZP, D -> elmt);
    tmp = new AuxField (work + NCOM * nTot, nZP, D -> elmt);
  }

  for (i = 0; i < NCOM; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
    U[i] = Us[i];
    N[i] = Uf[i];
    *N[i] = 0.0;
    (*Uphys[i] = *U[i]) . transform (INVERSE);
  }

  if (Geometry::cylindrical()) {

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i < NCOM; i++) {

	// -- Terms involving azimuthal derivatives and frame components.

	if (NCOM == 3) {
	  if (i == 1) N[1] -> times      (*Uphys[2], *Uphys[2]);
	  if (i == 2) N[2] -> timesMinus (*Uphys[2], *Uphys[1]);

	  if (nZ > 2) {
	    (*tmp = *U[i]) . gradient (2) . transform (INVERSE);
	    N[i] -> timesMinus (*Uphys[2], *tmp);
	  }
	}

	if (i == 2) N[i] -> divY ();

	// -- 2D convective derivatives.

	for (j = 0; j < 2; j++) {
	  (*tmp = *Uphys[i]) . gradient (j);
	  if (i <  2) tmp -> mulY ();
	  N[i] -> timesMinus (*Uphys[j], *tmp);
	}

	// -- Transform to Fourier space, smooth, add forcing.

	N[i] -> transform (FORWARD);
	master -> smooth (N[i]);

	ROOTONLY if (fabs (ff[i]) > EPSDP) {
	  Veclib::fill (nP, ff[i], work + NCOM * nTot, 1);
	  if (i < 2) master -> mulY (1, work + NCOM * nTot);
	  N[i] -> addToPlane (0, work + NCOM * nTot);
	}
      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NCOM; i++) {

	// -- Terms involving azimuthal derivatives and frame components.

	if (i == 0) N[0] -> timesMinus (*Uphys[0], *Uphys[1]);
	if (i == 1) N[1] -> timesMinus (*Uphys[1], *Uphys[1]);

	if (NCOM == 3) {

	  if (i == 1) N[1] -> timesPlus (*Uphys[2], *Uphys[2]);
	  if (i == 2) {
	    tmp -> times (*Uphys[2], *Uphys[1]);
	    N[2] -> axpy (-2., *tmp);
	  }

	  if (nZ > 2) {
	    tmp -> times (*Uphys[i], *Uphys[2]);
	    (*tmp) . transform (FORWARD). gradient (2). transform (INVERSE);
	    *N[i] -= *tmp;
	  }
	}

	if (i == 2) N[i] -> divY ();

	// -- 2D conservative derivatives.
     
	for (j = 0; j < 2; j++) {
	  (*tmp). times (*Uphys[j], *Uphys[i]) . gradient (j);
	  if (i <  2) tmp -> mulY ();
	  *N[i] -= *tmp;
	}

	// -- Transform to Fourier space, smooth, add forcing.

	N[i] -> transform (FORWARD);
	master -> smooth (N[i]);

	ROOTONLY if (fabs (ff[i]) > EPSDP) {
	  Veclib::fill (nP, ff[i], work + NCOM * nTot, 1);
	  if (i < 2) master -> mulY (1, work + NCOM * nTot);
	  N[i] -> addToPlane (0, work + NCOM * nTot);
	}
      }
    }

  } else {			// -- Cartesian coordinates.

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i < NCOM; i++) {
	for (j = 0; j < NDIM; j++) {
      
	  // -- Perform n_i -= u_j d(u_i) / dx_j.

	  if (j == 2) (*tmp = *U[i]) . gradient (j) . transform (INVERSE);
	  else    (*tmp = *Uphys[i]) . gradient (j);
	  N[i] -> timesMinus (*Uphys[j], *tmp);
	}

	master -> smooth (N[i]);
	N[i] -> transform (FORWARD);
	ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, ff[i]);
      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NCOM; i++) {
	for (j = 0; j < NDIM; j++) {

	  // -- Perform n_i -= d(u_i u_j) / dx_j.

	  tmp -> times (*Uphys[i], *Uphys[j]);
	  if (j == 2) tmp -> transform (FORWARD);
	  tmp -> gradient (j);
	  if (j == 2) tmp -> transform (INVERSE);
	  *N[i] -= *tmp;
	}

	master -> smooth (N[i]);
	N[i] -> transform (FORWARD);
	ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, ff[i]);
      }

    }
  }
 
  toggle = 1 - toggle;
}
