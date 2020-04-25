///////////////////////////////////////////////////////////////////////////////
// nonlinear.cpp
//
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u) + f.
//
// Here N(u) represents the nonlinear advection terms in the N--S
// equations transposed to the RHS and f is a vector of body force
// per unit mass (with possible space-time dependency).
//
// Copyright (c) 1994 <--> $Date: 2019/05/30 06:36:03 $, Hugh Blackburn
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

static char RCS[] = "$Id: nonlinear.cpp,v 9.1 2019/05/30 06:36:03 hmb Exp $";

#include <dns.h>


void skewSymmetric (Domain*     D ,
		    BCmgr*      B ,
		    AuxField**  Us,
		    AuxField**  Uf,
		    FieldForce* FF)
// ---------------------------------------------------------------------------
// Velocity field data areas of D (which on entry contain velocity
// data from the previous timestep, in Fourier space) and first level
// of Us are swapped (so that subsequently Us stores the old velocity
// data), then the next stage of nonlinear forcing terms N(u) are
// computed from velocity fields and left in the first (i.e. the
// supplied) level of Uf.
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
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).
//
// If a scalar, c, is present, also compute the advection term
// -0.5*[u.grad(c)+div(uc)].  This is a straightforward extension to
// the above for Cartesian coordinates (the loop over i is extended by
// 1 and the last component of 'u_i' is c), whereas in cylindrical
// coordinates we have
//
//           Nc = -0.5 {ud(c)/dx + vd(c)/dy + d(uc)/dx + d(vc)/dy +
//                 1/y [wd(c)/dz + d(wc)/dz + cv ]}
//
// Data are transformed to physical space for most of the operations.
// For gradients in the Fourier direction however, the data must be
// transferred back to Fourier space prior to differention in z, then
// trsnsformed back.  The final form of N delivered at the end of the
// routine is in Fourier space.  Note that there is no longer any
// provision for dealiasing in the Fourier direction and that all
// storage areas of D->u (including pressure) are overwritten here.
//  
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).  
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp    = D -> u[NADV]; // -- Pressure is used for scratch.
  Field*            master = D -> u[0];	   // -- For smoothing operations.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (master, Uphys, NCOM);

  if (Geometry::cylindrical()) {

    for (i = 0; i < NADV; i++) {

      if (i == 0) N[0] -> timesMinus (*Uphys[0], *Uphys[1]);
      if (i == 1) N[1] -> timesMinus (*Uphys[1], *Uphys[1]);

      // -- Terms involving azimuthal derivatives and frame components.

      if (NCOM == 3) {

        if (i == 1) {
          tmp -> times (*Uphys[2], *Uphys[2]);
          N[1] -> axpy ( 2.0, *tmp);
        }

        if (i == 2) {
          tmp -> times (*Uphys[2], *Uphys[1]);
          N[2] -> axpy (-3.0, *tmp);
        }

	if (i == 3)
	  N[3] -> timesMinus (*Uphys[3], *Uphys[1]);
	  
        if (NDIM == 3) {
          (*tmp = *U[i]) . gradient (2) . transform (INVERSE);
          N[i] -> timesMinus (*Uphys[2], *tmp);

          tmp -> times (*Uphys[i], *Uphys[2]);
          (*tmp) . transform (FORWARD). gradient (2). transform (INVERSE);
          *N[i] -= *tmp;
        }
      }

      if (i >= 2) N[i] -> divY ();
      
      // -- 2D convective derivatives.

      for (j = 0; j < 2; j++) {
        (*tmp = *Uphys[i]) . gradient (j);
        if (i < 2) tmp -> mulY ();
        N[i] -> timesMinus (*Uphys[j], *tmp);
      }

      // -- 2D conservative derivatives.

      for (j = 0; j < 2; j++) {
        (*tmp). times (*Uphys[j], *Uphys[i]) . gradient (j);
        if (i < 2) tmp -> mulY ();
        *N[i] -= *tmp;
      }

      *N[i] *= 0.5;      // -- Average the two forms to get skew-symmetric.

      if (i < NCOM) {
        FF     -> addPhysical (N[i], tmp, i, Uphys);
        N[i]   -> transform   (FORWARD);
        FF     -> addFourier  (N[i], tmp, i, U);
        master -> smooth      (N[i]);
      } else {
        master -> smooth    (N[i]);
        N[i]   -> transform (FORWARD);
      }
    }

  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NADV; i++) {
      for (j = 0; j < NDIM; j++) {

        // -- Perform n_i -= u_j d(u_i) / dx_j.

        if (j == 2) (*tmp = *U[i]) . gradient (j) . transform (INVERSE);
        else    (*tmp = *Uphys[i]) . gradient (j);
        N[i] -> timesMinus (*Uphys[j], *tmp);

        // -- Perform n_i -= d(u_i u_j) / dx_j.

        tmp -> times (*Uphys[i], *Uphys[j]);
        if (j == 2) tmp -> transform (FORWARD);
        tmp -> gradient (j);
        if (j == 2) tmp -> transform (INVERSE);
        *N[i] -= *tmp;

      }

      *N[i] *= 0.5;      // -- Average the two forms to get skew-symmetric.
      
      if (i < NCOM) {    // -- Body forces act only on momentum equations.
        FF     -> addPhysical (N[i], tmp, i, Uphys);
        N[i]   -> transform   (FORWARD);
        FF     -> addFourier  (N[i], tmp, i, U);
        master -> smooth     (N[i]);
      } else {
        master -> smooth    (N[i]);
        N[i]   -> transform (FORWARD);
      }
    }
  }
}


void altSkewSymmetric (Domain*     D ,
		       BCmgr*      B ,
		       AuxField**  Us,
		       AuxField**  Uf,
		       FieldForce* FF)
// ---------------------------------------------------------------------------
// This is a cheaper approximation to the full skew-symmetric
// computation (see previous routine).  The formulation of nonlinear
// terms used here is so-called "alternating skew symmetric" method
// (first documented by Bob Kerr) which uses the non-conservative and
// conservative forms of the nonlinear terms on alternating
// timesteps. This has shown in testing to be as robust as the full
// skew symmetric method but costs half as much.
//
// NB: this is now used as the default advection scheme.  As noted, it
// is about as robust but computationally cheaper than full
// skew-symmetric.  However, it has a slight defect, too: owing to
// alternation, it may produce small high-frequency temporal
// oscillations in the solution, of period 2*D_T.  If that is a
// problem, try full skew-symmetric (previous routine) or the
// non-conservative form of the nonlinear terms (next routine).
//
// Here are the two components of skew symmetric form in cylindrical
// coords (cf. what is written above for the full skew symmetric form):
//
// Non-conservative form:
//
//           Nx = -{ud(u)/dx + vd(u)/dy + 1/y [wd(u)/dz     ]}
//           Ny = -{ud(v)/dx + vd(v)/dy + 1/y [wd(v)/dz - ww]}
//           Nz = -{ud(w)/dx + vd(w)/dy + 1/y [wd(w)/dz + wv]}
//           Nc = -{ud(c)/dx + vd(c)/dy + 1/y [wd(c)/dz     ]}
//
// Conservative form:
//
//           Nx = -{d(uu)/dx + d(vu)/dy + 1/y [d(uw)/dz + vu     ]}
//           Ny = -{d(uv)/dx + d(vv)/dy + 1/y [d(vw)/dz + vv - ww]}
//           Nz = -{d(uw)/dx + d(vw)/dy + 1/y [d(ww)/dz + 2wv    ]}
//           Nc = -{d(uc)/dx + d(vc)/dy + 1/y [d(wc)/dz + cv     ]}
//
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).  
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();
  
  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp    = D -> u[NADV]; // -- Pressure is used for scratch.
  Field*            master = D -> u[0];	   // -- For smoothing operations.
  int_t             i, j;
  static int        toggle = 1;            // -- Switch u.grad(u) or div(uu).
  
  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (master, Uphys, NCOM);

  if (Geometry::cylindrical()) {

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i < NADV; i++) {

        // -- Terms involving azimuthal derivatives and frame components.

        if (NCOM == 3) {
          if (i == 1) N[1] -> times      (*Uphys[2], *Uphys[2]);
          if (i == 2) N[2] -> timesMinus (*Uphys[2], *Uphys[1]);

          if (NDIM == 3) {
            (*tmp = *U[i]) . gradient (2) . transform (INVERSE);
            N[i] -> timesMinus (*Uphys[2], *tmp);
          }
        }

        if (i >= 2) N[i] -> divY ();

        // -- 2D convective derivatives.

        for (j = 0; j < 2; j++) {
          (*tmp = *Uphys[i]) . gradient (j);
          if (i <  2) tmp -> mulY ();
          N[i] -> timesMinus (*Uphys[j], *tmp);
        }

        if (i < NCOM) {     // -- Body forces act only on momentum equations.
          FF     -> addPhysical (N[i], tmp, i, Uphys);
          N[i]   -> transform   (FORWARD);
          FF     -> addFourier  (N[i], tmp, i, U);
          master -> smooth      (N[i]);
        } else {
          N[i]   -> transform   (FORWARD);
          master -> smooth      (N[i]);
        }
      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NADV; i++) {

        // -- Terms involving azimuthal derivatives and frame components.

        if (i == 0) N[0] -> timesMinus (*Uphys[0], *Uphys[1]);
        if (i == 1) N[1] -> timesMinus (*Uphys[1], *Uphys[1]);

        if (NCOM == 3) {
          if (i == 1) N[1] -> timesPlus (*Uphys[2], *Uphys[2]);
          if (i == 2) {
            tmp -> times (*Uphys[2], *Uphys[1]);
            N[2] -> axpy (-2., *tmp);
          }
	  if (i == 3) N[3] -> timesMinus (*Uphys[3], *Uphys[1]);
          if (NDIM == 3) {
            tmp -> times (*Uphys[i], *Uphys[2]);
            (*tmp) . transform (FORWARD). gradient (2). transform (INVERSE);
            *N[i] -= *tmp;
          }
        }

        if (i >= 2) N[i] -> divY ();

        // -- 2D conservative derivatives.

        for (j = 0; j < 2; j++) {
          (*tmp). times (*Uphys[j], *Uphys[i]) . gradient (j);
          if (i <  2) tmp -> mulY ();
          *N[i] -= *tmp;
        }

        if (i < NCOM) {     // -- Body forces act only on momentum equations.
          FF     -> addPhysical (N[i], tmp, i, Uphys);
          N[i]   -> transform   (FORWARD);
          FF     -> addFourier  (N[i], tmp, i, U);
          master -> smooth      (N[i]);
        } else {
          N[i]   -> transform   (FORWARD);
          master -> smooth      (N[i]);
        }
      }
    }

  } else {			// -- Cartesian coordinates.

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i < NADV; i++) {
        for (j = 0; j < NDIM; j++) { // -- Perform n_i -= u_j d(u_i) / dx_j.

          if (j == 2) (*tmp = *U[i]) . gradient (j) . transform (INVERSE);
          else    (*tmp = *Uphys[i]) . gradient (j);
          N[i] -> timesMinus (*Uphys[j], *tmp);
        }

        if (i < NCOM) {
          FF     -> addPhysical (N[i], tmp, i, Uphys);
          N[i]   -> transform   (FORWARD);
          FF     -> addFourier  (N[i], tmp, i, U);
          master -> smooth      (N[i]);
        } else {
          master -> smooth      (N[i]);
          N[i]   -> transform   (FORWARD);
        }
      }

    } else { // -- Conservative component div(uu).

      for (i = 0; i < NADV; i++) {
        for (j = 0; j < NDIM; j++) { // -- Perform n_i -= d(u_i u_j) / dx_j.

          tmp -> times (*Uphys[i], *Uphys[j]);
          if (j == 2) tmp -> transform (FORWARD);
          tmp -> gradient (j);
          if (j == 2) tmp -> transform (INVERSE);
          *N[i] -= *tmp;
        }

        if (i < NCOM) {     // -- Body forces act only on momentum equations.
          FF     -> addPhysical (N[i], tmp, i, Uphys);
          N[i]   -> transform   (FORWARD);
          FF     -> addFourier  (N[i], tmp, i, U);
          master -> smooth      (N[i]);
        } else {
          N[i]   -> transform   (FORWARD);
          master -> smooth      (N[i]);
        }
      }
    }
  }

  toggle = 1 - toggle;
}


void convective (Domain*     D ,
		 BCmgr*      B ,
		 AuxField**  Us,
		 AuxField**  Uf,
		 FieldForce* FF)
// ---------------------------------------------------------------------------
// Nonlinear terms N(u) in (non-conservative) convective form are
//                 ~ ~
//           N  = - u . grad u
//           ~      ~        ~
//
// This has a similar operation count to altSkewSymmetric() but seems
// somewhat less stable/robust at high Reynolds numbers.  However, it
// *is* smooth in time, so if stability isn't an issue, but temporal
// smoothness and computational cost are, give it a try if you don't
// want to pay for full skew-symmetric.
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
//           Nc = -{ud(c)/dx + vd(c)/dy + 1/y [wd(c)/dz     ]}  
//
// NB: for the cylindrical coordinate formulation we actually here 
// compute y*Nx, y*Ny, Nz, as outlined in Blackburn & Sherwin (2004).
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp    = D -> u[NADV]; // -- Pressure is used for scratch.
  Field*            master = D -> u[0];	   // -- For smoothing operations.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (master, Uphys, NCOM);

  if (Geometry::cylindrical()) {

    for (i = 0; i < NADV; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (NCOM == 3) {
        if (i == 1) N[1] -> times      (*Uphys[2], *Uphys[2]);
        if (i == 2) N[2] -> timesMinus (*Uphys[2], *Uphys[1]);

        if (NDIM == 3) {
          (*tmp = *U[i]) . gradient (2) . transform (INVERSE);
          N[i] -> timesMinus (*Uphys[2], *tmp);
        }
      }

      if (i >= 2) N[i] -> divY ();

      // -- 2D convective derivatives.

      for (j = 0; j < 2; j++) {
        (*tmp = *Uphys[i]) . gradient (j);
        if (i < 2) tmp -> mulY ();
        N[i] -> timesMinus (*Uphys[j], *tmp);
      }

      if (i < NCOM) {     // -- Body forces act only on momentum equations.
        FF     -> addPhysical (N[i], tmp, i, Uphys);
        N[i]   -> transform   (FORWARD);
        FF     -> addFourier  (N[i], tmp, i, U);
        master -> smooth      (N[i]);
      } else {
        N[i]   -> transform   (FORWARD);
        master -> smooth      (N[i]);
      }
    }

  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NADV; i++) {
      for (j = 0; j < NDIM; j++) { // -- Perform n_i -= u_j d(u_i) / dx_j.

        if (j == 2) (*tmp = *U[i]) . gradient (j) . transform (INVERSE);
        else    (*tmp = *Uphys[i]) . gradient (j);
        N[i] -> timesMinus (*Uphys[j], *tmp);
      }

      if (i < NCOM) {     // -- Body forces act only on momentum equations.
        FF     -> addPhysical (N[i], tmp, i, Uphys);
        N[i]   -> transform   (FORWARD);
        FF     -> addFourier  (N[i], tmp, i, U);
        master -> smooth      (N[i]);
      } else {
        N[i]   -> transform   (FORWARD);
        master -> smooth      (N[i]);
      }
    }
  }
}


void Stokes (Domain*     D ,
	     BCmgr*      B ,
	     AuxField**  Us,
	     AuxField**  Uf,
	     FieldForce* FF)
// ---------------------------------------------------------------------------
// Stokes flow by definition does not have nonlinear terms but we
// still need to zero storage areas and add in body force ff.
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();
  const int_t NADV = D -> nAdvect();
  const int_t NCOM = D -> nVelCmpt();

  vector<AuxField*> U (NADV), N (NADV), Uphys (NADV);
  AuxField*         tmp    = D -> u[NADV]; // -- Pressure is used for scratch.
  Field*            master = D -> u[0];	   // -- For smoothing operations.
  int_t             i, j;

  for (i = 0; i < NADV; i++) {
    Uphys[i] = D -> u[i];
    N[i]     = Uf[i];
    U[i]     = Us[i];
    *N[i]    = 0.0;
    *U[i]    = *Uphys[i];
    Uphys[i] -> transform (INVERSE);
  }

  B -> maintainPhysical (master, Uphys, NCOM);
 
  for (i = 0; i < NCOM; i++) {
    FF     -> addPhysical (N[i], tmp, i, Uphys);
    N[i]   -> transform   (FORWARD);
    FF     -> addFourier  (N[i], tmp, i, U);
    master -> smooth      (N[i]);
  }

  if (D -> hasScalar()) {
    N[NCOM] -> transform  (FORWARD);
    master  -> smooth     (N[NCOM]);
  }
}
