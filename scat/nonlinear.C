///////////////////////////////////////////////////////////////////////////////
// nonlinear.C: 
//
// Copyright (C) 2001 <--> $Date$, Hugh Blackburn.
//
// This version of basic semtex DNS includes transport of a passive
// scalar 'c', carried as the NCOMth Field; pressure is
// NCOM+1. Cartesian and cylindrical coordinates.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <scat.h>


void nonlinear (Domain*         D ,
		AuxField**      Us,
		AuxField**      Uf,
		vector<real_t>& ff)
// ---------------------------------------------------------------------------
// Compute advection terms in Navier--Stokes equations: N(u).
//
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 ~ ~
//                   N = -0.5 ( u . grad u + div uu )
//                   ~          ~        ~       ~~
//
// i.e., in Cartesian component form
//
//              N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//               i           j    i      j      i j      j
//
// Alternating skew symmetric: alternately u.grad(u) and div(uu) are
// made on successive timesteps.
//
// The temperature transport term u . grad c is also built this way here.
// ---------------------------------------------------------------------------
{
  const int_t NDIM = Geometry::nDim();	// -- Number of space dimensions.
  const int_t NCOM = D -> nField() - 2;	// -- Number of velocity components.
  const int_t nP   = Geometry::planeSize();
  const int_t nZ   = Geometry::nZ();
  const int_t nZP  = Geometry::nZProc();
  const int_t nTot = Geometry::nTotProc();

  vector<AuxField*> U (NCOM + 1), N (NCOM + 1); 
  Field*            master = D -> u[0];
  int_t             i, j;

  static int               toggle = 1; // -- Switch u.grad(u) or div(uu).
  static real_t*           work;
  static vector<AuxField*> Uphys (NCOM + 1);
  static AuxField*         tmp;

  if (!work) { 			// -- First time through.
    work = new real_t [static_cast<size_t>((NCOM + 2) * nTot)];
    for (i = 0; i <= NCOM; i++)
      Uphys[i] = new AuxField (work + i * nTot, nZP, D -> elmt);
    tmp = new AuxField (work + (NCOM + 1) * nTot, nZP, D -> elmt);
  }

  for (i = 0; i <= NCOM; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
    U[i] = Us[i];
    N[i] = Uf[i];
    *N[i] = 0.0;
    (*Uphys[i] = *U[i]) . transform (INVERSE);
  }

#if 1				// -- New version.

  if (Geometry::cylindrical()) {

    if (toggle) { // -- Convective component u.grad(u).

      for (i = 0; i <= NCOM; i++) {

	// -- Terms involving azimuthal derivatives and frame components.

	if (NCOM > 2) {
	  if (i == 1) N[1] -> times      (*Uphys[2], *Uphys[2]);
	  if (i == 2) N[2] -> timesMinus (*Uphys[2], *Uphys[1]);

	  if (nZ > 2) {
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

      for (i = 0; i <= NCOM; i++) {

	// -- Terms involving azimuthal derivatives and frame components.

	if (i == 0) N[0] -> timesMinus (*Uphys[0], *Uphys[1]);
	if (i == 1) N[1] -> timesMinus (*Uphys[1], *Uphys[1]);

	if (NCOM > 2) {

	  if (i == 1) N[1] -> timesPlus (*Uphys[2], *Uphys[2]);

	  if (i == 2) {
	    tmp -> times (*Uphys[2], *Uphys[1]);
	    N[2] -> axpy (-2., *tmp);
	  }

	  if (i == 3) N[3] -> timesMinus (*Uphys[1], *Uphys[3]);

	  if (nZ > 2) {
	    tmp -> times (*Uphys[i], *Uphys[2]);
	    (*tmp) . transform (FORWARD). gradient (2). transform (INVERSE);
	    *N[i] -= *tmp;
	  }
	}

	if (i >= 2) N[i] -> divY ();

	// -- 2D conservative derivatives.
     
	for (j = 0; j < 2; j++) {
	  (*tmp). times (*Uphys[j], *Uphys[i]) . gradient (j);
	  if (i < 2) tmp -> mulY ();
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

      for (i = 0; i <= NCOM; i++) {
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

      for (i = 0; i <= NCOM; i++) {
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

#else
  if (Geometry::cylindrical()) {		// -- Cylindrical coordinates.

    for (i = 0; i <= NCOM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);

      if (NCOM > 2) {
	if (i == 1)		// -- radial compt.
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)		// -- azimuthal compt.
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);
	if (i == 3)		// -- scalar.
	  Veclib::vmul    (nTot32,       u32[3], 1, u32[1], 1,      n32[3], 1);

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

      if (i >= 2) master -> divY (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);

	if (i < 2) master -> mulY (nZ32, tmp);

	Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul       (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);

	if (i < 2) master -> mulY (nZ32, tmp);

	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) {
	Veclib::fill (nP, -2.0*ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	N[i] -> addToPlane (0, tmp);
      }

      *N[i] *= -0.5;		// -- Skew-symmetric NL.
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i <= NCOM; i++) {
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
      
      N[i]   -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -2.0*ff[i]);

      *N[i] *= -0.5;
    }
  }
#endif
}
