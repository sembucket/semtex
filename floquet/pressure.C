///////////////////////////////////////////////////////////////////////////////
// pressure.C: routines to deal with pressure field boundary conditions.
//
// Copyright (C) 1994, 2001 Hugh Blackburn
//
// Class variables Pn & Un provide storage for the mode equivalents of
//   Pn:  normal gradient of the pressure field,
//   Un: normal component of velocity,
// and are used to construct explicit extrapolative estimates of the natural
// BCs for the pressure field at the next time level.
//
// Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting
// methods for the incompressible Navier--Stokes equations", JCP 9(2).
//
// Pn & Un are indexed by time level, boundary, data plane, and location in
// that order (e.g. Pn[time][boundary][plane][i]).
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>

real**** PBCmgr::Pnx = 0;
real**** PBCmgr::Pny = 0;
real**** PBCmgr::Unx = 0;
real**** PBCmgr::Uny = 0;


void PBCmgr::build (const Field* P)
// ---------------------------------------------------------------------------
// Build class-scope storage structures required for high order PBCs.
//
// All BCs of kind HOPBC have their dP/dn values re-evaluated at each step.
// Here the storage to enable this is allocated and organized.
// There is some wastage as memory is also allocated for essential BCs.
// ---------------------------------------------------------------------------
{
  const integer np    = Geometry::nP();
  const integer nTime = (integer) Femlib::value ("N_TIME");
  const integer nEdge = P -> _nbound;
  const integer nZ    = P -> _nz;
  integer       i, j, k;

  Pnx = new real*** [(size_t) nTime];
  Pny = new real*** [(size_t) nTime];
  Unx = new real*** [(size_t) nTime];
  Uny = new real*** [(size_t) nTime];

  for (i = 0; i < nTime; i++) {
    Pnx[i] = new real** [(size_t) (4 * nEdge)];
    Pny[i] = Pnx[i] + nEdge;
    Unx[i] = Pny[i] + nEdge;
    Uny[i] = Unx[i] + nEdge;

    for (j = 0; j < nEdge; j++) {
      Pnx[i][j] = new real* [(size_t) (4 * nZ)];
      Pny[i][j] = Pnx[i][j] + nZ;
      Unx[i][j] = Pny[i][j] + nZ;
      Uny[i][j] = Unx[i][j] + nZ;

      for (k = 0; k < nZ; k++) {
	Pnx[i][j][k] = new real [(size_t) (4 * np)];
	Pny[i][j][k] = Pnx[i][j][k] + np;
	Unx[i][j][k] = Pny[i][j][k] + np;
	Uny[i][j][k] = Unx[i][j][k] + np;

	Veclib::zero (4 * np, Pnx[i][j][k], 1);
      }
    }
  }
}


void PBCmgr::maintain (const integer    step   ,
		       const Field*     P      ,
		       const AuxField** Us     ,
		       const AuxField** Uf     ,
		       const integer    timedep)
// ---------------------------------------------------------------------------
// Update storage for evaluation of high-order pressure boundary condition.
// Storage order for each edge represents a CCW traverse of element boundaries.
//
// If the velocity field varies in time on HOPB field boundaries (e.g. due
// to time-varying BCs) the local fluid acceleration will be estimated
// from input velocity fields by explicit extrapolation if timedep is true.
// This correction cannot be carried out at the first timestep, since the
// required extrapolation cannot be done.  If the acceleration is known,
// (for example, a known reference frame acceleration) it is probably better
// to leave timedep unset, and to use PBCmgr::accelerate() to add in the
// accelerative term.  Note also that since grad P is dotted with n, the
// unit outward normal, at a later stage, timedep only needs to be set if
// there are wall-normal accelerative terms.
//
// Field* master gives a list of pressure boundary conditions with which to
// traverse storage areas (note this assumes equal-order interpolations).
//
// No smoothing is done to high-order spatial derivatives computed here.
// ---------------------------------------------------------------------------
{
  const real      nu    =           Femlib::value ("KINVIS");
  const real      invDt = 1.0     / Femlib::value ("D_T");
  const integer   nTime = (integer) Femlib::value ("N_TIME");
  const integer   nEdge = P -> _nbound;
  const integer   nZ    = P -> _nz;
  const integer   nP    =  Geometry::nP();

  const AuxField* Ux = Us[0];
  const AuxField* Uy = Us[1];
  const AuxField* Uz = (Geometry::nPert() == 3) ? Us[2] : 0;
  real            *UxRe, *UxIm, *UyRe, *UyIm, *UzRe, *UzIm, *tmp;

  const AuxField* Nx = Uf[0];
  const AuxField* Ny = Uf[1];

  const vector<Boundary*>& BC = P -> _bsys -> BCs (0);
  register Boundary*       B;
  register integer         i, k, q;
  integer                  offset, skip, Je;

  vector<real> work (4 * nP + Integration::OrderMax + 1);

  // -- Roll grad P storage area up, load new level of nonlinear terms Uf.

  rollv (Pnx, nTime);
  rollv (Pny, nTime);

  for (i = 0; i < nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();

    for (k = 0; k < nZ; k++) {
      Veclib::copy (nP, Nx -> _plane[k] + offset, skip, Pnx[0][i][k], 1);
      Veclib::copy (nP, Ny -> _plane[k] + offset, skip, Pny[0][i][k], 1);
    }
  }

  // -- Add in -nu * curl curl u. There are 3 cases to deal with:
  //    perturbation is real, half-complex or full-complex.

  real* xr    = work();
  real* xi    = xr + nP;
  real* yr    = xi + nP;
  real* yi    = yr + nP;
  real* alpha = yi + nP;

  for (i = 0; i < nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();

    if (Geometry::nZ() == 1) {

      UxRe = Ux -> _plane[0];
      UyRe = Uy -> _plane[0];

      if (Geometry::nPert() == 2) { // -- Real perturbation.
	B -> curlCurl (0, UxRe, 0, UyRe, 0, 0, 0,    xr, 0, yr, 0);
      } else {			    // -- Half-complex perturbation.
	UzIm = Uz -> _plane[0];
	B -> curlCurl (1, UxRe, 0, UyRe, 0, 0, UzIm, xr, 0, yr, 0);
      }
      Blas::axpy (nP, -nu, xr, 1, Pnx[0][i][0], 1);
      Blas::axpy (nP, -nu, yr, 1, Pny[0][i][0], 1);
  
    } else {			    // -- Full complex peturbation.
      UxRe = Ux -> _plane[0];
      UxIm = Ux -> _plane[1];
      UyRe = Uy -> _plane[0];
      UyIm = Uy -> _plane[1];
      UzRe = Uz -> _plane[0];
      UzIm = Uz -> _plane[1];

      B -> curlCurl (1, UxRe,UxIm, UyRe,UyIm, UzRe,UzIm, xr,xi, yr,yi);

      Blas::axpy (nP, -nu, xr, 1, Pnx[0][i][0], 1);
      Blas::axpy (nP, -nu, xi, 1, Pnx[0][i][1], 1);
      Blas::axpy (nP, -nu, yr, 1, Pny[0][i][0], 1);
      Blas::axpy (nP, -nu, yi, 1, Pny[0][i][1], 1);
    }
  }

  if (timedep) {

    // -- Estimate -du / dt by backwards differentiation and add in.
    
    if (step > 1) {
      Je  = min (step - 1, nTime);
      tmp = xr;
      Integration::StifflyStable (Je, alpha);
      
      for (i = 0; i < nEdge; i++) {
	B      = BC[i];
	offset = B -> dOff ();
	skip   = B -> dSkip();

	for (k = 0; k < nZ; k++) {
	  ROOTONLY if (k == 1) continue;

	  Veclib::copy (nP, Ux -> _plane[k] + offset, skip, tmp, 1);
	  Blas::scal   (nP, alpha[0], tmp, 1);
	  for (q = 0; q < Je; q++)
	    Blas::axpy (nP, alpha[q + 1], Unx[q][i][k], 1, tmp, 1);
	  Blas::axpy (nP, -invDt, tmp, 1, Pnx[0][i][k], 1);
	  
	  Veclib::copy (nP, Uy -> _plane[k] + offset, skip, tmp, 1);
	  Blas::scal   (nP, alpha[0], tmp, 1);
	  for (q = 0; q < Je; q++)
	    Blas::axpy (nP, alpha[q + 1], Uny[q][i][k], 1, tmp, 1);
	  Blas::axpy (nP, -invDt, tmp, 1, Pny[0][i][k], 1);
	}
      }
    }

    // -- Roll velocity storage area up, load new level.

    rollv (Unx, nTime);
    rollv (Uny, nTime);
      
    for (i = 0; i < nEdge; i++) {
      B      = BC[i];
      offset = B -> dOff ();
      skip   = B -> dSkip();
    
      for (k = 0; k < nZ; k++) {
	Veclib::copy (nP, Ux -> _plane[k] + offset, skip, Unx[0][i][k], 1);
	Veclib::copy (nP, Uy -> _plane[k] + offset, skip, Uny[0][i][k], 1);
      }
    }
  }
}


void PBCmgr::evaluate (const integer id   ,
		       const integer np   ,
		       const integer plane,
		       const integer step ,
		       const real*   nx   ,
		       const real*   ny   ,
		       real*         tgt  )
// ---------------------------------------------------------------------------
// Load PBC value with values obtained from HOBC multi-level storage.
//
// The boundary condition for evaluation is
//
//   dP       /                           du  \
//   -- = n . | N(u) - a + f + \nu*L(u) - --  |  =  n . grad P.
//   dn   ~   \ ~ ~    ~   ~       ~ ~    dt  /     ~
//
// Grad P is estimated at the end of the current timestep using explicit
// extrapolation, then dotted into n.
// ---------------------------------------------------------------------------
{
  if (step < 1) return;

  register integer q, Je = (integer) Femlib::value ("N_TIME");
  vector<real>     work (Integration::OrderMax + 2 * np);
  real*            beta  = work();
  real*            tmpX  = beta + Integration::OrderMax;
  real*            tmpY  = tmpX + np;

  Je = min (step, Je);
  Integration::Extrapolation (Je, beta);
  Veclib::zero (2 * np, tmpX, 1);
  
  for (q = 0; q < Je; q++) {
    Blas::axpy (np, beta[q], Pnx[q][id][plane], 1, tmpX, 1);
    Blas::axpy (np, beta[q], Pny[q][id][plane], 1, tmpY, 1);
  }
    
  Veclib::vmul  (np, nx, 1, tmpX, 1, tgt, 1);
  Veclib::vvtvp (np, ny, 1, tmpY, 1, tgt, 1, tgt, 1);
}


void PBCmgr::accelerate (const Vector& a,
			 const Field*  u)
// ---------------------------------------------------------------------------
// Add in frame acceleration term on boundaries that correspond to
// essential velocity BCs (a = - du/dt).  Note that the acceleration
// time level should correspond to the time level in the most recently
// updated pressure gradient storage.  Work only takes place on zeroth
// Fourier mode.
//
// Yes, this is a HACK!
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = u -> _bsys -> BCs (0);
  register Boundary*       B;
  register integer         i;

  for (i = 0; i < u -> _nbound; i++) {
    B = BC[i];

    B -> addForGroup ("velocity", a.x, Pnx[0][i][0]);
    B -> addForGroup ("velocity", a.y, Pny[0][i][0]);
  }
}
