///////////////////////////////////////////////////////////////////////////////
// pressure.C: routines to deal with pressure field boundary conditions.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// Class variables _Pn & _Un provide storage for the mode equivalents of
//   _Pn:  normal gradient of the pressure field,
//   _Un: normal component of velocity,
// and are used to construct explicit extrapolative estimates of the
// natural BCs for the pressure field at the next time level.
//
// Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting
// methods for the incompressible Navier--Stokes equations", JCP 9(2).
//
// _Pn & _Un are indexed by time level, boundary, data plane, and
// location in that order (e.g. _Pn[time][boundary][plane][i]).
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>

real_t**** PBCmgr::_Pnx = 0;
real_t**** PBCmgr::_Pny = 0;
real_t**** PBCmgr::_Unx = 0;
real_t**** PBCmgr::_Uny = 0;


void PBCmgr::build (const Field* P)
// ---------------------------------------------------------------------------
// Build class-scope storage structures required for high order PBCs.
//
// All BCs of kind HOPBC have their dP/dn values re-evaluated at each step.
// Here the storage to enable this is allocated and organized.
// There is some wastage as memory is also allocated for essential BCs.
// ---------------------------------------------------------------------------
{
  const int_t np    = Geometry::nP();
  const int_t nTime = Femlib::ivalue ("N_TIME");
  const int_t nEdge = P -> _nbound;
  const int_t nZ    = P -> _nz;
  int_t       i, j, k;

  _Pnx = new real_t*** [static_cast<size_t>(nTime)];
  _Pny = new real_t*** [static_cast<size_t>(nTime)];
  _Unx = new real_t*** [static_cast<size_t>(nTime)];
  _Uny = new real_t*** [static_cast<size_t>(nTime)];

  for (i = 0; i < nTime; i++) {
    _Pnx[i] = new real_t** [static_cast<size_t>(4 * nEdge)];
    _Pny[i] = _Pnx[i] + nEdge;
    _Unx[i] = _Pny[i] + nEdge;
    _Uny[i] = _Unx[i] + nEdge;

    for (j = 0; j < nEdge; j++) {
      _Pnx[i][j] = new real_t* [static_cast<size_t>(4 * nZ)];
      _Pny[i][j] = _Pnx[i][j] + nZ;
      _Unx[i][j] = _Pny[i][j] + nZ;
      _Uny[i][j] = _Unx[i][j] + nZ;

      for (k = 0; k < nZ; k++) {
	_Pnx[i][j][k] = new real_t [static_cast<size_t>(4 * np)];
	_Pny[i][j][k] = _Pnx[i][j][k] + np;
	_Unx[i][j][k] = _Pny[i][j][k] + np;
	_Uny[i][j][k] = _Unx[i][j][k] + np;

	Veclib::zero (4 * np, _Pnx[i][j][k], 1);
      }
    }
  }
}


void PBCmgr::maintain (const int_t    step   ,
		       const Field*     P      ,
		       const AuxField** Us     ,
		       const AuxField** Uf     ,
		       const bool       timedep)
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
  const real_t nu    = Femlib::value ("KINVIS");
  const real_t invDt = 1.0 / Femlib::value ("D_T");
  const int_t  nTime = Femlib::ivalue ("N_TIME");
  const int_t  nEdge = P -> _nbound;
  const int_t  nZ    = P -> _nz;
  const int_t  nP    =  Geometry::nP();

  const AuxField* Ux = Us[0];
  const AuxField* Uy = Us[1];
  const AuxField* Uz = (Geometry::nPert() == 3) ? Us[2] : 0;

  const AuxField* Nx = Uf[0];
  const AuxField* Ny = Uf[1];

  const vector<Boundary*>& BC = P -> _bsys -> BCs (0);
  register Boundary*       B;
  register int_t           i, k, q;
  int_t                    offset, skip, Je;

  // -- Roll grad P storage area up, load new level of nonlinear terms Uf.

  rollv (_Pnx, nTime);
  rollv (_Pny, nTime);

  for (i = 0; i < nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();

    for (k = 0; k < nZ; k++) {
      Veclib::copy (nP, Nx -> _plane[k] + offset, skip, _Pnx[0][i][k], 1);
      Veclib::copy (nP, Ny -> _plane[k] + offset, skip, _Pny[0][i][k], 1);

      // -- For cylindrical coordinates, N_ are radius-premultiplied. Cancel.

      if (Geometry::cylindrical()) {
	B -> divY (_Pnx[0][i][k]);
	B -> divY (_Pny[0][i][k]);
      }
    }
  }

  // -- Add in -nu * curl curl u. There are 3 cases to deal with:
  //    perturbation is real_t, half-complex or full-complex.

  vector<real_t> work (5 * sqr(nP) + 7 * nP + Integration::OrderMax + 1);
  real_t         *UxRe, *UxIm, *UyRe, *UyIm, *UzRe, *UzIm, *tmp;
  real_t*        wrk   = &work[0];
  real_t*        xr    = wrk + 5*sqr(nP) + 3*nP;
  real_t*        xi    = xr + nP;
  real_t*        yr    = xi + nP;
  real_t*        yi    = yr + nP;
  real_t*        alpha = yi + nP;

  for (i = 0; i < nEdge; i++) {
    B = BC[i];

    if (Geometry::nZ() == 1) {

      UxRe = Ux -> _plane[0];
      UyRe = Uy -> _plane[0];

      if (Geometry::problem() == Geometry::O2_2D ||
	  Geometry::problem() == Geometry::SO2_2D ) { // -- Real perturbation.
	B->curlCurl(0,UxRe,0,UyRe,0,0,0,xr,0,yr,0,wrk);
      } else {			    // -- Half-complex perturbation.
	UzIm = Uz -> _plane[0];
	B->curlCurl(1,UxRe,0,UyRe,0,0,UzIm,xr,0,yr,0,wrk);
      }
      Blas::axpy (nP, -nu, xr, 1, _Pnx[0][i][0], 1);
      Blas::axpy (nP, -nu, yr, 1, _Pny[0][i][0], 1);
  
    } else {			    // -- Full complex perturbation.
      UxRe = Ux -> _plane[0];
      UxIm = Ux -> _plane[1];
      UyRe = Uy -> _plane[0];
      UyIm = Uy -> _plane[1];
      UzRe = Uz -> _plane[0];
      UzIm = Uz -> _plane[1];

      B->curlCurl(1,UxRe,UxIm,UyRe,UyIm,UzRe,UzIm,xr,xi,yr,yi,wrk);

      Blas::axpy (nP, -nu, xr, 1, _Pnx[0][i][0], 1);
      Blas::axpy (nP, -nu, xi, 1, _Pnx[0][i][1], 1);
      Blas::axpy (nP, -nu, yr, 1, _Pny[0][i][0], 1);
      Blas::axpy (nP, -nu, yi, 1, _Pny[0][i][1], 1);
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
	    Blas::axpy (nP, alpha[q + 1], _Unx[q][i][k], 1, tmp, 1);
	  Blas::axpy (nP, -invDt, tmp, 1, _Pnx[0][i][k], 1);
	  
	  Veclib::copy (nP, Uy -> _plane[k] + offset, skip, tmp, 1);
	  Blas::scal   (nP, alpha[0], tmp, 1);
	  for (q = 0; q < Je; q++)
	    Blas::axpy (nP, alpha[q + 1], _Uny[q][i][k], 1, tmp, 1);
	  Blas::axpy (nP, -invDt, tmp, 1, _Pny[0][i][k], 1);
	}
      }
    }

    // -- Roll velocity storage area up, load new level.

    rollv (_Unx, nTime);
    rollv (_Uny, nTime);
      
    for (i = 0; i < nEdge; i++) {
      B      = BC[i];
      offset = B -> dOff ();
      skip   = B -> dSkip();
    
      for (k = 0; k < nZ; k++) {
	Veclib::copy (nP, Ux -> _plane[k] + offset, skip, _Unx[0][i][k], 1);
	Veclib::copy (nP, Uy -> _plane[k] + offset, skip, _Uny[0][i][k], 1);
      }
    }
  }
}


void PBCmgr::evaluate (const int_t   id   ,
		       const int_t   np   ,
		       const int_t   plane,
		       const int_t   step ,
		       const real_t* nx   ,
		       const real_t* ny   ,
		       real_t*       tgt  )
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

  register int_t q, Je = Femlib::ivalue ("N_TIME");
  vector<real_t> work (Integration::OrderMax + 2 * np);
  real_t*        beta  = &work[0];
  real_t*        tmpX  = beta + Integration::OrderMax;
  real_t*        tmpY  = tmpX + np;

  Je = min (step, Je);
  Integration::Extrapolation (Je, beta);
  Veclib::zero (2 * np, tmpX, 1);
  
  for (q = 0; q < Je; q++) {
    Blas::axpy (np, beta[q], _Pnx[q][id][plane], 1, tmpX, 1);
    Blas::axpy (np, beta[q], _Pny[q][id][plane], 1, tmpY, 1);
  }
    
  Veclib::vvtvvtp (np, nx, 1, tmpX, 1, ny, 1, tmpY, 1, tgt, 1);
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
  register int_t           i;

  for (i = 0; i < u -> _nbound; i++) {
    B = BC[i];

    B -> addForGroup ("velocity", a.x, _Pnx[0][i][0]);
    B -> addForGroup ("velocity", a.y, _Pny[0][i][0]);
  }
}
