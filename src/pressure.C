/*****************************************************************************
 * pressure.C: routines to deal with pressure field and boundary conditions.
 *****************************************************************************/

// $Id$


#include "Fem.h"


#ifdef __DECCXX
  #pragma define_template roll<double*>
#endif


// -- Static PBCmanager class variables: 

HOBC*  PBCmanager::store     = 0;
BC*    PBCmanager::essential = 0;
BC*    PBCmanager::hopbc     = 0;


void  PBCmanager::build (Field& f)
// ---------------------------------------------------------------------------
// Build class-scope storage structures required for high order PBCs.
// Reset Field-boundary BCs to have type either HOPBC or ESSENTIAL (0.0).
//
// Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting
// methods for the incompressible Navier--Stokes equations", JCP 9(2).
//
// All BCs of kind HOPBC have their dP/dn values re-evaluated at each step.
// Here the storage to enable this is allocated and organized.
// There is some wastage as memory is also allocated for ESSENTIAL BCs.
// ---------------------------------------------------------------------------
{
  // -- Initialize static storage.

  essential = new BC;
  essential -> id    = 0;
  essential -> kind  = ESSENTIAL;
  essential -> value = 0.0;

  hopbc = new BC;
  hopbc -> id        = 0;
  hopbc -> kind      = HOPBC;
  hopbc -> value     = 0.0;

  int  nOrder = iparam   ("N_TIME");
  int  nEdge  = f.nBound ();
  int  ntot, i, q;

  ntot = f.switchPressureBCs (hopbc, essential);

  store = new HOBC [nEdge];

  store[0].Px = rmatrix (nOrder, ntot);
  store[0].Py = rmatrix (nOrder, ntot);
  store[0].Ux = rmatrix (nOrder, ntot);
  store[0].Uy = rmatrix (nOrder, ntot);

  Veclib::zero (nOrder*ntot, *store[0].Px, 1);
  Veclib::zero (nOrder*ntot, *store[0].Py, 1);
  Veclib::zero (nOrder*ntot, *store[0].Ux, 1);
  Veclib::zero (nOrder*ntot, *store[0].Uy, 1);

  i = 0;
  ntot = 0;
  for (ListIterator<Boundary*> j(f.boundary_list);
       j.more(); j.next(), i++) {
    if (i) {
      store[i].Px = new real* [nOrder];
      store[i].Py = new real* [nOrder];
      store[i].Ux = new real* [nOrder];
      store[i].Uy = new real* [nOrder];
      for (q = 0; q < nOrder; q++) {
	store[i].Px[q]  = store[0].Px[q] + ntot;
	store[i].Py[q]  = store[0].Py[q] + ntot;
	store[i].Ux[q]  = store[0].Ux[q] + ntot;
	store[i].Uy[q]  = store[0].Uy[q] + ntot;
      }
    }
    ntot += j.current() -> nKnot();
  }
}


void  PBCmanager::maintain (int  step, Field***  Us, Field***  Uf)
// ---------------------------------------------------------------------------
// Update storage for evaluation of high-order pressure boundary condition.
//
// No smoothing is done to high-order spatial derivatives computed here.
//
// We add estimates of velocity local acceleration to the linear terms;
// since these must be constructed from velocity fields, there must be the
// same amount of storage as for the time order of the scheme, (since, e.g.
// for a first order scheme we need two levels of velocities to estimate a
// time derivative: these come from the new one passed in, and the boundary
// store).  Note also that this term cannot be estimated on the first step.
// ---------------------------------------------------------------------------
{
  int         Je     = iparam ("N_TIME");
  int         nOrder = Je;
  real        alpha[TIME_ORDER_MAX], gamma;
  const real  nu    = dparam ("KINVIS");
  const real  invDt = 1.0 / dparam ("DELTAT");

  const Field* Ux = Us[0][0];
  const Field* Uy = Us[1][0];
  const Field* Nx = Uf[0][0];
  const Field* Ny = Uf[1][0];

  int   i, q, np, offset, skip;

  // -- Roll grad P storage area up, load new level of nonlinear terms.

  Boundary *B1;
  
  ListIterator<Boundary*> j(Nx -> boundary_list);

  for (i = 0; j.more(); j.next(), i++) {
    B1 = j.current ();

    np     = B1 -> nKnot ();
    offset = B1 -> nOff  ();
    skip   = B1 -> nSkip ();

    roll (store[i].Px, nOrder);
    Veclib::copy (np, Nx -> data + offset, skip, store[i].Px[0], 1);

    roll (store[i].Py, nOrder);
    Veclib::copy (np, Ny -> data + offset, skip, store[i].Py[0], 1);
  }

  // -- Add in -nu * curl curl u and -du/dt.

  real *tmp, *wx, *wy;

  for (i = 0, j = Ux -> boundary_list; j.more(); j.next(), i++) {
    B1 = j.current ();

    np     = B1 -> nKnot ();
    offset = B1 -> nOff  ();
    skip   = B1 -> nSkip ();

    wx = rvector (np);
    wy = rvector (np);

    B1 -> curlCurl (Ux -> data, Uy -> data, wx, wy);

    Blas::axpy (np, -nu, wy, 1, store[i].Px[0], 1);
    Blas::axpy (np,  nu, wx, 1, store[i].Py[0], 1);

    // -- Estimate -du/dt by backwards differentiation, add in.

    if (step > 1 && !B1 -> isEssential()) {

      Je    = min (step - 1, Je);
      tmp   = wx;
      gamma = Icoef[Je - 1].gamma;
      Veclib::copy (Je, Icoef[Je - 1].alpha, 1, alpha, 1);
      
      Veclib::copy (np, Ux -> data + offset, skip, tmp, 1);
      Blas::scal   (np, gamma, tmp, 1);
      for (q = 0; q < Je; q++)
	Blas::axpy (np, -alpha[q], store[i].Ux[q], 1, tmp, 1);
      Blas::axpy (np, -invDt, tmp, 1, store[i].Px[0], 1);

      Veclib::copy (np, Uy -> data + offset, skip, tmp, 1);
      Blas::scal   (np, gamma, tmp, 1);
      for (q = 0; q < Je; q++)
	Blas::axpy (np, -alpha[q], store[i].Uy[q], 1, tmp, 1);
      Blas::axpy (np, -invDt, tmp, 1, store[i].Py[0], 1);
    }

    // -- Roll velocity storage area up, load new level.

    roll (store[i].Ux, nOrder);
    Veclib::copy (np, Ux -> data + offset, skip, store[i].Ux[0], 1);

    roll (store[i].Uy, nOrder);
    Veclib::copy (np, Uy -> data + offset, skip, store[i].Uy[0], 1);
    
    freeVector (wx);
    freeVector (wy);
  }
}


void  PBCmanager::evaluate (int   id,     int   np,  int   step,
			    real* value,  real* nx,  real* ny  )
// ---------------------------------------------------------------------------
// Load PBC storage with values obtained from HOBC multi-level storage.
//
// The boundary condition for evaluation is
//
//   dP       /                           du  \
//   -- = n . | N(u) - a + f + \nu*L(u) - --  |  =  n . grad P.
//   dn   ~   \ ~ ~    ~   ~       ~ ~    dt  /     ~
//
// Grad P is estimated at the end of the current timestep using explicit
// extrapolation.
// ---------------------------------------------------------------------------
{
  int    q, Je = iparam ("N_TIME");
  real   beta[TIME_ORDER_MAX];
  
  real  *tmpX = rvector (np);
  real  *tmpY = rvector (np);

  Je = min (step, Je);
  Veclib::copy (Je, Icoef[Je - 1].beta, 1, beta, 1);

  Veclib::zero (np, tmpX, 1);
  Veclib::zero (np, tmpY, 1);
  
  --id;

  for (q = 0; q < Je; q++) {
    Blas::axpy (np, beta[q], store[id].Px[q], 1, tmpX, 1);
    Blas::axpy (np, beta[q], store[id].Py[q], 1, tmpY, 1);
  }
    
  Veclib::vmul  (np, nx, 1, tmpX, 1, value, 1);
  Veclib::vvtvp (np, ny, 1, tmpY, 1, value, 1, value, 1);

  freeVector (tmpX);
  freeVector (tmpY);
}
