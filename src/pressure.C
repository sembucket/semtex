/*****************************************************************************
 * pressure.C: routines to deal with pressure field and boundary conditions.
 *****************************************************************************/

static char 
RCSid[]="$Id$";

#include "Fem.h"

#ifdef __DECCXX
  #pragma define_template roll<double*>
#endif


// -- Static PBCmanager class variables: 

HOBC*  PBCmanager::store     = 0;
BC*    PBCmanager::essential = 0;
BC*    PBCmanager::hopbc     = 0;


void  PBCmanager::build (Field& P)
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
  essential -> kind  = BC::essential;
  essential -> value = 0.0;

  hopbc = new BC;
  hopbc -> id        = 0;
  hopbc -> kind      = BC::hopbc;
  hopbc -> value     = 0.0;

  int  nTime = iparam ("N_TIME");
  int  nEdge = P.nBound    ();
  int  ntot  = P.resetPBCs (hopbc, essential);
  int  i, q;

  store = new HOBC [nEdge];

  store[0].Px = rmatrix (nTime, ntot);
  store[0].Py = rmatrix (nTime, ntot);
  store[0].Ux = rmatrix (nTime, ntot);
  store[0].Uy = rmatrix (nTime, ntot);

  Veclib::zero (nTime*ntot, *store[0].Px, 1);
  Veclib::zero (nTime*ntot, *store[0].Py, 1);
  Veclib::zero (nTime*ntot, *store[0].Ux, 1);
  Veclib::zero (nTime*ntot, *store[0].Uy, 1);

  ntot = 0;
  i    = 0;
  for (ListIterator<Boundary*> j(P.boundary_list); j.more(); j.next()) {
    if (ntot) {
      store[i].Px = new real* [nTime];
      store[i].Py = new real* [nTime];
      store[i].Ux = new real* [nTime];
      store[i].Uy = new real* [nTime];
      for (q = 0; q < nTime; q++) {
	store[i].Px[q]  = store[0].Px[q] + ntot;
	store[i].Py[q]  = store[0].Py[q] + ntot;
	store[i].Ux[q]  = store[0].Ux[q] + ntot;
	store[i].Uy[q]  = store[0].Uy[q] + ntot;
      }
    }
    ntot += j.current() -> nKnot();
    i    += 1;
  }
}


void  PBCmanager::maintain (int            step  ,
			    const Field*   master,
			    const Field*** Us    ,
			    const Field*** Uf    )
// ---------------------------------------------------------------------------
// Update storage for evaluation of high-order pressure boundary condition.
//
// No smoothing is done to high-order spatial derivatives computed here.
//
// We add estimates of local rate of change of velocity to the linear terms;
// since these must be extrapolated from velocity fields, there must be the
// same amount of storage as for the time order of the scheme, (since, e.g.
// for a first order scheme we need two levels of velocities to estimate a
// time derivative: these come from the new one passed in, and the boundary
// store).  Note also that this term cannot be estimated on the first step.
//
// Field* master gives a list of pressure boundary conditions with which to
// traverse storage areas (note this assumes equal-order interpolations).
// ---------------------------------------------------------------------------
{
  const int   nTime = iparam ("N_TIME");
  const real  nu    = dparam ("KINVIS");
  const real  invDt = 1.0 / dparam ("DELTAT");

  const Field* Ux = Us[0][0];
  const Field* Uy = Us[1][0];
  const Field* Nx = Uf[0][0];
  const Field* Ny = Uf[1][0];

  int   i, np, offset, skip;

  // -- Roll grad P storage area up, load new level of nonlinear terms.

  register     Boundary*   B;
  ListIterator<Boundary*>  j (master -> boundary_list);

  for (i = 0; j.more (); j.next (), i++) {
    B = j.current ();

    np     = B -> nKnot ();
    offset = B -> nOff  ();
    skip   = B -> nSkip ();

    roll (store[i].Px, nTime);
    Veclib::copy (np, Nx -> data + offset, skip, store[i].Px[0], 1);

    roll (store[i].Py, nTime);
    Veclib::copy (np, Ny -> data + offset, skip, store[i].Py[0], 1);
  }

  // -- Add in -nu * curl curl u and -du/dt.

  int    Je, q;
  real  *tmp, *wx, *wy;
  real*  alpha = rvector (Integration::OrderMax + 1);

  for (i = 0, j.reset(); j.more(); j.next(), i++) {
    B = j.current ();

    np     = B -> nKnot ();
    offset = B -> nOff  ();
    skip   = B -> nSkip ();

    wx = rvector (np);
    wy = rvector (np);

    B -> curlCurl (Ux -> data, Uy -> data, wx, wy);

    Blas::axpy (np, -nu, wy, 1, store[i].Px[0], 1);
    Blas::axpy (np,  nu, wx, 1, store[i].Py[0], 1);

    // -- Estimate -du/dt by backwards differentiation, add in on HOPBs.

    if (step > 1 && !B -> isEssential()) {

      Je    = min (step - 1, nTime);
      tmp   = wx;
      Integration::StifflyStable (Je, alpha);
      
      Veclib::copy (np, Ux -> data + offset, skip, tmp, 1);
      Blas::scal   (np, alpha[0], tmp, 1);
      for (q = 0; q < Je; q++)
	Blas::axpy (np, alpha[q + 1], store[i].Ux[q], 1, tmp, 1);
      Blas::axpy (np, -invDt, tmp, 1, store[i].Px[0], 1);

      Veclib::copy (np, Uy -> data + offset, skip, tmp, 1);
      Blas::scal   (np, alpha[0], tmp, 1);
      for (q = 0; q < Je; q++)
	Blas::axpy (np, alpha[q + 1], store[i].Uy[q], 1, tmp, 1);
      Blas::axpy (np, -invDt, tmp, 1, store[i].Py[0], 1);
    }

    // -- Roll velocity storage area up, load new level.

    roll (store[i].Ux, nTime);
    Veclib::copy (np, Ux -> data + offset, skip, store[i].Ux[0], 1);

    roll (store[i].Uy, nTime);
    Veclib::copy (np, Uy -> data + offset, skip, store[i].Uy[0], 1);
    
    freeVector (wx);
    freeVector (wy);
  }

  freeVector (alpha);
}


void  PBCmanager::evaluate (int   id,     int   np,  int   step,
			    real* value,  real* nx,  real* ny  )
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
  register int  q, Je = iparam ("N_TIME");

  real*  beta  = rvector (Integration::OrderMax);
  real*  tmpX  = rvector (np);
  real*  tmpY  = rvector (np);

  Je = min (step, Je);
  Integration::Extrapolation (Je, beta);

  Veclib::zero (np, tmpX, 1);
  Veclib::zero (np, tmpY, 1);
  
  --id;

  for (q = 0; q < Je; q++) {
    Blas::axpy (np, beta[q], store[id].Px[q], 1, tmpX, 1);
    Blas::axpy (np, beta[q], store[id].Py[q], 1, tmpY, 1);
  }
    
  Veclib::vmul  (np, nx, 1, tmpX, 1, value, 1);
  Veclib::vvtvp (np, ny, 1, tmpY, 1, value, 1, value, 1);

  freeVector (beta);
  freeVector (tmpX);
  freeVector (tmpY);
}
