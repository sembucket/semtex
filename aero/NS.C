///////////////////////////////////////////////////////////////////////////////
// NS.C:  Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// Copyright (c) Hugh Blackburn 1997--2001.
//
// This version includes body coupling and solves Navier--Stokes in an
// accelerating reference frame.  In addition, it allows Smagorinsky
// based LES to be used: for this, define LES during compilation.
//
// It is assumed that any velocity boundary conditions are not functions
// of "z", i.e. that all the information is in the zeroth Fourier mode.
// Supplied velocity boundary conditions must belong to GROUP "velocity".
// Body-wall boundaries must belong to another group.
//
// References:
// [1] Karniadakis, Israeli & Orszag 1991.  "High-order splitting
//     methods for the incompressible Navier--Stokes equations", JCP 9(2).
// [2] Blackburn & Henderson 1996.  "Lock-in behaviour in simulated
//     vortex-induced vibration", Exptl Thermal & Fluid Sci., 12(2), 184--189.
// [3] Blackburn & Henderson 1999.  "A study of two-dimensional flow past
//     an oscillating cylinder", JFM 385, 255--286.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "aero.h"

typedef ModalMatrixSys Msys;
static  integer        NDIM, NORD;

static void   nonLinear (Domain*, AuxField**, AuxField**, AuxField*, Vector&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const integer, AuxField*, Msys*);


void NavierStokes (Domain*       D,
		   Body*         B,
		   AeroAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains initialized storage for velocity Fields
// 'u', 'v' ('w' for 3D) and constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NDIM = Geometry::nDim();
  NORD = (integer) Femlib::value ("N_TIME");

  integer       i, j, k;
  const real    dt     =           Femlib::value ("D_T");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();
  const integer ntot   = Geometry::nTotProc();
  real*         alloc  = new real [(size_t) 2 * NDIM * NORD * ntot];

  // -- Initialize body motion coupling terms.

  Vector a = {0.0, 0.0, 0.0};   // -- Body/frame acceleration.
  Vector v = {0.0, 0.0, 0.0};	// -- Body velocity, to adjust velocity BCs.

  // -- Create & initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [(size_t) 2 * NORD];
  AuxField*** Uf = Us + NORD;

  for (k = 0, i = 0; i < NORD; i++) {
    Us[i] = new AuxField* [(size_t) 2 * NDIM];
    Uf[i] = Us[i] + NDIM;
    for (j = 0; j < NDIM; j++) {
      *(Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
      *(Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
    }
  }

#if defined (LES)

#if defined (NOMODEL)
  Femlib::value ("REFVIS", Femlib::value ("2.0 * KINVIS"));
#endif

  if (Femlib::value ("REFVIS") > 0.0) {
    real kinVis = Femlib::value ("REFVIS");
    real refVis = Femlib::value ("KINVIS");
    Femlib::value ("KINVIS", kinVis);
    Femlib::value ("REFVIS", refVis);
  } 

  AuxField* EV = new AuxField (new real [(size_t) ntot], nZ, D -> elmt, 'e');
  *EV = 0.0;
  ROOTONLY EV -> addToPlane (0, Femlib::value ("REFVIS - KINVIS"));

#else

  AuxField* EV = 0;
#endif

  // -- Create global matrix systems.

  Msys** MMS = preSolve (D);

  // -- Create multi-level storage for pressure BCS.

  Field* Pressure = D -> u[NDIM];
  PBCmgr::build (Pressure);

  // -- Timestepping loop.

  while (D -> step < nStep) {

    // -- Set a to old value of frame acceleration.

    ROOTONLY a = B -> acceleration();
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

#if defined (LES)
    // -- Compute spatially-varying kinematic eddy viscosity.
    eddyViscosity (D, Us[0], Uf[0], EV);
#endif

    // -- Unconstrained forcing substep.

    nonLinear (D, Us[0], Uf[0], EV, a);

    waveProp  (D, (const AuxField***)Us, (const AuxField***)Uf);

    // -- Body motion, get velocity at new time level.

    ROOTONLY {
      B -> move (D -> step);
      v = B -> velocity();
    }

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
		      (const AuxField**)Us[0],
		      (const AuxField**)Uf[0], 0);
    ROOTONLY PBCmgr::accelerate (a, D -> u[0]);

    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NDIM);
    setPForce ((const AuxField**)Us[0], Uf[0]);
    Solve     (D, NDIM,  Uf[0][0], MMS[NDIM]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NDIM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NDIM);

    // -- Viscous correction substep, adjust velocity BCs for frame motion.

    for (i = 0; i < NDIM; i++) {
      ROOTONLY {
	D -> u[i] -> evaluateM0Boundaries (D -> step);
	if      (i == 0) D -> u[0] -> addToM0Boundaries (-v.x, "velocity");
	else if (i == 1) D -> u[1] -> addToM0Boundaries (-v.y, "velocity");
      }
      Solve (D, i, Uf[0][i], MMS[i]);
    }

    // -- Find forces with new domain fields.

    ROOTONLY B -> force (D);

    // -- Process results of this step.

    A -> analyse (Us[0]);
  }

#if defined (LES)

  // -- Dump ratio eddy/molecular viscosity to file visco.fld.

  ofstream          evfl;
  vector<AuxField*> visco (1);

  visco[0] = EV;

  ROOTONLY {
    evfl.open ("visco.fld", ios::out);
    EV -> addToPlane (0, Femlib::value ("KINVIS - REFVIS"));
  }
  (*EV /= Femlib::value ("REFVIS")) . transform (INVERSE);

  writeField (evfl, D -> name, D -> step, D -> time, visco);

  ROOTONLY evfl.close();

#endif
}


static void nonLinear (Domain*    D ,
		       AuxField** Us,
		       AuxField** Uf,
		       AuxField*  EV,
		       Vector&    a )
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations
//                   N(u) + div(2*EV*S) - a.
//
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 ~ ~
//                   N = -0.5 ( u . grad u + div uu )
//                   ~          ~        ~       ~~
//
// i.e., in Cartesian component form
//
//              N  = -0.5 ( u  d(u u ) / dx  + d(u u ) / dx ).
//               i           j    i j      j      i j      j
//
// All product terms are evaluated pseudospectrally, in physical space.
// ---------------------------------------------------------------------------
{
  integer           i, j;
  vector<real>      A (NDIM);
  const real        EPS    = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;
  const integer     nZ     = Geometry::nZ();
  const integer     nZP    = Geometry::nZProc();
  const integer     nP     = Geometry::planeSize();
  const integer     nPP    = Geometry::nBlock();
  const integer     nPR    = Geometry::nProc();
  const integer     nTot   = Geometry::nTotProc();
  const integer     nZ32   = Geometry::nZ32();
  const integer     nTot32 = nZ32 * nP;
  vector<real>      work ((2 * NDIM + 1) * nTot32);
  vector<real*>     u32 (NDIM);
  vector<real*>     n32 (NDIM);
  vector<AuxField*> U   (NDIM);
  vector<AuxField*> N   (NDIM);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * NDIM * nTot32;

  ROOTONLY A[0] = a.x; A[1] = a.y; if (NDIM == 3) A[2] = a.z;

  Veclib::zero ((2 * NDIM + 1) * nTot32, work(), 1); // -- A catch-all cleanup.

  for (i = 0; i < NDIM; i++) {
    u32[i] = work() +  i        * nTot32;
    n32[i] = work() + (i + NDIM) * nTot32;
  }

#if defined (LES)

  // -- Start with contribution from divergence of SGS.

  EV -> transform32 (INVERSE, tmp);
  Blas::scal (nTot32, -4.0, tmp, 1); // -- 2 = -0.5 * -4 ... see below.

  // -- Diagonal stress-divergence terms.

  for (i = 0; i < NDIM; i++) {

    Us[i] -> transform32 (INVERSE, n32[i]);
    Veclib::vmul (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

    if (i == 2) {
      Femlib::exchange   (n32[2], nZ32,        nP, FORWARD);      
      Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, FORWARD);
      Veclib::zero       (nTot32 - nTot, n32[2] + nTot, 1);
      master -> gradient (nZ, nPP, n32[2], 2);
      Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, INVERSE);
      Femlib::exchange   (n32[2], nZ32,        nP, INVERSE);      
    } else
      master -> gradient (nZ32, nP, n32[i], i);
  }

  // -- Off-diagonal stress-divergence terms.

  for (i = 0; i < NDIM; i++)
    for (j = i + 1; j < NDIM; j++) {

      Uf[i + j - 1] -> transform32 (INVERSE, u32[0]);
      Veclib::vmul                 (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy                 (nTot32,         u32[0], 1, u32[1], 1);

      // -- Super-diagonal.
      
      if (j == 2) {
	Femlib::exchange   (u32[0], nZ32,        nP, FORWARD);
	Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, FORWARD);
	Veclib::zero       (nTot32 - nTot, u32[0] + nTot, 1);
	master -> gradient (nZ, nPP, u32[0], 2);
	Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, INVERSE);
	Femlib::exchange   (u32[0], nZ32,        nP, INVERSE);
      } else
	master -> gradient (nZ32, nP, u32[0], j);
      Veclib::vadd (nTot32, u32[0], 1, n32[i], 1, n32[i], 1);

      // -- Sub-diagonal.

      master -> gradient (nZ32, nP, u32[1], i);
      Veclib::vadd       (nTot32, u32[1], 1, n32[j], 1, n32[j], 1);
    }
#endif

  for (i = 0; i < NDIM; i++) {
    U[i] = Us[i];
    N[i] = Uf[i];
    AuxField::swapData  (D -> u[i], U[i]);
    U[i] -> transform32 (INVERSE, u32[i]);
  }

  for (i = 0; i < NDIM; i++) {
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
    ROOTONLY if (fabs (A[i]) > EPS) N[i] -> addToPlane (0, 2.0*A[i]);
    *N[i] *= -0.5;
  }
}


static void waveProp (Domain*           D ,
		      const AuxField*** Us,
		      const AuxField*** Uf)
// ---------------------------------------------------------------------------
// Compute the first substep of stiffly-stable timestepping scheme.
//
// On entry, the most recent velocity fields are in Us, and the most
// recent nonlinear terms in Uf.  The intermediate velocity field u^ is
// computed and left in D's velocity areas. 
// ---------------------------------------------------------------------------
{
  integer           i, q;
  vector<AuxField*> H (NDIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NDIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const integer Je = min (D -> step, NORD);
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

  for (i = 0; i < NDIM; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[q][i]);
      H[i] -> axpy ( beta [q]    , *Uf[q][i]);
    }
}


static void setPForce (const AuxField** Us,
		       AuxField**       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in Us.  Create div u^ / D_T
// in the first dimension of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) (*Uf[i]  = *Us[i]) . gradient (i);
  for (i = 1; i < NDIM; i++)  *Uf[0] += *Uf[i];

  *Uf[0] /= dt;
}


static void project (const Domain* D ,
		     AuxField**    Us,
		     AuxField**    Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in Us.  Constrain velocity field:
//                    u^^ = u^ - D_T * grad P;
// u^^ is left in Uf.
//
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to
// create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt    = Femlib::value ("D_T");
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < NDIM; i++) {

    (*Uf[i] = *D -> u[NDIM]) . gradient (i);

    Us[i] -> axpy (-dt, *Uf[i]);
    Field::swapData (Us[i], Uf[i]);

    *Uf[i] *= alpha;
  }
}


static Msys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE == 1 selects iterative solvers for velocity components,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  const integer           nmodes = Geometry::nModeProc();
  const integer           base   = Geometry::baseMode();
  const integer           itLev  = (integer) Femlib::value ("ITERATIVE");
  const real              beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [(size_t) (NDIM + 1)];
  integer                 i;
  vector<real>            alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, alpha());
  const real              lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NDIM; i++)
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i],    (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NDIM] = new Msys
      (0.0,     beta, base, nmodes, E, D -> b[NDIM], (itLev<2)?DIRECT:JACPCG);

  return M;
}


static void Solve (Domain*       D,
		   const integer i,
		   AuxField*     F,
		   Msys*         M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D->u[i], using F as a forcing Field.
// Iterative or direct solver selected on basis of field type, step,
// time order and command-line arguments.
// ---------------------------------------------------------------------------
{
  const integer step = D -> step;

  if (i < NDIM && step < NORD) { // -- We need a temporary matrix system.
    const integer Je      = min (step, NORD);    
    const integer base    = Geometry::baseMode();
    const integer nmodes  = Geometry::nModeProc();
    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
    const real    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real    beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
