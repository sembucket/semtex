//////////////////////////////////////////////////////////////////////////////
// NS.C:  Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// This version includes body coupling and solves Navier--Stokes in an
// accelerating reference frame.
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
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <aero.h>

#ifdef __DECCXX
  #pragma define_template roll<AuxField*>
#endif

typedef ModalMatrixSystem ModeSys;
static  int               DIM;

static void  nonLinear (Domain*, AuxField***, AuxField***, Vector&);
static void  waveProp  (Domain*, const AuxField***, const AuxField***);
static void  setPForce (const AuxField***, AuxField***);
static void  project   (const Domain*, AuxField***, AuxField***);

static ModeSys** preSolve (const Domain*);
static void      Solve    (Field*, AuxField*, ModeSys*, const int, const int);


void NavierStokes (Domain*   D,
		   Body*     B,
		   Analyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains initialized storage for velocity Fields
// 'u', 'v' ('w' for 3D) and constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  DIM = D -> nField() - 1;

  int        i, j;
  const real dt       =       Femlib::value ("D_T");
  const int  nOrder   = (int) Femlib::value ("N_TIME");
  const int  nStep    = (int) Femlib::value ("N_STEP");
  const int  nZ       = (int) Femlib::value ("N_Z");

  Field*     Pressure = D -> u[DIM];
  ModeSys**  MMS      = preSolve (D);

  Vector  a = {0.0, 0.0, 0.0};   // -- Frame acceleration for N--S.
  Vector  v = {0.0, 0.0, 0.0};	 // -- Body velocity, to adjust velocity BCs.

  // -- Set up multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [DIM];
  AuxField*** Uf = new AuxField** [DIM];

  for (i = 0; i < DIM; i++) {
    Us[i] = new AuxField* [nOrder];
    Uf[i] = new AuxField* [nOrder];
    for (j = 0; j < nOrder; j++) {
      Us[i][j] = new AuxField (D -> Esys, nZ);
      Uf[i][j] = new AuxField (D -> Esys, nZ);
    }
  }

  // -- Set up multi-level storage for pressure BCS.

  PBCmgr::build (Pressure);

  // -- Timestepping loop.

  while (D -> step < nStep) {

    // -- Set a to old value of frame acceleration.

    a = B -> acceleration();
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us, Uf, a);
    waveProp  (D, (const AuxField***) Us, (const AuxField***) Uf);

    // -- Body motion, get velocity at new time level.

    B -> move (D -> step);
    v = B -> velocity();

    // -- Pressure projection substep.

    PBCmgr::maintain   (D -> step, Pressure,
		        (const AuxField***) Us,
			(const AuxField***) Uf,
		        0);	// -- *NOT* time dependent!
    PBCmgr::accelerate (a, D -> u[0]);

    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < DIM; i++) {
      AuxField::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce ((const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Update multilevel velocity storage.

    for (i = 0; i < DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll (Us[i], nOrder);
    }

    // -- Viscous correction substep, adjust velocity BCs for frame motion.

    for (i = 0; i < DIM; i++) {
      D -> u[i] -> evaluateM0Boundaries (D -> step);
      if      (i == 0) D -> u[0] -> addToM0Boundaries (-v.x, "velocity");
      else if (i == 1) D -> u[1] -> addToM0Boundaries (-v.y, "velocity");
      Solve (D -> u[i], Uf[i][0], MMS[i], D -> step, nOrder);
    }

    // -- Find forces with new domain fields.

    B -> force (*D);

    // -- Process results of this step.

    A -> analyse (Us);
  }
}


static void nonLinear (Domain*     D ,
		       AuxField*** Us,
		       AuxField*** Uf,
		       Vector&     a )
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u) - a.
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
// ---------------------------------------------------------------------------
{
  int               i, j;
  const int         nZ     = Geometry::nZ();
  const int         nP     = Geometry::planeSize();
  const int         nTot   = nZ * nP;
  const int         nZ32   = (3 * nZ) >> 1;
  const int         nTot32 = nZ32 * nP;
  const real        EPS    = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;
  vector<real>      work ((2 * DIM + 1) * nTot32);
  vector<real*>     u32 (DIM);
  vector<real*>     n32 (DIM);
  vector<AuxField*> U (DIM);
  vector<AuxField*> N (DIM);
  vector<real>      A (DIM);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * DIM * nTot32;

  A[0] = a.x; A[1] = a.y; if (DIM == 3) A[2] = a.z;

  for (i = 0; i < DIM; i++) {
    u32[i] = work() +  i        * nTot32;
    n32[i] = work() + (i + DIM) * nTot32;

    AuxField::swapData (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    U[i] -> transform32  (u32[i], -1);

    N[i] = Uf[i][0];
    Veclib::zero (nTot32, n32[i],  1);
  }
  
  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      
      // -- Perform n_i = u_j d(u_i) / dx_j.

      Veclib::copy (nTot32, u32[i], 1, tmp,  1);
      if (j == 2) {
	Femlib::DFTr (tmp, nZ32, nP, +1);
	Veclib::zero (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, tmp, j);
	Femlib::DFTr (tmp, nZ32, nP, -1);
      } else {
	master -> gradient (nZ32, tmp, j);
      }
      Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

      // -- Perform n_i += d(u_i u_j) / dx_j.

      Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
      if (j == 2) {
	Femlib::DFTr (tmp, nZ32, nP, +1);
	Veclib::zero (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, tmp, j);
	Femlib::DFTr (tmp, nZ32, nP, -1);
      } else {
	master -> gradient (nZ32, tmp, j);
      }
      Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

    }
    N[i]   -> transform32 (n32[i], +1);
    master -> smooth (N[i]);
    if (fabs (A[i]) > EPS) N[i] -> addToPlane (0, 2.0*A[i]);
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
  int               i, q, Je;
  vector<AuxField*> H (DIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  Je = (int) Femlib::value ("N_TIME");
  Je = min (D -> step, Je);

  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(), 1);

  for (i = 0; i < DIM; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[i][q]);
      H[i] -> axpy ( beta [q]    , *Uf[i][q]);
    }
}


static void setPForce (const AuxField*** Us,
		       AuxField***       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  int        i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {
   *Uf[i][0] = *Us[i][0];
    Uf[i][0] -> gradient (i);
  }
  
  for (i = 1; i < DIM; i++) *Uf[0][0] += *Uf[i][0];

  *Uf[0][0] /= dt;
}


static void project (const Domain* D ,
		     AuxField***   Us,
		     AuxField***   Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in lowest level of Us.  Constrain velocity field:
//                    u^^ = u^ - D_T * grad P;
// u^^ is left in lowest level of Uf.
//
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to create
// forcing for viscous step.
// ---------------------------------------------------------------------------
{
  int        i;
  const real dt    = Femlib::value ("D_T");
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM];
    Uf[i][0] -> gradient (i);
  
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Field::swapData (Us[i][0], Uf[i][0]);

    *Uf[i][0] *= alpha;
  }
}


static ModeSys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE == 1 selects iterative solvers for velocity components,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  char                 name;
  int                  i;
  const int            nSys   = D -> Nsys.getSize();
  const int            nZ     = Geometry::nZ();
  const int            nModes = (nZ + 1) >> 1;
  const int            base   = 0;
  const int            itLev  = (int) Femlib::value ("ITERATIVE");
  const int            nOrder = (int) Femlib::value ("N_TIME");
  const real           beta   = Femlib::value ("BETA");
  ModeSys**            M      = new ModeSys* [DIM + 1];
  vector<Element*>&    E      = ((Domain*) D) -> Esys;
  const NumberSystem** N      = new const NumberSystem* [3];

  // -- Velocity systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real   lambda2 = alpha[0] / Femlib::value ("KINVIS*D_T");

    for (i = 0; i < DIM; i++) {
      name = D -> u[i] -> name();
      D -> setNumber (name, N);
      M[i] = new ModalMatrixSystem (lambda2, beta, name, base, nModes, E, N);
    }
  } else
    for (i = 0; i < DIM; i++) M[i] = 0;

  // -- Pressure system.

  if (itLev < 2) {
    name = D -> u[DIM] -> name();
    D -> setNumber (name, N);
    M[DIM] = new ModalMatrixSystem (0.0, beta, name, base, nModes, E, N);
  } else
    M[DIM] = 0;

  return M;
}


static void Solve (Field*    U     ,
		   AuxField* Force ,
		   ModeSys*  M     ,
		   const int step  ,
		   const int nOrder)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  char routine[] = "Solve";

  const int  iterative = M == 0;
  const char name      = U -> name();
  const int  velocity  = name == 'u' || name == 'v' || name == 'w';
  const int  pressure  = name == 'p';

  if (!(velocity || pressure))
    message (routine, "input field type not recognized", ERROR);

  if (pressure) {
    if   (iterative) U -> solve (Force, 0.0);
    else             U -> solve (Force, M);
    return;
  }

  if (velocity) {
    if (iterative || step < nOrder) {
      const int    Je = min (step, nOrder);
      vector<real> alpha (Je + 1);
      Integration::StifflyStable (Je, alpha());
      const real   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
      
      U -> solve (Force, lambda2);

    } else
      U -> solve (Force, M);

    return;
  }
}
