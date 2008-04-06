///////////////////////////////////////////////////////////////////////////////
// integrate.C: integrate unsteady linearised Navier--Stokes problem
// forward in time (or its adjoint, backwards in time).
//
// Copyright (c) 2000 <--> $Date$, Hugh Blackburn
//
// This version implements linearised forward and adjoint advection
// terms and evolves a single Fourier mode.  Both the number of
// velocity components in the perturbation velocity field
// (Geometry::nPert) and their complex scalar structure is variable,
// partly dependent on the number of velocity components used for the
// 2D base flow (Geometry::nBase).
//
// For cylindrical coordinates:
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <stab.h>

static int_t              NORD, NPERT, NBASE, NZ, CYL, PROB;
static vector<MatrixSys*> MS;

static void        waveProp   (Domain*, const AuxField***, const AuxField***);
static void        setPForce  (const AuxField**, AuxField**);
static void        project    (const Domain*, AuxField**, AuxField**);
static MatrixSys** preSolve   (const Domain*);
static void        Solve      (Domain*, const int_t, AuxField*, MatrixSys*);


void integrate (void            (*Advection)(Domain*, AuxField**, AuxField**),
		Domain*         D   ,
		StabAnalyser*   A   )
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// For stability code, also contains base velocity fields 'U', 'V'
// (and 'W').  This routine can be called repeatedly, so initialisation
// only occurs once.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NPERT = Geometry::nPert();
  NBASE = Geometry::nBase();
  NZ    = Geometry::nZ();
  PROB  = Geometry::problem();
  NORD  = Femlib::ivalue ("N_TIME");

  int_t        i, j, k;
  const real_t dt       = Femlib::value  ("D_T");
  const real_t period   = Femlib::value  ("BASE_PERIOD");
  const real_t tstart   = Femlib::value  ("T_OFFSET");
  const int_t  nStep    = Femlib::ivalue ("N_STEP");
  const bool   forwards = (Advection == linAdvect) ? true : false;

  static MatrixSys** MS;
  static AuxField*** Us;
  static AuxField*** Uf;
  static Field*      Pressure = D -> u[NPERT];

  if (!MS) {	      // -- Initialise static data (enable call-back).
    
    // -- Create global matrix systems

    MS = preSolve (D);
    
    // -- Create multi-level storage for velocities and forcing.

    const int_t ntot  = Geometry::nTotProc();
    real_t*     alloc = new real_t [static_cast<size_t>(2*NPERT*NORD*ntot)];

    Us = new AuxField** [static_cast<size_t>(NORD)];
    Uf = new AuxField** [static_cast<size_t>(NORD)];
    
    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [static_cast<size_t>(NPERT)];
      Uf[i] = new AuxField* [static_cast<size_t>(NPERT)];
      for (j = 0; j < NPERT; j++) {
	Us[i][j] = new AuxField (alloc + k++ * ntot, NZ, D -> elmt);
	Uf[i][j] = new AuxField (alloc + k++ * ntot, NZ, D -> elmt);
      }
    }

    // -- Create multi-level storage for pressure BCS.

    PBCmgr::build (Pressure);

    // -- Apply coupling to radial & azimuthal velocity BCs.

    if (Geometry::cylindrical() && 
	(PROB != Geometry::O2_2D || PROB != Geometry::SO2_2D))
      Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);
  }
    
  // -- Because we have to restart from scratch on each call, zero these:

  for (i = 0; i < NORD; i++)
    for (j = 0; j < NPERT; j++) {
      *Us[i][j] = 0.0;
      *Uf[i][j] = 0.0;
    }

  // -- Timestepping loop.

  D -> step = 0;
  D -> time = 0.0;

  // -- If base flow is periodic in time, start at the correct phase point.

  if (period > EPSDP) D -> time = (forwards) ? tstart : tstart + dt*nStep;

  while (D -> step < nStep) {

    // -- Reconstruct base velocity fields if appropriate.

    D -> updateBase();

    // -- Set the time to end of timestep.

    D -> step += 1; 
    D -> time += (forwards) ? dt : -dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.
    
    Advection (D, Us[0], Uf[0]);

    // -- Pressure substep.

    PBCmgr::maintain (D -> step, Pressure, 
		      const_cast<const AuxField**>(Us[0]),
		      const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (D -> step);

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp (D, const_cast<const AuxField***>(Us),
	         const_cast<const AuxField***>(Uf));

    for (i = 0; i < NPERT; i++) AuxField::swapData (D -> u[i], Us[0][i]);

    rollm     (Uf, NORD, NPERT);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NPERT,  Uf[0][0], MS[NPERT]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NPERT; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NPERT);

    // -- Viscous correction substep.

    if (Geometry::cylindrical() && 
	(PROB != Geometry::O2_2D || PROB != Geometry::SO2_2D)) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NPERT; i++) {
      ROOTONLY D -> u[i] -> evaluateM0Boundaries (D -> step);
      Solve (D, i, Uf[0][i], MS[i]);
    }
    if (Geometry::cylindrical() &&
	(PROB != Geometry::O2_2D || PROB != Geometry::SO2_2D))
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.
    
    A -> analyse (Us[0]);
  }
}


void linAdvect (Domain*    D ,
		AuxField** Us,
		AuxField** Uf)
// ---------------------------------------------------------------------------
// Compute linearised (forcing) terms in Navier--Stokes equations: N(u).
//
// Here N(u) represents the linearised advection terms in the N--S equations
// transposed to the RHS.
//
// Velocity field data areas of D and first level of Us are swapped,
// then the next stage of nonlinear forcing terms N(u) are computed
// from velocity fields and left in the first level of Uf.
//
// Linearised terms N(u) are computed in convective form
//                 
//           N  = - ( U.grad u + u.grad U )
//           ~        ~      ~   ~      ~
// The data are taken as being in Fourier space, but, as there are
// only two modes involved (the base flow and the perturbation) the
// convolution sums end up being 2D operations.
// 
// Things are made a little more complicated by the fact that while
// the base flow U is always purely real, with only one plane of data,
// the perturbation field u can have either one or two planes of data.
// To deal with this, Auxfield operations times and timesPlus
// (actually convolutions) are assumed to have a purely real Auxfield
// as the second operand.  Assignment and Gradient operators are also
// modified. And this difference in structure is the reason behind
// using two lots of temporary storage, T, and U[NBASE].
// ---------------------------------------------------------------------------
{
  int_t             i, j;
  vector<AuxField*> U(NBASE + 1), u(NPERT), N(NPERT);
  Field*            T = D -> u[0];

  // -- Set up local aliases.

  for (i = 0; i < NBASE + 1; i++)
    U[i] = D -> U[i];

  for (i = 0; i < NPERT; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
     u[i] = Us[i];
     N[i] = Uf[i];
    *N[i] = 0.0;
  }

  // -- Centrifugal, Coriolis terms for cylindrical coords.

  if (Geometry::cylindrical()) {
    if (NPERT == 3)
     *N[2] += T -> times (*u[2], *U[1]) . divY();
    if (NBASE == 3) {
      N[1] -> axpy (-2.0, T -> times (*u[2], *U[2]));
     *N[2] += T -> times (*u[1], *U[2]) . divY();
    }
  }

  // -- N_i += U_j d(u_i) / dx_j.

  for (i = 0; i < NPERT; i++)
    for (j = 0; j < NBASE; j++) {
      (*T = *u[i]) . gradient (j);
      if      (Geometry::cylindrical() && i <  2 && j <  2) T -> mulY();
      else if (Geometry::cylindrical() && i == 2 && j == 2) T -> divY();
      N[i] -> timesPlus (*T, *U[j]);
    }

  // -- N_i += u_j d(U_i) / dx_j; dU_i/dz=0.

  for (i = 0; i < NBASE; i++)
    for (j = 0; j < 2; j++) {
      (*U[NBASE] = *U[i]) . gradient (j);
      if (Geometry::cylindrical() && i < 2) U[NBASE] -> mulY();
      N[i] -> timesPlus (*u[j], *U[NBASE]);
    }

  for (i = 0; i < NPERT; i++) {
    T -> smooth (N[i]);
    *N[i] *= -1.0;
  }
}


void linAdvectT (Domain*    D ,
		 AuxField** Us,
		 AuxField** Uf)
// ---------------------------------------------------------------------------
// Compute ADJOINT linearised terms in Navier--Stokes equations: M(u).
//
// Adjoint terms M(u) are (NB sign change and transpose)
//
//           M  = - ( grad U.u - U.grad u )
//           ~             ~ ~   ~      ~
// ---------------------------------------------------------------------------
{
  int_t             i, j;
  vector<AuxField*> U(NBASE + 1), u(NPERT), N(NPERT);
  Field*            T = D -> u[0];

  // -- Set up local aliases.

  for (i = 0; i < NBASE + 1; i++)
    U[i] = D -> U[i];

  for (i = 0; i < NPERT; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
     u[i] = Us[i];
     N[i] = Uf[i];
    *N[i] = 0.0;
  }

  // -- Centrifugal, Coriolis terms for cylindrical coords.

  if (Geometry::cylindrical()) {
    if (NPERT == 3)
     *N[2] += T -> times (*u[2], *U[1]) . divY();
    if (NBASE == 3) {    // -- If NBASE = 3 then also NPERT = 3.
      *N[1] += T -> times (*u[2], *U[2]);
      N[2] -> axpy (-2.0, T -> times (*u[1], *U[2]) . divY() );
    }
  }

  // -- N_i -= U_j d(u_i) / dx_j.

  for (i = 0; i < NPERT; i++)
    for (j = 0; j < NBASE; j++) {
      (*T = *u[i]) . gradient (j);
      if      (Geometry::cylindrical() && i <  2 && j <  2) T -> mulY();
      else if (Geometry::cylindrical() && i == 2 && j == 2) T -> divY();
      N[i] -> timesMinus (*T, *U[j]);
    }

  // -- N_i += u_j d(U_j) / dx_i; dU_j/dz=0.

  for (i = 0; i < 2; i++)
    for (j = 0; j < NBASE; j++) {
      (*U[NBASE] = *U[j]) . gradient (i);
      if (Geometry::cylindrical()) U[NBASE] -> mulY();
      N[i] -> timesPlus (*u[j], *U[NBASE]);
    }

  for (i = 0; i < NPERT; i++) {
    T -> smooth (N[i]);
    *N[i] *= -1.0;
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
  int_t             i, q;
  vector<AuxField*> H (NPERT);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NPERT; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t    Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);
  Blas::scal (Je, Femlib::value ("D_T"), &beta[0], 1);

  for (i = 0; i < NPERT; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[q][i]);
      H[i] -> axpy ( beta [q]    , *Uf[q][i]);
    }
}


static void setPForce (const AuxField** Us,
		       AuxField**       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level
// of Us.  Create div u^ / D_T in the first component, first level
// storage of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t dt = Femlib::value ("D_T");

  for (i = 0; i < NPERT; i++) (*Uf[i] = *Us[i]) . gradient (i);
  if  (PROB == Geometry::O2_3D_SYMM) *Uf[2] *= -1.0;
  for (i = 1; i < NPERT; i++) *Uf[0] += *Uf[i];

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
  int_t        i;
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");

  for (i = 0; i < NPERT; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i == 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
  }

  for (i = 0; i < NPERT; i++) {
    (*Us[0] = *D -> u[NPERT]) . gradient (i);
    if (Geometry::cylindrical() && i < 2) Us[0] -> mulY();
    Uf[i] -> axpy (beta, *Us[0]);
  }
}


static MatrixSys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for system with only 1 Fourier mode.
// ---------------------------------------------------------------------------
{
  const int_t  mode  = (PROB == Geometry::O2_2D ||
			PROB == Geometry::SO2_2D)  ? 0 : 1;
  const int_t  bmode = mode * Femlib::ivalue ("BETA");
  const real_t beta  = mode * Femlib::value  ("BETA");
  const int_t  itLev = Femlib::ivalue ("ITERATIVE");

  vector<real_t> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  const real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  MatrixSys**      system = new MatrixSys* [static_cast<size_t>(NPERT + 1)];
  MatrixSys*       M;
  bool             found;
  solver_kind      method;
  real_t           betak2;
  const NumberSys* N;

  cout << "-- Installing matrices     : " << flush;

  method = (itLev < 1) ? DIRECT : JACPCG;

  // -- Velocities, starting with u.

  N      = D -> b[0] -> Nsys (bmode);
  betak2 = sqr (Field::modeConstant (D -> u[0] -> name(), mode, beta));
  M = new MatrixSys (lambda2, betak2, bmode, D -> elmt, D -> b[0], method);
  MS.insert (MS.end(), M);
  system[0] = M;
  cout << ((method == DIRECT) ? '*' : '&') << flush;

  vector<MatrixSys*>::iterator m;

  // -- v.

  N      = D -> b[1] -> Nsys (bmode);
  betak2 = sqr (Field::modeConstant (D -> u[1] -> name(), mode, beta));

  for (found = false, m = MS.begin(); !found && m != MS.end(); m++) {
    M = *m; found = M -> match (lambda2, betak2, N, method);
  }
  if (found) {
    system[1] = M;
    cout << "." << flush;
  } else {
    M = new MatrixSys (lambda2, betak2, bmode, D -> elmt, D -> b[1], method);
    MS.insert (MS.end(), M);
    system[1] = M;
    cout << ((method == DIRECT) ? '*' : '&') << flush;
  }

  // -- w.

  if (NPERT == 3) {
    N      = D -> b[2] -> Nsys (bmode);
    betak2 = sqr (Field::modeConstant (D -> u[2] -> name(), mode, beta));

    for (found = false, m = MS.begin(); !found && m != MS.end(); m++) {
      M = *m; found = M -> match (lambda2, betak2, N, method);
    }
    if (found) {
      system[2] = M;
      cout << "." << flush;
    } else {
      M = new MatrixSys (lambda2, betak2, bmode, D -> elmt, D -> b[2], method);
      MS.insert (MS.end(), M);
      system[2] = M;
      cout <<  ((method == DIRECT) ? '*' : '&') << flush;
    }
  }
    
  // -- Pressure.

  method = (itLev < 2) ? DIRECT : JACPCG;

  betak2 = sqr (Field::modeConstant (D -> u[NPERT] -> name(), mode, beta));
  M      = new MatrixSys (0.0, betak2,bmode, D -> elmt, D -> b[NPERT], method);
  MS.insert (MS.end(), M);
  system[NPERT] = M;
  cout << ((method == DIRECT) ? '*' : '&') << endl;
  
  return system;
}


static void Solve (Domain*     D,
		   const int_t i,
		   AuxField*   F,
		   MatrixSys*  M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order.
// ---------------------------------------------------------------------------
{
  const int_t step = D -> step;

  if (i < NPERT && step < NORD) {

    // -- We need a temporary matrix system for a viscous solve.

    const int_t     mode  = (PROB == Geometry::O2_2D ||
			     PROB == Geometry::SO2_2D ) ? 0 : 1;
    const int_t     bmode = mode * Femlib::ivalue ("BETA");
    const real_t    beta  = mode * Femlib::value  ("BETA");
    const int_t     Je    = min (step, NORD);    
    vector<real_t>  alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t betak2  = sqr(Field::modeConstant(D->u[i]->name(),mode,beta));
    const real_t lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

    MatrixSys* tmp =
      new MatrixSys (lambda2, betak2, bmode, D -> elmt, D -> b[i], JACPCG);

    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
