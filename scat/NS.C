///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// This version includes the Boussinesq approximation for heat transport
// and buoyant convection, with temperature carried as the DIM+1th Field.
// Solved in Cartesian space only.
//
// Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting
// methods for the incompressible Navier--Stokes equations", JCP 9(2).
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <NS.h>

#ifdef __DECCXX
  #pragma define_template roll<AuxField*>
#endif

typedef ModalMatrixSystem ModeSys;
static  int               DIM;

static void  nonLinear (Domain*, AuxField***, AuxField***);
static void  buoyancy  (Domain*, AuxField***, AuxField***, const Vector&);
static void  waveProp  (Domain*, const AuxField***, const AuxField***);
static void  setPForce (const Domain*, const AuxField***, AuxField***);
static void  project   (const Domain*, AuxField***, AuxField***);
static void  setUForce (const Domain*, AuxField***);

static ModeSys** preSolve (const Domain*);
static void      Solve    (Field*, AuxField*, ModeSys*, const int, const int);


void NavierStokes (Domain*   D,
		   Analyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// temperature Field 'T', followed by constraint/pressure field 'p'.
//
// Us is multi-level auxillary Field storage for velocities+temperature and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  DIM = D -> nField() - 2;

  int        i, j;
  const real dt       =       Femlib::value ("D_T");
  const int  nOrder   = (int) Femlib::value ("N_TIME");
  const int  nStep    = (int) Femlib::value ("N_STEP");
  const int  nZ       = (int) Femlib::value ("N_Z");
  Vector     g        = {0.0, 0.0, 0.0};
  Field*     Pressure = D -> u[DIM + 1];
  ModeSys**  MMS      = preSolve (D);

  // -- Set up multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [DIM + 1];
  AuxField*** Uf = new AuxField** [DIM + 1];

  for (i = 0; i <= DIM; i++) {
    Us[i] = new AuxField* [nOrder];
    Uf[i] = new AuxField* [nOrder];
    for (j = 0; j < nOrder; j++) {
      Us[i][j] = new AuxField (D -> Esys, nZ);
      Uf[i][j] = new AuxField (D -> Esys, nZ);
    }
  }

  // -- Set up multi-level storage for pressure BCS.

  PBCmgr::build (Pressure);

  // -- Set up gravity vector.  Note directional components g_1, g_2, g_3.
  
  {
    const real norm = Femlib::value ("sqrt (g_1*g_1 + g_2*g_2 + g_3*g_3)");

    g.x = Femlib::value ("g_1 * GRAVITY") / norm;
    g.y = Femlib::value ("g_2 * GRAVITY") / norm;
    g.z = Femlib::value ("g_3 * GRAVITY") / norm;
  }

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us, Uf);
    buoyancy  (D, Us, Uf, g);
    waveProp  (D, (const AuxField***) Us, (const AuxField***) Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
		      (const AuxField***) Us,
		      (const AuxField***) Uf, 1);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i <= DIM; i++) {
      AuxField::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce (D, (const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM + 1], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Viscous and thermal diffusion substep.

    setUForce (D, Uf);
    for (i = 0; i <= DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll  (Us[i], nOrder);
      Solve (D -> u[i], Uf[i][0], MMS[i], D -> step, nOrder);
    }

    // -- Process results of this step.

    A -> analyse (Us);
  }
}


static void nonLinear (Domain*     D ,
		       AuxField*** Us,
		       AuxField*** Uf)
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u).
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
// If STOKES is defined for compilation, the nonlinear terms are set to zero.
// This means that the Navier--Stokes equations become the Stokes equations.
// Also, the convective terms in the temperature equation are zero,
// so that the transient heat diffucion equation is solved in place
// of the convection/diffusion equation (i.e. temperature is uncoupled
// from the velocity field).
//
// Note that all gradient operations are performed on T = D -> u[0], and that
// an even number are performed in the Fourier direction, so that the
// frames/planes return to their original places by the end of the subroutine.
//
// The temperature transport term u . grad T is also built this way here.
// ---------------------------------------------------------------------------
{
  int               i, j;
  vector<AuxField*> U (DIM + 1);	// -- Shorthand for velocity  fields.
  vector<AuxField*> N (DIM + 1);	// -- Shorthand for nonlinear fields.

  for (i = 0; i <= DIM; i++) {
    AuxField::swapData (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    N[i] = Uf[i][0];
  }

  for (i = 0; i <= DIM; i++) *N[i] = 0.0;

#ifndef STOKES

  // -- Build skew-symmetric nonlinear terms.

  AuxField* T      = D -> u[0];		// -- Workspace.
  Field*    master = D -> u[0];	        // -- Template too.

  for (i = 0; i <= DIM; i++) {
    for (j = 0; j < DIM; j++) {
      *T = *U[i];
      T    -> gradient (j);
      N[i] -> addprod  (*U[j], *T);

      T -> product  (*U[i], *U[j]);
      T -> gradient (j);
      *N[i] += *T;
    }
    master -> smooth (N[i]);
    *N[i] *= -0.5;
  }

#endif
}


static void buoyancy (Domain*       D ,
		      AuxField***   Us,
		      AuxField***   Uf,
		      const Vector& g )
// ---------------------------------------------------------------------------
// 
// The buoyancy term in the momentum equation is
//
//                      - BETA_T * (T - T_REF) g.
//                                             ~
// We add an explicit estimate of this to the nonlinear terms in Uf.
// The first level of Us has the last values of the data Fields, D is free.
// ---------------------------------------------------------------------------
{
  int          i;
  const real   EPS  = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  AuxField*    T    = Us [DIM][0];
  AuxField*    work = D -> u[DIM];
  vector<real> G (3);
  
  G[0] = g.x; G[1] = g.y; G[2] = g.z;

  for (i = 0; i < DIM; i++) {
    if (fabs (G[i]) > EPS) {
      *work      = *T;
      *work     -=  Femlib::value ("T_REF" );
      *work     *=  Femlib::value ("BETA_T");
      *work     *=  G[i];
      *Uf[i][0] -= *work;
    }
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
//
// The intermediate temperature Field is also created in the same loop.
// ---------------------------------------------------------------------------
{
  int               i, q;
  vector<AuxField*> H (DIM + 1);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i <= DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  int Je = (int) Femlib::value ("N_TIME");
  Je = min (D -> step, Je);

  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

  for (i = 0; i <= DIM; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[i][q]);
      H[i] -> axpy ( beta [q]    , *Uf[i][q]);
    }
}


static void setPForce (const Domain*     D ,
		       const AuxField*** Us,
		       AuxField***       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
//
// A plane-swapping operation takes place in the first time level of the
// Fourier direction of Uf as a result of gradient operation.
// ---------------------------------------------------------------------------
{
  int         i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {
   *Uf[i][0] = *Us[i][0];
    Uf[i][0] -> gradient (i);
  }
  
  for (i = 1; i < DIM; i++)
    *Uf[0][0] += *Uf[i][0];

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
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf.  This returns to original place the swapping done
// by setPForce.
//
// In addition, the intermediate temperature Field is transferred from the
// lowest level of Us to Uf for subsequent operations.
// ---------------------------------------------------------------------------
{
  int        i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM + 1];
    Uf[i][0] -> gradient (i);
  
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    AuxField::swapData (Us[i][0], Uf[i][0]);
  }

  AuxField::swapData (Us[DIM][0], Uf[DIM][0]);
}


static void setUForce (const Domain* D ,
		       AuxField***   Uf)
// ---------------------------------------------------------------------------
// On entry, intermediate velocity storage u^^ is in lowest levels of Uf.
// Multiply by -1.0 / (D_T * KINVIS) to create forcing for diffusion step.
//
// A similar operation occurs for the intermediate temperature Field, except
// that the coefficient of thermal diffusivity is used in place of kinematic
// viscosity.
// ---------------------------------------------------------------------------
{
  int  i;
  real alpha;

  alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  for (i = 0; i < DIM; i++)
    *Uf[i][0] *= alpha;

  alpha = -1.0 / Femlib::value ("D_T * KINVIS / PRANDTL");
  *Uf[DIM][0] *= alpha;
}


static ModeSys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE == 1 selects iterative solvers for velocities and temperature,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  int                 i, j, found;
  const int           nSys   = D -> Nsys.getSize();
  const int           nZ     = D -> u[0] -> nZ();
  const int           nModes = (nZ + 1) >> 1;
  const int           itLev  = (int) Femlib::value ("ITERATIVE");
  const int           nOrder = (int) Femlib::value ("N_TIME");
  const real          beta   = Femlib::value ("BETA");
  ModeSys**           M      = new ModeSys* [DIM + 2];
  vector<Element*>&   E      = ((Domain*) D) -> Esys;
  const NumberSystem* N;

  // -- Velocity and temperature systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());

    // -- Velocities.

    real lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    N            = D -> u[0] -> system();
    M[0]         = new ModalMatrixSystem (lambda2, beta, nModes, E, N);

    for (i = 1; i < DIM; i++) {
      for (found = 0, j = 0; j < i; j++)
	if (found = D -> u[j] -> system() == D -> u[i] -> system()) break;
      
      if (found)
	M[i] = M[j];
      else {
	N    = D -> u[i] -> system();
	M[i] = new ModalMatrixSystem (lambda2, beta, nModes, E, N);
      }
    }

    // -- Temperature.

    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");
    N       = D -> u[DIM] -> system();
    M[DIM]  = new ModalMatrixSystem (lambda2, beta, nModes, E, N);

  } else
    for (i = 0; i <= DIM; i++) M[i] = 0;

  // -- Pressure system.

  if (itLev < 2) {
    N          = D -> u[DIM + 1] -> system();
    M[DIM + 1] = new ModalMatrixSystem (0.0, beta, nModes, E, N);

  } else
    M[DIM + 1] = 0;

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

  const int  iterative   = M == 0;
  const char name        = U -> name();
  const int  velocity    = name == 'u' || name == 'v' || name == 'w';
  const int  temperature = name == 'T';
  const int  pressure    = name == 'p';
  real       lambda2;

  if (!(velocity || pressure || temperature))
    message (routine, "input field type not recognized", ERROR);

  if (pressure) {
    lambda2 = 0.0;
    if   (iterative) U -> solve (Force, lambda2);
    else             U -> solve (Force, M);
    return;
  }

  if (velocity || temperature) {
    if (iterative || step < nOrder) {
      const int    Je = min (step, nOrder);
      vector<real> alpha (Je + 1);
      Integration::StifflyStable (Je, alpha());

      lambda2 = alpha[0];
      if   (velocity) lambda2 /= Femlib::value ("D_T * KINVIS");
      else            lambda2 /= Femlib::value ("D_T * KINVIS / PRANDTL");
      
      U -> solve (Force, lambda2);

    } else
      U -> solve (Force, M);

    return;
  }
}
