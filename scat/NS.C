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

#include <buoy.h>

typedef ModalMatrixSystem ModeSys;
static  integer           DIM;

static void  nonLinear (Domain*, AuxField***, AuxField***);
static void  buoyancy  (Domain*, AuxField***, AuxField***, vector<real>&);
static void  waveProp  (Domain*, const AuxField***, const AuxField***);
static void  setPForce (const AuxField***, AuxField***);
static void  project   (const Domain*, AuxField***, AuxField***);

static ModeSys** preSolve (const Domain*);
static void      Solve    (Field*, AuxField*, ModeSys*,
			   const integer, const integer);


void NavierStokes (Domain*   D,
		   Analyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// temperature Field 'c', followed by constraint/pressure field 'p'.
//
// Us is multi-level auxillary Field storage for velocities+temperature and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  DIM = Geometry::nDim();

  integer       i, j;
  const real    dt     =           Femlib::value ("D_T");
  const integer nOrder = (integer) Femlib::value ("N_TIME");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();

  // -- Set up gravity vector.  Note directional components g_1, g_2, g_3.
  
  vector<real> g(3);
  {
    const real norm = Femlib::value ("sqrt (g_1*g_1 + g_2*g_2 + g_3*g_3)");

    g[0] = Femlib::value ("g_1 * GRAVITY") / norm;
    g[1] = Femlib::value ("g_2 * GRAVITY") / norm;
    g[2] = Femlib::value ("g_3 * GRAVITY") / norm;
  }

  // -- Create global matrix systems.

  ModeSys** MMS = preSolve (D);

  // -- Create & initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [(size_t) (DIM + 1)];
  AuxField*** Uf = new AuxField** [(size_t) (DIM + 1)];

  for (i = 0; i <= DIM; i++) {
    Us[i] = new AuxField* [(size_t) nOrder];
    Uf[i] = new AuxField* [(size_t) nOrder];
    for (j = 0; j < nOrder; j++) {
      *(Us[i][j] = new AuxField (D -> Esys, nZ)) = 0.0;
      *(Uf[i][j] = new AuxField (D -> Esys, nZ)) = 0.0;
    }
  }

  // -- Set up multi-level storage for pressure BCS.

  Field* Pressure = D -> u[DIM + 1];
  PBCmgr::build (Pressure);

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
    setPForce ((const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM + 1], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Viscous and thermal diffusion substep.

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
  integer i, j;

#if defined(STOKES)

  for (i = 0; i <= DIM; i++) *Uf[i][0] = 0.0;

#else

  const integer     nZ     = Geometry::nZ();
  const integer     nZP    = Geometry::nZProc();
  const integer     nP     = Geometry::planeSize();
  const integer     nPP    = Geometry::nBlock();
  const integer     nPR    = Geometry::nProc();
  const integer     nTot   = nZP * nP;
  const integer     nZ32   = (nPR > 1) ? nZP : (3 * nZP) >> 1;
  const integer     nTot32 = nZ32 * nP;
  vector<real>      work ((2 * (DIM + 1) + 1) * nTot32);
  vector<real*>     u32 (DIM + 1);
  vector<real*>     n32 (DIM + 1);
  vector<AuxField*> U   (DIM + 1);
  vector<AuxField*> N   (DIM + 1);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * (DIM + 1) * nTot32;

  for (i = 0; i <= DIM; i++) {
    u32[i] = work() +  i            * nTot32;
    n32[i] = work() + (i + DIM + 1) * nTot32;

    AuxField::swapData (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    U[i] -> transform32  (u32[i], -1);

    N[i] = Uf[i][0];
    Veclib::zero (nTot32, n32[i],  1);
  }

  for (i = 0; i <= DIM; i++) {
    for (j = 0; j < DIM; j++) {
      
      // -- Perform n_i += u_j d(u_i) / dx_j.

      Veclib::copy (nTot32, u32[i], 1, tmp,  1);
      if (j == 2) {
	Femlib::transpose  (tmp, nZ32,        nP, +1);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ,  nPP, tmp, j);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	Femlib::transpose  (tmp, nZ32,        nP, -1);
      } else {
	master -> gradient (nZ32, nP, tmp, j);
      }
      Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

      // -- Perform n_i += d(u_i u_j) / dx_j.

      Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
      if (j == 2) {
	Femlib::transpose  (tmp, nZ32,        nP, +1);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ,  nPP, tmp, j);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	Femlib::transpose  (tmp, nZ32,        nP, -1);
      } else {
	master -> gradient (nZ32, nP, tmp, j);
      }
      Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
	
    }
    N[i]   -> transform32 (n32[i], +1);
    master -> smooth (N[i]);
    *N[i] *= -0.5;
  }

#endif
}


static void buoyancy (Domain*       D ,
		      AuxField***   Us,
		      AuxField***   Uf,
		      vector<real>& g )
// ---------------------------------------------------------------------------
// The buoyancy term in the momentum equation is
//
//                      - BETA_T * (T - T_REF) g.
//                                             ~
// We add an explicit estimate of this to the nonlinear terms in Uf.
// The first level of Us has the last values of the data Fields, D is free.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real EPS  = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  AuxField*  T    = Us [DIM][0];
  AuxField*  work = D -> u[DIM];

  for (i = 0; i < DIM; i++) {
    if (fabs (g[i]) > EPS) {
      *work      = *T;
      *work     -=  Femlib::value ("T_REF" );
      *work     *=  Femlib::value ("BETA_T");
      *work     *=  g[i];
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
  integer           i, q;
  vector<AuxField*> H (DIM + 1);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i <= DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  integer Je = (integer) Femlib::value ("N_TIME");
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


static void setPForce (const AuxField*** Us,
		       AuxField***       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  integer    i;
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
// In addition, the intermediate temperature Field is transferred from the
// lowest level of Us to Uf for subsequent operations.
//
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to create
// forcing for viscous step.
//
// A similar operation occurs for the intermediate temperature Field, except
// that the coefficient of thermal diffusivity is used in place of kinematic
// viscosity.
// ---------------------------------------------------------------------------
{
  int        i;
  const real dt = Femlib::value ("D_T");
  real       alpha;

  alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM + 1];
    Uf[i][0] -> gradient (i);
  
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    AuxField::swapData (Us[i][0], Uf[i][0]);

    *Uf[i][0] *= alpha;
  }

  alpha = -1.0 / Femlib::value ("D_T * KINVIS / PRANDTL");
  AuxField::swapData (Us[DIM][0], Uf[DIM][0]);
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
  char                 name;
  integer              i;
  const integer        nSys   = D -> Nsys.getSize();
  const integer        nmodes = Geometry::nModeProc();
  const integer        base   = Geometry::baseMode();
  const integer        itLev  = (integer) Femlib::value ("ITERATIVE");
  const integer        nOrder = (integer) Femlib::value ("N_TIME");
  const real           beta   = Femlib::value ("BETA");
  ModeSys**            M      = new ModeSys* [(size_t) (DIM + 2)];
  vector<Element*>&    E      = ((Domain*) D) -> Esys;
  const NumberSystem** N      = new const NumberSystem* [(size_t) 3];

  // -- Velocity and temperature systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());

    // -- Velocities.

    real lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

    for (i = 0; i < DIM; i++) {
      name = D -> u[i] -> name();
      D -> setNumber (name, N);
      M[i] = new ModalMatrixSystem (lambda2, beta, name, base, nmodes, E, N);
    }

    // -- Temperature.

    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");

    name    = D -> u[DIM] -> name();
    D -> setNumber (name, N);
    M[DIM]  = new ModalMatrixSystem (lambda2, beta, name, base, nmodes, E, N);

  } else
    for (i = 0; i <= DIM; i++) M[i] = 0;

  // -- Pressure system.

  if (itLev < 2) {
    name    = D -> u[DIM + 1] -> name();
    D -> setNumber (name, N);
    M[DIM + 1] = new ModalMatrixSystem (0.0, beta, name, base, nmodes, E, N);
  } else
    M[DIM + 1] = 0;

  return M;
}


static void Solve (Field*        U     ,
		   AuxField*     Force ,
		   ModeSys*      M     ,
		   const integer step  ,
		   const integer nOrder)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  char routine[] = "Solve";

  const integer  iterative   = M == 0;
  const char name            = U -> name();
  const integer  velocity    = name == 'u' || name == 'v' || name == 'w';
  const integer  temperature = name == 'c';
  const integer  pressure    = name == 'p';
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
      const integer    Je = min (step, nOrder);
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
