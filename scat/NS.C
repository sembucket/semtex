//////////////////////////////////////////////////////////////////////////////
// NS.C:  Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting
// methods for the incompressible Navier--Stokes equations", JCP 9(2).
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <NS.h>

#ifdef __DECCXX
  #pragma define_template roll<Field*>
  #pragma define_template min<int>
#endif

static void  nonLinear (Domain*,  Field***, Field***, Vector);
static void  waveProp  (Domain*,  Field***, Field***);
static void  setPForce (Field***, Field***);
static void  project   (const Domain*,  Field***, Field***);
static void  setUForce (Domain*,  Field***);

static void  preSolve  (SystemField*, const real&);
static void  Solve     (SystemField*, Field*, const int&, const int&);


void NavierStokes (Domain* D, Analyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity SystemFields 'u', 'v' and
// constraint SystemField 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  const real  dt     = dparam ("DELTAT");
  const int   nOrder = iparam ("N_TIME");
  const int   nStep  = iparam ("N_STEP");
  const int   DIM    = iparam ("N_VAR" );

  Vector  a  = {0.0, 0.0, 0.0};   // -- Frame acceleration for N--S.

  // -- Set up multi-level storage for velocities and forcing.
  
  Field***  Us = new Field** [DIM];
  Field***  Uf = new Field** [DIM];

  for (int i = 0; i < DIM; i++) {
    Us[i] = new Field* [nOrder];
    Uf[i] = new Field* [nOrder];
     for (int j = 0; j < nOrder; j++) {
      Us[i][j] = new Field (*D -> u[i]);
      Uf[i][j] = new Field (*D -> u[i]);
    }
  }

  // -- Set up to solve velocity viscous step.

  real* alpha = rvector (Integration::OrderMax + 1);
  Integration::StifflyStable (nOrder, alpha);
  const real lambda2 = alpha[0] / (dt * dparam ("KINVIS"));
  freeVector (alpha);  

  preSolve (D -> u[0], lambda2);

  for (i = 1; i < DIM; i++)
    D -> u[i] -> assemble (*D -> u[0]);

  // -- Set up pressure solution matrices.

  SystemField* Pressure = D -> u[DIM];
  preSolve (Pressure, 0.0);

  // -- Timestepping loop.

  while (D -> step () < nStep) {
 
    D -> step () += 1; 
    D -> time () += dt;
    setDparam ("t", D -> time ());

    // -- Unconstrained forcing substep.

    nonLinear (D, Us, Uf, a );
    waveProp  (D, Us, Uf);

    // -- Pressure projection substep.

    PBCmanager::maintain (D -> step (), *Pressure, Us, Uf);
    Pressure -> evaluateBoundaries (D -> step ());
    for (i = 0; i < DIM; i++) {
      Field::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce (Us, Uf);
    Solve (Pressure, Uf[0][0], D -> step (), nOrder);
    project (D, Us, Uf);

    // -- Viscous correction substep.

    setUForce (D, Uf);
    for (i = 0; i < DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll (Us[i], nOrder);
      D -> u[i] -> evaluateBoundaries (D -> step ());
      Solve (D -> u[i], Uf[i][0], D -> step (), nOrder);
    }

    // -- Process results of this step.

    A -> analyse ();
  }
}


static void nonLinear (Domain*   D ,
		       Field***  Us,
		       Field***  Uf,
		       Vector    a )
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u) - a.
//
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms are computed in skew-symmetric form (Zang 1991)
//                    -0.5 ( u . grad u + div uu ).
//                           ~        ~       ~~
// ---------------------------------------------------------------------------
{
  Field::swapData (D -> u[0], Us[0][0]);
  Field::swapData (D -> u[1], Us[1][0]);

  // -- Build nonlinear terms.

  const Field& Ux = *Us[0][0];
  const Field& Uy = *Us[1][0];

  Field& Tp = *D -> u[0];
  
  Field& Nx = *Uf[0][0];
  Field& Ny = *Uf[1][0];

#ifdef STOKES  /* -- No nonlinear terms. */
  
  Nx = 0.0;
  Ny = 0.0;

#else          /* -- Build skew-symmetric nonlinear terms. */

  // -- Conservative NL terms.

  Ny = Nx.prod (Ux, Uy);
  Nx.grad (0, 1);
  Ny.grad (1, 0);

  Tp . prod (Ux, Ux) . grad (1, 0);
  Nx += Tp;

  Tp . prod (Uy, Uy) . grad (0, 1);
  Ny += Tp;

  // -- Nonconservative NL terms.

  Tp = Ux;
  Nx . addprod (Ux, Tp . grad (1, 0));

  Tp = Ux;
  Nx . addprod (Uy, Tp . grad (0, 1));

  Tp = Uy;
  Ny . addprod (Ux, Tp . grad (1, 0));

  Tp = Uy;
  Ny . addprod (Uy, Tp . grad (0, 1));

  // -- Smooth result on domain velocity boundary system & scale.

  D -> u[0] -> smooth (&Nx);
  D -> u[0] -> smooth (&Ny);

  Nx *= -0.5;
  Ny *= -0.5;
  
#endif

  // -- Add in distributed forcing to complete construction of nonlinear terms.

  if (fabs (a.x) > EPSDP) Nx -= a.x;
  if (fabs (a.y) > EPSDP) Ny -= a.y;
}


static void  waveProp (Domain*   D ,
		       Field***  Us,
		       Field***  Uf)
// ---------------------------------------------------------------------------
// Compute the first substep of stiffly-stable timestepping scheme.
//
// On entry, the most recent velocity fields are in Us, and the most
// recent nonlinear terms in Uf.  The intermediate velocity field u^ is
// computed and left in D's velocity areas. 
// ---------------------------------------------------------------------------
{
  // -- Construct u^ in Hx & Hy.

  Field& Hx = *D -> u[0] = 0.0;
  Field& Hy = *D -> u[1] = 0.0;

  int  Je = iparam ("N_TIME");
  Je = min (D -> step (), Je);
  
  real* alpha = rvector (Integration::OrderMax + 1);
  real* beta  = rvector (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha);
  Integration::Extrapolation (Je, beta );
  Blas::scal (Je, dparam ("DELTAT"), beta,  1);

  for (int q = 0; q < Je; q++) {
    Hx . axpy (-alpha[q + 1], *Us[0][q]) . axpy (beta[q], *Uf[0][q]);
    Hy . axpy (-alpha[q + 1], *Us[1][q]) . axpy (beta[q], *Uf[1][q]);
  }

  freeVector (alpha);
  freeVector (beta );
}


static void setPForce (Field***  Us, Field***  Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first level of Us.
// Create div u^ / DELTAT in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  *Uf[0][0] = *Us[0][0];
  *Uf[1][0] = *Us[1][0];

  
  Uf[0][0] -> grad (1, 0);
  Uf[1][0] -> grad (0, 1);

  *Uf[0][0] += *Uf[1][0];
  *Uf[0][0] /= dparam ("DELTAT");
}


static void project (const Domain* D, Field*** Us, Field*** Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in lowest level of Us.  Constrain velocity field:
//                    u^^ = u^ - DELTAT * grad P;
// u^^ is left in lowest level of Uf.
// ---------------------------------------------------------------------------
{
  const int   DIM = D -> nField () - 1;
  const real  dt  = dparam ("DELTAT");

  for (int i = 0; i < DIM; i++)
    *Uf[i][0] = *D -> u[DIM];
  
  Uf[0][0] -> grad (1, 0);
  Uf[1][0] -> grad (0, 1);

  for (i = 0; i < DIM; i++) {
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Field::swapData (Us[i][0], Uf[i][0]);
  }
}


static void setUForce (Domain* D, Field*** Uf)
// ---------------------------------------------------------------------------
// On entry, intermediate velocity storage u^^ is in lowest levels of Us.
// Multiply by -1.0 / (DELTAT * KINVIS) to create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  const int   DIM   = D -> nField () - 1;
  const real  alpha = -1.0 / (dparam ("DELTAT") * dparam ("KINVIS"));

  for (int i = 0; i < DIM; i++)
    *Uf[i][0] *= alpha;
}


static void preSolve (SystemField* F, const real& lambda2)
// ---------------------------------------------------------------------------
// If F is selected for direct solution, form & factor matrices.
// ---------------------------------------------------------------------------
{
  char  routine[] = "preSolve";

  const char name      = F -> getName ();
  const int  iterative = option ("ITERATIVE");
  const int  velocity  = name == 'u' || name == 'v' || name == 'w';
  const int  pressure  = name == 'p';

  if (!(velocity || pressure))
    message (routine, "input field type not recognized", ERROR);

  if (velocity && !iterative) {
    message (routine, ": -- Building velocity matrices", REMARK);
    F -> assemble (lambda2);
    return;
  }

  if (pressure && iterative < 2) {
    message (routine, ": -- Building pressure matrices", REMARK);
    F -> assemble (lambda2);
    return;
  }
}


static void Solve (SystemField*  U     ,
		   Field*        F     ,
		   const int&    step  ,
		   const int&    nOrder)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  char routine[] = "Solve";

  const char name      = U -> getName ();
  const int  iterative = option ("ITERATIVE");
  const int  velocity  = name == 'u' || name == 'v' || name == 'w';
  const int  pressure  = name == 'p';

  if (!(velocity || pressure))
    message (routine, "input field type not recognized", ERROR);

  if (pressure) {
    if   (iterative > 1) U -> solve (F, 0.0);
    else                 U -> solve (F);
    return;
  }

  if (velocity) {
    if (iterative || step < nOrder) {
      const int   Je     = min (step, nOrder);
      real*       alpha  = rvector (Je + 1);
      const real  dt     = dparam  ("DELTAT");
      const real  KinVis = dparam  ("KINVIS");

      Integration::StifflyStable (Je, alpha);
      const real lambda2 = alpha[0] / (dt * KinVis);
      freeVector (alpha);
      
      U -> solve (F, lambda2);

    } else
      U -> solve (F);

    return;
  }
}
