///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting
// methods for the incompressible Navier--Stokes equations", JCP 9(2).
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <NS.h>


#ifdef __DECCXX
  #pramga define_template vector<real>
  #pragma define_template vector<Field*>
  #pragma define_template matrix<Field*>
  #pragma define_template roll<Field*>
  #pragma define_template min<int>
#endif


static void  nonLinear (Domain*, matrix<Field*>&, matrix<Field*>&, Vector&);
static void  waveProp  (Domain*, const matrix<Field*>&, const matrix<Field*>&);
static void  setPForce (const Domain*, const matrix<Field*>&, matrix<Field*>&);
static void  project   (const Domain*, matrix<Field*>&, matrix<Field*>&);
static void  setUForce (const Domain*, matrix<Field*>&);

static void  preSolve  (SystemField*, const real&);
static void  Solve     (SystemField*, Field*, const int&, const int&);


void NavierStokes (Domain*   D,
		   Analyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity SystemFields 'u', 'v' ('w') and
// constraint SystemField 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  const real dt     = Femlib::parameter ("DELTAT");
  const int  nOrder = Femlib::integer   ("N_TIME");
  const int  nStep  = Femlib::integer   ("N_STEP");
  const int  DIM    = D -> nField () - 1;

  Vector  a  = {0.0, 0.0, 0.0};   // -- Frame acceleration for N--S.

  // -- Set up multi-level storage for velocities and forcing.

  matrix<Field*> Us (DIM, nOrder);
  matrix<Field*> Uf (DIM, nOrder);

  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < nOrder; j++) {
      Us (i, j) = new Field (*D -> u[i]);
      Uf (i, j) = new Field (*D -> u[i]);
    }
  }

  // -- Set up to solve velocity viscous step.

  vector<real> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (nOrder, alpha ());
  const real lambda2 = alpha[0] / (dt * Femlib::parameter ("KINVIS"));
  delete (alpha);  

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

    nonLinear (D, Us, Uf, a);
    waveProp  (D, Us, Uf);

    // -- Pressure projection substep.

    PBCmanager::maintain (D -> step (), *Pressure, Us, Uf);
    Pressure -> evaluateBoundaries (D -> step ());
    for (i = 0; i < DIM; i++) {
      Field::swapData (D -> u[i], Us (i, 0));
      roll (Uf (i), nOrder);
    }
    setPForce (D, Us, Uf);
    Solve     (Pressure, Uf (0, 0), D -> step (), nOrder);
    project   (D, Us, Uf);

    // -- Viscous correction substep.

    setUForce (D, Uf);
    for (i = 0; i < DIM; i++) {
      *Us (i, 0) = *D -> u[i];
      roll (Us (i), nOrder);
      D -> u[i] -> evaluateBoundaries (D -> step ());
      Solve (D -> u[i], Uf (i, 0), D -> step (), nOrder);
    }

    // -- Process results of this step.

    A -> analyse ();
  }
}


static void nonLinear (Domain*         D ,
		       matrix<Field*>& Us,
		       matrix<Field*>& Uf,
		       Vector&         a )
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
//
// If STOKES is defined for compilation, the nonlinear terms are set to zero.
//
// Note that all gradient operations are performed on T = D -> u[0], and that
// an even number are performed in the Fourier direction, so that the
// frames return to their original places by the end of the subroutine.
// ---------------------------------------------------------------------------
{
  int        i, j;
  const int  DIM = D -> nField () - 1;
  const real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  vector<Field*> U (DIM);	// -- Shorthand for velocity  fields.
  vector<Field*> N (DIM);	// -- Shorthand for nonlinear fields.
  vector<real>   A (DIM);	// -- Vector form of frame acceleration.

  A[0] = a.x; A[1] = a.y; if (DIM == 3) A[2] = a.z;

  for (i = 0; i < DIM; i++) {
    Field::swapData (D -> u[i], Us (i, 0));
    U[i] = Us (i, 0);
    N[i] = Uf (i, 0);
  }

  for (i = 0; i < DIM; i++) *N[i] = 0.0;

#ifndef STOKES

  // -- Build skew-symmetric nonlinear terms.

  Field* T = D -> u[0];		// -- Workspace.

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      *T = *U[i];
      T    -> gradient (j);
      N[i] -> addprod  (*U[j], *T);

      T -> product  (*U[i], *U[j]);
      T -> gradient (j);
      *N[i] += *T;
    }
    T -> smooth (N[i]);
    *N[i] *= -0.5;
  }

#endif

  for (i = 0; i < DIM; i++) if (fabs (A[i]) > EPS) *N[i] -= A[i];
}


static void waveProp (Domain*               D ,
		      const matrix<Field*>& Us,
		      const matrix<Field*>& Uf)
// ---------------------------------------------------------------------------
// Compute the first substep of stiffly-stable timestepping scheme.
//
// On entry, the most recent velocity fields are in Us, and the most
// recent nonlinear terms in Uf.  The intermediate velocity field u^ is
// computed and left in D's velocity areas. 
// ---------------------------------------------------------------------------
{
  int            i;
  const int      DIM = D -> nField () - 1;
  vector<Field*> H (DIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  int  Je = Femlib::integer ("N_TIME");
  Je = min (D -> step (), Je);

  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha ());
  Integration::Extrapolation (Je, beta  ());
  Blas::scal (Je, Femlib::parameter ("DELTAT"), beta (),  1);

  for (int i = 0; i < DIM; i++)
    for (int q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us (i, q));
      H[i] -> axpy ( beta [q]    , *Uf (i, q));
    }
}


static void setPForce (const Domain*         D ,
		       const matrix<Field*>& Us,
		       matrix<Field*>&       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / DELTAT in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
//
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf as a result of gradient operation.
// ---------------------------------------------------------------------------
{
  int       i;
  const int DIM = D -> nField () - 1;

  for (i = 0; i < DIM; i++) {
   *Uf (i, 0) = *Us (i, 0);
    Uf (i, 0) -> gradient (i);
  }
  
  for (i = 1; i < DIM; i++) *Uf (0, 0) += *Uf (i, 0);

  *Uf (0, 0) /= Femlib::parameter ("DELTAT");
}


static void project (const Domain*   D ,
		     matrix<Field*>& Us,
		     matrix<Field*>& Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in lowest level of Us.  Constrain velocity field:
//                    u^^ = u^ - DELTAT * grad P;
// u^^ is left in lowest level of Uf.
//
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf.  This returns to original place the swapping done
// by setPForce.
// ---------------------------------------------------------------------------
{
  int        i;
  const int  DIM = D -> nField () - 1;
  const real dt  = Femlib::parameter ("DELTAT");

  for (int i = 0; i < DIM; i++) {

   *Uf (i, 0) = *D -> u[DIM];
    Uf (i, 0) -> gradient (i);
  
    Us (i, 0) -> axpy (-dt, *Uf (i, 0));
    Field::swapData (Us (i, 0), Uf (i, 0));
  }
}


static void setUForce (const Domain*   D ,
		       matrix<Field*>& Uf)
// ---------------------------------------------------------------------------
// On entry, intermediate velocity storage u^^ is in lowest levels of Uf.
// Multiply by -1.0 / (DELTAT * KINVIS) to create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  const int  DIM   = D -> nField () - 1;
  const real alpha = -1.0 /
    (Femlib::parameter ("DELTAT") * Femlib::parameter ("KINVIS"));

  for (int i = 0; i < DIM; i++) *Uf (i, 0) *= alpha;
}


static void preSolve (SystemField* F      ,
		      const real&  lambda2)
// ---------------------------------------------------------------------------
// If F is selected for direct solution, form & factor matrices.
// ---------------------------------------------------------------------------
{
  char  routine[] = "preSolve";

  const char name      = F -> getName ();
  const int  iterative = Femlib::option ("ITERATIVE");
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
  const int  iterative = Femlib::option ("ITERATIVE");
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
      const int    Je     = min (step, nOrder);
      vector<real> alpha (Je + 1);
      const real   dt     = Femlib::parameter ("DELTAT");
      const real   KinVis = Femlib::parameter ("KINVIS");

      Integration::StifflyStable (Je, alpha);
      const real lambda2 = alpha[0] / (dt * KinVis);
      
      U -> solve (F, lambda2);

    } else
      U -> solve (F);

    return;
  }
}
