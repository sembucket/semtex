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
  #pragma define_template roll<AuxField*>
#endif

typedef ModalMatrixSystem ModeSys;
static  int               DIM;

static void  nonLinear (Domain*, AuxField***, AuxField***, Vector&);
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
// constraint Field 'p'.
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
  Vector     a        = {0.0, 0.0, 0.0};
  Field*     Pressure = D -> u[DIM];
  ModeSys**  MMS      = preSolve (D);

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
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us, Uf, a);
    waveProp  (D, (const AuxField***) Us, (const AuxField***) Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
		      (const AuxField***) Us, (const AuxField***) Uf, 1);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < DIM; i++) {
      AuxField::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce (D, (const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Viscous correction substep.

    setUForce (D, Uf);
    for (i = 0; i < DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll (Us[i], nOrder);
      Solve (D -> u[i], Uf[i][0], MMS[i], D -> step, nOrder);
    }

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
//              N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//               i           j    i      j      i j      j
//
// If STOKES is defined for compilation, the nonlinear terms are set to zero.
//
// Note that all gradient operations are performed on T = D -> u[0], and that
// an even number are performed in the Fourier direction, so that the
// frames return to their original places by the end of the subroutine.
// ---------------------------------------------------------------------------
{
  int               i, j;
  const real        EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  vector<AuxField*> U (DIM);	// -- Shorthand for velocity  fields.
  vector<AuxField*> N (DIM);	// -- Shorthand for nonlinear fields.
  vector<real>      A (DIM);	// -- Vector form of frame acceleration.

  A[0] = a.x; A[1] = a.y; if (DIM == 3) A[2] = a.z;

  for (i = 0; i < DIM; i++) {
    AuxField::swapData (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    N[i] = Uf[i][0];
  }

  for (i = 0; i < DIM; i++) *N[i] = 0.0;

#ifndef STOKES

  // -- Build skew-symmetric nonlinear terms.

  AuxField* T      = D -> u[0];		// -- Workspace.
  Field*    master = D -> u[0];	        // -- Template too.

  for (i = 0; i < DIM; i++) {
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

  for (i = 0; i < DIM; i++) if (fabs (A[i]) > EPS) *N[i] -= A[i];
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
  int               i, q;
  vector<AuxField*> H (DIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  int  Je = (int) Femlib::value ("N_TIME");
  Je = min (D -> step, Je);

  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

  for (i = 0; i < DIM; i++)
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
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf as a result of gradient operation.
// ---------------------------------------------------------------------------
{
  int         i;
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
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf.  This returns to original place the swapping done
// by setPForce.
// ---------------------------------------------------------------------------
{
  int        i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM];
    Uf[i][0] -> gradient (i);
  
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Field::swapData (Us[i][0], Uf[i][0]);
  }
}


static void setUForce (const Domain* D ,
		       AuxField***   Uf)
// ---------------------------------------------------------------------------
// On entry, intermediate velocity storage u^^ is in lowest levels of Uf.
// Multiply by -1.0 / (D_T * KINVIS) to create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  int        i;
  const real alpha = -1.0 / Femlib::value ("D_T*KINVIS");

  for (i = 0; i < DIM; i++) *Uf[i][0] *= alpha;
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
  int                 i, j, found;
  const int           nSys   = D -> Nsys.getSize();
  const int           nZ     = D -> u[0] -> nZ();
  const int           nModes = (nZ + 1) >> 1;
  const int           itLev  = (int) Femlib::value ("ITERATIVE");
  const int           nOrder = (int) Femlib::value ("N_TIME");
  const real          beta   = Femlib::value ("BETA");
  ModeSys**           M      = new ModeSys* [DIM + 1];
  vector<Element*>&   E      = ((Domain*) D) -> Esys;
  const NumberSystem* N;

  // -- Velocity systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real   lambda2 = alpha[0] / Femlib::value ("KINVIS*D_T");

    N = D -> u[0] -> system();
    M[0] = new ModalMatrixSystem (lambda2, beta, nModes, E, N);

    for (i = 1; i < DIM; i++) {
      for (found = 0, j = 0; j < i; j++)
	if (found = D -> u[j] -> system() == D -> u[i] -> system()) break;
      
      if (found)
	M[i] = M[j];
      else {
	N = D -> u[i] -> system();
	M[i] = new ModalMatrixSystem (lambda2, beta, nModes, E, N);
      }
    }

  } else
    for (i = 0; i < DIM; i++) M[i] = 0;

  // -- Pressure system.

  if (itLev < 2) {
    N      = D -> u[DIM] -> system();
    M[DIM] = new ModalMatrixSystem (0.0, beta, nModes, E, N);

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
      const real   lambda2 = alpha[0] / Femlib::value ("D_T*KINVIS");
      
      U -> solve (Force, lambda2);

    } else
      U -> solve (Force, M);

    return;
  }
}
