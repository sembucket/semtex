///////////////////////////////////////////////////////////////////////////////
// transport.C: Unsteady advection--diffusion solver, Cartesian cords.
//
// Here there is only a single velocity component, and we solve
// for the transport of the scalar field 'c'.  The single component
// of fluid velocity, 'w', is in the Fourier/z direction; this is a 
// linear transport problem:
//
//   dc     dc
//   -- + w -- = D div grad c.
//   dt     dz
//
// We use incompressibility of the velocity field and implement 
// the advection terms in a "skew-symmetric" formulation:
//
//   dc         /   dc   dwc \
//   -- = - 0.5 | w -- + --- | + D div grad c.
//   dt         \   dz    dz /
//
// Time split is discussed in: Karniadakis, Israeli & Orszag,
// "High-order splitting methods for the incompressible Navier--Stokes
// equations", JCP 9(2). 1991.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <chroma.h>

typedef ModalMatrixSystem ModeSys;

static void     advection (AuxField*, AuxField*, const real*);
static void     extrap    (AuxField*, AuxField**, AuxField**, const int);
static ModeSys* preSolve  (const Domain*);
static void     Solve     (Field*, AuxField*, ModeSys*, const int, const int);

static real dt;
static int  nZ;


void transport (Domain*     D,
		RunInfo*    I,
		const real* w)
// ---------------------------------------------------------------------------
// On entry, D contains storage for scalar Field 'c', pre-transformed in
// the z direction.  Input w contains the z-component velocity field.
//
// Us is multi-level auxillary Field storage for scalar history and 
// Uf is multi-level auxillary Field storage for advection history terms.
// ---------------------------------------------------------------------------
{
  int       i;
  const int Je     = (int) Femlib::value ("N_TIME");
  const int nStep  = (int) Femlib::value ("N_STEP");
  Field*    scalar = D -> u[0];

  // -- Set file-scope globals;

  dt = Femlib::value ("D_T");
  nZ = Geometry::nZProc();

  // -- Create global matrix systems.

  ModeSys* MMS = preSolve (D);

  // -- Create & initialise multi-level storage for old fields and advection.

  AuxField** Us = new AuxField* [(size_t) Je];
  AuxField** Uf = new AuxField* [(size_t) Je];

  for (i = 0; i < Je; i++) {
    *(Us[i] = new AuxField (D -> Esys, nZ)) = 0.0;
    *(Uf[i] = new AuxField (D -> Esys, nZ)) = 0.0;
  }

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Explicit update for advection and old time levels.

    AuxField::swapData (scalar, Us[0]);
    advection (Us[0],  Uf[0], w);
    extrap    (scalar, Us, Uf, D -> step);
    roll      (Us, Je);
    roll      (Uf, Je);
    AuxField::swapData (scalar, Us[0]);

    // -- Diffusion substep.

    Solve (scalar, Us[0], MMS, D -> step, Je);

    // -- Process results of this step.

    I -> report (Us[0]);

  }
}


static void advection (AuxField*   C,
		       AuxField*   A,
		       const real* w)
// ---------------------------------------------------------------------------
// Compute the Fourier transform of the advection terms
// 
//       ^
//      dc
//    w --  = i beta k w c
//      dz                 k
//         k
//
// and leave in A.  This is a linear operation, no dealiasing needed.
// Old time level of scalar is passed in as C.
// ---------------------------------------------------------------------------
{
  int i;

  for (i = 0; i < nZ; i++)
    A -> setPlane (i, w);

  (A -> times (*A, *C)) . gradient (2);
}


static void extrap (AuxField*  C  ,
		    AuxField** Us ,
		    AuxField** Uf ,
		    const int  step)
// ---------------------------------------------------------------------------
// Compute the RHS explicit approximation of forcing terms.
//
// On entry, the most recent scalar fields are in Us, and the most
// recent advection terms in Uf.  The forcing field is
// computed and left in C's storage areas.
// ---------------------------------------------------------------------------
{
  int          q, Je = (int) Femlib::value ("N_TIME");
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);

  Je = min (step, Je);  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je + 1, 1.0 / dt, alpha(),  1);

  *C = 0.0;

  for (q = 0; q < Je; q++) {
    C -> axpy (alpha[q + 1], *Us[q]);
    C -> axpy (beta [q]    , *Uf[q]);
  }
}


static ModeSys* preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystem for D -> u[0].  If iterative solution
// (ITERATIVE == 1) is selected, return 0.
// ---------------------------------------------------------------------------
{
  const char           name   = D -> u[0] -> name();
  const int            nSys   = D -> Nsys.getSize();
  const int            nmodes = Geometry::nModeProc();
  const int            base   = Geometry::baseMode();
  const int            itLev  = (int) Femlib::value ("ITERATIVE");
  const int            nOrder = (int) Femlib::value ("N_TIME");
  const real           beta   = Femlib::value ("BETA");
  ModeSys*             M;
  vector<Element*>&    E      = ((Domain*) D) -> Esys;
  const NumberSystem** N      = new const NumberSystem* [(size_t) 3];

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real lambda2 = alpha[0] / dt;

    D -> setNumber (name, N);
    M = new ModalMatrixSystem (lambda2, beta, name, base, nmodes, E, N);

  } else M = 0;

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
  if (M == 0 || step < nOrder) {

    const int    Je = min (step, nOrder);
    vector<real> alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
    const real   lambda2 = alpha[0] / dt;

    U -> solve (Force, lambda2);

  } else
    U -> solve (Force, M);

}
