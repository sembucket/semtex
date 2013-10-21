//////////////////////////////////////////////////////////////////////////////
// integrate.C: Unsteady Navier--Stokes solver, using "stiffly-stable"
// time integration.  Geometries may be 2- or 3-dimensional, Cartesian
// or cylindrical.  Fourier expansions are used in the homogeneous
// direction.  This file now provides integrate as a call-back
// routine; after initialisation, integrate may be called repeatedly
// without reinitialising internal storage.
//
// The form of the advection terms is selected by the input scheme.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "newt.h"

typedef ModalMatrixSys Msys;

static int_t  NDIM, NCOM, NORD, C3D;

static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);


void integrate (Domain*        D,
		Advection scheme)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NDIM = Geometry::nDim();	// -- Number of space dimensions.
  NCOM = D -> nField() - 1;	// -- Number of velocity components.
  NORD = Femlib::ivalue ("N_TIME");
  C3D  = Geometry::cylindrical() && NDIM == 3;

  int_t        i, j, k;
  const real_t dt    = Femlib:: value ("D_T");
  const int_t  nStep = Femlib::ivalue ("N_STEP");
  const int_t  nZ    = Geometry::nZProc();
  Msys**       MMS;

  vector<real_t> ff (3);

  if (scheme == linear ) {
    ff[0] = 0.0;
    ff[1] = 0.0;
    ff[2] = 0.0;
  } else {
    ff[0] = Femlib::value ("FFX");
    ff[1] = Femlib::value ("FFY");
    ff[2] = Femlib::value ("FFZ");
  }

  static Msys**      MMSL;
  static Msys**      MMSN;
  static AuxField*** Us;
  static AuxField*** Uf;
  static Field*      Pressure = D -> u[NCOM];

  // -- Initialise static storage.

  if (!MMSN && scheme == nonlinear) MMSN = preSolve (D);
  if (!MMSL && scheme ==    linear) MMSL = preSolve (D);

  if (!Us) {
    
    // -- Create multi-level storage for velocities and forcing.
    
    const int_t ntot  = Geometry::nTotProc();
    real_t*     alloc = new real_t [static_cast<size_t> (2*NCOM*NORD*ntot)];
    Us                = new AuxField** [static_cast<size_t> (2*NORD)];
    Uf                = Us + NORD;

    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [static_cast<size_t> (2*NCOM)];
      Uf[i] = Us[i] + NCOM;
      for (j = 0; j < NCOM; j++) {
	Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
	Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
      }
    }

    // -- Create multi-level storage for pressure BCS.

    PBCmgr::build (Pressure);

    // -- Apply coupling to radial & azimuthal velocity BCs.
    
    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);
  }
    
  // -- Because we may restart from scratch on each call, zero these:

  for (i = 0; i < NORD; i++)
    for (j = 0; j < NCOM; j++) {
      *Us[i][j] = 0.0;
      *Uf[i][j] = 0.0;
    }

  // -- Select modal matrix systems according to advection scheme.

  MMS = (scheme == linear) ? MMSL : MMSN;

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Compute advection terms from previous velocity field.

    scheme (D, Us[0], Uf[0], ff);

    // -- Update high-order pressure BC storage.

    PBCmgr::maintain (D -> step, Pressure, 
		      const_cast<const AuxField**>(Us[0]),
		      const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (D -> step);
    
    // -- Complete unconstrained advective substep.

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp (D, const_cast<const AuxField***>(Us), 
                 const_cast<const AuxField***>(Uf));

    // -- Compute pressure.

    for (i = 0; i < NCOM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NCOM);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NCOM,  Uf[0][0], MMS[NCOM]);

    // -- Correct velocities for pressure gradient. 

    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NCOM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NCOM);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NCOM; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);
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
// This is the only routine that makes explicit use of the multi time
// level structure of Us & Uf.
// ---------------------------------------------------------------------------
{
  int_t             i, q;
  vector<AuxField*> H (NCOM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NCOM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t    Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);
  Blas::scal (Je, Femlib::value ("D_T"), &beta[0],  1);

  for (i = 0; i < NCOM; i++)
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
  int_t        i;
  const real_t dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) (*Uf[i] = *Us[i]) . gradient (i);

  for (i = 1; i < NDIM; i++) *Uf[0] += *Uf[i];

  *Uf[0] /= dt;
}


static void project (const Domain* D ,
		     AuxField**    Us,
		     AuxField**    Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in Us.  Constrain velocity field:
//
//                    u^^ = u^ - D_T * grad P,
//
// then scale by -1.0 / (D_T * KINVIS) to create forcing for viscous step.
//
// u^^ is left in Uf.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");

  for (i = 0; i < NCOM; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i == 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
  }

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[NCOM]) . gradient (i);
    if (Geometry::cylindrical() && i <  2) Us[0] -> mulY();
    Uf[i] -> axpy (beta, *Us[0]);
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
  const int_t             nmodes = Geometry::nModeProc();
  const int_t             base   = Geometry::baseMode();
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const real_t            beta   = Femlib:: value ("BETA");
  const vector<Element*>& E = D -> elmt;
  Msys**                  M = new Msys* [static_cast<size_t>(NCOM + 1)];
  int_t                   i;

  vector<real_t> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  const real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NCOM; i++)
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i],    (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NCOM] = new Msys
      (0.0,     beta, base, nmodes, E, D -> b[NCOM], (itLev<2)?DIRECT:JACPCG);

  return M;
}


static void Solve (Domain*     D,
		   const int_t i,
		   AuxField*   F,
		   Msys*       M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D->u[i], using F as a forcing Field.
// Iterative or direct solver selected on basis of field type, step,
// time order and command-line arguments.
// ---------------------------------------------------------------------------
{
  const int_t step = D -> step;

  if (i < NCOM && step < NORD) { // -- We need a temporary matrix system.
    const int_t    Je      = min (step, NORD);    
    const int_t    base    = Geometry::baseMode();
    const int_t    nmodes  = Geometry::nModeProc();

    vector<real_t> alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real_t   beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}