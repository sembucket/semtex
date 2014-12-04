///////////////////////////////////////////////////////////////////////////////
// integrate.C: unsteady Navier--Stokes DNS, using "stiffly-stable"
// time integration.
//
// Copyright (C) 2001 <--> $Date$, Hugh Blackburn.
//
// This version of basic semtex DNS includes transport of a passive
// scalar 'c', carried as the NCOMth Field; pressure is
// NCOM+1. Cartesian and cylindrical coordinates.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <scat.h>

typedef ModalMatrixSys Msys;
static  int_t          NCOM, NDIM, NORD;
static  bool           C3D;

void   nonlinear (Domain*, BCmgr*, AuxField**, AuxField**, vector<real_t>&);

static void   buoyancy  (Domain*, AuxField**, AuxField**, vector<real_t>&);
static void   tempGrad  (AuxField**, AuxField**);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);


void NavierStokes (Domain*       D,
		   BCmgr*        B,
		   ScatAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// temperature Field 'c', followed by constraint/pressure field 'p'.
//
// Us is multi-level auxillary Field storage for velocities+temperature and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NDIM = Geometry::nDim();	// -- Number of space dimensions.
  NCOM = D -> nField() - 2;	// -- Number of velocity components.
  NORD = Femlib::ivalue ("N_TIME");
  C3D  = Geometry::cylindrical() && NDIM == 3;

  int_t        i, j, k;
  const real_t dt     = Femlib:: value ("D_T");
  const int_t  nStep  = Femlib::ivalue ("N_STEP");
  const int_t  TBCS   = Femlib::ivalue ("TBCS");
  const int_t  nZ     = Geometry::nZProc();
  const int_t  ntot   = Geometry::nTotProc();
  real_t*      alloc  = new real_t [static_cast<size_t>(2*(NCOM+1)*NORD*ntot)];

  // -- Create global matrix systems.

  Msys** MMS = preSolve (D);

  // -- Create & initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [static_cast<size_t>(2 * NORD)];
  AuxField*** Uf = Us + NORD;

  for (k = 0, i = 0; i < NORD; i++) {
    Us[i] = new AuxField* [static_cast<size_t>(2*(NCOM+1))];
    Uf[i] = Us[i] + NCOM + 1;
    for (j = 0; j <= NCOM; j++) {
      *(Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
      *(Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
    }
  }

  // -- Set up multi-level storage for pressure BCS.

  Field* Pressure = D -> u[NCOM + 1];
  B -> buildComputedBCs (Pressure);

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Create spatially-constant forcing terms.

  vector<real_t> ff (4);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = (NCOM == 3) ? Femlib::value ("FFZ") : Femlib::value ("FFC");
  ff[3] = (NCOM == 3) ? Femlib::value ("FFC") : 0.0;
  
  // -- Set up gravity vector.  Note directional components g_1, g_2, g_3.
  
  vector<real_t> g(3);
  {
    const real_t norm = Femlib::value ("sqrt (g_1*g_1 + g_2*g_2 + g_3*g_3)");

    g[0] = Femlib::value ("g_1 * GRAVITY") / norm;
    g[1] = Femlib::value ("g_2 * GRAVITY") / norm;
    g[2] = Femlib::value ("g_3 * GRAVITY") / norm;
  }

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substeps.

    nonlinear (D, B, Us[0], Uf[0], ff);

    buoyancy  (D, Us[0], Uf[0], g);

    ROOTONLY tempGrad (Us[0], Uf[0]);

    // -- Update high-order pressure BC storage.

    B -> maintainFourier (D -> step, Pressure, 
			  const_cast<const AuxField**>(Us[0]),
			  const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (Pressure, D -> step);

    // -- Complete unconstrained advective substep and compute pressure.

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp  (D, const_cast<const AuxField***>(Us),
	          const_cast<const AuxField***>(Uf));
    for (i = 0; i <= NCOM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NCOM+1);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NCOM+1,  Uf[0][0], MMS[NCOM+1]);

    // -- Correct velocities for pressure. 

    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity+scalar storage.

    for (i = 0; i <= NCOM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NCOM+1);

    // -- Re-evaluate (possibly time-dependent) BCs.

    for (i = 0; i < NCOM; i++)  {
      D -> u[i] -> evaluateBoundaries (0, D -> step, false);
      D -> u[i] -> bTransform         (FORWARD);
      D -> u[i] -> evaluateBoundaries (Pressure, D -> step, true);
    }
    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

    // -- Viscous and thermal diffusion substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i <= NCOM; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0], Uf[0]);
  }
}


static void buoyancy (Domain*         D ,
		      AuxField**      Us,
		      AuxField**      Uf,
		      vector<real_t>& g )
// ---------------------------------------------------------------------------
// The Boussinesq buoyancy term in the momentum equation is
//
//                      - BETA_T * (T - T_REF) g.
//                                             ~
// We add an explicit estimate of this to the nonlinear terms in Uf.
// The first level of Us has the last values of the data Fields, D is free.
//
// Cylindrical coordinates: for simplicity, we only presently allow
// axial buoyancy. To allow off-axis gravity would require the
// buoyancy to be done pseudospectrally - not so hard.  Note that the
// axial (and radial, if we computed them) terms need to be multiplied
// by y here to match what we do with the formation of the nonlinear
// terms.
// ---------------------------------------------------------------------------
{
  int_t     i;
  AuxField* T    = Us [NCOM];
  AuxField* work = D -> u[NCOM];

  if (Geometry::cylindrical()) {
    if (fabs (g[0]) > EPSDP) {
      *work      = *T;
      ROOTONLY work -> addToPlane (0, -Femlib::value ("T_REF"));
      *work     *=  Femlib::value ("BETA_T");
      *work     *=  g[0];
       work     -> mulY();
      *Uf[0] -= *work;
    }
  } else
    for (i = 0; i < NCOM; i++)
      if (fabs (g[i]) > EPSDP) {
	*work      = *T;
	ROOTONLY work -> addToPlane (0, -Femlib::value ("T_REF"));
	*work     *=  Femlib::value ("BETA_T");
	*work     *=  g[i];
	*Uf[i] -= *work;
      }
}


static void tempGrad (AuxField** Us,
		      AuxField** Uf)
// ---------------------------------------------------------------------------
// If there is heat input we may wish to take off a linear correction
// for bulk temperature gradient DTBDX, only in defined for the first
// coordinate direction. This procedure only drives mode 0, we work in
// Fourier. 
// ---------------------------------------------------------------------------
{
  const real_t dtbdx = Femlib::value ("DTBDX");

  if (fabs (dtbdx) < EPSDP) return;

  const int_t     nP = Geometry::planeSize();
  AuxField*       U  = Us[0];
  AuxField*       N  = Uf[NCOM]; // -- Advection term for scalar.
  vector<real_t>  work(nP);
  real_t*         u = &work[0];

  U -> getPlane  (0, u);
  Blas::scal     (nP, -dtbdx, u, 1);

  N -> addToPlane (0, u);
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
  int_t             i, q;
  vector<AuxField*> H (NCOM + 1);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i <= NCOM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t    Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);
  Blas::scal (Je, Femlib::value ("D_T"), &beta[0],  1);

  for (i = 0; i <= NCOM; i++)
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
// A similar operation occurs for the intermediate scalar Field,
// except that the coefficient of diffusivity is used in place of
// kinematic viscosity.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");
  const real_t Pr    =        Femlib::value ("PRANDTL");

  for (i = 0; i <= NCOM; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i >= 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
  }
  *Uf[NCOM] *= Pr;

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[NCOM+1]) . gradient (i);
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
// ITERATIVE == 1 selects iterative solvers for velocities and temperature,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const int_t             nmodes = Geometry::nModeProc();
  const int_t             base   = Geometry::baseMode();
  const real_t            beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [static_cast<size_t>(NCOM + 2)];
  vector<real_t>          alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  int_t                   i;
  real_t                  lambda2;

  // -- Velocity systems.

  lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
  for (i = 0; i < NCOM; i++)
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i],     (itLev<1)?DIRECT:JACPCG);

  // -- Scalar system.

  lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");
  M[NCOM] = new Msys
    (lambda2, beta, base, nmodes, E, D -> b[NCOM],    (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NCOM+1] = new Msys
    (0.0,     beta, base, nmodes, E, D -> b[NCOM+1],  (itLev<2)?DIRECT:JACPCG);

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

  if (i <= NCOM && step < NORD) { // -- We need a temporary matrix system.
    const int_t Je      = min (step, NORD);    
    const int_t base    = Geometry::baseMode();
    const int_t nmodes  = Geometry::nModeProc();

    vector<real_t>  alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t    beta    = Femlib::value ("BETA");
    real_t          lambda2;

    if (i < NCOM) lambda2 = alpha[0] / Femlib::value("D_T * KINVIS");
    else          lambda2 = alpha[0] / Femlib::value("D_T * KINVIS / PRANDTL");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
