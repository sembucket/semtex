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

static void   advection (Domain*, AuxField**, AuxField**, vector<real_t>&);
static void   buoyancy  (Domain*, AuxField**, AuxField**, vector<real_t>&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);


void NavierStokes (Domain*       D,
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
  PBCmgr::build (Pressure);

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Create spatially-constant forcing terms.

  vector<real_t> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");
  
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

    // -- Unconstrained forcing substep.

    advection (D, Us[0], Uf[0], ff);
    buoyancy  (D, Us[0], Uf[0], g);

    // -- Update high-order pressure BC storage.

    PBCmgr::maintain (D -> step, Pressure, 
		      const_cast<const AuxField**>(Us[0]),
		      const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (D -> step);

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

    // -- Re-evaluate (time-dependent) BCs?

    if (TBCS == 1)
      // -- 2D/mode0 base BCs (only).
      for (i = 0; i <= NCOM; i++)
	ROOTONLY D -> u[i] -> evaluateM0Boundaries (D -> step);
    else if (TBCS == 2) {
      // -- All modes.
      for (i = 0; i <= NCOM; i++) {
	D -> u[i] -> evaluateBoundaries (0, false);
	D -> u[i] -> bTransform (FORWARD);
      }
      if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);
    }

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


static void advection (Domain*         D ,
		       AuxField**      Us,
		       AuxField**      Uf,
		       vector<real_t>& ff)
// ---------------------------------------------------------------------------
// Compute advection terms in Navier--Stokes equations: N(u).
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
// If STOKES is defined for compilation, the nonlinear terms are set
// to zero.  This means that the Navier--Stokes equations become the
// Stokes equations.  Also, the convective terms in the scalar
// equation are zero, so that the transient diffusion equation is
// solved in place of the convection/diffusion equation (i.e. scalar
// is uncoupled from the velocity field).
//
// The temperature transport term u . grad c is also built this way here.
//
// NB: no dealiasing for concurrent execution (or if ALIAS is defined).
// ---------------------------------------------------------------------------
{
  int_t i, j;

#if defined(STOKES)

  for (i = 0; i <= NCOM; i++) {
    *Uf[i] = 0.0;
    AuxField::swapData (D -> u[i], Us[i]);
  }

#else

  const int_t       nZ     = Geometry::nZ();
  const int_t       nZP    = Geometry::nZProc();
  const int_t       nP     = Geometry::planeSize();
  const int_t       nPP    = Geometry::nBlock();
  const int_t       nPR    = Geometry::nProc();
  const int_t       nTot   = Geometry::nTotProc();
#if defined (ALIAS)
  const int_t       nZ32   = Geometry::nZProc();
#else
  const int_t       nZ32   = Geometry::nZ32();
#endif 
  const int_t       nTot32 = nZ32 * nP;
  vector<real_t>    work ((2 * (NCOM+1) + 1) * nTot32);
  vector<real_t*>   u32 (NCOM+1);
  vector<real_t*>   n32 (NCOM+1);
  vector<AuxField*> U   (NCOM+1);
  vector<AuxField*> N   (NCOM+1);
  Field*            master = D -> u[0];
  real_t*           tmp    = &work[0] + 2 * (NCOM+1) * nTot32;

  Veclib::zero ((2*(NCOM+1)+1) * nTot32, &work[0], 1); // -- Catch-all cleanup.

  for (i = 0; i <= NCOM; i++) {
    u32[i] = &work[0] +  i             * nTot32;
    n32[i] = &work[0] + (i + NCOM + 1) * nTot32;

    AuxField::swapData (D -> u[i], Us[i]);

    U[i] = Us[i];
    U[i] -> transform32 (INVERSE, u32[i]);
    N[i] = Uf[i];
  }

  if (Geometry::cylindrical()) {		// -- Cylindrical coordinates.

    for (i = 0; i <= NCOM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);

      if (NCOM > 2) {
	if (i == 1)		// -- radial compt.
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)		// -- azimuthal compt.
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);
	if (i == 3)		// -- scalar.
	  Veclib::vmul    (nTot32,       u32[3], 1, u32[1], 1,      n32[3], 1);

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);
	
	  Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
	}
      }

      if (i >= 2) master -> divY (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);

	if (i < 2) master -> mulY (nZ32, tmp);

	Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul       (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);

	if (i < 2) master -> mulY (nZ32, tmp);

	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) {
	Veclib::fill (nP, -2.0*ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	N[i] -> addToPlane (0, tmp);
      }

      *N[i] *= -0.5;		// -- Skew-symmetric NL.
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i <= NCOM; i++) {
      for (j = 0; j < NDIM; j++) {
      
	// -- Perform n_i += u_j d(u_i) / dx_j.

	Veclib::copy (nTot32, u32[i], 1, tmp,  1);
	if (j == 2) {
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

	// -- Perform n_i += d(u_i u_j) / dx_j.

	Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
	if (j == 2) {
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
	
      }

      // -- Transform to Fourier space, smooth, add forcing.
      
      N[i]   -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);

      ROOTONLY if (fabs (ff[i]) > EPSDP) N[i] -> addToPlane (0, -2.0*ff[i]);

      *N[i] *= -0.5;
    }
  }

#endif
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
// Cylindrical coordinates: for simplicity, we only currently allow
// axial buoyancy. To allow off-axis gravity would require the
// buoyancy to be done pseudospectrally - not so hard.  Note that the
// axial and radial terms need to be multiplied by y here to match
// what we do with the formation of the nonlinear terms.
// ---------------------------------------------------------------------------
{
  int_t     i;
  AuxField* T    = Us [NCOM];
  AuxField* work = D -> u[NCOM];

  if (Geometry.cylindrical()) {
    if (fabs (g[0]) > EPSDP) {
      *work      = *T;
      *work     -=  Femlib::value ("T_REF" );
      *work     *=  Femlib::value ("BETA_T");
      *work     *=  g[0];
      *work     -> mulY();
      *Uf[0] -= *work;
    }
  } else
    for (i = 0; i < NCOM; i++)
      if (fabs (g[i]) > EPSDP) {
	*work      = *T;
	*work     -=  Femlib::value ("T_REF" );
	*work     *=  Femlib::value ("BETA_T");
	*work     *=  g[i];
	*Uf[i] -= *work;
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
