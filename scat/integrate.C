///////////////////////////////////////////////////////////////////////////////
// integrate.C: unsteady Navier--Stokes DNS, using "stiffly-stable"
// time integration.
//
// Copyright (C) 2001 Hugh Blackburn.
//
// This version of basic semtex DNS includes transport of a passive
// scalar 'c', carried as the DIMth Field; pressure is
// DIM+1. Cartesian and cylindrical coordinates.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <scat.h>

typedef ModalMatrixSys Msys;
static  integer        NDIM, NORD, CYL, C3D;

static void   advection (Domain*, AuxField**, AuxField**);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const integer, AuxField*, Msys*);


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
  NDIM = Geometry::nDim();
  NORD = static_cast<integer>(Femlib::value ("N_TIME"));
  CYL  = Geometry::system() == Geometry::Cylindrical;
  C3D  = CYL && NDIM == 3;

  integer       i, j, k;
  const real    dt     =                      Femlib::value ("D_T");
  const integer nStep  = static_cast<integer>(Femlib::value ("N_STEP"));
  const integer nZ     = Geometry::nZProc();
  const integer ntot   = Geometry::nTotProc();
  real*         alloc  = new real [static_cast<size_t>(2*(NDIM+1)*NORD*ntot)];

  // -- Create global matrix systems.

  Msys** MMS = preSolve (D);

  // -- Create & initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [static_cast<size_t>(2 * NORD)];
  AuxField*** Uf = Us + NORD;

  for (k = 0, i = 0; i < NORD; i++) {
    Us[i] = new AuxField* [static_cast<size_t>(2*(NDIM+1))];
    Uf[i] = Us[i] + NDIM + 1;
    for (j = 0; j <= NDIM; j++) {
      *(Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
      *(Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
    }
  }

  // -- Set up multi-level storage for pressure BCS.

  Field* Pressure = D -> u[NDIM + 1];
  PBCmgr::build (Pressure);

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    advection (D, Us[0], Uf[0]);
    waveProp  (D, const_cast<const AuxField***>(Us),
	          const_cast<const AuxField***>(Uf));

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
		      const_cast<const AuxField**>(Us[0]),
		      const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i <= NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NDIM+1);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NDIM+1,  Uf[0][0], MMS[NDIM+1]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i <= NDIM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NDIM+1);

#if defined (UNSTEADY)		// -- Unsteady velocity BCs.
    for (i = 0; i < NDIM; i++) {
      D -> u[i] -> evaluateBoundaries (D -> step);
      D -> u[i] -> bTransform (FORWARD);
    }
#endif

    // -- Viscous and thermal diffusion substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i <= NDIM; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0]);
  }
}


static void advection (Domain*    D ,
		       AuxField** Us,
		       AuxField** Uf)
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
  integer i, j;

#if defined(STOKES)

  for (i = 0; i <= NDIM; i++) {
    *Uf[i] = 0.0;
    AuxField::swapData (D -> u[i], Us[i]);
  }

#else

  const integer     nZ     = Geometry::nZ();
  const integer     nZP    = Geometry::nZProc();
  const integer     nP     = Geometry::planeSize();
  const integer     nPP    = Geometry::nBlock();
  const integer     nPR    = Geometry::nProc();
  const integer     nTot   = Geometry::nTotProc();
#if defined (ALIAS)
  const integer     nZ32   = Geometry::nZProc();
#else
  const integer     nZ32   = Geometry::nZ32();
#endif 
  const integer     nTot32 = nZ32 * nP;
  vector<real>      work ((2 * (NDIM+1) + 1) * nTot32);
  vector<real*>     u32 (NDIM+1);
  vector<real*>     n32 (NDIM+1);
  vector<AuxField*> U   (NDIM+1);
  vector<AuxField*> N   (NDIM+1);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * (NDIM+1) * nTot32;

  Veclib::zero ((2*(NDIM+1)+1) * nTot32, work(), 1); // -- A catch-all cleanup.

  for (i = 0; i <= NDIM; i++) {
    u32[i] = work() +  i             * nTot32;
    n32[i] = work() + (i + NDIM + 1) * nTot32;

    AuxField::swapData (D -> u[i], Us[i]);

    U[i] = Us[i];
    U[i] -> transform32 (INVERSE, u32[i]);
    N[i] = Uf[i];
  }

  if (CYL) {			// -- Cylindrical coordinates.

    for (i = 0; i <= NDIM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);

      if (NDIM > 2) {
	if (i == 1)
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);
	if (i == 3)
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

      master -> divR (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vvtvp      (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul       (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);
      *N[i] *= -0.5;
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i <= NDIM; i++) {
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
      *N[i] *= -0.5;
    }
  }

#endif
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
  vector<AuxField*> H (NDIM + 1);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i <= NDIM; i++) {
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

  for (i = 0; i <= NDIM; i++)
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
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) (*Uf[i] = *Us[i]) . gradient (i);

  if (C3D) Uf[2] -> divR();

  for (i = 1; i < NDIM; i++) *Uf[0] += *Uf[i];

  if (CYL) *Uf[0] += (*Uf[1] = *Us[1]) . divR();

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
  int        i;
  const real dt = Femlib::value ("D_T");
  real       alpha;

  alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  for (i = 0; i < NDIM; i++) {

    (*Uf[i] = *D -> u[NDIM+1]) . gradient (i);

    if (C3D && i == 2) Uf[2] -> divR();

    Us[i] -> axpy (-dt, *Uf[i]);
    Field::swapData (Us[i], Uf[i]);

    *Uf[i] *= alpha;
  }

  alpha = -1.0 / Femlib::value ("D_T * KINVIS / PRANDTL");
  AuxField::swapData (Us[NDIM], Uf[NDIM]);
  *Uf[NDIM] *= alpha;
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
  const integer itLev = static_cast<integer>(Femlib::value ("ITERATIVE"));
  const integer           nmodes = Geometry::nModeProc();
  const integer           base   = Geometry::baseMode();
  const real              beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [static_cast<size_t>(NDIM + 2)];
  vector<real>            alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, alpha());
  integer                 i;
  real                    lambda2;

  // -- Velocity systems.

  lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
  for (i = 0; i < NDIM; i++)
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i],     (itLev<1)?DIRECT:JACPCG);

  // -- Scalar system.

  lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");
  M[NDIM] = new Msys
    (lambda2, beta, base, nmodes, E, D -> b[NDIM],    (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NDIM+1] = new Msys
    (0.0,     beta, base, nmodes, E, D -> b[NDIM+1],  (itLev<2)?DIRECT:JACPCG);

  return M;
}


static void Solve (Domain*       D,
		   const integer i,
		   AuxField*     F,
		   Msys*         M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D->u[i], using F as a forcing Field.
// Iterative or direct solver selected on basis of field type, step,
// time order and command-line arguments.
// ---------------------------------------------------------------------------
{
  const integer step = D -> step;

  if (i <= NDIM && step < NORD) { // -- We need a temporary matrix system.
    const integer Je      = min (step, NORD);    
    const integer base    = Geometry::baseMode();
    const integer nmodes  = Geometry::nModeProc();

    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
    const real    beta    = Femlib::value ("BETA");
    real          lambda2;

    if (i < NDIM) lambda2 = alpha[0] / Femlib::value("D_T * KINVIS");
    else          lambda2 = alpha[0] / Femlib::value("D_T * KINVIS / PRANDTL");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
