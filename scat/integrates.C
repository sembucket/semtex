///////////////////////////////////////////////////////////////////////////////
// integrates.C: (unsteady) linear scalar advection-diffusion in a
// prescribed/fixed velocity field. Essentially the same as scat but
// with no pressure field required, and no evolution of velocity. The
// session file needs to have included fields and BCs for velocity but
// not pressure. Obviously there is also no provision for Boussinesq
// buoyancy.
//
// Copyright (C) 2010 <--> $Date$, Hugh Blackburn.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <scat.h>

typedef ModalMatrixSys Msys;
static  int_t          NCOM, NDIM, NORD;
static  bool           C3D;

static void  advect   (Domain*, AuxField*, AuxField*);
static void  extrap   (AuxField*,const int_t,const AuxField**,const AuxField**);
static Msys* preSolve (const Domain*);
static void  Solve    (Domain*, AuxField*, Msys*);


void AdvectDiffuse (Domain*       D,
		    ScatAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w'),
// scalar Field 'c', and pressure 'p'. Integrate advection-diffusion
// for c, keeping velocity fields frozen and ignoring 'p'.
//
// Us is multi-level auxillary Field storage for scalar and 
// Uf is multi-level auxillary Field storage for advection terms.
// ---------------------------------------------------------------------------
{
  NDIM = Geometry::nDim();	// -- Number of space dimensions.
  NCOM = D -> nField() - 2;	// -- Number of velocity components. 
  NORD = Femlib::ivalue ("N_TIME");
  C3D  = Geometry::cylindrical() && NDIM == 3;

  int_t        i;
  const real_t dt      = Femlib:: value ("D_T");
  const int_t  nStep   = Femlib::ivalue ("N_STEP");
  const int_t  TBCS    = Femlib::ivalue ("TBCS");
  const int_t  nZ      = Geometry::nZProc();
  const int_t  ntot    = Geometry::nTotProc();
  real_t*      alloc   = new real_t [static_cast<size_t>(2*NORD*ntot)];
  Field        *scalar = D -> u[NCOM];
  AuxField     **Us,  **Uf;

  // -- Create global matrix systems.

  Msys* MMS = preSolve (D);

  // -- Create & initialize multi-level storage for scalar and advection terms.

  Us = new AuxField* [NORD + NORD]; Uf = Us + NORD;
  for (i = 0; i < NORD; i++) {
    *(Us[i] = new AuxField (alloc + i        * ntot, nZ, D -> elmt)) = 0.0;
    *(Uf[i] = new AuxField (alloc + (NORD+i) * ntot, nZ, D -> elmt)) = 0.0;
  }

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; D -> time += dt; Femlib::value ("t", D -> time);

    // -- Explicit update for advection and old time levels.

    advect (D, Us[0], Uf[0]);

    // -- Construct Helmholz forcing in D -> u[NCOM].

    extrap (scalar, D -> step,
	    const_cast<const AuxField**> (Us),
	    const_cast<const AuxField**> (Uf));

    // -- Rearrange storage.

    rollv (Us, NORD);
    rollv (Uf, NORD);
    AuxField::swapData (scalar, Uf[0]);

    // -- Re-evaluate (time-dependent) BCs?

    if (TBCS == 1)        // -- 2D/mode0 BCs (only).
      ROOTONLY scalar -> evaluateM0Boundaries (D -> step);
    else if (TBCS == 2) { // -- All modes.
      scalar -> evaluateBoundaries (0, false);
      scalar -> bTransform (FORWARD);
    }

    // -- Diffusion substep.

    Solve (D, Uf[0], MMS);

    // -- Process results of this step. Presently this is omitted in
    //    this "frozen velocity" code as there are insufficient
    //    storage locations in Us & Uf.
    //    A -> analyse (Us, Uf);

    // -- For now we just check if a field dump is required.

    D -> dump();
  }
}


static void advect (Domain*   D ,
		    AuxField* Us,
		    AuxField* Uf)
// ---------------------------------------------------------------------------
// Compute advection terms N = u . grad c, built in alternating
// skew-symmetric form, toggling between convective and conservative
// formulations of advection terms.
//
// On odd step numbers, make (convective)           u . grad c
//                                                  ~
// while on even step numbers, make (conservative)  div uc
//                                                      ~
// This seems to be efficient and about as robust as full skew-symmetric.
//
// Scalar field data area of D Us are swapped, then the next stage of
// nonlinear forcing terms N(u) are computed from velocity and scalar
// fields and left in Uf.
//
// If STOKES is defined for compilation, the convective terms in the
// scalar equation are zero, so that the transient diffusion equation
// is solved in place of the convection/diffusion equation
// (i.e. scalar is uncoupled from the velocity field).
//
// NB: no dealiasing for concurrent execution (or if ALIAS is defined).
// ---------------------------------------------------------------------------
{
  int_t i, j;

#if defined(STOKES)

  *Uf = 0.0;
  AuxField::swapData (D -> u[NCOM], Us);

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

  Field*            master = D -> u[NCOM];
  AuxField*         N;

  static vector<real_t>    work ((NCOM + 3) * nTot32);
  static vector<real_t*>   u32   (NCOM + 1);
  static vector<AuxField*> U     (NCOM + 1);
  static real_t*           tmp = NULL; // -- First-time flag.
  static real_t*           n32 = NULL;

  static int               toggle = 1; // -- Alternation flag.

  if (!n32) {
    for (i = 0; i <= NCOM; i++) u32[i] = &work[i * nTot32];
    for (i = 0; i <  NCOM; i++) {
      U[i] = D -> u[i];
      U[i] -> transform32 (INVERSE, u32[i]);
    }
    n32 = &work[(NCOM + 1) * nTot32];
    tmp = &work[(NCOM + 2) * nTot32];
  }

  N = Uf;
  Veclib::zero (3 * nTot32, &work[NCOM * nTot32], 1);

  AuxField::swapData (D -> u[NCOM], Us);
  (U[NCOM] = Us) -> transform32 (INVERSE, u32[NCOM]);

  if (Geometry::cylindrical()) { // -- Cylindrical coordinates.

    if (toggle) {		// -- convective/non-conservative.

      if (nZ > 2) {
	Veclib::copy       (nTot32, u32[NCOM], 1, tmp, 1);
	Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, nPP, tmp, 2);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32, 1, n32, 1);
	
	master -> divY     (nZ32, n32);
      }

      // -- 2D derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, u32[NCOM], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vvtvp      (nTot32, u32[j], 1, tmp, 1, n32, 1, n32, 1);
      }

    } else {			// -- conservative.

      if (nZ > 2) {
	Veclib::vmul       (nTot32, u32[NCOM], 1, u32[2], 1, tmp, 1);
	Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, nPP, tmp, 2);
	Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	Veclib::vadd       (nTot32, tmp, 1, n32, 1, n32, 1);

	master -> divY     (nZ32, n32);
      }

      // -- 2D derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul       (nTot32, u32[j], 1, u32[NCOM], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vadd       (nTot32, tmp, 1, n32, 1, n32, 1);
      }
    }
  
  } else {			// -- Cartesian coordinates.

    if (toggle) {	       // -- Perform n_i += u_j d(u_i) / dx_j.

      for (j = 0; j < NDIM; j++) {
	Veclib::copy (nTot32, u32[NCOM], 1, tmp,  1);
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
	Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32, 1, n32, 1);
      }

    } else {		       // -- Perform n_i += d(u_i u_j) / dx_j.

      for (j = 0; j < NDIM; j++) {
	Veclib::vmul  (nTot32, u32[NCOM], 1, u32[j], 1, tmp,  1);
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
	Veclib::vadd (nTot32, tmp, 1, n32, 1, n32, 1);
      }
    }
  }

  // -- Transform to Fourier space and smooth.
      
  N -> transform32 (FORWARD, n32);
  master -> smooth (N);

  toggle = 1 - toggle;

#endif
}


static void extrap (AuxField*        C   ,
		    const int_t      step,
		    const AuxField** Us  ,
		    const AuxField** Uf  )
// ---------------------------------------------------------------------------
// On entry, the most recent scalar field is in Us[0], and the most
// recent advection term in Uf[0].  The intermediate scalar field
// (forcing for Helmholtz equation) is computed and left in D's
// scalar storage (passed in as C).
// ---------------------------------------------------------------------------
{
  int_t          q;
  const int_t    Je = min (step, NORD);
  const real_t   dt = Femlib::value ("D_T");
  const real_t   DiffusionCoeff = Femlib::value ("KINVIS / PRANDTL");
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);

  Blas::scal (Je + 1, 1.0 / (dt * DiffusionCoeff), &alpha[0], 1);
  Blas::scal (Je,     1.0 / DiffusionCoeff       , &beta [0], 1);

  *C = 0.0;

  for (q = 0; q < Je; q++) {
    C -> axpy (alpha[q + 1], *Us[q]);
    C -> axpy (beta [q]    , *Uf[q]);
  }
}


static Msys* preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystem for scalar c.
// ---------------------------------------------------------------------------
{
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const int_t             nmodes = Geometry::nModeProc();
  const int_t             base   = Geometry::baseMode();
  const real_t            beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  vector<real_t>          alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  real_t                  lambda2 = 
                            alpha[0]/Femlib::value ("D_T * KINVIS / PRANDTL");
  return new Msys
    (lambda2, beta, base, nmodes, E, D -> b[NCOM], (itLev<1) ? DIRECT:JACPCG);
}


static void Solve (Domain*   D,
		   AuxField* F,
		   Msys*     M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D->u[NCOM], using F as a forcing Field.
// Iterative or direct solver selected on basis of step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  const int_t step = D -> step;

  if (step < NORD) { // -- We need a temporary matrix system.
    const int_t Je     = min (step, NORD);    
    const int_t base   = Geometry::baseMode();
    const int_t nmodes = Geometry::nModeProc();

    vector<real_t> alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t   beta    = Femlib::value ("BETA");
    real_t         lambda2 = alpha[0]/Femlib::value ("D_T * KINVIS / PRANDTL");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[NCOM], JACPCG);
    D -> u[NCOM] -> solve (F, tmp);
    delete tmp;

  } else D -> u[NCOM] -> solve (F, M);
}
