///////////////////////////////////////////////////////////////////////////////
// integrateNS.C: Unsteady Navier--Stokes solver, using "stiffly-stable"
// time integration.  Geometries may be 2- or 3-dimensional, Cartesian
// or cylindrical.  Fourier expansions are used in the homogeneous
// direction.  This file now provides integrateNS as a call-back
// routine; after initialisation, integrate may be called repeatedly
// without reinitialising internal storage.
//
// Copyright (C) 1994,2003 Hugh Blackburn
//
//
// For cylindrical coordinates:
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <dns.h>

// -- Dealiased product routine for testing.

void daprod (const int, const real*, const real*, real*);

typedef ModalMatrixSys Msys;

static integer NDIM, NCOM, NORD, CYL, C3D;

static void   nonLinear (Domain*, AuxField**, AuxField**, vector<real>&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const integer, AuxField*, Msys*);


void integrateNS (Domain*      D,
		  DNSAnalyser* A)
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
  NORD = (integer) Femlib::value ("N_TIME");
  CYL  = Geometry::system() == Geometry::Cylindrical;
  C3D  = CYL && NDIM == 3;

  integer       i, j, k;
  const real    dt     =           Femlib::value ("D_T");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();

  static Msys**      MMS;
  static AuxField*** Us;
  static AuxField*** Uf;
  static Field*      Pressure = D -> u[NCOM];

  if (!MMS) {			// -- Initialise static storage.

    // -- Create global matrix systems.

    MMS = preSolve (D);

    // -- Create, multi-level storage for velocities and forcing.
    
    const integer ntot  = Geometry::nTotProc();
    real*         alloc = new real [(size_t) 2 * NCOM * NORD * ntot];
    Us                  = new AuxField** [(size_t) 2 * NORD];
    Uf                  = Us + NORD;

    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [(size_t) 2 * NCOM];
      Uf[i] = Us[i] + NCOM;
      for (j = 0; j < NCOM; j++) {
	*(Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
	*(Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
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

  // -- Create spatially-constant forcing terms.

  vector<real> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us[0], Uf[0], ff);
    waveProp  (D, (const AuxField***)Us, (const AuxField***)Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure, 
		      (const AuxField**)Us[0],
		      (const AuxField**)Uf[0]);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < NCOM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NCOM);
    setPForce ((const AuxField**)Us[0], Uf[0]);
    Solve     (D, NCOM,  Uf[0][0], MMS[NCOM]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NCOM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NCOM);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NCOM; i++) {
#if defined (TBCS) // -- Re-evaluate the (time-varying) 2D base BCs (only).
      ROOTONLY D -> u[i] -> evaluateM0Boundaries (D -> step);
#endif
      Solve (D, i, Uf[0][i], MMS[i]);
    }
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0], Uf[0]);
  }

}


static void nonLinear (Domain*       D ,
		       AuxField**    Us,
		       AuxField**    Uf,
		       vector<real>& ff)
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u) + ff.
//
// Here N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS and ff is a vector of body force per unit mass.
//
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 ~ ~
//           N  = -0.5 ( u . grad u + div uu )
//           ~           ~        ~       ~~
//
// i.e., in Cartesian component form
//
//           N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//            i           j    i      j      i j      j
//
// in cylindrical coordinates
//
//           Nx = -0.5 {ud(u)/dx + vd(u)/dy +  d(uu)/dx + d(vu)/dy +
//                 1/y [wd(u)/dz + d(uw)/dz + vu      ]}
//           Ny = -0.5 {ud(v)/dx + vd(v)/dy +  d(uv)/dx + d(vv)/dy +
//                 1/y [wd(v)/dz + d(vw)/dz + vv - 2ww]}
//           Nz = -0.5 {ud(w)/dx + vd(w)/dy +  d(uw)/dx + d(vw)/dy +
//                 1/y [wd(w)/dz + d(ww)/dz + 3wv     ]}
//
// If STOKES is defined for compilation, the nonlinear terms are set to zero.
//
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// NB: no dealiasing for concurrent execution (or if ALIAS is defined).
// ---------------------------------------------------------------------------
{
  integer    i, j;
  const real EPS = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;

#if defined(STOKES)

  for (i = 0; i < NCOM; i++) {
    *Uf[i] = 0.0;
    AuxField::swapData (D -> u[i], Us[i]);
    ROOTONLY if (fabs (ff[i]) > EPS) Uf[i] -> addToPlane (0, -ff[i]);
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
  vector<real>      work ((3 * NCOM + 1) * nTot32);
  vector<real*>     u32 (NCOM);
  vector<real*>     n32 (NCOM);
  vector<AuxField*> U   (NCOM);
  vector<AuxField*> N   (NCOM);
  Field*            master = D -> u[0];
  real*             tmp    = &work[0] + 2 * NCOM * nTot32;
  real*             tmp1   = &work[0] + 3 * NCOM * nTot32;

  Veclib::zero ((2 * NCOM + 1) * nTot32, &work[0], 1); // -- Catch-all cleanup.

  for (i = 0; i < NCOM; i++) {
    u32[i] = &work[0] +  i         * nTot32;
    n32[i] = &work[0] + (i + NCOM) * nTot32;

    AuxField::swapData (D -> u[i], Us[i]);

    U[i] = Us[i];
    U[i] -> transform32 (INVERSE, u32[i]);
    N[i] = Uf[i];
  }

  if (CYL) {			// -- Cylindrical coordinates.

    for (i = 0; i < NCOM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

#ifndef NCONSERVE 
      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);
#endif

      if (NCOM == 3) {

#ifndef NCONSERVE
	if (i == 1)
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);
#else
	if (i == 1)
	  Veclib::svvttvp (nTot32, -1.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  1.0, u32[2], 1, u32[1], 1,      n32[2], 1);
#endif

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);

#ifndef NCONSERVE	
	  Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	  Femlib::exchange   (tmp, nZ32,        nP, FORWARD);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (tmp, nZ32,        nP, INVERSE);
	  Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
#endif
	}
      }

      master -> divY (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vvtvp      (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

#ifndef NCONSERVE
      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul       (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }
#endif

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (FORWARD, n32[i]);
      master -> smooth (N[i]);
#ifndef NCONSERVE
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
#else
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -ff[i]);
      *N[i] *= -1.0;
#endif
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NCOM; i++) {
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
#if 0
	Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);
#else
	daprod (nZ32, u32[j], tmp, tmp1);
	Veclib::vadd (nTot32, tmp1, 1, n32[i], 1, n32[i], 1);
#endif

#ifndef NCONSERVE
	// -- Perform n_i += d(u_i u_j) / dx_j.
#if 0
	Veclib::vmul (nTot32, u32[i], 1, u32[j], 1, tmp, 1);
#else
	daprod (nZ32, u32[i], u32[j], tmp);
#endif
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
#endif	
      }

      // -- Transform to Fourier space, smooth, add forcing.
      
      N[i]   -> transform32 (FORWARD, n32[i]);
//      master -> smooth (N[i]);
#if 0 				// -- do projection-stabilization.
      N[i] -> projStab (0.05, *master); 
#endif

#ifndef NCONSERVE
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
#else
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -ff[i]);
      *N[i] *= -1.0;
#endif
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
// This is the only routine that makes explicit use of the multi time
// level structure of Us & Uf.
// ---------------------------------------------------------------------------
{
  integer           i, q;
  vector<AuxField*> H (NCOM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NCOM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const integer Je = min (D -> step, NORD);
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
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
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) (*Uf[i] = *Us[i]) . gradient (i);

  if (C3D) Uf[2] -> divY();

  for (i = 1; i < NDIM; i++) *Uf[0] += *Uf[i];

  if (CYL) *Uf[0] += (*Uf[1] = *Us[1]) . divY();

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
  integer    i;
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real beta  =  1.0 / Femlib::value ("KINVIS");

  for (i = 0; i < NCOM; i++) {
    Field::swapData (Us[i], Uf[i]);
    *Uf[i] *= alpha;
  }

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[NCOM]) . gradient (i);
    if (C3D && i == 2) Us[0] -> divY();
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
  const integer           nmodes = Geometry::nModeProc();
  const integer           base   = Geometry::baseMode();
  const integer           itLev  = (integer) Femlib::value ("ITERATIVE");
  const real              beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [(size_t) (NCOM + 1)];
  integer                 i;
  vector<real>            alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  const real              lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NCOM; i++)
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i],    (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NCOM] = new Msys
      (0.0,     beta, base, nmodes, E, D -> b[NCOM], (itLev<2)?DIRECT:JACPCG);

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

  if (i < NCOM && step < NORD) { // -- We need a temporary matrix system.
    const integer Je      = min (step, NORD);    
    const integer base    = Geometry::baseMode();
    const integer nmodes  = Geometry::nModeProc();
    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real    beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
