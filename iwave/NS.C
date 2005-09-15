///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
// Geometries may be 2- or 3-dimensional, Cartesian or cylindrical.
// Fourier expansions are used in the homogeneous direction.
//
// References:
// 1.  Karniadakis, Israeli & Orszag 1991.  "High-order splitting methods
//     for the incompressible Navier--Stokes equations", JCP 9(2).
// 2.  Tomboulides, Orszag & Karniadakis 1993.  "Direct and
//     large-eddy simulation of axisymmetric wakes",  AIAA-93-0546.
//
// For cylindrical coordinates:
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <dns.h>

typedef ModalMatrixSys Msys;

static int_t NDIM, NORD, CYL, C3D;

static void   nonLinear (Domain*, AuxField**, AuxField**);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);


void NavierStokes (Domain*      D,
		   DNSAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NDIM = Geometry::nDim();
  NORD = Femlib::ivalue ("N_TIME");
  CYL  = Geometry::system() == Geometry::Cylindrical;
  C3D  = CYL && NDIM == 3;

  int_t        i, j, k;
  const real_t dt     = Femlib::value  ("D_T");
  const int_t  nStep  = Femlib::ivalue ("N_STEP");
  const int_t  nZ     = Geometry::nZProc();
  const int_t  nP     = Geometry::planeSize();
  const int_t  ntot   = Geometry::nTotProc();
  real_t*      alloc  = new real_t [(size_t) 2 * NDIM * NORD * ntot];

  if (!C3D || nZ < 4)
    message ("iwave", "must be Cylindrical, 3D, nZ/Proc > 6", ERROR);

  // -- Create global matrix systems.

  Msys** MMS = preSolve (D);

  // -- Create, initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [(size_t) 2 * NORD];
  AuxField*** Uf = Us + NORD;

  for (k = 0, i = 0; i < NORD; i++) {
    Us[i] = new AuxField* [(size_t) 2 * NDIM];
    Uf[i] = Us[i] + NDIM;
    for (j = 0; j < NDIM; j++) {
      *(Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
      *(Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
    }
  }

  // -- Create multi-level storage for pressure BCS.

  Field* Pressure = D -> u[NDIM];
  PBCmgr::build (Pressure);

  // -- Apply coupling to radial & azimuthal velocity BCs.

  Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us[0], Uf[0]);
    waveProp  (D, (const AuxField***)Us, (const AuxField***)Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure, 
		      (const AuxField**)Us[0],
		      (const AuxField**)Uf[0]);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NDIM);
    setPForce ((const AuxField**)Us[0], Uf[0]);
    Solve     (D, NDIM, Uf[0][0], MMS[NDIM]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NDIM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NDIM);

    // -- Viscous correction substep.

    AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
    AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    for (i = 0; i < NDIM; i++) Solve (D, i, Uf[0][i], MMS[i]);
    AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0], Uf[0]);
  }
}


static void nonLinear (Domain*    D ,
		       AuxField** Us,
		       AuxField** Uf)
// ---------------------------------------------------------------------------
// Compute nonlinear (forcing) terms in Navier--Stokes equations: N(u).
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
// Forcing for inertia-wave problem: see below.
//
// NB: no dealiasing for concurrent execution.
// ---------------------------------------------------------------------------
{
  int_t i, j;

#if defined(STOKES)

  for (i = 0; i < NDIM; i++) {
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
  const int_t       nZ32   = Geometry::nZ32();
  const int_t       nTot32 = nZ32 * nP;
  vector<real_t>    work ((2 * NDIM + 1) * nTot32);
  vector<real_t*>   u32 (NDIM);
  vector<real_t*>   n32 (NDIM);
  vector<AuxField*> U   (NDIM);
  vector<AuxField*> N   (NDIM);
  Field*            master = D -> u[0];
  real_t*           tmp    = &work[0] + 2 * NDIM * nTot32;

  Veclib::zero ((2*NDIM + 1) * nTot32, &work[0], 1); // -- A catch-all cleanup.

  for (i = 0; i < NDIM; i++) {
    u32[i] = &work[0] +  i         * nTot32;
    n32[i] = &work[0] + (i + NDIM) * nTot32;

    AuxField::swapData (D -> u[i], Us[i]);

    U[i] = Us[i];
    U[i] -> transform32 (INVERSE, u32[i]);
    N[i] = Uf[i];
  }

  for (i = 0; i < NDIM; i++) {

    // -- Terms involving azimuthal derivatives and frame components.

    if (i == 0)
      Veclib::vmul    (nTot32,       u32[0],1, u32[1],1,            n32[0], 1);
    if (i == 1) {
      Veclib::vmul    (nTot32,       u32[1],1, u32[1],1,            n32[1], 1);
      Veclib::svvttvp (nTot32, -2.0, u32[2],1, u32[2],1, n32[1], 1, n32[1], 1);
    }
    if (i == 2)
      Veclib::svvtt   (nTot32,  3.0, u32[2],1, u32[1],1,            n32[2], 1);

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

    master -> divY (nZ32, n32[i]);
    
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

    // -- Average NL terms.

    Blas::scal (nTot32, -0.5, n32[i], 1);
  }

#endif      
  
  // -- Add forcing for inertia-wave problem (transient version):
  //    (NB: T here is actually OMEGA_1*t.)
  //
  // Fx: +[2*w2/w1*sin(THETA)+THETA_DDOT] *                y*sin(T-z)
  //     +[w2*THETA_DOT/(w1*w1)*cos(THETA)-THETA_DOT/w1] * y*cos(T-z)
  //     +2*{w2/w1*sin(THETA)*cos(T)-THETA_DOT/w1*sin(T)}*{v*cos(z)-w*sin(z)}
  //     +2*{w2/w1*sin(THETA)*sin(T)-THETA_DOT/w1*cos(T)}*{v*sin(z)+w*cos(z)}
  //
  // Fy: +2*[w2/w1*cos(THETA)+1]*w
  //     -2*[w2/w1*sin(THETA)*cos(T-z)+THETA_DOT/w1*sin(T-z)]*u
  //     -{[w2*THETA_DOT/(w1*w1)*cos(THETA)-THETA_DOT/w1]*cos(T-z)
  //                                       +THETA_DDOT   *sin(T-z)}*x
  //
  // Fz: -2*[w2/w1*cos(THETA)+1]*v
  //     -2*[w2/w1*sin(THETA)*sin(T-z)-THETA_DOT/w1*cos(T-z)]*u
  //     -{[w2*THETA_DOT/(w1*w1)*cos(THETA)-THETA_DOT/w1]*sin(T-z)
  //                                       -THETA_DDOT   *cos(T-z)}*x
  //
  // The model adopted for angular motion of spin axis is harmonic:
  // if t < T_RISE:
  //
  // THETA      = 0.5*      ANG_MAX*(1-cos(w3*t)), otherwise THETA = ANG_MAX
  // THETA_DOT  = 0.5*   w3*ANG_MAX*   sin(w3*t),  otherwise THETA_DOT  = 0
  // THETA_DDOT = 0.5*w3*w3*ANG_MAX*   cos(w3*t),  otherwise THETA_DDOT = 0
  //
  // where w3 = PI/T_RISE.

  static const real_t Tr    = Femlib::value ("T_RISE");
  static const real_t Th    = Femlib::value ("ANG_MAX");
  static const real_t w1    = Femlib::value ("OMEGA_1");
  static const real_t w2    = Femlib::value ("OMEGA_2");
  static const real_t w3    = Femlib::value ("PI/T_RISE");
  static const real_t dz    = Femlib::value ("TWOPI/BETA") / (nZ32 * nPR);
  static const real_t dt    = Femlib::value ("D_T");
  static const real_t w2ow1 = w2/w1;

  const real_t t            = Femlib::value ("t") - dt;
  const real_t T            = w1 * t;
  const int  transient    = t < Tr;

  real_t theta, thetaD, thetaDD;
  
  if (transient) {
    theta   = 0.5*Th*(1.0-cos(w3*t));
    thetaD  = 0.5*w3*Th*sin(w3*t);
    thetaDD = 0.5*w3*w3*Th*cos(w3*t);
  } else {
    theta   = Th;
    thetaD  = thetaDD = 0.0;
  }

  // -- V & w cross terms.

  const real_t cross = 2.0*w2ow1*cos(theta)+1.0;
 
  Blas::axpy (nTot32,  cross, u32[2], 1, n32[1], 1);
  Blas::axpy (nTot32, -cross, u32[1], 1, n32[2], 1);

  // -- Everything else requires z position.

  const real_t cost = cos (T);
  const real_t sint = sin (T);
  const real_t x1   =  2.0*w2ow1*sin(theta)+thetaDD;
  const real_t x2   =  w2*thetaD/(w1*w1)*cos(theta)-thetaD/w1;
  const real_t x3   =  2.0*(w2ow1*sin(theta)*cost-thetaD/w1*sint);
  const real_t x4   =  2.0*(w2ow1*sin(theta)*sint-thetaD/w1*cost);

  real_t z, costmz, sintmz;

  for (i = 0; i < nZ32; i++) {

    z      = dz * (i + nZ32 * Geometry::procID());
    costmz = cos (T - z);
    sintmz = sin (T - z);

    Veclib::fill (nP, 1.0, tmp, 1);
    master -> mulY (1, tmp);
  
    // -- Axial momentum.

    Blas::axpy (nP, x1*sintmz+x2*costmz, tmp,         1, n32[0]+i*nP, 1);
    Blas::axpy (nP, x3*cos(z)+x4*sin(z), u32[1]+i*nP, 1, n32[0]+i*nP, 1);
    Blas::axpy (nP, x4*cos(z)-x3*sin(z), u32[2]+i*nP, 1, n32[0]+i*nP, 1);

    Veclib::fill (nP, 1.0, tmp, 1);
    master -> mulX (1, tmp);

    // -- Radial momentum.

    Blas::axpy (nP, -2.0*(w2ow1*sin(theta)*costmz+thetaD/w1*sintmz),
		u32[0]+i*nP, 1, n32[1]+i*nP, 1);
    if (transient)
      Blas::axpy (nP, -x2*costmz-thetaDD*sintmz, tmp, 1, n32[1]+i*nP, 1);

    // -- Angular momentum.

    Blas::axpy (nP, -2.0*(w2ow1*sin(theta)*sintmz+thetaD/w1*costmz),
		u32[0]+i*nP, 1, n32[2]+i*nP, 1);
    if (transient)
      Blas::axpy (nP, -x2*sintmz+thetaDD*costmz, tmp, 1, n32[2]+i*nP, 1);
  }

  for (i = 0; i < NDIM; i++) {
    N[i]   -> transform32 (FORWARD, n32[i]);
    master -> smooth (N[i]);
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
  vector<AuxField*> H (NDIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NDIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);
  Blas::scal (Je, Femlib::value ("D_T"), &beta[0],  1);

  for (i = 0; i < NDIM; i++)
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
//                    u^^ = u^ - D_T * grad P;
// u^^ is left in Uf.
//
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to
// create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t dt    = Femlib::value ("D_T");
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < NDIM; i++) {

    (*Uf[i] = *D -> u[NDIM]) . gradient (i);

    if (C3D) Uf[2] -> divY();

    Us[i] -> axpy (-dt, *Uf[i]);
    Field::swapData (Us[i], Uf[i]);

    *Uf[i] *= alpha;
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
  const real_t            beta   = Femlib::value  ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [(size_t) (NDIM + 1)];
  int_t                   i;
  vector<real_t>          alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  const real_t            lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NDIM; i++)
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i],    (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NDIM] = new Msys
      (0.0,     beta, base, nmodes, E, D -> b[NDIM], (itLev<2)?DIRECT:JACPCG);

  return M;
}


static void Solve (Domain*       D,
		   const int_t i,
		   AuxField*     F,
		   Msys*         M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D->u[i], using F as a forcing Field.
// Iterative or direct solver selected on basis of field type, step,
// time order and command-line arguments.
// ---------------------------------------------------------------------------
{
  const int_t step = D -> step;

  if (i < NDIM && step < NORD) { // -- We need a temporary matrix system.
    const int_t Je      = min (step, NORD);    
    const int_t base    = Geometry::baseMode();
    const int_t nmodes  = Geometry::nModeProc();
    vector<real_t>  alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real_t    beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
