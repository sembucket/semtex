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

static integer NDIM, NORD, CYL, C3D;

static void   nonLinear (Domain*, AuxField**, AuxField**);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const integer, AuxField*, Msys*);


void NavierStokes (Domain*      D      ,
		   DNSAnalyser* A      ,
		   const char*  masklag)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NDIM = Geometry::nDim();
  NORD = (integer) Femlib::value ("N_TIME");
  CYL  = Geometry::system() == Geometry::Cylindrical;
  C3D  = CYL && NDIM == 3;

  integer       i, j, k;
  const real    dt     =           Femlib::value ("D_T");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();
  const integer nP     = Geometry::planeSize();
  const integer ntot   = Geometry::nTotProc();
  real*         alloc  = new real [(size_t) 2 * NDIM * NORD * ntot];

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

  // -- Create filter mask.

  AuxField* mask = 0;
  if (masklag) {
    ROOTONLY cout << "-- Setting up filter mask." << endl;
    mask = new AuxField ((alloc = new real [(size_t) ntot]), nZ, D->elmt, 'm');
    mask -> buildMask (masklag);
    ROOTONLY {
      Veclib::fill (4*nP, 1.0, alloc,      1);
      Veclib::zero (  nP,      alloc + nP, 1);
    }
  }

  // -- Apply coupling to radial & azimuthal velocity BCs.

  Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us[0], Uf[0]);
#if 1
    // -- Apply masking.

    if (mask) for (i = 0; i < NDIM; i++) *Uf[0][i] *= *mask;
#endif
    waveProp  (D, Us, Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure, Us[0], Uf[0]);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NDIM);
    setPForce (Us[0], Uf[0]);
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

    A -> analyse (Us[0]);
  }
}


static void nonLinear (Domain*       D ,
		       AuxField**    Us,
		       AuxField**    Uf)
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
// Forcing for inertia-wave problem:
//
// Fx: -2*w2/w1*sin(THETA)*{ y*cos(t+z)                      + 
//                           [sin(z)*cos(t)+cos(z)*sin(t)]*v +
//                           [cos(z)*cos(t)-sin(z)*sin(t)]*w }
// Fy: +2*[1+w2/w1*cos(THETA)]*w + 2*w2/w1*sin(THETA)*sin(t+z)*u
// Fz: -2*[1+w2/w1*cos(THETA)]*v + 2*w2/w1*sin(THETA)*cos(t+z)*u
//
// NB: no dealiasing for concurrent execution.
// ---------------------------------------------------------------------------
{
  integer i, j;

#if defined(STOKES)

  for (i = 0; i < NDIM; i++) {
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
  const integer     nZ32   = Geometry::nZ32();
  const integer     nTot32 = nZ32 * nP;
  vector<real>      work ((2 * NDIM + 1) * nTot32);
  vector<real*>     u32 (NDIM);
  vector<real*>     n32 (NDIM);
  vector<AuxField*> U   (NDIM);
  vector<AuxField*> N   (NDIM);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * NDIM * nTot32;

  Veclib::zero ((2 * NDIM + 1) * nTot32, work(), 1); // -- A catch-all cleanup.

  for (i = 0; i < NDIM; i++) {
    u32[i] = work() +  i         * nTot32;
    n32[i] = work() + (i + NDIM) * nTot32;

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

    // -- Average NL terms.

    Blas::scal (-0.5, nTot32, n32[i], 1);
  }

#endif      
  
  // -- Add forcing for inertia-wave problem:
  //
  // Fx: -2*w2/w1*sin(THETA)*{ y*cos(t+z)                      + 
  //                           [sin(z)*cos(t)+cos(z)*sin(t)]*v +
  //                           [cos(z)*cos(t)-sin(z)*sin(t)]*w }
  // Fy: +2*[1+w2/w1*cos(THETA)]*w + 2*w2/w1*sin(THETA)*sin(t+z)*u
  // Fz: -2*[1+w2/w1*cos(THETA)]*v + 2*w2/w1*sin(THETA)*cos(t+z)*u

  real cos0, sin0, z, A1ctcz, A1stcz, A1ctsz, A1stsz, A1ctpz, A1stpz;

  if ((integer) Femlib::value ("LINEAR")) {
    cos0 = 1.0;
    sin0 = Femlib::value ("THETA");
  } else {
    cos0 = Femlib::value ("cos(THETA)");
    sin0 = Femlib::value ("sin(THETA)");
  }

  const real t     = Femlib::value ("t");
  const real cost  = cos (t);
  const real sint  = sin (t);
  const real omega = Femlib::value ("OMEGA");
  const real w2ow1 = 2.0 * omega - 1.0;
  const real A1    = 2.0 * w2ow1 * sin0;
  const real A2    = 2.0 * (1.0 + w2ow1 * cos0);
  const real dz    = Femlib::value ("TWOPI/BETA") / (nZ32 * nPR);

  // -- V & w cross terms.
 
  Blas::axpy (nTot32,  A2, u32[2], 1, n32[1], 1);
  Blas::axpy (nTot32, -A2, u32[1], 1, n32[2], 1);

  // -- Now the things that require us to know z position.

  for (i = 0; i < nZ32; i++) {

    z      = dz * (i + nZ32 * Geometry::procID());
    A1ctcz = A1 * cost * cos (z);
    A1stcz = A1 * sint * cos (z);
    A1ctsz = A1 * cost * sin (z);
    A1stsz = A1 * sint * sin (z);
    A1ctpz = A1 * cos (t + z);
    A1stpz = A1 * sin (t + z);

    Veclib::fill (nP, 1.0, tmp, 1);
    master -> mulR (1, tmp);
    Blas::axpy (nP, -A1ctpz,        tmp,    1, n32[0] + i * nP, 1);
    Blas::axpy (nP, -A1stcz-A1ctsz, u32[1], 1, n32[0] + i * nP, 1);
    Blas::axpy (nP, -A1ctcz+A1stsz, u32[2], 1, n32[0] + i * nP, 1);
							
    Blas::axpy (nP,  A1stpz,        u32[0], 1, n32[1] + i * nP, 1);
					
    Blas::axpy (nP,  A1ctpz,        u32[0], 1, n32[2] + i * nP, 1);
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
  integer           i, q;
  vector<AuxField*> H (NDIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NDIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const integer Je = min (D -> step, NORD);
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

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
// level u^ is stored in Us.  Constrain velocity field:
//                    u^^ = u^ - D_T * grad P;
// u^^ is left in Uf.
//
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to
// create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt    = Femlib::value ("D_T");
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < NDIM; i++) {

    (*Uf[i] = *D -> u[NDIM]) . gradient (i);

    if (C3D) Uf[2] -> divR();

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
  const integer           nmodes = Geometry::nModeProc();
  const integer           base   = Geometry::baseMode();
  const integer           itLev  = (integer) Femlib::value ("ITERATIVE");
  const real              beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [(size_t) (NDIM + 1)];
  integer                 i;
  vector<real>            alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, alpha());
  const real              lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

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

  if (i < NDIM && step < NORD) { // -- We need a temporary matrix system.
    const integer Je      = min (step, NORD);    
    const integer base    = Geometry::baseMode();
    const integer nmodes  = Geometry::nModeProc();
    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
    const real    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real    beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
