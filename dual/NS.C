///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
// Geometries may be 2- or 3-dimensional, Cartesian or cylindrical.
// Fourier expansions are used in the homogeneous direction.
//
// NB: Modified for use with dual.  Nonlinear terms are computed using
// convolution sums.
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

#include <dual.h>

typedef ModalMatrixSys Msys;

static integer NDIM, NORD, CYL, C3D;

static void   nonLinear (Domain*, AuxField**, AuxField**, vector<real>&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const integer, AuxField*, Msys*);


void NavierStokes (Domain*       D,
		   DualAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NDIM  = Geometry::nDim();
  NORD  = (integer) Femlib::value ("N_TIME");
  CYL   = Geometry::system() == Geometry::Cylindrical;
  C3D   = CYL && NDIM == 3;

  integer       i, j, k;
  const real    dt     =           Femlib::value ("D_T");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();
  const integer ntot   = Geometry::nTotProc();
  real*         alloc  = new real [(size_t) 2 * NDIM * NORD * ntot];

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

  // -- Create spatially-constant forcing terms.

  vector<real> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us[0], Uf[0], ff);
    waveProp  (D, Us, Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure, Us[0], Uf[0]);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NDIM);
    setPForce (Us[0], Uf[0]);
    Solve     (D, NDIM,  Uf[0][0], MMS[NDIM]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NDIM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NDIM);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NDIM; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0]);
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
// For use with dual, velocity components are always held in the
// Fourier-transformed state, and non-linear terms are computed using
// convoultion sums.  There are only two modes (0 and the
// "Fundamental").  The nonlinear terms are dealiased.  We can assume
// that we are always operating in 3 space dimensions.
// ---------------------------------------------------------------------------
{
  integer    i, j;
  const real EPS = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;

#if defined(STOKES)

  for (i = 0; i < NDIM; i++) {
    *Uf[i] = 0.0;
    AuxField::swapData (D -> u[i], Us[i]);
    if (fabs (ff[i]) > EPS) Uf[i] -> addToPlane (0, -ff[i]);
  }

#else

  const integer     nZ = Geometry::nZ();
  const integer     nP = Geometry::planeSize();
  vector<AuxField*> U (NDIM), N(NDIM), T(NDIM);
  Field*            master = D -> u[0];

  for (i = 0; i < NDIM; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
    U[i] = Us[i];
    N[i] = Uf[i];
    T[i] = D -> u[i];
  }

  if (CYL) {			// -- Cylindrical coordinates.

     N[0] -> convolve (*U[0], *U[1]);
     N[1] -> convolve (*U[1], *U[1]);
     N[2] -> convolve (*U[1], *U[2]);
     T[0] -> convolve (*U[2], *U[2]);
    *N[2] *=  3.0;
    *T[0] *= -2.0;
    *N[1] += *T[0];

    for (i = 0; i < NDIM; i++) {

      *T[0]  = *U[i];
       T[0] -> gradient (2);
       T[0] -> convolve (*U[2], *T[0]);
       T[1] -> convolve (*U[2], *U[i]);
       T[1] -> gradient (2);
      *N[i] += *T[0];
      *N[i] += *T[1];
       N[i] -> divR();

      for (j = 0; j < 2; j++) {
      
	// -- Perform n_i += u_j d(u_i) / dx_j.

	*T[0]  = *U[i];
	 T[0] -> gradient (j);
	 T[1] -> convolve (*U[j], *T[0]);
	*N[i] += *T[1];

	// -- Perform n_i += d(u_i u_j) / dx_j.
	
	 T[0] -> convolve (*U[i], *U[j]);
	 T[0] -> gradient (j);
	*N[i] += *T[0];

      }

      // -- Smooth, add forcing.
      
      master -> smooth (N[i]);
      if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NDIM; i++) {
      for (j = 0; j < NDIM; j++) {
      
	// -- Perform n_i += u_j d(u_i) / dx_j.

	*T[0] = *U[i];
	 T[0] -> gradient (j);
	 T[1] -> convolve (*U[j], *T[0]);
	*N[i] = *T[1];

	// -- Perform n_i += d(u_i u_j) / dx_j.

	 T[0] -> convolve (*U[i], *U[j]);
	 T[0] -> gradient (j);
	*N[i] += *T[0];

      }

      // -- Smooth, add forcing.
      
      master -> smooth (N[i]);
      if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
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
      (lambda2, beta, 0, 2, E, D -> b[i],    (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NDIM] = new Msys
      (0.0,     beta, 0, 2, E, D -> b[NDIM], (itLev<2)?DIRECT:JACPCG);

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
    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
    const real    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real    beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, 0, 2, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
