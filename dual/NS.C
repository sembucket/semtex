///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
// Geometries may be 2- or 3-dimensional, Cartesian or cylindrical.
// Fourier expansions are used in the homogeneous direction.
//
// Copyright (c) 2000 <--> $Date$, Hugh Blackburn
//
// NB: Modified for use with dual.  Nonlinear terms are computed using
// convolution sums since fields are in Fourier space always.
//
// References:
// 1.  Karniadakis, Israeli & Orszag 1991.  "High-order splitting methods
//     for the incompressible Navier--Stokes equations", JCP 9(2).
// 2. Guermond & Shen (2003) "Velocity correction projection methods for 
//    incompressible flows", SIAM J Numer Anal 41(1).
// 3. Blackburn & Sherwin (2004) "Formulation of a Galerkin spectral
//    element--Fourier method for three-dimensional incompressible flows
//    in cylindrical geometries", JCP.
//
// For cylindrical coordinates:
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "dual.h"

typedef ModalMatrixSys Msys;

static int_t NDIM, NORD, C3D;

static void   nonLinear (Domain*, AuxField**, AuxField**, vector<real_t>&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);


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
  NDIM = Geometry::nDim();
  NORD = Femlib::ivalue ("N_TIME");
  C3D  = Geometry::cylindrical() && NDIM == 3;

  int_t        i, j, k;
  const real_t dt     = Femlib:: value ("D_T");
  const int_t  nStep  = Femlib::ivalue ("N_STEP");
  const int_t  nZ     = Geometry::nZProc();
  const int_t  ntot   = Geometry::nTotProc();
  real_t*      alloc  = new real_t [static_cast<size_t> (2*NDIM*NORD*ntot)];

  // -- Create global matrix systems.

  Msys** MMS = preSolve (D);

  // -- Create, initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [static_cast<size_t> (2*NORD)];
  AuxField*** Uf = Us + NORD;

  for (k = 0, i = 0; i < NORD; i++) {
    Us[i] = new AuxField* [static_cast<size_t> (2*NDIM)];
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

  vector<real_t> ff (3);

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

    // -- Compute nonlinear terms from previous velocity field.

    nonLinear (D, Us[0], Uf[0], ff);

    // -- Update high-order pressure BC storage.

    PBCmgr::maintain (D -> step, Pressure,
		      const_cast<const AuxField**>(Us[0]),
		      const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (D -> step);

    // -- Complete unconstrained advective substep.

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }
    waveProp  (D, const_cast<const AuxField***>(Us),
	          const_cast<const AuxField***>(Uf));
    for (i = 0; i < NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm (Uf, NORD, NDIM);

    // -- Compute pressure.

    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NDIM, Uf[0][0], MMS[NDIM]);

    // -- Correct velocities for pressure.

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

    A -> analyse (Us[0], Uf[0]);
  }
}


static void nonLinear (Domain*         D ,
		       AuxField**      Us,
		       AuxField**      Uf,
		       vector<real_t>& ff)
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
// As for current semtex version, in cylindrical coordinates we
// actually compute y*Nx, y*Ny, Nz, as outlined in reference [3].   
//
// If STOKES is defined for compilation, the nonlinear terms are set to zero.
//
// For use with dual, velocity components are always held in the
// Fourier-transformed state, and non-linear terms are computed using
// convoultion sums.  There are only two modes (0 and the
// "Fundamental").  The nonlinear terms are dealiased.  We can assume
// that we are always operating in 3 space dimensions (and have 3
// velocity components).
// ---------------------------------------------------------------------------
{
  int_t        i, j;
  const real_t EPS = (sizeof(real_t) == sizeof(double)) ? EPSDP : EPSSP;

#if defined (STOKES)

  for (i = 0; i < NDIM; i++) {
    *Uf[i] = 0.0;
    AuxField::swapData (D -> u[i], Us[i]);
    if (fabs (ff[i]) > EPS) Uf[i] -> addToPlane (0, -ff[i]);
  }

#else

  vector<AuxField*> U (NDIM), N (NDIM), T (NDIM);
  Field*            master = D -> u[0];

  for (i = 0; i < NDIM; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
    U[i] = Us[i];
    N[i] = Uf[i];
    T[i] = D -> u[i];
  }

  if (Geometry::cylindrical()) {

#if defined (SKEW)
     N[0] -> convolve (*U[0], *U[1]);
     N[1] -> convolve (*U[1], *U[1]);
     N[2] -> convolve (*U[1], *U[2]);
     T[0] -> convolve (*U[2], *U[2]);
    *N[2] *=  3.0;
    *T[0] *= -2.0;
    *N[1] += *T[0];
#else
    *N[0]  = 0.0;
     N[1] -> convolve (*U[2], *U[2]);
    *N[1] *= -1.0;
     N[2] -> convolve (*U[1], *U[2]);
#endif

    for (i = 0; i < NDIM; i++) {

#if defined (SKEW)
      *T[0]  = *U[i];
       T[0] -> gradient (2);
       T[0] -> convolve (*U[2], *T[0]);
       T[1] -> convolve (*U[2], *U[i]);
       T[1] -> gradient (2);
      *N[i] += *T[0];
      *N[i] += *T[1];
#else
      *T[0]  = *U[i];
       T[0] -> gradient (2);
       T[0] -> convolve (*U[2], *T[0]);
      *N[i] += *T[0];
#endif

      if (i == 2) N[i] -> divY();

      for (j = 0; j < 2; j++) {
      
	// -- Perform n_i += u_j d(u_i) / dx_j.

	*T[0]  = *U[i];
	 T[0] -> gradient (j);
	 if (i < 2) T[0] -> mulY();
	 T[1] -> convolve (*U[j], *T[0]);
	*N[i] += *T[1];

	// -- Perform n_i += d(u_i u_j) / dx_j.
#if defined (SKEW)	
	 T[0] -> convolve (*U[i], *U[j]);
	 T[0] -> gradient (j);
	 if (i < 2) T[0] -> mulY();
	*N[i] += *T[0];
#endif
      }

      // -- Smooth, add forcing.
      
      master -> smooth (N[i]);
#if defined (SKEW) 		//  -- Check forcing - does it need mult by Y?
      if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
#else
      if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -ff[i]);
      *N[i] *= -1.0;
#endif
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NDIM; i++) {
      *N[i] = 0.0;

      for (j = 0; j < NDIM; j++) {
      
	// -- Perform n_i += u_j d(u_i) / dx_j.

	*T[0] = *U[i];
	 T[0] -> gradient (j);
	 T[1] -> convolve (*U[j], *T[0]);
	*N[i] += *T[1];

	// -- Perform n_i += d(u_i u_j) / dx_j.

#if defined (SKEW)
	 T[0] -> convolve (*U[i], *U[j]);
	 T[0] -> gradient (j);
	*N[i] += *T[0];
#endif
      }

      // -- Smooth, add forcing.
      
      master -> smooth (N[i]);
#if defined (SKEW)
      if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
#else
      if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -ff[i]);
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
  for (i = 1; i < NDIM; i++) *Uf[0] += *Uf[i];

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
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");

  for (i = 0; i < NDIM; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i == 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
    (*Us[0] = *D -> u[NDIM]) . gradient (i);
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
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const real_t            beta   = Femlib:: value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [static_cast<size_t> (NDIM + 1)];
  int_t                   i;

  vector<real_t>          alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  const real_t            lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NDIM; i++)
    M[i] = new Msys (lambda2, beta, 0, 2, E, D->b[i], (itLev<1)?DIRECT:JACPCG);

  // -- Pressure system.

  M[NDIM] = new Msys (0.0, beta, 0, 2, E, D->b[NDIM], (itLev<2)?DIRECT:JACPCG);

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

  if (i < NDIM && step < NORD) { // -- We need a temporary matrix system.
    const int_t    Je = min (step, NORD);    
    vector<real_t> alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real_t   beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys (lambda2, beta, 0, 2, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
