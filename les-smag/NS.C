///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// Copyright (c) Hugh Blackburn 1998--1999.
//
// This version incorporates LES.  This is handled within the
// framework of the DNS solver by splitting the SGS stress into two
// parts: one which can be dealt with implicitly as diffusion of
// momentum with a spatially-constant kinematic viscosity KINVIS, the
// other dealt with explicitly as the divergence of a stress.  This
// way the velocity BCs are still set correctly in the viscous
// substep.  The splitting is convenient for eddy-viscosity approaches
// to computing SGS stresses, but may need re-evaluation for other
// techniques.
//
// The splitting technique has been discussed by
//
//  @InCollection{koy93,
//  author = 	 "George E. Karnidakis and Steven A. Orszag and Victor Yakhot",
//  title = 	 "Renormalization Group Theory Simulation of Transitional
//                and Turbulent Flow over a Backward-Facing Step",
//  booktitle =	 "Large Eddy Simulation of Complex Engineering and
//		  Geophysical Flows",
//  publisher =	 "Cambridge",
//  year =	 1993,
//  editor =	 "Boris Galperin and Steven A. Orszag",
//  chapter =	 8,
//  pages =	 "159--177"
//  }
//
// For testing, define DEBUG during compilation.  This will set the
// kinematic viscosity equal to twice the input value and formulate
// the SGS stresses as those supplied by the product of the strain
// rate tensor with a spatially-constant NEGATIVE viscosity -KINVIS:
// the results should be the same as for a DNS with total viscosity
// KINVIS.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <les.h>

typedef ModalMatrixSys Msys;
static  integer        NORD, NDIM, CYL, C3D;

static void   nonLinear (Domain*, AuxField**, AuxField**,
			 AuxField*, vector<real>&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const integer, AuxField*, Msys*);


void NavierStokes (Domain*      D,
		   LESAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NORD = (integer) Femlib::value ("N_TIME");
  NDIM = Geometry::nDim();
  CYL  = Geometry::system() == Geometry::Cylindrical;
  C3D  = CYL && NDIM == 3;

  integer       i, j, k;
  const real    dt     =           Femlib::value ("D_T"   );
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();
  const integer ntot   = Geometry::nTotProc();
  real*         alloc  = new real [(size_t) (2 * NDIM + 1) * NORD * ntot];

  // -- Create global matrix systems: rename viscosities beforehand.

#if defined(DEBUG)
  Femlib::value ("REFVIS", Femlib::value ("2.0 * KINVIS"));
#endif

  if (Femlib::value ("REFVIS") > 0.0) {
    real kinVis = Femlib::value ("REFVIS");
    real refVis = Femlib::value ("KINVIS");
    Femlib::value ("KINVIS", kinVis);
    Femlib::value ("REFVIS", refVis);
  } 

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

  // -- Create and initialize eddy-viscosity storage.

  AuxField* EV = new AuxField (alloc + k * ntot, nZ, D -> elmt, 'e');
  *EV = 0.0;
  ROOTONLY EV -> addToPlane (0, Femlib::value ("REFVIS - KINVIS"));

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Compute spatially-varying kinematic eddy viscosity \epsilon.

    eddyViscosity (D, Us[0], Uf[0], EV);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us[0], Uf[0], EV, ff);
    waveProp  (D, Us, Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure, Us[0], Uf[0]);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm (Uf, NORD, NDIM);
  
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

  // -- Dump ratio eddy/molecular viscosity to file visco.fld.

  ofstream          evfl;
  vector<AuxField*> visco (1);

  visco[0] = EV;

  ROOTONLY {
    evfl.open ("visco.fld", ios::out);
    EV -> addToPlane (0, Femlib::value ("KINVIS"));
  }

  (*EV /= Femlib::value ("REFVIS")) . transform (INVERSE);

  writeField (evfl, D -> name, D -> step, D -> time, visco);

  ROOTONLY evfl.close();
}


static void nonLinear (Domain*       D ,
		       AuxField**    Us,
		       AuxField**    Uf,
		       AuxField*     EV,
		       vector<real>& ff)
// ---------------------------------------------------------------------------
// Compute nonlinear + forcing terms in Navier--Stokes equations
//                       N(u) + div(2*EV*S) + ff.
//
// Here N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS, div(2*EV*S) is the divergence of the non-constant
// component of the eddy-viscosity SGS terms and ff is a vector of body
// force per unit mass.
//
// On entry, D contains the old velocity (and pressure) fields, the
// lowest levels of Us & Uf contain the components of the old
// strain-rate tensor and vV contains the spatially-varying viscosity.
// On exit, the velocity field storage areas of D are free, the zeroth
// level of Us contains the old velocities and the zeroth level of Uf
// contains the most recent explicit forcing terms.  Velocity field
// data areas of D and first level of Us are swapped, then the next
// stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 ~ ~
//           N  = -0.5 ( u . grad u + div uu )
//           ~           ~        ~       ~~
//
// i.e. in Cartesian component form
//
//           N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ),
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
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// NB: dealiasing is not used for multiprocessor operation.
// ---------------------------------------------------------------------------
{
  integer           i, j;
  const real        EPS    = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;
  const integer     nZ     = Geometry::nZ();
  const integer     nZP    = Geometry::nZProc();
  const integer     nP     = Geometry::planeSize();
  const integer     nPP    = Geometry::nBlock();
  const integer     nPR    = Geometry::nProc();
  const integer     nTot   = nZP * nP;
  const integer     nZ32   = (nPR > 1) ? nZP : (3 * nZ) >> 1;
  const integer     nTot32 = nZ32 * nP;
  vector<real>      work ((2 * NDIM + 1) * nTot32);
  vector<real*>     u32 (NDIM);
  vector<real*>     n32 (NDIM);
  vector<AuxField*> U   (NDIM);
  vector<AuxField*> N   (NDIM);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * NDIM * nTot32;

  for (i = 0; i < NDIM; i++) {
    u32[i] = work() +  i         * nTot32;
    n32[i] = work() + (i + NDIM) * nTot32;
  }

  // -- Start with contribution from divergence of SGS.

  EV -> transform32 (INVERSE, tmp);
  Blas::scal (nTot32, -4.0, tmp, 1); // -- 2 = -0.5 * -4 ... see below.

  if (CYL) {			// -- Cylindrical coordinates.

    // -- 2D stress-divergence terms.

    Us[0] -> transform32 (INVERSE, n32[0]);
    Veclib::vmul         (nTot32,  tmp, 1, n32[0], 1, n32[0], 1);
    master   -> gradient (nZ32, nP, n32[0], 0);

    Us[1] -> transform32 (INVERSE, n32[1]);
    Veclib::vmul         (nTot32,  tmp, 1, n32[1], 1, n32[1], 1);
#if 1
// -- Testing new formulation.
    Veclib::copy         (nTot32, n32[1], 1, u32[0], 1);
    master   -> gradient (nZ32, nP, n32[1], 1);
    master   -> divR     (nZ32,   u32[0]);
    Veclib::vadd         (nTot32, n32[1], 1, u32[0], 1, n32[1], 1);
#else
    master   -> mulR     (nZ32, n32[1]);
    master   -> gradient (nZ32, nP, n32[1], 1);
    master   -> divR     (nZ32, n32[1]);
#endif

    Uf[0] -> transform32 (INVERSE, u32[0]);
    Veclib::vmul         (nTot32,  tmp, 1, u32[0], 1, u32[0], 1);
    Veclib::copy         (nTot32,          u32[0], 1, u32[1], 1);
#if 1
// -- Testing.
    master   -> gradient (nZ32, nP, u32[1], 0);
    Veclib::vadd         (nTot32,   u32[1], 1, n32[1], 1, n32[1], 1);
    Veclib::copy         (nTot32,              u32[0], 1, u32[1], 1);
    master   -> gradient (nZ32, nP, u32[0], 1);
    Veclib::vadd         (nTot32,   u32[0], 1, n32[0], 1, n32[0], 1);
    master   -> divR     (nZ32,     u32[1]);
    Veclib::vadd         (nTot32,   u32[1], 1, n32[0], 1, n32[0], 1);
#else
    master   -> mulR     (nZ32, u32[0]);
    master   -> gradient (nZ32, nP, u32[0], 1);
    master   -> divR     (nZ32, u32[0]);
    master   -> gradient (nZ32, nP, u32[1], 0);
    Veclib::vadd         (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
    Veclib::vadd         (nTot32, u32[1], 1, n32[1], 1, n32[1], 1);
#endif

    // -- 3D stress-divergence terms.

    if (C3D) {		

      Us[2] -> transform32 (INVERSE, n32[2]);
      Veclib::vmul         (nTot32, tmp, 1, n32[2], 1, n32[2], 1);
      Veclib::copy         (nTot32,         n32[2], 1, u32[0], 1);
      master -> divR       (nZ32, u32[0]);
      Veclib::vsub         (nTot32, n32[1], 1, u32[0], 1, n32[1], 1);
      Femlib::exchange     (n32[2], nZ32,        nP, FORWARD);
      Femlib::DFTr         (n32[2], nZ32 * nPR, nPP, FORWARD);
      Veclib::zero         (nTot32 - nTot, n32[2] + nTot, 1);
      master -> gradient   (nZ, nPP, n32[2], 2);
      Femlib::DFTr         (n32[2], nZ32 * nPR, nPP, INVERSE);
      Femlib::exchange     (n32[2], nZ32,        nP, INVERSE);
      master -> divR       (nZ32, n32[2]);

      Uf[1] -> transform32 (INVERSE, u32[0]);
      Veclib::vmul         (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy         (nTot32,         u32[0], 1, u32[1], 1);
      Femlib::exchange     (u32[0], nZ32,        nP, FORWARD);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, FORWARD);
      Veclib::zero         (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient   (nZ, nPP, u32[0], 2);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, INVERSE);
      Femlib::exchange     (u32[0], nZ32,        nP, INVERSE);
      master -> divR       (nZ32, u32[0]);
      master -> gradient   (nZ32, nP, u32[1], 0);
      Veclib::vadd         (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
      Veclib::vadd         (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);

      Uf[2] -> transform32 (INVERSE, u32[0]);
      Veclib::vmul         (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy         (nTot32,         u32[0], 1, u32[1], 1);
      Veclib::copy         (nTot32,         u32[0], 1, u32[2], 1);
      Femlib::exchange     (u32[0], nZ32,        nP, FORWARD);      
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, FORWARD);
      Veclib::zero         (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient   (nZ, nPP, u32[0], 2);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, INVERSE);
      Femlib::exchange     (u32[0], nZ32,        nP, INVERSE);      
      master -> divR       (nZ32, u32[0]);
      master -> gradient   (nZ32, nP, u32[1], 1);
      master -> divR       (nZ32, u32[2]);
      Veclib::vadd         (nTot32, u32[0], 1, n32[1], 1, n32[1], 1);
      Veclib::vadd         (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);
      Blas::axpy           (nTot32, 2.0, u32[2], 1,       n32[2], 1);
    }

  } else {			// -- Cartesian coordinates.

    // -- Diagonal stress-divergence terms.

    for (i = 0; i < NDIM; i++) {

      Us[i] -> transform32 (INVERSE, n32[i]);
      Veclib::vmul (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

      if (i == 2) {
	Femlib::exchange   (n32[2], nZ32,        nP, FORWARD);      
	Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, FORWARD);
	Veclib::zero       (nTot32 - nTot, n32[2] + nTot, 1);
	master -> gradient (nZ, nPP, n32[2], 2);
	Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, INVERSE);
	Femlib::exchange   (n32[2], nZ32,        nP, INVERSE);      
      } else
	master -> gradient (nZ32, nP, n32[i], i);
    }

    // -- Off-diagonal stress-divergence terms.

    for (i = 0; i < NDIM; i++)
      for (j = i + 1; j < NDIM; j++) {

	Uf[i + j - 1] -> transform32 (INVERSE, u32[0]);
	Veclib::vmul                 (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
	Veclib::copy                 (nTot32,         u32[0], 1, u32[1], 1);

	// -- Super-diagonal.

	if (j == 2) {
	  Femlib::exchange   (u32[0], nZ32,        nP, FORWARD);
	  Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, FORWARD);
	  Veclib::zero       (nTot32 - nTot, u32[0] + nTot, 1);
	  master -> gradient (nZ, nPP, u32[0], 2);
	  Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, INVERSE);
	  Femlib::exchange   (u32[0], nZ32,        nP, INVERSE);
	} else
	  master -> gradient (nZ32, nP, u32[0], j);
	Veclib::vadd (nTot32, u32[0], 1, n32[i], 1, n32[i], 1);

	// -- Sub-diagonal.

	master -> gradient (nZ32, nP, u32[1], i);
	Veclib::vadd       (nTot32, u32[1], 1, n32[j], 1, n32[j], 1);
      }
  }

  // -- Nonlinear terms.

  for (i = 0; i < NDIM; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
    U[i] = Us[i];
    N[i] = Uf[i];
    U[i] -> transform32 (INVERSE, u32[i]);
    if (CYL) N[i] -> mulR (nZ32, n32[i]);
  }
  
  if (CYL) {			// -- Cylindrical coordinates.

    for (i = 0; i < NDIM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (i == 0)
	Veclib::vvtvp (nTot32, u32[0], 1, u32[1], 1, n32[0], 1, n32[0], 1);
      if (i == 1)
	Veclib::vvtvp (nTot32, u32[1], 1, u32[1], 1, n32[1], 1, n32[1], 1);

      if (NDIM == 3) {
	if (i == 1)
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  3.0, u32[2], 1, u32[1], 1,      n32[2], 1);

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
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
    }
 
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NDIM; i++) {
      for (j = 0; j < NDIM; j++) {
      
	// -- Nonconservative contribution,  n_i += u_j d(u_i) / dx_j.

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

	// -- Conservative contribution,  n_i += d(u_i u_j) / dx_j.

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
      ROOTONLY if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
    }
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
