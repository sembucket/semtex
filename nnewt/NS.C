///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// Copyright (c) 1998 <--> $Date$,
//   Murray Rudman, Hugh Blackburn
//
// This version incorporates spatially varying generalized Newtonian
// viscosity.  This is handled within the framework of the DNS solver
// by splitting the viscous stress into two parts: one which can be
// dealt with implicitly as diffusion of momentum with a
// spatially-constant kinematic viscosity KINVIS, the other dealt with
// explicitly as the divergence of a stress.  This way the velocity
// BCs are still set correctly in the viscous substep.
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
// For testing, define NOMODEL during compilation.  This will set the
// kinematic viscosity equal to twice the input value and formulate
// the SGS stresses as those supplied by the product of the strain
// rate tensor with a spatially-constant NEGATIVE viscosity -KINVIS:
// the results should be the same as for a DNS with total viscosity
// KINVIS.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <nnewt.h>

typedef ModalMatrixSys Msys;
static  int_t          NORD, NDIM, NCOM;
static  bool           C3D;

static void   nonLinear (Domain*, AuxField**, AuxField**,
			 AuxField*, vector<real_t>&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);


void NavierStokes (Domain*        D,
		   nnewtAnalyser* A)
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
  NORD = Femlib::ivalue ("N_TIME");
  C3D  = Geometry::cylindrical() && NCOM == 3;

  int_t        i, j, k;
  const real_t dt    = Femlib:: value ("D_T");
  const int_t  nStep = Femlib::ivalue ("N_STEP");
  const int_t  nZ    = Geometry::nZProc();
  const int_t  nP    = Geometry::planeSize();
  const int_t  ntot  = Geometry::nTotProc();
  real_t*      alloc = new real_t [static_cast<size_t>(2*NCOM + 1)*NORD*ntot];

  // -- Create global matrix systems: rename viscosities beforehand.

  if (Femlib::value ("REFVIS") > 0.0) {
    real_t kinVis = Femlib::value ("REFVIS");
    real_t refVis = Femlib::value ("KINVIS");
    Femlib::value ("KINVIS", kinVis);
    Femlib::value ("REFVIS", refVis);
  } 

  Msys** MMS = preSolve (D);

  // -- Create, initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [static_cast<size_t> (2*NORD)];
  AuxField*** Uf = Us + NORD;

  for (k = 0, i = 0; i < NORD; i++) {
    Us[i] = new AuxField* [static_cast<size_t> (2*NCOM)];
    Uf[i] = Us[i] + NCOM;
    for (j = 0; j < NCOM; j++) {
      *(Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
      *(Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt)) = 0.0;
    }
  }

  // -- Create multi-level storage for pressure BCS.

  Field* Pressure = D -> u[NCOM];
  PBCmgr::build (Pressure);

  // -- Create spatially-constant forcing terms.

  vector<real_t> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Create and initialize non-Newtonian viscosity storage.
  
  AuxField* NNV = new AuxField (alloc + k * ntot, nZ, D -> elmt, 'n');
  *NNV = 0.0;
  ROOTONLY NNV -> addToPlane (0, Femlib::value ("-KINVIS"));

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Compute spatially-varying kinematic eddy viscosity \epsilon.

    viscosity (D, Us[0], Uf[0], NNV);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us[0], Uf[0], NNV, ff);

    // -- Update high-order pressure BC storage.

    PBCmgr::maintain (D -> step, Pressure, 
		      const_cast<const AuxField**> (Us[0]),
		      const_cast<const AuxField**> (Uf[0]));
    Pressure -> evaluateBoundaries (D -> step);

    // -- Complete unconstrained advective substep.

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }
    waveProp (D, const_cast<const AuxField***> (Us),
	         const_cast<const AuxField***> (Uf));
    for (i = 0; i < NCOM; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NCOM);

    // -- Compute pressure.

    setPForce (const_cast<const AuxField**> (Us[0]), Uf[0]);
    Solve     (D, NCOM,  Uf[0][0], MMS[NCOM]);
    
    // -- Correct velocities for pressure gradient.

    project (D, Us[0], Uf[0]);

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

    A -> analyse (Us[0], Uf[0], NNV);
  }

  // -- Dump non-Newtonian viscosity to file visco.fld.

  ofstream          nnvfl;
  vector<AuxField*> visco (1);

  visco[0] = NNV;

  ROOTONLY {
    nnvfl.open ("visco.fld", ios::out);
    NNV -> addToPlane (0, Femlib::value ("KINVIS"));
  }
  NNV -> transform (INVERSE);

  writeField (nnvfl, D -> name, D -> step, D -> time, visco);

  ROOTONLY nnvfl.close();
}


static void nonLinear (Domain*         D ,
		       AuxField**      Us,
		       AuxField**      Uf,
		       AuxField*       EV,
		       vector<real_t>& ff)
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
// Refer to Bird, Stewart & Lightfoot (1960) Appendix A for
// stress-divergence operations in cylindrical & Cartesian
// coordinates.
//
// NB: dealiasing is not used for multiprocessor operation.
// ---------------------------------------------------------------------------
{
  int_t           i, j;
  const real_t    EPS    = (sizeof(real_t) == sizeof(double)) ? EPSDP : EPSSP;
  const int_t     nZ     = Geometry::nZ();
  const int_t     nZP    = Geometry::nZProc();
  const int_t     nP     = Geometry::planeSize();
  const int_t     nPP    = Geometry::nBlock();
  const int_t     nPR    = Geometry::nProc();
  const int_t     nTot   = nZP * nP;
#if defined (ALIAS)
  const int_t     nZ32   = nZP; assert (Geometry::nProc() == 1);
#else
  const int_t     nZ32   = (nPR > 1) ? nZP : (3 * nZ) >> 1;
#endif
  const int_t     nTot32 = nZ32 * nP;
  vector<real_t>  work ((2 * NCOM + 1) * nTot32);
  vector<real_t*> u32(NCOM), n32(NCOM);
  Field*          master = D -> u[0];
  real_t*         tmp    = &work[0] + 2 * NCOM * nTot32;

  for (i = 0; i < NCOM; i++) {
    u32[i] = &work[0] +  i         * nTot32;
    n32[i] = &work[0] + (i + NCOM) * nTot32;
  }

  Veclib::zero (3*nTot32, n32[0], 1);

  // -- Nonlinear terms proper.
 
  for (i = 0; i < NCOM; i++) D -> u[i] -> transform32 (INVERSE, u32[i]);

  if (Geometry::cylindrical()) {	// -- Cylindrical coordinates.

    for (i = 0; i < NCOM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (i == 0)
	Veclib::vmul (nTot32, u32[0], 1, u32[1], 1, n32[0], 1);
      if (i == 1)
	Veclib::vmul (nTot32, u32[1], 1, u32[1], 1, n32[1], 1);

      if (NCOM == 3) {
	if (i == 1)
	  Veclib::svvttvp (nTot32, -2.0, u32[2],1,u32[2],1,n32[1],1,n32[1], 1);
	if (i == 2)
	  Veclib::svvtt   (nTot32,  3.0, u32[2],1,u32[1],1,         n32[2], 1);

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

      if (i == 2) master -> divY (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	if (i < 2) master -> mulY (nZ32, tmp);
	Veclib::vvtvp (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- 2D conservative derivatives.
     
      for (j = 0; j < 2; j++) {
	Veclib::vmul (nTot32, u32[j], 1, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	if (i < 2) master -> mulY (nZ32, tmp);
	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }
    }

  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NCOM; i++) {
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
    }
  }

  // -- Add contribution from divergence of SGS.

  EV -> transform32 (INVERSE, tmp);
  Blas::scal (nTot32, -4.0, tmp, 1); // -- 2 = -0.5 * -4 ... see below.

  if (Geometry::cylindrical()) {	// -- Cylindrical coordinates.

    // -- 2D stress-divergence terms.

    Us[0] -> transform32 (INVERSE, u32[0]);
    Veclib::vmul         (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
    master -> gradient   (nZ32, nP, u32[0], 0);
    master -> mulY       (nZ32, u32[0]);
    Veclib::vadd         (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);

    Us[1] -> transform32 (INVERSE, u32[0]);
    Veclib::vmul         (nTot32, tmp,    1, u32[0], 1, u32[0], 1);
    Veclib::vadd         (nTot32, n32[1], 1, u32[0], 1, n32[1], 1);
    master -> gradient   (nZ32, nP, u32[0], 1);
    master -> mulY       (nZ32,   u32[0]);
    Veclib::vadd         (nTot32, n32[1], 1, u32[0], 1, n32[1], 1);

    Uf[0] -> transform32 (INVERSE, u32[0]);
    Veclib::vmul         (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
    Veclib::copy         (nTot32, u32[0], 1, u32[1], 1);
    Veclib::vadd         (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
    master -> gradient   (nZ32, nP, u32[0], 0);
    master -> gradient   (nZ32, nP, u32[1], 1);
    master -> mulY       (nZ32, u32[0]);
    master -> mulY       (nZ32, u32[1]);
    Veclib::vadd         (nTot32, u32[0], 1, n32[1], 1, n32[1], 1);
    Veclib::vadd         (nTot32, u32[1], 1, n32[0], 1, n32[0], 1);

    // -- 3D stress-divergence terms.

    if (C3D) {		

      Us[2] -> transform32 (INVERSE, u32[0]);
      Veclib::vmul         (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::vsub         (nTot32, n32[1], 1, u32[0], 1, n32[1], 1);
      Femlib::exchange     (u32[0], nZ32,        nP, FORWARD);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, FORWARD);
      Veclib::zero         (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient   (nZ, nPP, u32[0], 2);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, INVERSE);
      Femlib::exchange     (u32[0], nZ32,        nP, INVERSE);
      Veclib::vadd         (nTot32, n32[2], 1, u32[0], 1, n32[2], 1);

      Uf[1] -> transform32 (INVERSE, u32[0]);
      Veclib::vmul         (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy         (nTot32,         u32[0], 1, u32[1], 1);
      Femlib::exchange     (u32[0], nZ32,        nP, FORWARD);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, FORWARD);
      Veclib::zero         (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient   (nZ, nPP, u32[0], 2);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, INVERSE);
      Femlib::exchange     (u32[0], nZ32,        nP, INVERSE);
      Veclib::vadd         (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
      master -> gradient   (nZ32, nP, u32[1], 0);
      master -> mulY       (nZ32, u32[1]);
      Veclib::vadd         (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);

      Uf[2] -> transform32 (INVERSE, u32[0]);
      Veclib::vmul         (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy         (nTot32,         u32[0], 1, u32[1], 1);
      Blas::axpy           (nTot32, 2.0,    u32[0], 1, n32[2], 1);
      Femlib::exchange     (u32[0], nZ32,        nP, FORWARD);      
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, FORWARD);
      Veclib::zero         (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient   (nZ, nPP, u32[0], 2);
      Femlib::DFTr         (u32[0], nZ32 * nPR, nPP, INVERSE);
      Femlib::exchange     (u32[0], nZ32,        nP, INVERSE);
      Veclib::vadd         (nTot32, u32[0], 1, n32[1], 1, n32[1], 1);
      master -> gradient   (nZ32, nP, u32[1], 1);
      master -> mulY       (nZ32, u32[1]);
      Veclib::vadd         (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);

    }

  } else {			// -- Cartesian coordinates.

    // -- Diagonal stress-divergence terms.

    for (i = 0; i < NCOM; i++) {

      Us[i] -> transform32 (INVERSE, u32[i]);
      Veclib::vmul (nTot32, tmp, 1, u32[i], 1, u32[i], 1);

      if (i == 2) {
	Femlib::exchange   (u32[2], nZ32,        nP, FORWARD);      
	Femlib::DFTr       (u32[2], nZ32 * nPR, nPP, FORWARD);
	Veclib::zero       (nTot32 - nTot, u32[2] + nTot, 1);
	master -> gradient (nZ, nPP, u32[2], 2);
	Femlib::DFTr       (u32[2], nZ32 * nPR, nPP, INVERSE);
	Femlib::exchange   (u32[2], nZ32,        nP, INVERSE);      
      } else
	master -> gradient (nZ32, nP, u32[i], i);

      Veclib::vadd (nTot32, u32[i], 1, n32[i], 1, n32[i], 1);
    }

    // -- Off-diagonal stress-divergence terms.

    for (i = 0; i < NCOM; i++)
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

  for (i = 0; i < NCOM; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
    Uf[i]  -> transform32 (FORWARD, n32[i]);
    master -> smooth (Uf[i]);
    ROOTONLY if (fabs (ff[i]) > EPS)
      if (Geometry::cylindrical()) {
	Veclib::fill (nP, -2.0*ff[i], tmp, 1);
	if (i < 2) master -> mulY (1, tmp);
	Uf[i] -> addToPlane (0, tmp);
      } else
	Uf[i] -> addToPlane (0, -2.0*ff[i]);
    *Uf[i] *= -0.5;
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
  vector<AuxField*> H (NCOM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NCOM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t    Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);
  
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
//
//                    u^^ = u^ - D_T * grad P,
//
// then scale by -1.0 / (D_T * KINVIS) to create forcing for viscous step.
//
// u^^ is left in Uf.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");

  for (i = 0; i < NCOM; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i == 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
  }

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[NCOM]) . gradient (i);
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
  const int_t             nmodes = Geometry::nModeProc();
  const int_t             base   = Geometry::baseMode();
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const real_t            beta   = Femlib:: value ("BETA");
  const vector<Element*>& E = D -> elmt;
  Msys**                  M = new Msys* [static_cast<size_t>(NCOM + 1)];
  int_t                   i;

  vector<real_t> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  const real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

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

  if (i < NCOM && step < NORD) { // -- We need a temporary matrix system.
    const int_t    Je     = min (step, NORD);    
    const int_t    base   = Geometry::baseMode();
    const int_t    nmodes = Geometry::nModeProc();

    vector<real_t> alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real_t   beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
