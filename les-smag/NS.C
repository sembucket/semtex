///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
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
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>

void eddyViscosity (const Domain*, AuxField***, AuxField***, AuxField*);

typedef ModalMatrixSystem  ModeSys;
static  integer            DIM, CYL, C3D;

static void      nonLinear (Domain*, AuxField***, AuxField***,
			    AuxField*, vector<real>&);
static void      waveProp  (Domain*, const AuxField***, const AuxField***);
static void      setPForce (const AuxField***, AuxField***);
static void      project   (const Domain*, AuxField***, AuxField***);
static ModeSys** preSolve  (const Domain*);
static void      Solve     (Field*, AuxField*, ModeSys*,
			    const integer, const integer);


void NavierStokes (Domain*   D,
		   Analyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  DIM = Geometry::nDim();
  CYL = Geometry::system() == Geometry::Cylindrical;
  C3D = CYL && DIM == 3;

  integer       i, j;
  const real    dt     =           Femlib::value ("D_T"   );
  const integer nOrder = (integer) Femlib::value ("N_TIME");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();

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

  ModeSys** MMS = preSolve (D);

  // -- Create & initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [(size_t) DIM];
  AuxField*** Uf = new AuxField** [(size_t) DIM];

  for (i = 0; i < DIM; i++) {
    Us[i] = new AuxField* [(size_t) nOrder];
    Uf[i] = new AuxField* [(size_t) nOrder];
    for (j = 0; j < nOrder; j++) {
      *(Us[i][j] = new AuxField (D -> Esys)) = 0.0;
      *(Uf[i][j] = new AuxField (D -> Esys)) = 0.0;
    }
  }

  // -- Create multi-level storage for pressure BCS.

  Field* Pressure = D -> u[DIM];
  PBCmgr::build (Pressure);

  // -- Create spatially-constant forcing terms.

  vector<real> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], +1);

  // -- Create and initialize eddy-viscosity storage.

  AuxField* EV = new AuxField (D -> Esys, 'e');
  *EV = 0.0;
  ROOTONLY EV -> addToPlane (0, Femlib::value ("REFVIS - KINVIS"));

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Compute spatially-varying kinematic eddy viscosity \epsilon.

    eddyViscosity (D, Us, Uf, EV);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us, Uf, EV, ff);
    waveProp  (D, (const AuxField***) Us, (const AuxField***) Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
		      (const AuxField***) Us, (const AuxField***) Uf, 1);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < DIM; i++) {
      AuxField::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce ((const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Update multilevel velocity storage.

    for (i = 0; i < DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll (Us[i], nOrder);
    }

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [1][0], Uf [2][0], +1);
      AuxField::couple (D -> u[1], D -> u[2], +1);
    }
    for (i = 0; i < DIM; i++)
      Solve (D -> u[i], Uf[i][0], MMS[i], D -> step, nOrder);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], -1);

    // -- Process results of this step.

    A -> analyse (Us);
  }

  // -- Dump ratio eddy/molecular viscosity to file visco.fld.

  ofstream          evfl;
  vector<AuxField*> visco (1);

  visco[0] = EV;

  ROOTONLY {
    evfl.open ("visco.fld", ios::out);
    EV -> addToPlane (0, Femlib::value ("KINVIS"));
  }

  (*EV /= Femlib::value ("REFVIS")) . transform (-1);

  writeField (evfl, D -> name, D -> step, D -> time, visco);

  ROOTONLY evfl.close();
}


static void nonLinear (Domain*       D ,
		       AuxField***   Us,
		       AuxField***   Uf,
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
// while for cylindrical coordinates
//
//           Nx = -0.5 {ud(u)/dx + vd(u)/dy +  d(uu)/dx +
//                 1/y [wd(u)/dz + d(uw)/dz + d(yvu)/dy      ]}
//           Ny = -0.5 {ud(v)/dx + vd(v)/dy +  d(uv)/dx +
//                 1/y [wd(v)/dz + d(vw)/dz + d(yvv)/dy - 2ww]}
//           Nz = -0.5 {ud(w)/dx + vd(w)/dy +  d(uw)/dx +
//                 1/y [wd(w)/dz + d(ww)/dz + d(yvw)/dy + 2wv]}.
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
  vector<real>      work ((2 * DIM + 1) * nTot32);
  vector<real*>     u32 (DIM);
  vector<real*>     n32 (DIM);
  vector<AuxField*> U   (DIM);
  vector<AuxField*> N   (DIM);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * DIM * nTot32;

  for (i = 0; i < DIM; i++) {
    u32[i] = work() +  i        * nTot32;
    n32[i] = work() + (i + DIM) * nTot32;
  }

  // -- Start with contribution from divergence of SGS.

  EV -> transform32 (tmp, -1);
  Blas::scal (nTot32, -4.0, tmp, 1); // -- 2 = -0.5 * -4 ... see below.

  if (CYL) {			// -- Cylindrical coordinates.

    // -- 2D stress-divergence terms.

    Us[0][0] -> transform32 (n32[0], -1);
    Veclib::vmul            (nTot32,  tmp, 1, n32[0], 1, n32[0], 1);
    master   -> gradient    (nZ32, nP, n32[0], 0);

    Us[1][0] -> transform32 (n32[1], -1);
    Veclib::vmul            (nTot32,  tmp, 1, n32[1], 1, n32[1], 1);
    master   -> mulR        (nZ32, n32[1]);
    master   -> gradient    (nZ32, nP, n32[1], 1);
    master   -> divR        (nZ32, n32[1]);

    Uf[0][0] -> transform32 (u32[0], -1);
    Veclib::vmul            (nTot32,  tmp, 1, u32[0], 1, u32[0], 1);
    Veclib::copy            (nTot32,          u32[0], 1, u32[1], 1);
    master   -> mulR        (nZ32, u32[0]);
    master   -> gradient    (nZ32, nP, u32[0], 1);
    master   -> divR        (nZ32, u32[0]);
    master   -> gradient    (nZ32, nP, u32[1], 0);
    Veclib::vadd            (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
    Veclib::vadd            (nTot32, u32[1], 1, n32[1], 1, n32[1], 1);

    // -- 3D stress-divergence terms.

    if (C3D) {		

      Us[2][0] -> transform32 (n32[2], -1);
      Veclib::vmul            (nTot32, tmp, 1, n32[2], 1, n32[2], 1);
      Veclib::copy            (nTot32,         n32[2], 1, u32[0], 1);
      master -> divR          (nZ32, u32[0]);
      Veclib::vsub            (nTot32, n32[1], 1, u32[0], 1, n32[1], 1);
      Femlib::transpose       (n32[2], nZ32,        nP, +1);
      Femlib::DFTr            (n32[2], nZ32 * nPR, nPP, +1);
      Veclib::zero            (nTot32 - nTot, n32[2] + nTot, 1);
      master -> gradient      (nZ, nPP, n32[2], 2);
      Femlib::DFTr            (n32[2], nZ32 * nPR, nPP, -1);
      Femlib::transpose       (n32[2], nZ32,        nP, -1);
      master -> divR          (nZ32, n32[2]);

      Uf[1][0] -> transform32 (u32[0], -1);
      Veclib::vmul            (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy            (nTot32,         u32[0], 1, u32[1], 1);
      Femlib::transpose       (u32[0], nZ32,        nP, +1);
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, +1);
      Veclib::zero            (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient      (nZ, nPP, u32[0], 2);
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, -1);
      Femlib::transpose       (u32[0], nZ32,        nP, -1);
      master -> divR          (nZ32, u32[0]);
      master -> gradient      (nZ32, nP, u32[1], 0);
      Veclib::vadd            (nTot32, u32[0], 1, n32[0], 1, n32[0], 1);
      Veclib::vadd            (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);

      Uf[2][0] -> transform32 (u32[0], -1);
      Veclib::vmul            (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy            (nTot32,         u32[0], 1, u32[1], 1);
      Veclib::copy            (nTot32,         u32[0], 1, u32[2], 1);
      Femlib::transpose       (u32[0], nZ32,        nP, +1);      
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, +1);
      Veclib::zero            (nTot32 - nTot, u32[0] + nTot, 1);
      master -> gradient      (nZ, nPP, u32[0], 2);
      Femlib::DFTr            (u32[0], nZ32 * nPR, nPP, -1);
      Femlib::transpose       (u32[0], nZ32,        nP, -1);      
      master -> divR          (nZ32, u32[0]);
      master -> gradient      (nZ32, nP, u32[1], 1);
      master -> divR          (nZ32, u32[2]);
      Veclib::vadd            (nTot32, u32[0], 1, n32[1], 1, n32[1], 1);
      Veclib::vadd            (nTot32, u32[1], 1, n32[2], 1, n32[2], 1);
      Blas::axpy              (nTot32, 2.0, u32[2], 1,       n32[2], 1);
    }

  } else {			// -- Cartesian coordinates.

    // -- Diagonal stress-divergence terms.

    for (i = 0; i < DIM; i++) {

      Us[i][0] -> transform32 (n32[i], -1);
      Veclib::vmul (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

      if (i == 2) {
	Femlib::transpose  (n32[2], nZ32,        nP, +1);      
	Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, +1);
	Veclib::zero       (nTot32 - nTot, n32[2] + nTot, 1);
	master -> gradient (nZ, nPP, n32[2], 2);
	Femlib::DFTr       (n32[2], nZ32 * nPR, nPP, -1);
	Femlib::transpose  (n32[2], nZ32,        nP, -1);      
      } else
	master -> gradient (nZ32, nP, n32[i], i);
    }

    // -- Off-diagonal stress-divergence terms.

    for (i = 0; i < DIM; i++)
      for (j = i + 1; j < DIM; j++) {

	Uf[i + j - 1][0] -> transform32 (u32[0], -1);
	Veclib::vmul                    (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
	Veclib::copy                    (nTot32,         u32[0], 1, u32[1], 1);

	// -- Super-diagonal.

	if (j == 2) {
	  Femlib::transpose  (u32[0], nZ32,        nP, +1);
	  Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, u32[0] + nTot, 1);
	  master -> gradient (nZ, nPP, u32[0], 2);
	  Femlib::DFTr       (u32[0], nZ32 * nPR, nPP, -1);
	  Femlib::transpose  (u32[0], nZ32,        nP, -1);
	} else
	  master -> gradient (nZ32, nP, u32[0], j);
	Veclib::vadd (nTot32, u32[0], 1, n32[i], 1, n32[i], 1);

	// -- Sub-diagonal.

	master -> gradient (nZ32, nP, u32[1], i);
	Veclib::vadd       (nTot32, u32[1], 1, n32[j], 1, n32[j], 1);
      }
  }

  // -- Nonlinear terms.

  for (i = 0; i < DIM; i++) {
    AuxField::swapData (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    N[i] = Uf[i][0];
    U[i] -> transform32 (u32[i], -1);
    if (CYL) N[i] -> mulR (nZ32, n32[i]);
  }
  
  if (CYL) {			// -- Cylindrical coordinates.

    for (i = 0; i < DIM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (DIM == 3) {
	if      (i == 1)	// -- Centripetal.
	  Veclib::svvtt (nTot32, -2.0, u32[2], 1, u32[2], 1, n32[1], 1);
	else if (i == 2)	// -- Coriolis.
	  Veclib::svvtt (nTot32,  2.0, u32[2], 1, u32[1], 1, n32[2], 1);

	if (nZ > 2) {
	  Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	  Femlib::transpose  (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::transpose  (tmp, nZ32,        nP, -1);
	  Veclib::vvtvp      (nTot32, u32[2], 1, tmp, 1, n32[i], 1, n32[i], 1);
	
	  Veclib::vmul       (nTot32, u32[i], 1, u32[2], 1, tmp, 1);
	  Femlib::transpose  (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ, nPP, tmp, 2);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::transpose  (tmp, nZ32,        nP, -1);
	  Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
	}
      }

      // -- Terms made with product of radius.
      
      Veclib::vmul       (nTot32, u32[1], 1, u32[i], 1, tmp,  1);
      master -> mulR     (nZ32, tmp);
      master -> gradient (nZ32, nP, tmp, 1);
      Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      master -> divR     (nZ32, n32[i]);

      // -- 2D non-conservative derivatives.

      for (j = 0; j < 2; j++) {
	Veclib::copy       (nTot32, u32[i], 1, tmp, 1);
	master -> gradient (nZ32, nP, tmp, j);
	Veclib::vvtvp      (nTot32, u32[j], 1, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Remaining conservative derivative.
      
      Veclib::vmul       (nTot32, u32[0], 1, u32[i], 1, tmp, 1);
      master -> gradient (nZ32, nP, tmp, 0);
      Veclib::vadd       (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

      // -- Transform to Fourier space, smooth, add forcing.

      N[i]   -> transform32 (n32[i], +1);
      master -> smooth (N[i]);
      if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
      *N[i] *= -0.5;
    }
  
  } else {			// -- Cartesian coordinates.

    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
      
	// -- Nonconservative contribution,  n_i += u_j d(u_i) / dx_j.

	Veclib::copy (nTot32, u32[i], 1, tmp,  1);
	if (j == 2) {
	  Femlib::transpose  (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::transpose  (tmp, nZ32,        nP, -1);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

	// -- Conservative contribution,  n_i += d(u_i u_j) / dx_j.

	Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
	if (j == 2) {
	  Femlib::transpose  (tmp, nZ32,        nP, +1);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, +1);
	  Veclib::zero       (nTot32 - nTot, tmp + nTot, 1);
	  master -> gradient (nZ,  nPP, tmp, j);
	  Femlib::DFTr       (tmp, nZ32 * nPR, nPP, -1);
	  Femlib::transpose  (tmp, nZ32,        nP, -1);
	} else {
	  master -> gradient (nZ32, nP, tmp, j);
	}
	Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      }

      // -- Transform to Fourier space, smooth, add forcing.
      
      N[i]   -> transform32 (n32[i], +1);
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
  vector<AuxField*> H (DIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const integer Je = min (D -> step, (integer) Femlib::value ("N_TIME"));
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

  for (i = 0; i < DIM; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[i][q]);
      H[i] -> axpy ( beta [q]    , *Uf[i][q]);
    }
}


static void setPForce (const AuxField*** Us,
		       AuxField***       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {
   *Uf[i][0] = *Us[i][0];
    Uf[i][0] -> gradient (i);
  }

  if (C3D)
    Uf[2][0] -> divR();

  for (i = 1; i < DIM; i++) *Uf[0][0] += *Uf[i][0];

  if (CYL) {
    *Uf[1][0]  = *Us[1][0];
    Uf[1][0] ->  divR();
    *Uf[0][0] += *Uf[1][0];
  }

  *Uf[0][0] /= dt;
}


static void project (const Domain* D ,
		     AuxField***   Us,
		     AuxField***   Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in lowest level of Us.  Constrain velocity field:
//                    u^^ = u^ - D_T * grad P;
// u^^ is left in lowest level of Uf.
//
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to create
// forcing for viscous step.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt    = Femlib::value ("D_T");
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM];
    Uf[i][0] -> gradient (i);

    if (C3D) Uf[2][0] -> divR();

    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Field::swapData (Us[i][0], Uf[i][0]);

    *Uf[i][0] *= alpha;
  }
}


static ModeSys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE == 1 selects iterative solvers for velocity components,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  char                 name;
  integer              i;
  const integer        nSys   = D -> Nsys.getSize();
  const integer        nmodes = Geometry::nModeProc();
  const integer        base   = Geometry::baseMode();
  const integer        itLev  = (integer) Femlib::value ("ITERATIVE");
  const integer        nOrder = (integer) Femlib::value ("N_TIME");
  const real           beta   = Femlib::value ("BETA");
  ModeSys**            M      = new ModeSys* [(size_t) (DIM + 1)];
  vector<Element*>&    E      = ((Domain*) D) -> Esys;
  const NumberSystem** N      = new const NumberSystem* [3];

  // -- Velocity systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

    for (i = 0; i < DIM; i++) {
      name = D -> u[i] -> name();
      D -> setNumber (name, N);
      M[i] = new ModalMatrixSystem (lambda2, beta, name, base, nmodes, E, N);
    }
  } else
    for (i = 0; i < DIM; i++) M[i] = 0;

  // -- Pressure system.

  if (itLev < 2) {
    name = D -> u[DIM] -> name();
    D -> setNumber (name, N);
    M[DIM] = new ModalMatrixSystem (0.0, beta, name, base, nmodes, E, N);
  } else
    M[DIM] = 0;

  return M;
}


static void Solve (Field*        U     ,
		   AuxField*     Force ,
		   ModeSys*      M     ,
		   const integer step  ,
		   const integer nOrder)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  const char    routine[] = "Solve";
  const integer iterative = M == 0;
  const char    name      = U -> name();
  const integer velocity  = name == 'u' || name == 'v' || name == 'w';
  const integer pressure  = name == 'p';

  if (!(velocity || pressure))
    message (routine, "input field type not recognized", ERROR);

  if (pressure) {
    if   (iterative) U -> solve (Force, 0.0);
    else             U -> solve (Force, M);
    return;
  }

  if (velocity) {
    if (iterative || step < nOrder) {
      const integer    Je = min (step, nOrder);
      vector<real> alpha (Je + 1);
      Integration::StifflyStable (Je, alpha());
      const real   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
      
      U -> solve (Force, lambda2);

    } else
      U -> solve (Force, M);

    return;
  }
}
