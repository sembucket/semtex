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
// the total effect should be the same as for a DNS with total
// viscosity KINVIS.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <les.h>

#ifdef __DECCXX
  #pragma define_template roll<AuxField*>
#endif

typedef ModalMatrixSystem ModeSys;
static  int               DIM, NZ;

static void nonLinear (Domain*, AuxField***, AuxField***,
		       AuxField*, vector<real>&);
static void waveProp  (Domain*, const AuxField***, const AuxField***);
static void setPForce (const Domain*, const AuxField***, AuxField***);
static void project   (const Domain*, AuxField***, AuxField***);

static ModeSys** preSolve (const Domain*);
static void      Solve    (Field*, AuxField*, ModeSys*, const int, const int);


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
#if defined (DEBUG)
  Femlib::value ("KINVIS", 2.0 * Femlib::value ("KINVIS"));
#endif

  DIM = Geometry::nDim();
  NZ  = Geometry::nZ();

  int        i, j;
  const real dt     =       Femlib::value ("D_T"   );
  const int  nOrder = (int) Femlib::value ("N_TIME");
  const int  nStep  = (int) Femlib::value ("N_STEP");

  // -- Create global matrix systems.

  ModeSys** MMS = preSolve (D);

  // -- Create multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [DIM];
  AuxField*** Uf = new AuxField** [DIM];

  for (i = 0; i < DIM; i++) {
    Us[i] = new AuxField* [nOrder];
    Uf[i] = new AuxField* [nOrder];
    for (j = 0; j < nOrder; j++) {
      Us[i][j] = new AuxField (D -> Esys, NZ);
      Uf[i][j] = new AuxField (D -> Esys, NZ);
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

  // -- Create eddy-viscosity storage.

  AuxField* EV = new AuxField (D -> Esys, NZ);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Find old SGS stress.

    eddyViscosity (D, Us, Uf, EV);

    // -- Unconstrained forcing substep.

    nonLinear (D, Us, Uf, EV, ff);
    waveProp  (D, (const AuxField***) Us, (const AuxField***) Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain(D -> step, Pressure,
		     (const AuxField***)Us, (const AuxField***)Uf, 1);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < DIM; i++) {
      AuxField::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce (D, (const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Viscous correction substep.

    for (i = 0; i < DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll (Us[i], nOrder);
      Solve (D -> u[i], Uf[i][0], MMS[i], D -> step, nOrder);
    }

    // -- Process results of this step.

    A -> analyse (Us);
  }
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
// lowest levels of Us & Uf contain the components of the old strain-rate
// tensor and vV contains the spatially-varying viscosity.  On exit, the
// velocity field storage areas of D are free, the zeroth level of Us contains
// the old velocities and the zeroth level of Uf contains the most recent
// explcit forcing terms.
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
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
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
// ---------------------------------------------------------------------------
{
  int               i, j;
  const int         nP     = Geometry::planeSize();
  const int         nTot   = NZ * nP;
  const int         nZ32   = (3 * NZ) >> 1;
  const int         nTot32 = nZ32 * nP;
  const real        EPS    = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;
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

  EV -> transform32 (tmp, -1);
  Blas::scal (nTot32, -4.0, tmp, 1); // -- 2 = -0.5 * -4 ... see below.

  // -- Diagonal stress-divergence terms.

  for (i = 0; i < DIM; i++) {

    Us[i][0] -> transform32 (n32[i], -1);
    Veclib::vmul (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

    if (i == 2) {
      Femlib::DFTr (n32[2], nZ32, nP, +1);
      Veclib::zero (nTot32 - nTot, n32[2] + nTot, 1);
      master -> gradient (NZ, n32[2], 2);
      Femlib::DFTr (n32[2], nZ32, nP, -1);
    } else
      master -> gradient (nZ32, n32[i], i);
  }

  // -- Off-diagonal stress-divergence terms.

  for (i = 0; i < DIM; i++) {
    for (j = i + 1; j < DIM; j++) {

      Uf[i + j - 1][0] -> transform32 (u32[0], -1);
      Veclib::vmul (nTot32, tmp, 1, u32[0], 1, u32[0], 1);
      Veclib::copy (nTot32,         u32[0], 1, u32[1], 1);

      // -- Super-diagonal.

      if (j == 2) {
	Femlib::DFTr (u32[0], nZ32, nP, +1);
	Veclib::zero (nTot32 - nTot, u32[0] + nTot, 1);
	master -> gradient (NZ, u32[0], 2);
	Femlib::DFTr (u32[2], nZ32, nP, -1);
      } else
	master -> gradient (nZ32, u32[0], j);
      Veclib::vadd (nTot32, u32[0], 1, n32[i], 1, n32[i], 1);

      // -- Sub-diagonal.

      master -> gradient (nZ32, u32[1], i);
      Veclib::vadd (nTot32, u32[1], 1, n32[i], 1, n32[i], 1);
    }
  }

  // -- Nonlinear terms.

  for (i = 0; i < DIM; i++) {
    AuxField::swapData  (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    N[i] = Uf[i][0];
    U[i] -> transform32 (u32[i], -1);
  }
  
  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      
      // -- Perform n_i = u_j d(u_i) / dx_j.

      Veclib::copy (nTot32, u32[i], 1, tmp,  1);
      if (j == 2) {
	Femlib::DFTr (tmp, nZ32, nP, +1);
	Veclib::zero (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (NZ, tmp, j);
	Femlib::DFTr (tmp, nZ32, nP, -1);
      } else {
	master -> gradient (nZ32, tmp, j);
      }
      Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

      // -- Perform n_i += d(u_i u_j) / dx_j.

      Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
      if (j == 2) {
	Femlib::DFTr (tmp, nZ32, nP, +1);
	Veclib::zero (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (NZ, tmp, j);
	Femlib::DFTr (tmp, nZ32, nP, -1);
      } else {
	master -> gradient (nZ32, tmp, j);
      }
      Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

    }
    N[i]   -> transform32 (n32[i], +1);
    master -> smooth (N[i]);

    if (fabs (ff[i]) > EPS) N[i] -> addToPlane (0, -2.0*ff[i]);
    *N[i] *= -0.5;
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
  int               i, q;
  vector<AuxField*> H (DIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  int  Je = (int) Femlib::value ("N_TIME");
  Je = min (D -> step, Je);

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


static void setPForce (const Domain*     D ,
		       const AuxField*** Us,
		       AuxField***       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
//
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf as a result of gradient operation.
// ---------------------------------------------------------------------------
{
  int         i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {
   *Uf[i][0] = *Us[i][0];
    Uf[i][0] -> gradient (i);
  }
  
  for (i = 1; i < DIM; i++) *Uf[0][0] += *Uf[i][0];

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
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to compute
// forcing for viscous step.
// ---------------------------------------------------------------------------
{
  int        i;
  const real dt    = Femlib::value ("D_T");
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM];
    Uf[i][0] -> gradient (i);
  
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Field::swapData (Us[i][0], Uf[i][0]);

    *Uf[i][0] *= alpha;
  }
}


static ModeSys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Create ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE == 1 selects iterative solvers for velocity components,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  char                 name;
  int                  i;
  const int            nSys   = D -> Nsys.getSize();
  const int            nModes = (NZ + 1) >> 1;
  const int            base   = 0;
  const int            itLev  = (int) Femlib::value ("ITERATIVE");
  const int            nOrder = (int) Femlib::value ("N_TIME");
  const real           beta   = Femlib::value ("BETA");
  ModeSys**            M      = new ModeSys* [DIM + 1];
  vector<Element*>&    E      = ((Domain*) D) -> Esys;
  const NumberSystem** N      = new const NumberSystem* [3];

  // -- Velocity systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real   lambda2 = alpha[0] / Femlib::value ("KINVIS*D_T");

    for (i = 0; i < DIM; i++) {
      name = D -> u[i] -> name();
      D -> setNumber (name, N);
      M[i] = new ModalMatrixSystem (lambda2, beta, name, base, nModes, E, N);
    }
  } else
    for (i = 0; i < DIM; i++) M[i] = 0;

  // -- Pressure system.

  if (itLev < 2) {
    name = D -> u[DIM] -> name();
    D -> setNumber (name, N);
    M[DIM] = new ModalMatrixSystem (0.0, beta, name, base, nModes, E, N);
  } else
    M[DIM] = 0;

  return M;
}


static void Solve (Field*    U     ,
		   AuxField* Force ,
		   ModeSys*  M     ,
		   const int step  ,
		   const int nOrder)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  char routine[] = "Solve";

  const int  iterative = M == 0;
  const char name      = U -> name();
  const int  velocity  = name == 'u' || name == 'v' || name == 'w';
  const int  pressure  = name == 'p';

  if (!(velocity || pressure))
    message (routine, "input field type not recognized", ERROR);

  if (pressure) {
    if   (iterative) U -> solve (Force, 0.0);
    else             U -> solve (Force, M);
    return;
  }

  if (velocity) {
    if (iterative || step < nOrder) {
      const int    Je = min (step, nOrder);
      vector<real> alpha (Je + 1);
      Integration::StifflyStable (Je, alpha());
      const real   lambda2 = alpha[0] / Femlib::value ("D_T*KINVIS");
      
      U -> solve (Force, lambda2);

    } else
      U -> solve (Force, M);

    return;
  }
}
