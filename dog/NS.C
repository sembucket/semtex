///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//
// This version implements linearised advection terms and evolves a
// single Fourier mode.
//
// For cylindrical coordinates:
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
//
// "$Id$";
///////////////////////////////////////////////////////////////////////////////

#include "stab.h"

static integer NORD, CYL, C3D;

integer NVEC = 0;  // -- Number of components in perturbation velocity.

static void        Linearised (Domain*, AuxField**, AuxField**);
static void        waveProp   (Domain*, const AuxField***, const AuxField***);
static void        setPForce  (const AuxField**, AuxField**);
static void        project    (const Domain*, AuxField**, AuxField**);
static MatrixSys** preSolve   (const Domain*);
static void        Solve      (Domain*, const integer, AuxField*, MatrixSys*);


void NavierStokes (Domain*       D,
		   STABAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Now also contains base velocity fields 'U' and 'V' .
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NVEC = D -> nField() - 1;
  NORD = (integer) Femlib::value ("N_TIME");
  CYL  = Geometry::system() == Geometry::Cylindrical;
  C3D  = CYL && Geometry::nDim() == 3;

  integer       i, j, k;
  const real    dt    =           Femlib::value ("D_T");
  const integer nStep = (integer) Femlib::value ("N_STEP");

  static MatrixSys** MS;
  static AuxField*** Us;
  static AuxField*** Uf;
  static Field*      Pressure;

  if (D -> step == 0) {			// -- Initialise static data.
    
    // -- Create global matrix systems

    MS = preSolve (D);
    
    // -- Create, initialise multi-level storage for velocities and forcing.

    const integer ntot = Geometry::nTotProc();
    real* alloc        = new real [(size_t) 2 * NVEC * NORD * ntot];    

    Us = new AuxField** [(size_t) NORD];
    Uf = new AuxField** [(size_t) NORD];
    
    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [(size_t) NVEC];
      Uf[i] = new AuxField* [(size_t) NVEC];

      for (j = 0; j < NVEC; j++) {
	*(Us[i][j] = new AuxField (alloc + k++ * ntot, 1, D -> elmt)) = 0.0;
	*(Uf[i][j] = new AuxField (alloc + k++ * ntot, 1, D -> elmt)) = 0.0;
      }
    }

    // -- Create multi-level storage for pressure BCS.

    PBCmgr::build (Pressure = D -> u[NVEC]);

    // -- Apply coupling to radial & azimuthal velocity BCs.

    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);
  }

  // -- Timestepping loop.

  do {
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Update domain velocity fields if periodic base fields exist.

    if (D -> n_basefiles > 1) D -> Base_update();
    
    // -- Unconstrained forcing substep.

    Linearised (D, Us[0], Uf[0]);
    waveProp   (D, (const AuxField***)Us, (const AuxField***)Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain
      (D -> step, Pressure, (const AuxField**)Us[0], (const AuxField**)Uf[0]);

    Pressure -> evaluateBoundaries (D -> step);

    for (i = 0; i < NVEC; i++) AuxField::swapData (D -> u[i], Us[0][i]);
    rollm     (Uf, NORD, NVEC);
    setPForce ((const AuxField**)Us[0], Uf[0]);
    Solve     (D, NVEC, Uf[0][0], MS[NVEC]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NVEC; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NVEC);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NVEC; i++) Solve (D, i, Uf[0][i], MS[i]);
    if (C3D) AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0]);

  } while (D -> step % nStep);
}


static void Linearised (Domain*    D ,
			AuxField** Us,
			AuxField** Uf)
// ---------------------------------------------------------------------------
// Compute linearised (forcing) terms in Navier--Stokes equations: N(u) + ff.
//
// Here N(u) represents the linearised advection terms in the N--S equations
// transposed to the RHS and ff is a vector of body force per unit mass.
//
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// linearised terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 
//           N  = -1/2 ( U.grad u + div Uu )
//           ~           ~      ~       ~~
//                -1/2 ( u.grad U + div uU )
//                       ~      ~       ~~
// The data are taken as being in Fourier space, but, as there are
// only two modes involved (the base flow and the perturbation) the
// convolution sums end up being 2D operations.
// ---------------------------------------------------------------------------
{
  integer           i, j;
  vector<AuxField*> U (2), u (NVEC), N (NVEC);
  Field*            T = D -> u[0];

  // -- Set up local aliases.

  U[0] = D -> U[0]; U[1] = D -> U[1];

  for (i = 0; i < NVEC; i++) {
    AuxField::swapData (D -> u[i], Us[i]);
     u[i] = Us[i];
     N[i] = Uf[i];
    *N[i] = 0.0;
  }

  if (CYL) {			// -- Cylindrical coordinates.

    message ("Linearised", "cylindrical coordinates not implemented", ERROR);

  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NVEC; i++) {
      for (j = 0; j < 2; j++) {
      
	// -- N_i += U_j d(u_i) / dx_j.

	N[i] -> timesPlus(*U[j],(*T=*u[i]).gradient(j));

	// -- N_i += d(U_j u_i) / dx_j.

	*N[i] += T->times(*U[j],*u[i]).gradient(j);

	if (i < 2) {		// -- Since U[2] = 0.

	  // -- N_i += u_j d(U_i) / dx_j.

	  N[i]->timesPlus(*u[j],(*T=*U[i]).gradient(j));

	  // -- N_i += d(u_j U_i) / dx_j.

	  *N[i] += T->times(*u[j],*U[i]).gradient(j);
	
	  if (NVEC == 3)	// -- N_i += d(u_2 U_i) / dx_2.

	    *N[i] -= T->times(*u[2],*U[i]).gradient(2); // -- u[2] is Imag.
	}
      }

      T -> smooth (N[i]);
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
  vector<AuxField*> H (NVEC);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NVEC; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const integer Je = min (D -> step, NORD);
  vector<real>  alpha (Integration::OrderMax + 1);
  vector<real>  beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

  for (i = 0; i < NVEC; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[q][i]);
      H[i] -> axpy ( beta [q]    , *Uf[q][i]);
    }
}


static void setPForce (const AuxField** Us,
		       AuxField**       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first NVECension, first level storage
// of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < NVEC; i++) (*Uf[i] = *Us[i]) . gradient(i);

  if (NVEC == 3) *Uf[2] *= -1.0; // -- Since u[2] is Imaginary.

  if (C3D) Uf[2] -> divR();

  for (i = 1; i < NVEC; i++) *Uf[0] += *Uf[i];  

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

  for (i = 0; i < NVEC; i++) {

    (*Uf[i] = *D -> u[NVEC]) . gradient (i);

    if (C3D) Uf[2] -> divR();

    Us[i] -> axpy (-dt, *Uf[i]);
    Field::swapData (Us[i], Uf[i]);
   
    *Uf[i] *= alpha;
  }
}


static MatrixSys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for system with only 1 Fourier mode.
// ---------------------------------------------------------------------------
{
  const real              betak2 = (NVEC<3) ? 0.0 : Femlib::value("BETA*BETA");
  const vector<Element*>& E      = D -> elmt;
  MatrixSys**             M      = new MatrixSys* [(size_t) (NVEC + 1)];
  const integer           itLev  = (integer) Femlib::value ("ITERATIVE");

  vector<real> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, alpha());
  const real   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.
  
  cout << "-- Installing matrices for fields [u, v,";

  M[0] = new MatrixSys 
    (lambda2, betak2, 0, E, D -> b[0],     (itLev < 1) ? DIRECT : JACPCG);

  M[1] = M[0]; if (NVEC == 3) { cout << " w,"; M[2] = M[0]; }

  // -- Pressure system.

  cout << " p]" << endl;

  M[NVEC] = new MatrixSys
    (0.0,     betak2, 0, E, D -> b[NVEC], (itLev < 2) ? DIRECT : JACPCG);

  return M;
}


static void Solve (Domain*       D,
		   const integer i,
		   AuxField*     F,
		   MatrixSys*    M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order.
// ---------------------------------------------------------------------------
{
  const integer step = D -> step;

  if (i < NVEC && step < NORD) {

    // -- We need a temporary matrix system for a viscous solve.

    const integer Je = min (step, NORD);    
    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
    const real    betak2  = (NVEC < 3) ? 0.0 : Femlib::value ("BETA * BETA");
    const real    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

    MatrixSys* tmp = new MatrixSys
      (lambda2, betak2, 0, D -> elmt, D -> b[0], JACPCG);

    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
