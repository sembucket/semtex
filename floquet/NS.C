///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
// Geometries may be 2- or 3-dimensional, Cartesian or cylindrical.
// Fourier expansions are used in the homogeneous direction.
//
// References:
// 1.  Karniadakis, Israeli & Orszag 1991.  "High-order splitting methods
//     for the incompressible Navier--Stokes equations", JCP 9(2).
// 2.  Tomboulides, Orszag & Karniadakis 1993.  "Direct and
//     large-eddy simulation of axisymmetric wake",  AIAA-93-0546.
//
// For cylindrical coordinates:
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
//
// "$Id$";
///////////////////////////////////////////////////////////////////////////////

#include <stab.h>
#include <Sem.h>

typedef MatrixSys Mx_sys;

static integer verbose = 0;

static integer NORD, CYL, C3D;
integer NDIM = 0;

static integer firstTime = 1;   // is this the first use (for memory alloc)

static void     Linearised (Domain*, AuxField**, AuxField**, vector<real>&);
static void     waveProp   (Domain*, const AuxField***, const AuxField***);
static void     setPForce  (const AuxField**, AuxField**);
static void     project    (const Domain*, AuxField**, AuxField**);
static Mx_sys** preSolve   (const Domain*);
static void     Solve      (Domain*, const integer, AuxField*, Mx_sys*);


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
  NDIM = D->nField()-1;  // can't use nDim() cos' based on nZ planes.
  NORD = (integer) Femlib::value ("N_TIME");
  CYL  = Geometry::system() == Geometry::Cylindrical;
  C3D  = CYL && Geometry::nDim()== 3;
  verbose = (integer) Femlib::value ("VERBOSE");

  integer       i, j, k;
  const real    dt     =           Femlib::value ("D_T");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer ntot   = Geometry::nTotProc();
  real*         alloc  = new real [(size_t) 2 * NDIM * NORD * ntot];

  static Mx_sys**     MS;
  static AuxField***  Us;
  static AuxField***  Uf;
  static Field*       Pressure;

  // -- Create global matrix systems.
  if (firstTime) {
    MS = preSolve (D);
    firstTime = 0;
    
    // -- Create, initialize multi-level storage for velocities and forcing.
    
    Us = new AuxField** [(size_t) NORD];
    Uf = new AuxField** [(size_t) NORD];
    
    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [(size_t) NDIM];
      Uf[i] = new AuxField* [(size_t) NDIM];

      for (j = 0; j < NDIM; j++) {
	*(Us[i][j] = new AuxField (alloc + k++ * ntot, 1, D -> elmt)) = 0.0;
	*(Uf[i][j] = new AuxField (alloc + k++ * ntot, 1, D -> elmt)) = 0.0;
      }
    }
    // -- Create multi-level storage for pressure BCS.
    Pressure = D -> u[NDIM];
    PBCmgr::build(Pressure);

  }
  else {
  // initialize Us & Uf
    for (i = 0; i < NORD; i++) {
      for (j = 0; j < NDIM; j++) {
	*Us[i][j] = 0.0;
	*Uf[i][j] = 0.0;
      }
    }
  } 

  // -- Create spatially-constant forcing terms.
  vector<real> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");
  VERBOSE cout << "-- FF(X, Y, Z) = " << ff[0] <<", "
	       << ff[1] << ", "<<  ff[2] <<  endl;

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Timestepping loop.

  while (D -> step < nStep) {    

    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // update domain velocity fields if periodic base fields exist.

    if ( D->n_basefiles > 1 ) D -> Base_update();
    
    // -- Unconstrained forcing substep.

    Linearised (D, Us[0], Uf[0], ff);
    waveProp   (D, (const AuxField***)Us, (const AuxField***)Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
    	      (const AuxField**)Us[0],
    	      (const AuxField**)Uf[0],1);

    Pressure -> evaluateBoundaries (D -> step);

    for (i = 0; i < NDIM; i++) AuxField::swapData (D -> u[i], Us[0][i]);

    rollm (Uf, NORD, NDIM);

    setPForce ((const AuxField**)Us[0], Uf[0]);
    
    Solve (D, NDIM, Uf[0][0], MS[NDIM]);
   
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.
    for (i = 0; i < NDIM; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NDIM);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf[0][1], Uf[0][2], +1);
      AuxField::couple (D -> u[1], D -> u[2], +1);
    }
    for (i = 0; i < NDIM; i++) Solve (D, i, Uf[0][i], MS[i]);
    if (C3D) AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.

    A -> analyse (Us[0]);
  }
}


static void Linearised (Domain*       D ,
			AuxField**    Us,
			AuxField**    Uf,
			vector<real>& ff)
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
//                       ~      ~       ~~
//                -1/2 ( u.grad U + div uU )
//                       ~      ~       ~~
//
// i.e., in Cartesian component form
//// If STOKES is defined for compilation, the nonlinear terms are set to zero.

//           N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//            i           j    i      j      i j      j
//
// in cylindrical coordinates
//
//              This has yet to be evaluated !!!!!
//
//
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
//
// ---------------------------------------------------------------------------
{
  integer i, j;
  const real    beta   = Femlib::value ("BETA");
  const integer     nP     = Geometry::planeSize();
  vector<real>      work ((2 * NDIM + 3) * nP);
  vector<real*>     u (NDIM);
  vector<real*>     U (NDIM);
  vector<real*>     n (NDIM);
  vector<AuxField*> N (NDIM);
  vector<AuxField*> Ut   (NDIM);

  Field*            master = D -> u[0];
  real*             tmp    = work() + (2 * NDIM + 2) * nP;

  Veclib::zero ((2 * NDIM + 3) * nP, work(), 1); // -- A catch-all cleanup.


  for (i = 0; i < NDIM; i++) {
    u[i] = work() +  i             * nP;
    n[i] = work() + (i + NDIM)     * nP;

    AuxField::swapData (D -> u[i], Us[i]);

    if (i < 2) {
      U[i] = work() + (i + 2 * NDIM) * nP;
      D -> U[i] -> getPlane (0, U[i]);
    }
    //D -> u[i] -> getPlane (0, u[i]);

    Ut[i] = Us[i];
    Ut[i] -> transform32 (INVERSE, u[i]);
    N[i] = Uf[i];
  }

  if (CYL) {			// -- Cylindrical coordinates.

    cerr << " Cylindrical coordinates are not yet calculated" << endl;
    
  } else {			// -- Cartesian coordinates.
    for (i = 0; i < NDIM; i++) {
      for (j = 0; j < 2; j++) { // 3rd component later if required
      
	// -- Perform N_i += U_j d(u_i) / dx_j.
	Veclib::copy (nP, u[i], 1, tmp,  1);
	master -> gradient (1, nP, tmp, j);
	Veclib::vvtvp (nP, U[j], 1, tmp,  1, n[i], 1, n[i], 1);

	// -- Perform N_i += d(U_j u_i) / dx_j.
	Veclib::vmul  (nP, U[j], 1, u[i], 1, tmp,  1);
	master -> gradient (1, nP, tmp, j);
	Veclib::vadd (nP, tmp, 1, n[i], 1, n[i], 1);
      
	// only if i not z-dir'n => U_2 (W) = 0
	
	if ( i < 2 ) {

	  // -- Perform N_i += u_j d(U_i) / dx_j
	  Veclib::copy (nP, U[i], 1, tmp, 1);
	  master -> gradient (1, nP, tmp, j);
	  Veclib::vvtvp (nP, u[j], 1, tmp,  1, n[i], 1, n[i], 1);

	  // -- Perform N_i += d(u_j U_i) / dx_j.
	  Veclib::vmul  (nP, u[j], 1, U[i], 1, tmp,  1);
	  master -> gradient (1, nP, tmp, j);
	  Veclib::vadd (nP, tmp, 1, n[i], 1, n[i], 1);
	
	  // -- If 3D Perform N_i += d(u_2 U_i) / dx_2 = - w B U_i

	  if ( NDIM == 3 ) {
	    Veclib::vmul  (nP, u[2], 1, U[i], 1, tmp,  1);
	    Blas::scal(nP, beta, tmp, 1); 
	    Veclib::vadd (nP, tmp, 1, n[i], 1, n[i], 1);
	  }
	}
      }

     
      N[i] -> setPlane (0, n[i]);
      master -> smooth (N[i]);
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
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first NDIMension, first level storage
// of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) {
    *Uf[i] = *Us[i];
    Uf[i]->gradient (i);
  }

  if (NDIM == 3) *Uf[2] *= -1.0;

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

    *Uf[i] = *D -> u[NDIM];

    Uf[i] -> gradient(i);

    if (C3D) Uf[2] -> divR();

    Us[i] -> axpy (-dt, *Uf[i]);
    Field::swapData (Us[i], Uf[i]);
   
    *Uf[i] *= alpha;
  }
}

static Mx_sys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for system with only 1 fourier mode
// ---------------------------------------------------------------------------
{
  const real              betak2 = (D->nField()==3) ? 0.0 : Femlib::value ("BETA*BETA");
  const vector<Element*>& E      = D -> elmt;
  Mx_sys**                M      = new Mx_sys* [(size_t) ( NDIM + 1)];
  const integer           itLev  = (integer) Femlib::value ("ITERATIVE");

  vector<real> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, alpha());
  const real   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.
  
  VERBOSE cout << "-- Installing matrices for field [u *, v .,";

  M[0] = new Mx_sys 
    (lambda2, betak2, 0, E, D -> b[0], ((itLev<1)?DIRECT:JACPCG));

  M[1] = M[0];

  // z matrix system - if applicable
  if ( NDIM == 3 ) {
    VERBOSE cout << " w .,";
    M[2] = M[0];
  }
  // -- Pressure system.

  VERBOSE cout << " p *]" << endl;
  M[NDIM] = new Mx_sys 
    (0.0, betak2, 0, E, D -> b[NDIM], ((itLev<1)?DIRECT:JACPCG));

  return M;
}

static void Solve (Domain*       D,
		   const integer i,
		   AuxField*     F,
		   Mx_sys*       M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  const integer step = D -> step;

  if (i < NDIM && step < NORD) { // -- We need a temporary matrix system.
    const integer Je      = min (step, NORD);    
    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
  const real              betak2 = (D->nField()==3) ? 0.0 : Femlib::value ("BETA*BETA");
    const real    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

    Mx_sys* tmp = new Mx_sys
      (lambda2, betak2, 0, D -> elmt, D -> b[0], JACPCG);
    D -> u[i] -> solve_MS (F, tmp, 0, 0);
    delete tmp;

  } else D -> u[i] -> solve_MS (F, M, 0, 0);

}
