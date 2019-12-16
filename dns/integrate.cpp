///////////////////////////////////////////////////////////////////////////////
// integrate.cpp: Unsteady Navier--Stokes solver, using
// "stiffly-stable" time integration [1,2].  Geometries may be 2- or
// 3-dimensional, Cartesian or cylindrical [3].  Fourier expansions
// are used in the homogeneous (z) direction.  This file provides
// integrate as a call-back routine; after initialisation, integrate
// may be called repeatedly without reinitialising internal storage.
//
// For cylindrical coordinates (Fourier in azimuth):
//   u <==> axial     velocity,  x <==> axial     coordinate direction,
//   v <==> radial    velocity,  y <==> radial    coordinate direction,
//   w <==> azimuthal velocity,  z <==> azimuthal coordinate direction.
//
// For Cartesian coordinates (Fourier in z):
//   u <==> x-component  velocity
//   v <==> y-component  velocity
//   w <==> z-component  velocity
//
// In either system, the w velocity component is optional for 2D
// (N_Z=1) (i.e. can have 2D2C or 2D3C).  If 3D (N_Z > 1), w should
// appear in session.
//
// Optionally integrate concentration of advected scalar field c.
//
// Copyright (c) 1994 <--> $Date: 2019/06/21 13:22:31 $, Hugh Blackburn
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag (1991) "High-order splitting methods
//     for the incompressible Navier--Stokes equations", JCP 97:414--443
// [2] Guermond & Shen (2003) "Velocity correction projection methods for
//     incompressible flows", SIAM J Numer Anal 41:112-134
// [3] Blackburn & Sherwin (2004) "Formulation of a Galerkin spectral
//     element--Fourier method for three-dimensional incompressible flows
//     in cylindrical geometries", JCP 179:759-778
// [4] Dong, Karniakadis & Chryssostomides (2014) "A robust and
//     accurate outflow boundary condition for incompressible flow
//     simulations on severely-truncated unbounded domains", JCP 261:83-105.
// [5] Blackburn, Lee, Albrecht & Singh (2019) "Semtex: a spectral
//     element–Fourier solver for the incompressible Navier–Stokes
//     equations in cylindrical or Cartesian coordinates", CPC.
// --
// This file is part of Semtex.
//
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: integrate.cpp,v 9.2 2019/06/21 13:22:31 hmb Exp $";

#include <dns.h>
#include <mpi.h>

typedef ModalMatrixSys Msys;

// -- File-scope constants and routines:

static int_t NDIM, NCOM, NORD, NADV;
static bool  C3D;

static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const int_t, AuxField*, Msys*);

//#define ID_DIAGNOSTIC 1

#ifdef ID_DIAGNOSTIC
bool alloc_diagnostics = true;
AuxField** vort;
AuxField* enst;
AuxField* pres;

void diagnostics(Domain* domain) {
  double prod, diss, int_dudy;
  Vector du;
  ofstream file;

  // allocate if not already done
  if(alloc_diagnostics) {
    alloc_diagnostics = false;
    vort = new AuxField*[static_cast<size_t>(3)];
    for(int ii = 0; ii < 3; ii++) {
      vort[ii] = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    }
    enst = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    pres = new AuxField(new real_t[Geometry::nTotal()], Geometry::nZProc(), domain->elmt);
    return;
  }

  // pressure field was stomped on by the fieldforce, copy back from temporary variable
  pres->transform(INVERSE);
  // transform state into physical space in order to perform pointwise multiplications
  for(int ii = 0; ii < 3; ii++) domain->u[ii]->transform(INVERSE);

  // energy production, compute as: I = \int_{V} u.GRAD p dV
  *enst = 0.0;
  for(int ii = 0; ii < 3; ii++) {
    *vort[ii] = *pres;
    if(ii == 2) vort[ii]->transform(FORWARD);
    vort[ii]->gradient(ii);
    if(ii == 2) vort[ii]->transform(INVERSE);
    if(ii == 2) vort[ii]->divY();

    // constant pressure gradient forcing
    //if(ii == 0) *vort[ii] += 4.0*Femlib::value("KINVIS");
    // constant mass flux forcing
    if(ii == 0) {
      *domain->u[3] = *domain->u[0];
      domain->u[3]->transform(FORWARD);
      domain->u[3]->gradient(1);
      du = Field::normTraction(domain->u[3]);
      int_dudy = -2.0 * du.y * Femlib::value("KINVIS") / Femlib::value("XMAX");
      MPI_Bcast(&int_dudy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      *vort[ii] += int_dudy;
    }

    *vort[ii] *= *domain->u[ii];
    *enst += *vort[ii];
  }
  enst->transform(FORWARD);
  prod = enst->integral() / Femlib::value("KINVIS");

  // energy dissipation, compute as: D = \int_{V} [v^2 + w^2 - 2 w dv/\theta + 2 v dwd/theta]/y^2 + |GRAD u|^2 + |GRAD v|^2 + |GRAD w|^2 dV
  // first do the |GRAD U|^2, U = {u,v,w} terms
  *enst = 0.0;

  for(int ii = 0; ii < 3; ii++) {
    for(int jj = 0; jj < 3; jj++) {
      if(jj == ii) continue;

      *vort[0] = *domain->u[ii];
      if(jj == 2) vort[0]->transform(FORWARD);
      vort[0]->gradient(jj);
      if(jj == 2) vort[0]->transform(INVERSE);
      if(jj == 2) vort[0]->divY();

      *vort[0] *= *vort[0];
      *enst += *vort[0];
    }

    int jj = (ii+1)%3;
    int kk = (ii+2)%3;

    *vort[0] = *domain->u[jj];
    *vort[1] = *domain->u[kk];

    if(kk == 2) vort[0]->transform(FORWARD);
    vort[0]->gradient(kk);
    if(kk == 2) vort[0]->transform(INVERSE);
    if(kk == 2) vort[0]->divY();

    if(jj == 2) vort[1]->transform(FORWARD);
    vort[0]->gradient(jj);
    if(jj == 2) vort[1]->transform(INVERSE);
    if(jj == 2) vort[1]->divY();

    *vort[0] *= *vort[1];
    *vort[0] *= 2.0;
    *enst -= *vort[0];
  }

  *vort[0] = *domain->u[2];
  vort[0]->gradient(1);
  vort[0]->divY();
  *vort[0] *= *domain->u[2];
  *vort[0] *= 2.0;
  *enst += *vort[0];

  *vort[0] = *domain->u[1];
  vort[0]->transform(FORWARD);
  vort[0]->gradient(2);
  vort[0]->transform(INVERSE);
  vort[0]->divY();
  *vort[0] *= *domain->u[2];
  vort[0]->divY();
  *vort[0] *= 2.0;
  *enst -= *vort[0];

  // integeate in fourier space
  enst->transform(FORWARD);
  diss = enst->integral();

  // transform state back into fourier space
  for(int ii = 0; ii < 3; ii++) domain->u[ii]->transform(FORWARD);
  pres->transform(FORWARD);
  *domain->u[3] = *pres;

  if(!Geometry::procID()) {
    file.open("production_dissipation.txt", ios::app);
    file.precision(12);
    file << domain->step << "\t" << prod << "\t" << diss << "\n";
    file.close();
  }
}
#endif

void integrate (void (*advection) (Domain*    , 
                                   BCmgr*     ,
                                   AuxField** , 
                                   AuxField** ,
                                   FieldForce*),
                Domain*      D ,
                BCmgr*       B ,
                DNSAnalyser* A ,
                FieldForce*  FF)
// ---------------------------------------------------------------------------
// On entry, D contains storage (in the following order!) for:
// -- velocity Fields 'u', 'v' (and 'w' if 2D3C or 3D),
// -- optional scalar Field 'c',
// -- constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NCOM = D -> nVelCmpt();              // -- Number of velocity components.
  NADV = D -> nAdvect();               // -- Number of advected fields.
  NDIM = Geometry::nDim();	       // -- Number of space dimensions.
  NORD = Femlib::ivalue ("N_TIME");    // -- Time integration order.
  C3D  = Geometry::cylindrical() && NDIM == 3;
  
  int_t              i, j, k;
  const real_t       dt    = Femlib:: value ("D_T");
  const int_t        nStep = Femlib::ivalue ("N_STEP");
  const int_t        nZ    = Geometry::nZProc();
  static Msys**      MMS;
  static AuxField*** Us;
  static AuxField*** Uf;
  Field*             Pressure = D -> u[NADV];

  if (!MMS) {			// -- Initialise static storage.

    // -- Create multi-level storage for velocities and forcing.

    const int_t ntot  = Geometry::nTotProc();
    real_t*     alloc = new real_t [static_cast<size_t>(2 * NADV*NORD * ntot)];
    Us                = new AuxField** [static_cast<size_t>(2 * NORD)];
    Uf                = Us + NORD;

    for (k = 0, i = 0; i < NORD; i++) {
      Us[i] = new AuxField* [static_cast<size_t>(2 * NADV)];
      Uf[i] = Us[i] + NADV;
      for (j = 0; j < NADV; j++) {
        Us[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
        Uf[i][j] = new AuxField (alloc + k++ * ntot, nZ, D -> elmt);
      }
    }

    // -- Create global matrix systems.

    MMS = preSolve (D);

    // -- Create multi-level storage for pressure BCS.

    B -> buildComputedBCs (Pressure);

    // -- Apply coupling to radial & azimuthal velocity BCs.

    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);
  }

  // -- Because we may restart from scratch on each call, zero these:

  *Pressure = 0.0;

  for (i = 0; i < NORD; i++)
    for (j = 0; j < NADV; j++) {
      *Us[i][j] = 0.0;
      *Uf[i][j] = 0.0;
    }

  // -- Solve the Stokes flow problem with unit forcing
  if(!D->grn[0] && fabs(Femlib::value("Q_BAR")) > 1.0e-6) {
    real_t            L_x   = Femlib::value("XMAX");
    vector<AuxField*> tmp;
    tmp.resize(4);

    for (i = 0; i < 4; i++) {
      D->grn[i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), D->elmt, 'g'+i);
      tmp[i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), D->elmt, 'k'+i);
      *tmp[i] = *D->u[i];
      *D->u[i] = 0.0;
    }

    // -- Set the constant forcing
    for (i = 0; i < NCOM; i++) *Uf[0][i] = 0.0;
    ROOTONLY {
      Uf[0][0] -> addToPlane (0, 1.0);
      if (Geometry::cylindrical()) Uf[0][0] -> mulY ();
    }

    D->step += 1;

    // -- Update high-order pressure BC storage.
    B -> maintainFourier (D -> step, Pressure, const_cast<const AuxField**>(Us[0]), const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (Pressure, D -> step);

    // -- Complete unconstrained advective substep and compute pressure.
    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp (D, const_cast<const AuxField***>(Us), const_cast<const AuxField***>(Uf));
    for (i = 0; i < NADV; i++) AuxField::swapData (D -> u[i], Us[0][i]);

    rollm     (Uf, NORD, NADV);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NADV,  Uf[0][0], MMS[NADV]);

    // -- Correct velocities for pressure.
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.
    for (i = 0; i < NADV; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NADV);

    // -- Re-evaluate velocity (possibly time-dependent) BCs.
    for (i = 0; i < NADV; i++)  {
      D -> u[i] -> evaluateBoundaries (NULL,     D -> step, false);
      D -> u[i] -> bTransform         (FORWARD);
      D -> u[i] -> evaluateBoundaries (Pressure, D -> step, true);
    }
    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

    // -- Viscous correction substep.
    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NADV; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D) AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    D->step -= 1;


/*
    // read in constant forcing stokes solution from file
    char       buf[StrMax], file[StrMax];
    ifstream   grnfunc (strcat (strcpy (file, D -> name), ".grn"));
    Header     header;
    grnfunc >> header;
    for (i = 0; i < 4; i++) {
      grnfunc >> *D->u[i];
      if (header.swab()) D->u[i] -> reverse();
    }
    for (i = 0; i < 4; i++) D->u[i] -> transform (FORWARD);
*/

    // this will be broadcast to the other procs when applied
    D->Qg  = 2.0 * M_PI * D->u[0]->integral(0);
    D->Qg /= (M_PI * 1.0 * 1.0 * L_x);
    if(!Geometry::procID()) cout << "Stokes + unit forcing volumetric flux: " << D->Qg << endl;
    if(!Geometry::procID()) cout << "                          ux integral: " << 2.0 * M_PI * D->u[0]->integral(0) << endl;
    if(!Geometry::procID()) cout << "                          pipe length: " << L_x << endl;

    // -- Resetting fields
    for (i = 0; i < 4; i++) {
      *D->grn[i] = *D->u[i];
      *D->u[i]   = *tmp[i];
    }

    *Pressure = 0.0;
    for (i = 0; i < NORD; i++)
      for (j = 0; j < NADV; j++) {
        *Us[i][j] = 0.0;
        *Uf[i][j] = 0.0;
      }
    // do we really need to do this again??
    //B -> buildComputedBCs (Pressure);
  }

  // -- The following timestepping loop implements equations (15--18) in [5].

#ifdef ID_DIAGNOSTIC
  // setup only
  diagnostics(D);
#endif
  
  while (D -> step < nStep) {

    // -- Compute nonlinear terms from previous velocity field.
    //    Add physical space forcing, again at old time level.

    advection (D, B, Us[0], Uf[0], FF);
    
    // -- Now update the time (remainder including BCs at new time level).

    D -> step += 1;
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Update high-order pressure BC storage.

    B -> maintainFourier (D -> step, Pressure,
			  const_cast<const AuxField**>(Us[0]),
			  const_cast<const AuxField**>(Uf[0]));
    Pressure -> evaluateBoundaries (Pressure, D -> step);

    // -- Complete unconstrained advective substep and compute pressure.

    if (Geometry::cylindrical()) { Us[0][0] -> mulY(); Us[0][1] -> mulY(); }

    waveProp (D, const_cast<const AuxField***>(Us),
	         const_cast<const AuxField***>(Uf));
    for (i = 0; i < NADV; i++) AuxField::swapData (D -> u[i], Us[0][i]);

    rollm     (Uf, NORD, NADV);
    setPForce (const_cast<const AuxField**>(Us[0]), Uf[0]);
    Solve     (D, NADV,  Uf[0][0], MMS[NADV]);

#ifdef ID_DIAGNOSTIC
    // copy over before this gets stomped on by the fieldforce
    *pres = *D->u[NADV];
#endif

    // -- Correct velocities for pressure.

    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NADV; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NADV);

    // -- Re-evaluate velocity (possibly time-dependent) BCs.

    for (i = 0; i < NADV; i++)  {
      D -> u[i] -> evaluateBoundaries (NULL,     D -> step, false);
      D -> u[i] -> bTransform         (FORWARD);
      D -> u[i] -> evaluateBoundaries (Pressure, D -> step, true);
    }
    if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }

    if(D->step == 1) Femlib::ivalue("N_TIME", 1);
    if(D->step == 2) Femlib::ivalue("N_TIME", 2);
    for (i = 0; i < NADV; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if(D->step == 1) Femlib::ivalue("N_TIME", NORD);
    if(D->step == 2) Femlib::ivalue("N_TIME", NORD);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // constant flow rate
    if(fabs(Femlib::value("Q_BAR")) > 1.0e-6) {
      real_t L_x       = Femlib::value("XMAX");
      real_t _refQ     = Femlib::value("Q_BAR");
      real_t getQ, dP;
/*
      getQ  = 2.0 * M_PI * D->u[0]->integral(0);
      getQ /= (M_PI * 1.0 * 1.0 * L_x);

      ROOTONLY {
        dP = (_refQ - getQ) / D->Qg;

        //cout << "mass flux forcing: " << _refQ << "\t" << D->Qg << "\t" << getQ << "\t" << dP << endl;

        for(i = 0; i < Geometry::nProc(); i++) Femlib::send(&dP, 1, i);
      } else {
        Femlib::recv(&dP, 1, 0);
      }
      Femlib::synchronize();
*/
      if(!Geometry::procID()) {
        getQ  = 2.0 * M_PI * D->u[0]->integral(0);
        getQ /= (M_PI * 1.0 * 1.0 * L_x);
        dP    = (_refQ - getQ) / D->Qg;
      }
      //if(!Geometry::procID()) cout << D->step << ":\tmass flux forcing: " << _refQ << "\t" << D->Qg << "\t" << getQ << "\t" << dP << endl;
      MPI_Bcast(&dP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      for(i = 0; i < NADV; i++) D->u[i]->axpy(dP, *D->grn[i]);
    }

    // -- Process results of this step.
    
    //A -> analyse (Us[0], Uf[0]);
#ifdef ID_DIAGNOSTIC
    diagnostics(D);
#endif

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
  vector<AuxField*> H (NADV);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NADV; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const int_t    Je = min (D -> step, NORD);
  vector<real_t> alpha (Integration::OrderMax + 1);
  vector<real_t> beta  (Integration::OrderMax);

  Integration::StifflyStable (Je, &alpha[0]);
  Integration::Extrapolation (Je, &beta [0]);
  Blas::scal (Je, Femlib::value ("D_T"), &beta[0],  1);

  for (i = 0; i < NADV; i++)
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
// then scale by -1.0 / (D_T * KINVIS) to create forcing for viscous step
// (this is -1.0 / (D_T  * diffusivity) in the case of a scalar field).
//
// u^^ is left in Uf.
// ---------------------------------------------------------------------------
{
  int_t        i;
  const real_t alpha = -1.0 / Femlib::value ("D_T * KINVIS");
  const real_t beta  =  1.0 / Femlib::value ("KINVIS");
  const real_t Pr    =        Femlib::value ("PRANDTL");

  for (i = 0; i < NADV; i++) {
    Field::swapData (Us[i], Uf[i]);
    if (Geometry::cylindrical() && i >= 2) Uf[i] -> mulY();
    *Uf[i] *= alpha;
  }

  // -- For scalar, use diffusivity instead of viscosity.
  if (NADV > NCOM) *Uf[NCOM] *= Pr;

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[NADV]) . gradient (i);
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
// ITERATIVE >= 1 selects iterative solver for velocity components,
// ITERATIVE >= 2 selects iterative solver for non-zero pressure Fourier modes.
// ---------------------------------------------------------------------------
{
  const int_t             nmodes = Geometry::nModeProc();
  const int_t             base   = Geometry::baseMode(); const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const real_t            beta   = Femlib:: value ("BETA");
  const vector<Element*>& E = D -> elmt;
  Msys**                  M = new Msys* [static_cast<size_t>(NADV + 1)];
  int_t                   i;

  vector<real_t> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NCOM; i++) {
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i], (itLev) ? JACPCG : DIRECT);
  }

  // -- Scalar system.

  if (NADV != NCOM) {
    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");
    M[NCOM] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[NCOM],(itLev < 1)?DIRECT:JACPCG);
  }

  // -- Pressure system.

  if (itLev > 1)
    M[NADV] = new Msys
      (0.0, beta, base, nmodes, E, D -> b[NADV], MIXED);
  else
    M[NADV] = new Msys
      (0.0, beta, base, nmodes, E, D -> b[NADV], DIRECT);

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

  if (i < NADV && step < NORD) { // -- We need a temporary matrix system.
    const int_t Je     = min (step, NORD);
    const int_t base   = Geometry::baseMode();
    const int_t nmodes = Geometry::nModeProc();

    vector<real_t> alpha (Je + 1);
    Integration::StifflyStable (Je, &alpha[0]);
    const real_t   lambda2 = (i == NCOM) ? // -- True for scalar diffusion.
      alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL") :
      alpha[0] / Femlib::value ("D_T * KINVIS");
    const real_t   beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}
