//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element DNS for incompressible flows.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// USAGE:
// -----
// dns [options] session
//   options:
//   -h       ... print usage prompt
//   -i       ... use iterative solver for viscous [and pressure] steps
//   -t[t]    ... select time-varying BCs for mode 0 [or all modes]
//   -v[v...] ... increase verbosity level
//   -chk     ... turn off checkpoint field dumps [default: selected]
//   -S|C|N   ... regular skew-symm || convective || Stokes advection
//
// AUTHOR:
// ------
// Hugh Blackburn
// Department of Mechanical & Aerospace Engineering
// Monash University
// Vic 3800
// Australia
// hugh.blackburn@monash.edu
//
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
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <dns.h>
#include "rpo_integrate.h"

static char prog[] = "dns";

static int_t NDIM, NCOM, NORD, NADV;
static bool  C3D;

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
  const int_t             base   = Geometry::baseMode();
  const int_t             itLev  = Femlib::ivalue ("ITERATIVE");
  const real_t            beta   = Femlib:: value ("BETA");
  const vector<Element*>& E = D -> elmt;
  Msys**                  M = new Msys* [static_cast<size_t>(NADV + 1)];
  int_t                   i;

  vector<real_t> alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, &alpha[0]);
  real_t   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.
  for (i = 0; i < NCOM; i++)
    M[i] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[i], (itLev) ? JACPCG : DIRECT);

  // -- Scalar system.
  if(NADV != NCOM) {
    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS / PRANDTL");
    M[NCOM] = new Msys
      (lambda2, beta, base, nmodes, E, D -> b[NCOM],    (itLev<1)?DIRECT:JACPCG);
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
    const real_t   lambda2 = (i==NCOM) ? alpha[0] / Femlib::value("D_T * KINVIS / PRANDTL") :
                                         alpha[0] / Femlib::value("D_T * KINVIS");
    const real_t   beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
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
  vector<AuxField*> H (NADV);   // -- Mnemonic for u^{Hat}.

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
// then scale by -1.0 / (D_T * KINVIS) to create forcing for viscous step.
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
  if(NADV>NCOM) { // -- Rescale the temperature equation by the Prandtl number.
    *Uf[NCOM] *= Pr;
  }

  for (i = 0; i < NDIM; i++) {
    (*Us[0] = *D -> u[NADV]) . gradient (i);
    if (Geometry::cylindrical() && i <  2) Us[0] -> mulY();
    Uf[i] -> axpy (beta, *Us[0]);
  }
}

void integrate (void (*advection) (Domain*, 
                                   BCmgr*,
                                   AuxField**, 
                                   AuxField**,
                                   FieldForce*),
                Domain*      D ,
                BCmgr*       B ,
                FieldForce*  FF, int nStep)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  bool do_scat = strchr(D -> field, 'c'); // -- Do scalar transport.

  NDIM = Geometry::nDim();	// -- Number of space dimensions.
  NCOM = (do_scat) ? D -> nField() - 2 : D -> nField() - 1;	// -- Number of velocity components.
  NORD = Femlib::ivalue ("N_TIME");
  C3D  = Geometry::cylindrical() && NDIM == 3;
  NADV = D -> nField() - 1; // -- Number of advected fields.

  int_t              i, j, k;
  const real_t       dt    = Femlib:: value ("D_T");
  const int_t        TBCS  = Femlib::ivalue ("TBCS");
  const int_t        nZ    = Geometry::nZProc();

  static Msys**      MMS;
  static AuxField*** Us;
  static AuxField*** Uf;
  Field*             Pressure = D -> u[NADV];

  if (!MMS) {			// -- Initialise static storage.
    // -- Create multi-level storage for velocities and forcing.
    const int_t ntot  = Geometry::nTotProc();
    real_t*     alloc = new real_t [static_cast<size_t>(2*NADV*NORD*ntot)];
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
  for (i = 0; i < NORD; i++) for (j = 0; j < NADV; j++) {*Us[i][j] = 0.0; *Uf[i][j] = 0.0;}

  D -> step = 0;
  while (D -> step < nStep) {
    // -- Compute nonlinear terms from previous velocity field.
    //    Add physical space forcing, again at old time level.
    advection (D, B, Us[0], Uf[0], FF);
    //if(do_scat) ROOTONLY tempGrad (Us[0], Uf[0]);

    // -- Now update the time.
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

    // -- Correct velocities for pressure.
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.
    for (i = 0; i < NADV; i++) *Us[0][i] = *D -> u[i];
    rollm (Us, NORD, NADV);

    // -- Re-evaluate velocity (possibly time-dependent) BCs.
    for (i = 0; i < NADV; i++)  {
      D -> u[i] -> evaluateBoundaries (Pressure, D -> step, false);
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
  }
}

