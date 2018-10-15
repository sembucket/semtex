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
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

static char prog[] = "dns";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);

typedef ModalMatrixSys Msys;

static int_t NDIM, NCOM, NORD, NADV;
static bool  C3D;

struct Context {
    int              nSlice;
    int              nDofs;
    int              nDofsPlane;
    Mesh*            mesh;
    vector<Element*> elmt;
    Domain*          domain;
    BCmgr*           bman;
    FieldForce*      ff;
    vector<Field*>   ui;
    vector<Field*>   fi;
    real_t*          theta_i;
    real_t*          phi_i;
    real_t*          tau_i;
    // regular grid points in elements
    int_t*           el;
    real_t*          x;
    real_t*          y;
    real_t*          r;
    real_t*          s;
};

#define DOF 1
#define SLICE_DT 10.0
#define XMIN 0.0
#define XMAX 10.0
#define YMIN 0.0
#define YMAX 20.0
#define NELS_X 10
#define NELS_Y 12

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

void elements_to_logical(real_t* data_els, real_t* data_log) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int pt_x, pt_y;
  int index = -1;

  for(int el_y = 0; el_y < NELS_Y; el_y++) {
    for(int el_x = 0; el_x < NELS_X; el_x++) {
      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        index++;

        // skip over right and top edges for each element, as these are redundant
        if(pt_i%(elOrd+1) == elOrd || pt_i/(elOrd+1) == elOrd) continue;

        pt_x = el_x*elOrd + pt_i%(elOrd+1);
        pt_y = el_y*elOrd + pt_i/(elOrd+1);

        data_log[pt_y*NELS_X*elOrd + pt_x] = data_els[index];
      }
    }
  }
}

void logical_to_elements(real_t* data_log, real_t* data_els) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int shift_els, pt_r, pt_s, pt_x, pt_y;
  
  for(int el_y = 0; el_y < NELS_Y; el_y++) {
    for(int el_x = 0; el_x < NELS_X; el_x++) {
      shift_els = (el_y*NELS_X + el_x)*nodes_per_el;

      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        pt_r = pt_i%(elOrd+1);
        pt_s = pt_i/(elOrd+1);

        pt_x = el_x*elOrd + pt_r;
        pt_y = el_y*elOrd + pt_s;
        // asseume periodic in x
        if(pt_x == NELS_X*elOrd) pt_x = 0;

        data_els[shift_els+pt_i] = data_log[pt_y*NELS_X*elOrd + pt_x];
      }
    }
  }
}

void SEM_to_Fourier(int plane_k, Context* context, AuxField* us, real_t* data_f) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  int pt_i;
  Element* elmt;
  real_t* data = new real_t[NELS_X*elOrd];

  for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < nNodesX; pt_x++) {
      pt_i = pt_y*nNodesX + pt_x;
      elmt = context->elmt[context->el[pt_i]];
      data[pt_x] = us->probe(elmt, context->r[pt_i], context->s[pt_i], plane_k);
    }
    dDFTr(data, nNodesX, 1, +1);
    for(int pt_x = 0; pt_x < nModesX; pt_x++) {
      data_f[pt_y*nModesX + pt_x] = data[pt_x];
    }
  }

  delete[] data;
}

void Fourier_to_SEM(int plane_k, Context* context, AuxField* us, real_t* data_f) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  int pt_r;
  double dx, xr, theta;
  real_t* data = new real_t[nNodesX];
  const real_t *qx;

  Femlib::quadrature(&qx, 0, 0, 0  , elOrd+1, GLJ, 0.0, 0.0);

  dx = (XMAX - XMIN)/NELS_X;

  for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < nModesX; pt_x++) {
      data[pt_x] = data_f[pt_y*nNodesX + pt_x];
    }
    // fourier interpolation to GLL grid
    for(int pt_x = 0; pt_x < nModesX; pt_x++) {
      // coordinate in real space
      xr = XMIN + (pt_x/elOrd)*dx + (1.0 + qx[pt_x%(elOrd+1)])*dx;
      // coordinate in fourier space
      theta = 2.0*M_PI*xr/(XMAX - XMIN);

      data_f[pt_y*nNodesX + pt_x] = data[0];
      // ignore the nyquist frequency (entry [1])
      for(int mode_k = 1; mode_k < nModesX/2; mode_k++) {
        data_f[pt_y*nNodesX + pt_x] += data[2*mode_k+0]*cos(mode_k*theta);
        data_f[pt_y*nNodesX + pt_x] += data[2*mode_k+1]*sin(mode_k*theta);
      }
    }
  }

  logical_to_elements(data_f, us->plane(plane_k));

  delete[] data;
}

void UnpackX(Context* context, vector<Field*> fields, real_t* theta, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, slice_i, index;
  int nZ = Geometry::nZ();
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  AuxField* field;
  PetscScalar *xArray;
  real_t* data = new real_t[NELS_X*elOrd*NELS_Y*elOrd];

  VecGetArray(x, &xArray);

  index = 0;
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    field = fields[slice_i * context->domain->nField() + DOF];

    for(kk = 0; kk < nZ; kk++) {
      SEM_to_Fourier(kk, context, field, data);

      // skip over redundant real dofs
      for(jj = 0; jj < NELS_Y*elOrd; jj++) {
        for(ii = 0; ii < nModesX; ii++) {
          xArray[index++] = data[jj*nNodesX + ii];
        }
      }
    }
    theta[slice_i] = xArray[index++];
    tau[slice_i]   = xArray[index++];
  }

  VecRestoreArray(x, &xArray);

  delete[] data;
}

void RepackX(Context* context, vector<Field*> fields, real_t* theta, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, slice_i, index;
  int nZ = Geometry::nZ();
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  AuxField* field;
  PetscScalar *xArray;
  real_t* data = new real_t[NELS_X*elOrd*NELS_Y*elOrd];

  index = 0;
  VecGetArray(x, &xArray);
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    field = fields[slice_i * context->domain->nField() + DOF];

    for(kk = 0; kk < nZ; kk++) {
      // skip over redundant real dofs
      for(jj = 0; jj < NELS_Y*elOrd; jj++) {
        for(ii = 0; ii < nModesX; ii++) {
          data[jj*nNodesX + ii] = xArray[index++];
        }
      }

      Fourier_to_SEM(kk, context, field, data);
    }
    xArray[index++] = theta[slice_i];
    xArray[index++] = tau[slice_i];
  }

  VecRestoreArray(x, &xArray);

  delete[] data;
}

PetscErrorCode _snes_jacobian(SNES snes, Vec x, Mat J, Mat P, void* ctx) {
  cout << "jacobian assembly function - shouldn't be here... aborting.\n";

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);

  abort();
  return 0;
}

PetscErrorCode _snes_function(SNES snes, Vec x, Vec f, void* ctx) {
  Context* context = (Context*)ctx;
  const real_t dt = Femlib::value ("D_T");
  int nField = context->domain->nField();
  int slice_i, field_i, mode_i, dof_i, nStep;
  real_t ckt, skt, rTmp, cTmp;
  register real_t* data_r;
  register real_t* data_c;

  UnpackX(context, context->ui, context->theta_i, context->tau_i, x);

  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    // initialise the flow map fields with the solution fields
    for(field_i = 0; field_i < nField; field_i++) {
        *context->fi[slice_i * nField + field_i] = *context->ui[slice_i * nField + field_i];
        *context->domain->u[field_i] = *context->fi[slice_i * nField + field_i];
    }

    // solve the flow map (single step)
    nStep = (int)(context->tau_i[slice_i]/dt);
    integrate(skewSymmetric, context->domain, context->bman, context->ff, nStep);

    // phase shift
    for(mode_i = 1; mode_i < Geometry::nZ()/2; mode_i++) {
      ckt = cos(mode_i*context->theta_i[mode_i]);
      skt = sin(mode_i*context->theta_i[mode_i]);

      data_r = context->domain->u[DOF]->plane(2*mode_i+0);
      data_c = context->domain->u[DOF]->plane(2*mode_i+1);

      for(dof_i = 0; dof_i < context->nDofs; dof_i++) {
        rTmp = ckt*data_r[dof_i] - skt*data_c[dof_i];
        cTmp = skt*data_r[dof_i] + ckt*data_c[dof_i];
        data_r[dof_i] = rTmp;
        data_c[dof_i] = cTmp;
      }
    }

    // set f
    *context->fi[slice_i * nField + DOF]  = *context->domain->u[DOF];
    *context->fi[slice_i * nField + DOF] *= -1.0;
    *context->fi[slice_i * nField + DOF] += *context->ui[slice_i * nField + DOF];
  }

  RepackX(context, context->fi, context->theta_i, context->tau_i, f);

  return 0;
}

void rpo_solve(int nSlice, Mesh* mesh, vector<Element*> elmt, BCmgr* bman, Domain* domain, FieldForce* FF, vector<Field*> ui, vector<Field*> fi) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  BoundarySys* bsys;
  const NumberSys* nsys;
  real_t dx, dy, er, es;
  const real_t* qx;
  int_t pt_x, pt_y;
  const bool guess = true;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  Vec x, f;
  Mat J, P;
  SNES snes;

  context->nSlice = nSlice;
  context->mesh   = mesh;
  context->elmt   = elmt;
  context->domain = domain;
  context->bman   = bman;
  context->ff     = FF;
  context->ui     = ui;
  context->fi     = fi;

  bsys = ui[DOF]->bsys();
  nsys = bsys->Nsys(0);
  //context->nDofsPlane = nsys->nGlobal() + Geometry::nInode();
  context->nDofsPlane = (((NELS_X*elOrd)/2 + 2)*NELS_Y*elOrd);
  // add dofs for theta and tau for each time slice
  context->nDofs = Geometry::nZ() * context->nDofsPlane + 2;

  context->theta_i = new real_t[nSlice];
  context->phi_i   = new real_t[nSlice];
  context->tau_i   = new real_t[nSlice];
  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    context->theta_i[slice_i] = 0.0;
    context->phi_i[slice_i] = 0.0;
    context->tau_i[slice_i] = SLICE_DT;
  }

  // setup the fourier mapping data
  context->el = new int_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->x  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->y  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->r  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->s  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];

  dx = (XMAX - XMIN)/(NELS_X*elOrd);
  dy = (YMAX - YMIN)/(NELS_Y*elOrd);

  Femlib::quadrature(&qx, 0, 0, 0  , elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < NELS_X*elOrd*NELS_Y*elOrd; pt_i++) {
    pt_x = pt_i%(NELS_X*elOrd);
    pt_y = pt_i/(NELS_X*elOrd);
    context->x[pt_i] = XMIN + pt_x*dx;
    context->y[pt_i] = YMIN + 0.5*dy*(qx[pt_y] + 1.0); // y coordinates are still on the GLL grid
    for(int el_i = 0; el_i < mesh->nEl(); el_i++) {
      // pass er and es by reference?
      if(elmt[el_i]->locate(context->x[pt_i], context->y[pt_i], er, es, &work[0], guess)) {
        context->el[pt_i] = el_i;
        context->r[pt_i] = er;
        context->s[pt_i] = es;
        break;
      }
    }
  }

  VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nSlice*context->nDofs, &x);
  VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nSlice*context->nDofs, &f);

  MatCreate(MPI_COMM_WORLD, &J);
  MatSetType(J, MATMPIAIJ);
  MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, nSlice*context->nDofs, nSlice*context->nDofs);
  MatMPIAIJSetPreallocation(J, 4*nSlice, PETSC_NULL, 4*nSlice, PETSC_NULL);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, PETSC_DECIDE, PETSC_DECIDE, nSlice*context->nDofs, nSlice*context->nDofs);
  MatMPIAIJSetPreallocation(P, 4*nSlice, PETSC_NULL, 4*nSlice, PETSC_NULL);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESSetFunction(snes, f, _snes_function, (void*)context);
  SNESSetJacobian(snes, J, P, _snes_jacobian, (void*)context);
  SNESSetFromOptions(snes);

  RepackX(context, context->ui, context->theta_i, context->tau_i, x);
  SNESSolve(snes, NULL, x);
  RepackX(context, context->ui, context->theta_i, context->tau_i, x);

  VecDestroy(&x);
  VecDestroy(&f);
  MatDestroy(&J);
  MatDestroy(&P);
  delete[] context->el;
  delete[] context->x;
  delete[] context->y;
  delete[] context->r;
  delete[] context->s;
}

int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
#ifdef _GNU_SOURCE
  feenableexcept (FE_OVERFLOW);    // -- Force SIG8 crash on FP overflow.
#endif

  char*            session;
  bool             freeze = false;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  Domain*          domain;
  FieldForce*      FF;
  static char      help[] = "petsc";
  int              nSlice = 8;
  vector<Field*>   ui; // Solution fields for velocities, pressure at the i time slices
  vector<Field*>   fi; // Solution fields for flow maps at the i time slices
  vector<Field*>   uTmp;
  char             session_i[100];

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);

  preprocess (session, file, mesh, elmt, bman, domain, FF);

  //domain -> restart ();
  //ROOTONLY domain -> report ();
  
  // load in the time slices
  ui.resize(nSlice * domain->nField());
  fi.resize(nSlice * domain->nField());
  uTmp = domain->u;
  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    sprintf(session_i, "%s_%.4u", session, slice_i);
    FEML* file_i = new FEML(session_i);
    delete domain;
    domain = new Domain (file_i, elmt, bman);
    for(int field_i = 0; field_i < domain->nField(); field_i++) {
        ui[slice_i*domain->nField()+field_i] = domain->u[field_i];
    }
    domain->restart();
    delete file_i;
  }
  domain->u = uTmp;

  // solve the newton-rapheson problem
  rpo_solve(nSlice, mesh, elmt, bman, domain, FF, ui, fi);

  // dump the output

  Femlib::finalize ();

  PetscFinalize();

  return EXIT_SUCCESS;
}

static void getargs (int    argc   ,
		     char** argv   ,
		     bool&  freeze ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -f       ... freeze velocity field (scalar advection/diffusion only)\n"
    "  -i       ... use iterative solver for viscous steps\n"
    "  -v[v...] ... increase verbosity level\n"
    "  -chk     ... turn off checkpoint field dumps [default: selected]\n"
    "  -S|C|N   ... regular skew-symm || convective || Stokes advection\n";

  Femlib::ivalue ("ADVECTION", 1); // -- Default is alternating skew symmetric.

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'S': Femlib::ivalue ("ADVECTION", 0); break;
    case 'C': Femlib::ivalue ("ADVECTION", 2); break;
    case 'N': Femlib::ivalue ("ADVECTION", 3); break;
    case 'f':
      freeze = true;
      break;
    case 'i':
      do			// -- Only allowing ITERATIVE=1 (Viscous).
	Femlib::ivalue ("ITERATIVE", 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE",   Femlib::ivalue ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv))     Femlib::ivalue ("CHKPOINT",    0);
      else { fprintf (stdout, usage, prog); exit (EXIT_FAILURE); }
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc != 1) message (routine, "no session definition file", ERROR);
  else             session = *argv;

  Femlib::value ("DTBDX", 0.0);
}

static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			Domain*&          domain ,
			FieldForce*&      FF     )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed above in order of creation.
// ---------------------------------------------------------------------------
{

  const char routine[] = "preprocess";
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  int_t              i, np, nz, nel, procid, seed;

  // -- Initialise problem and set up mesh geometry.
  VERBOSE cout << "Building mesh ..." << endl;
  file = new FEML (session);
  mesh = new Mesh (file);
  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.
  VERBOSE cout << "Setting geometry ... ";
  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set (np, nz, nel, space);
  VERBOSE cout << "done" << endl;

  // -- If token RANSEED > 0 then initialize the random number
  //    generator based on wall clock time and process ID (i.e. a "truly"
  //    pseudo-random number).  NB: it is important to have done this
  //    before any other possible call to random number routines.

  if (Femlib::ivalue("RANSEED") > 0) {
    procid = Geometry::procID();
    seed   = -abs((procid + 1) * (char) time(NULL));
  } else seed = -1;
  Veclib::ranInit (seed);

  // -- Build all the elements.
  VERBOSE cout << "Building elements ... ";
  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);
  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.
  VERBOSE cout << "Building boundary condition manager ..." << endl;
  bman = new BCmgr (file, elmt);
  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.
  VERBOSE cout << "Building domain ..." << endl;
  domain = new Domain (file, elmt, bman);
  VERBOSE cout << "done" << endl;

  // -- Build field force.
  VERBOSE cout << "Building field force ..." << endl;
  FF = new FieldForce (domain, file);
  VERBOSE cout << "done" << endl;
}
