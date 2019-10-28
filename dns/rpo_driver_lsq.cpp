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

#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

#include <fftw3.h>

#include <dns.h>
#include "rpo_utils.h"
#include "rpo_preconditioner.h"

static char prog[] = "rpo";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

#define TESTING
#define X_FOURIER

#define NFIELD 3
#define NSLICE 1
#define THREE 3

static PetscErrorCode RPOVecNormL2_Hookstep(void* ctx,Vec v,PetscScalar* norm) {
  PetscInt ierr;
  ierr = VecNorm(v, NORM_2, norm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDot_Hookstep(void* ctx,Vec v1,Vec v2,PetscScalar* dot) {
  PetscInt ierr;
  ierr = VecDot(v1, v2, dot);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDiff_Hookstep(void* ctx,Vec y,Vec F,PetscScalar h) {
  PetscInt ierr;
  ierr = VecAXPY(y,-1.0,F);CHKERRQ(ierr);
  ierr = VecScale(y,1.0/h);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void __assign_scatter_semtex(Context* context) {
  int   nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int   nDofsCube_g = Geometry::nZ()     * context->nDofsPlane;
  int   ind_i       = 0;
  int*  inds;
  IS    isl, isg;
  Vec   vl, vg;

  inds = new int[context->localSize];

  context->lShift = new int*[context->nSlice];
  for(int slice_i = 0; slice_i < context->nSlice; slice_i++) {
    context->lShift[slice_i] = new int[context->nField];

    for(int field_i = 0; field_i < context->nField; field_i++) {
      context->lShift[slice_i][field_i] = slice_i * context->nDofsSlice + 
                                          field_i * nDofsCube_g + 
                                          Geometry::procID() * nDofsCube_l;

      for(int ind_j = 0; ind_j < nDofsCube_l; ind_j++) {
        inds[ind_i] = context->lShift[slice_i][field_i] + ind_j;
        ind_i++;
      }
    }
  }

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nSlice * context->nDofsSlice, &vg);

  ISCreateStride(MPI_COMM_SELF, context->localSize, 0, 1, &isl);
  ISCreateGeneral(MPI_COMM_WORLD, context->localSize, inds, PETSC_COPY_VALUES, &isg);

  VecScatterCreate(vg, isg, vl, isl, &context->global_to_semtex);
  
  VecDestroy(&vl);
  VecDestroy(&vg);
  ISDestroy(&isl);
  ISDestroy(&isg);
  delete[] inds;
}

void __UnpackX(Context* context, vector<AuxField*> fields, Vec x, bool add_ubar) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int mode_l, index;
  real_t* data_r;
  real_t* data_i;
  AuxField* field;
  const PetscScalar *xArray;
  Vec xl;
  double scale;

  data_r = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];
  data_i = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);

  VecScatterBegin(context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArrayRead(xl, &xArray);

  for(int field_i = 0; field_i < THREE; field_i++) {
    field = fields[field_i];

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      for(int point_y = 0; point_y < ney*elOrd; point_y++) {
        for(int point_x = 0; point_x < context->nModesX; point_x++) {
          mode_l = (point_x <= context->nModesX/2) ? point_x : point_x - context->nModesX; // fftw ordering of complex data
          scale  = GetScale(context, field_i, plane_i, point_x, point_y);

          index = LocalIndex(context, field_i, plane_i+0, point_x, point_y);
          data_r[point_y*context->nModesX+point_x] = xArray[index] / scale;
          // divergence free: mean component of radial velocity is 0
          // note that we are in \tilde{} variables, and the nyquist frequency is also 0
          if(field_i  > 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) data_r[point_y*context->nModesX+point_x] = 0.0;

          index = LocalIndex(context, field_i, plane_i+1, point_x, point_y);
          data_i[point_y*context->nModesX+point_x] = xArray[index] / scale;
          // don't include the nyquist frequency
          if(field_i == 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) data_i[point_y*context->nModesX+point_x] = 0.0;
        }
      }

      Fourier_to_SEM(plane_i, context, field, data_r, data_i, field_i);
    }
  }
  VecRestoreArrayRead(xl, &xArray);

  if(add_ubar) *fields[0] += *context->uBar;

  VecDestroy(&xl);
  delete[] data_r;
  delete[] data_i;
}

void __RepackX(Context* context, vector<AuxField*> fields, Vec x, bool rmv_ubar) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int mode_l, index;
  real_t* data_r;
  real_t* data_i;
  AuxField* field;
  PetscScalar *xArray;
  Vec xl;
  double scale;

  data_r = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];
  data_i = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);

  if(rmv_ubar) *fields[0] -= *context->uBar;

  VecGetArray(xl, &xArray);
  for(int field_i = 0; field_i < THREE; field_i++) {
    field = fields[field_i];

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      SEM_to_Fourier(plane_i, context, field, data_r, data_i);

      for(int point_y = 0; point_y < ney*elOrd; point_y++) {
        for(int point_x = 0; point_x < context->nModesX; point_x++) {
          mode_l = (point_x <= context->nModesX/2) ? point_x : point_x - context->nModesX; // fftw ordering of complex data
          scale  = GetScale(context, field_i, plane_i, point_x, point_y);

          index = LocalIndex(context, field_i, plane_i+0, point_x, point_y);
          xArray[index] = data_r[point_y*context->nModesX+point_x] * scale;
          // divergence free: mean component of radial velocity is 0
          // note that we are in \tilde{} variables, and the nyquist frequency is also 0
          if(field_i  > 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) xArray[index] = 0.0;

          index = LocalIndex(context, field_i, plane_i+1, point_x, point_y);
          xArray[index] = data_i[point_y*context->nModesX+point_x] * scale;
          // don't include the nyquist frequency
          if(field_i == 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) xArray[index] = 0.0;
        }
      }
    }
  }
  VecRestoreArray(xl, &xArray);

  if(rmv_ubar) *fields[0] += *context->uBar;

  VecScatterBegin(context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&xl);
  delete[] data_r;
  delete[] data_i;
}

bool OmitAxis(int field_i, int pt_y, int plane_i) {
  if(pt_y                                     ) return false;
  if(pt_y == 0 &&                 plane_i >= 4) return true; // k > 1
  if(pt_y == 0 && field_i == 0 && plane_i >= 2) return true; // k = 1
  if(pt_y == 0 && field_i == 1 && plane_i >= 2) return true; // k = 1
  if(pt_y == 0 && field_i == 1 && plane_i <= 1) return true; // k = 0
  if(pt_y == 0 && field_i == 2 && plane_i <= 1) return true; // k = 0
  return false;
}

void build_tangents(Context* context, Vec T_theta, Vec T_phi, Vec T_tau) {
  int          elOrd       = Geometry::nP() - 1;
  int          nNodesX     = context->nElsX * elOrd;
  int          nModesX     = context->nModesX;
  int          index;
  int          nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int          nl          = context->nField * nDofsCube_l;
  int          plane_j, mode_j;
  int          pt_j, el_j;
  double       k_x, k_z;
  double*      data_r      = new double[context->nDofsPlane];
  double*      data_i      = new double[context->nDofsPlane];
  PetscScalar  *thetaArray, *phiArray, *tauArray;
  Vec          theta_l, phi_l, tau_l;
  int          nStep;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &theta_l);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &phi_l  );
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &tau_l  );

  VecGetArray(theta_l, &thetaArray);
  VecGetArray(phi_l,   &phiArray  );
  VecGetArray(tau_l,   &tauArray  );

  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] = *context->u0[field_i];
    }

    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);
    context->domain->time = 0.0;
    context->domain->step = 0;
    Femlib::value("t", 0.0);

    //delete context->analyst;
    //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
    AuxField::couple(context->domain->u[1], context->domain->u[2], INVERSE);
    integrate(convective, context->domain, context->bman, context->analyst, context->ff);
    AuxField::couple(context->domain->u[1], context->domain->u[2], FORWARD);
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] -= *context->u0[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    } 
    Femlib::ivalue("N_STEP", nStep);
  }

  for(int dof_i = 0; dof_i < nl; dof_i++) { thetaArray[dof_i] = phiArray[dof_i] = tauArray[dof_i] = 0.0; }

  for(int field_i = 0; field_i < THREE; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      SEM_to_Fourier(plane_i, context, context->u0[field_i], data_r, data_i);
      for(int node_j = 0; node_j < context->nElsY * elOrd; node_j++) {
        for(int mode_i = 0; mode_i < nModesX; mode_i++) {
          mode_j = (mode_i <= nModesX/2) ? mode_i : mode_i - nModesX; // fftw ordering of complex data

          k_x = (2.0 * M_PI / (context->xmax /*-XMIN*/)) * (mode_j);

          index = LocalIndex(context, field_i, plane_i+0, mode_i, node_j);
          thetaArray[index] = -k_x * data_i[node_j * nModesX + mode_i];
          if(OmitAxis(field_i, node_j, plane_j)) thetaArray[index] = 0.0;

          index = LocalIndex(context, field_i, plane_i+1, mode_i, node_j);
          thetaArray[index] = +k_x * data_r[node_j * nModesX + mode_i];
          if(OmitAxis(field_i, node_j, plane_j)) thetaArray[index] = 0.0;
          // remove nyquist frequency
          if(Geometry::procID() == 0 && plane_i == 0 && mode_j == 0) thetaArray[index] = 0.0;
        }
      }
    }

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      SEM_to_Fourier(plane_i, context, context->u0[field_i], data_r, data_i);
      for(int dof_i = 0; dof_i < context->nElsY * elOrd * nModesX; dof_i++) {
        k_z = 1.0 * (plane_j / 2);

        index = LocalIndex(context, field_i, plane_i+0, dof_i%nModesX, dof_i/nModesX);
        phiArray[index] = -k_z * data_i[dof_i];
        if(OmitAxis(field_i, dof_i%nModesX, plane_j)) phiArray[index] = 0.0;
        index = LocalIndex(context, field_i, plane_i+1, dof_i%nModesX, dof_i/nModesX);
        phiArray[index] = +k_z * data_r[dof_i];
        if(OmitAxis(field_i, dof_i%nModesX, plane_j)) phiArray[index] = 0.0;
        // remove nyquist frequency
        if(Geometry::procID() == 0 && plane_i == 0 && mode_j == 0) phiArray[index] = 0.0;
      }
    }

    if(!context->travelling_wave) {
      for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
        plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
        SEM_to_Fourier(plane_i, context, context->domain->u[field_i], data_r, data_i);

        for(int dof_i = 0; dof_i < context->nDofsPlane; dof_i++) {
          index = LocalIndex(context, field_i, plane_i+0, dof_i%nModesX, dof_i/nModesX);
          tauArray[index] = data_r[dof_i];
          if(OmitAxis(field_i, dof_i%nModesX, plane_j)) tauArray[index] = 0.0;

          index = LocalIndex(context, field_i, plane_i+1, dof_i%nModesX, dof_i/nModesX);
          tauArray[index] = data_i[dof_i];
          if(OmitAxis(field_i, dof_i%nModesX, plane_j)) tauArray[index] = 0.0;
          // remove nyquist frequency
          if(Geometry::procID() == 0 && plane_i == 0 && dof_i%nModesX == 0) tauArray[index] = 0.0;
        }
      }
    }
  }

  VecRestoreArray(theta_l, &thetaArray);
  VecRestoreArray(phi_l, &phiArray);
  VecRestoreArray(tau_l, &tauArray);

  VecScatterBegin(context->global_to_semtex, theta_l, T_theta, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, theta_l, T_theta, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterBegin(context->global_to_semtex, phi_l,   T_phi,   INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, phi_l,   T_phi,   INSERT_VALUES, SCATTER_REVERSE);
  if(!context->travelling_wave) {
    VecScatterBegin(context->global_to_semtex, tau_l, T_tau,   INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(  context->global_to_semtex, tau_l, T_tau,   INSERT_VALUES, SCATTER_REVERSE);
  }

  delete[] data_r;
  delete[] data_i;
  VecDestroy(&theta_l);
  VecDestroy(&phi_l);
  VecDestroy(&tau_l);
}

PetscErrorCode _snes_jacobian(SNES snes, Vec x, Mat J, Mat P, void* ctx) {
  Context* context = (Context*)ctx;

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  J, MAT_FINAL_ASSEMBLY);

  return 0;
}

PetscErrorCode _snes_function(SNES snes, Vec x, Vec f, void* ctx) {
  Context* context = (Context*)ctx;
  int field_i, ksp_i, nSteps;
  double f_norm, x_norm, dt, time;
  char filename[100];
  PetscViewer viewer;
  Vec x_snes;
  KSP ksp;
  KSPConvergedReason reason;

  SNESGetKSP(snes, &ksp);
  KSPGetConvergedReason(ksp, &reason); // if reason != 0 then gmres solve has completed
  if(!Geometry::procID()) cout << "\tksp converged reason: " << reason << endl;

  // unpack the velocity field for use in the constraints assembly from the most recent acceptable guess
  SNESGetSolution(snes, &x_snes);
  VecNorm(x_snes, NORM_2, &x_norm);
  if(!Geometry::procID()) cout << "\tx_snes: " << x_snes << ", |x_snes|: " << x_norm << endl;
  if(!reason) {
    __UnpackX(context, context->u0, x_snes, true);
    if(!Geometry::procID()) cout << "\tunpacking constraints velocity from new_x\n";
    // get the current estimate of dx
    KSPGetIterationNumber(ksp, &ksp_i);
    if(!ksp_i) {
      if(!Geometry::procID())cout<<"\tcomputing the flow rate in integrate()...\n";
      context->domain->set_grn = true;
    } else {
      context->domain->set_grn = false;
    }
  } else {
    // write the current state vector
    sprintf(filename, "x_curr_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(x, viewer);
    PetscViewerDestroy(&viewer);
  }
  context->build_dx = true;

  __UnpackX(context, context->ui, x, true);
  KSPGetShifts_Hookstep(ksp, &context->theta_i[0], &context->phi_i[0], &time);
  if(!context->travelling_wave) context->tau_i[0] = time;

  // update the starting time for this slice
  context->domain->time = 0.0;
  context->domain->step = 0;
  Femlib::value("t", 0.0);

  // initialise the flow map fields with the solution fields
  for(field_i = 0; field_i < context->nField; field_i++) {
    *context->domain->u[field_i] = *context->ui[field_i];
  }

  // solve the flow map for time tau_i
  if(!context->travelling_wave) {
    nSteps = context->tau_i[0] / context->dt0;
    dt = context->tau_i[0] / nSteps;
    Femlib::value("D_T", dt);
    Femlib::ivalue("N_STEP", nSteps);
  } 

  if(!Geometry::procID()) {
    cout << "\trun time: " << Femlib::ivalue("N_STEP") * Femlib::value("D_T") << "\tnstep: " << Femlib::ivalue("N_STEP") << "\tdt: " << Femlib::value("D_T") << endl;
    cout << scientific << "\ttau:   " << context->tau_i[0] 
                       << "\ttheta: " << context->c_scale * context->theta_i[0] * (2.0*M_PI/context->xmax)
                       << "\tphi:   " << context->c_scale * context->phi_i[0] << endl;
  }

#ifdef TESTING
  if(!reason) {
    AuxField::couple(context->domain->u[1], context->domain->u[2], INVERSE);
    context->domain->dump();
    AuxField::couple(context->domain->u[1], context->domain->u[2], FORWARD);
  }
#endif
  // don't want to call the dns analysis, use custom integrate routine instead
  //delete context->analyst;
  //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
  AuxField::couple(context->domain->u[1], context->domain->u[2], INVERSE);
  integrate(convective, context->domain, context->bman, context->analyst, context->ff);
  AuxField::couple(context->domain->u[1], context->domain->u[2], FORWARD);
#ifdef TESTING
  if(!reason) {
    AuxField::couple(context->domain->u[1], context->domain->u[2], INVERSE);
    context->domain->dump();
    AuxField::couple(context->domain->u[1], context->domain->u[2], FORWARD);
  }
#endif

  if(!Geometry::procID()) cout << "\tapplying phase shifts: " << context->c_scale * context->theta_i[0] * (2.0 * M_PI / context->xmax) 
                               << ",\t" << context->c_scale * context->phi_i[0] << endl;
  // phase shift in theta (axial direction)
  phase_shift_x(context, context->c_scale * context->theta_i[0] * (2.0 * M_PI / context->xmax), -1.0, context->domain->u);
  // phase shift in phi (azimuthal direction)
  phase_shift_z(context, context->c_scale * context->phi_i[0], -1.0, context->domain->u);

  // set the residual vector
  for(field_i = 0; field_i < context->nField; field_i++) {
    *context->fi[field_i]  = *context->domain->u[field_i];
    /*if(context->travelling_wave)*/ *context->fi[field_i] -= *context->ui[field_i];
  }

  __RepackX(context, context->fi, f, false);

  VecNorm(x, NORM_2, &x_norm);
  VecNorm(f, NORM_2, &f_norm);
  context->iteration++;
  if(!Geometry::procID()) cout << "\t" << context->iteration << ":\tevaluating function, |x|: " << x_norm << "\t|f|: " << f_norm << endl;

  return 0;
}

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<AuxField*> ui, vector<AuxField*> fi, vector<AuxField*> u0, char* session) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X") * elOrd;
  real_t dx, dy, dy_sum, er, es, ex, ey;
  double norm;
  const real_t *qx, *wx;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j;
  bool found;
  PetscBool is_fgmres;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  Vec x, f;
  Mat P;
  KSP ksp;
  SNES snes;
  int np2 = Geometry::nP() * Geometry::nP();
  int nShifts;
  Vec T_theta, T_phi, T_tau;

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

  context->nSlice   = NSLICE;
  context->nField   = NFIELD;
  context->mesh     = mesh;
  context->elmt     = elmt;
  context->domain   = domain;
  context->bman     = bman;
  context->file     = file;
  context->analyst  = analyst;
  context->ff       = FF;
  context->ui       = ui;
  context->fi       = fi;
  context->u0       = u0;
  context->build_PC = true;
#ifdef X_FOURIER
  context->x_fourier = true;
#else
  context->x_fourier = false;
#endif
  context->travelling_wave = Femlib::ivalue("TRAV_WAVE");
  context->build_dx = false;
  context->nElsX    = Femlib::ivalue("NELS_X");
  context->nElsY    = Femlib::ivalue("NELS_Y");
  context->dt0      = Femlib::value("D_T");

  context->theta_i = new real_t[NSLICE];
  context->phi_i   = new real_t[NSLICE];
  context->tau_i   = new real_t[NSLICE];

  context->nModesX = nNodesX;
  context->xmax    = Femlib::value("XMAX");
  if(!Geometry::procID())cout<<"NELS_X: "<<context->nElsX<<", NELS_Y: "<<context->nElsY<<", XMAX: "<<context->xmax<<endl;
  if( context->travelling_wave && !Geometry::procID()) cout << "      travelling wave solution!!\n";
  if(!context->travelling_wave && !Geometry::procID()) cout << "NOT a travelling wave solution!!\n";
  context->iteration = Femlib::ivalue("RPO_LOAD_VEC");

  context->theta_i[0] = 0.0;
  context->phi_i[0]   = 0.0;
  context->tau_i[0]   = Femlib::ivalue("N_STEP") * Femlib::value("D_T");

  context->domain->set_grn = true;

  // setup the fourier mapping data
  context->el = new int_t[context->nElsX*elOrd*context->nElsY*elOrd];
  context->r  = new real_t[context->nElsX*elOrd*context->nElsY*elOrd];
  context->s  = new real_t[context->nElsX*elOrd*context->nElsY*elOrd];

  context->rad_weights = new double[context->nElsY*elOrd+1];
  context->rad_coords  = new double[context->nElsY*elOrd+1];

  dx = (context->xmax/* - XMIN*/)/nNodesX;
  Femlib::quadrature(&qx, &wx, 0, 0, elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < context->nElsX*elOrd*context->nElsY*elOrd; pt_i++) {
    pt_x = pt_i%(context->nElsX*elOrd);
    pt_y = pt_i/(context->nElsX*elOrd);
    el_x = pt_x/elOrd;
    el_y = pt_y/elOrd;
    el_i = el_y*context->nElsX + el_x;
    ex = /*XMIN +*/ pt_x*dx;
    // element y size increases with distance from the boundary
    ey = elmt[el_i]->_ymesh[(pt_y%elOrd)*(elOrd+1)];

    //ex += 1.0e-8; ey += 1.0e-8;
    found = false;  
    for(el_j = 0; el_j < mesh->nEl(); el_j++) {
      er = es = 0.0;
      if(elmt[el_j]->locate(ex, ey, er, es, &work[0], true) && !found) {
        if(fabs(er) < 1.0000000001 && fabs(es) < 1.0000000001) {
          context->el[pt_i] = el_j;

          if(er > +0.99999999) er = +0.99999999;
          if(er < -0.99999999) er = -0.99999999;
          if(es > +0.99999999) es = +0.99999999;
          if(es < -0.99999999) es = -0.99999999;

          context->r[pt_i] = er;
          context->s[pt_i] = es;
          found = true;
        }
      }
      if(found) break;
    }
    if(!found && !Geometry::procID()) {
      cout << Geometry::procID() << "\t:ERROR! element does not contain point: " << pt_i << "\tx: " << ex << "\ty: " << ey << endl;
      cout << "\tfound = " << found << endl;
      cout << "\tel i:   " << el_j << "\tnum els: " << mesh->nEl() << endl;
      cout << "\tpt x: " << pt_x << "\tpt y: " << pt_y << endl;
      abort();
    }
  }

  // compute the radial weights
  for(int pt_i = 0; pt_i <= context->nElsY*elOrd; pt_i++) context->rad_weights[pt_i] = 0.0;

  dy_sum = 0.0;
  for(int el_y = 0; el_y < context->nElsY; el_y++) {
    el_i = el_y*context->nElsX;
    dy = fabs(elmt[el_i]->_ymesh[elOrd*(elOrd+1)] - elmt[el_i]->_ymesh[0]);
    if(!Geometry::procID())printf("%d\tdy: %g\n",el_y,dy);
    for(int qp_i = 0; qp_i <= elOrd; qp_i++) {
      pt_y = el_y*elOrd + qp_i;
      context->rad_weights[pt_y] += 0.5 * dy * wx[qp_i];
      context->rad_coords[pt_y]   = dy_sum + 0.5 * dy * (qx[qp_i]+1.0);
    }
    dy_sum += dy;
  }
  if(!Geometry::procID())printf("\tdy_sum: %g\n",dy_sum);

  // add dofs for theta and tau for each time slice
#ifdef X_FOURIER
  context->nDofsPlane = context->nModesX*context->nElsY*elOrd;
#else
  context->nDofsPlane = context->nElsX*elOrd*context->nElsY*elOrd;
#endif
  context->nDofsSlice = context->nField * Geometry::nZ() * context->nDofsPlane;

  context->localSize  = NFIELD * Geometry::nZProc() * context->nDofsPlane;

  __assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &f);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, context->localSize, context->localSize, context->nDofsSlice, context->nDofsSlice);
  MatMPIAIJSetPreallocation(P, 1, PETSC_NULL, 1, PETSC_NULL);
  MatSetOptionsPrefix(P, "P_");
  MatSetFromOptions(P);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESGetKSP(snes, &ksp);
  SNESSetFunction(snes, f,    _snes_function, (void*)context);
  SNESSetJacobian(snes, P, P, _snes_jacobian, (void*)context);
  KSPSetType(ksp, KSPFGMRES);
  SNESSetType(snes, SNESNEWTONLS);
  SNESSetFromOptions(snes);

  KSPSetNorm_Hookstep(ksp,(void*)context,RPOVecNormL2_Hookstep);
  KSPSetDot_Hookstep(ksp,(void*)context,RPOVecDot_Hookstep);
  KSPSetDiff_Hookstep(ksp,(void*)context,RPOVecDiff_Hookstep);

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  // build the tangent vectors
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &T_theta);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &T_phi  );
  if(!context->travelling_wave) VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &T_tau);
  build_tangents(context, T_theta, T_phi, T_tau);
  nShifts = (context->travelling_wave) ? 2 : 3;
  KSPSetTangentVecs_Hookstep(ksp, nShifts, T_theta, T_phi, T_tau);
  KSPSetShifts_Hookstep(ksp, &context->theta_i[0], &context->phi_i[0], &context->tau_i[0]);

  context->uBar = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'p');
  base_profile(context, context->domain->u[0], Femlib::value("BASE_PROFILE_SCALE"), context->uBar);
  *context->ui[0] -= *context->uBar;
  velocity_scales(context);
  *context->ui[0] += *context->uBar;

  __RepackX(context, context->ui, x, true);
  VecNorm(x, NORM_2, &norm);
  if(!Geometry::procID()) cout << "|x_0|: " << norm << endl;
  if(norm < 1.0e-4) {
    if(!Geometry::procID()) cout << "ERROR: initial state vector norm is SMALL! "
                                 << "Are you sure you loaded the initial condition correctly??\n";
    abort();
  }

  context->c_scale = context->tau_i[0] / norm;
  if(!Geometry::procID()) printf("c scale: %g\n", context->c_scale);

  PetscObjectTypeCompare((PetscObject)ksp, KSPFGMRES, &is_fgmres);
  if(!is_fgmres) {
    if(!Geometry::procID()) cout << "ERROR: KSP type must be FGMRES for hookstep\n";
    abort();
  }

  // load state from petsc vectors
  if(Femlib::ivalue("RPO_LOAD_VEC")) {
    char filename[100];
    PetscViewer viewer;

    if(!Geometry::procID()) cout << "loading vectors at iteration: " << context->iteration << endl;
    sprintf(filename, "x_curr_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    VecZeroEntries(x);
    VecLoad(x, viewer);
    PetscViewerDestroy(&viewer);
  }

  SNESSolve(snes, NULL, x);
  __UnpackX(context, context->ui, x, true);

  if(!Geometry::procID()) cout << "rpo solve complete.\n";
  if(!Geometry::procID()) cout << "\tshift theta: " << context->c_scale * context->theta_i[0] * (2.0*M_PI/context->xmax) << endl;
  if(!Geometry::procID()) cout << "\tshift phi:   " << context->c_scale * context->phi_i[0] << endl;
  if(!Geometry::procID()) cout << "\tshift tau:   " << context->tau_i[0] << endl;

  VecDestroy(&x);
  VecDestroy(&f);
  VecDestroy(&T_theta);
  VecDestroy(&T_phi);
  if(!context->travelling_wave) VecDestroy(&T_tau);
  MatDestroy(&P);
  delete[] context->el;
  delete[] context->r;
  delete[] context->s;
  delete[] context->rad_weights;
  delete[] context->rad_coords;
}

int main (int argc, char** argv) {
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
  DNSAnalyser*     analyst;
  FieldForce*      FF;
  static char      help[] = "petsc";
  vector<AuxField*>   ui;  // Solution fields for velocities, pressure at the i time slices
  vector<AuxField*>   fi;  // Solution fields for flow maps at the i time slices
  vector<AuxField*>   u0;  // Initial guess for the solution velocities

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  domain -> restart ();
  //ROOTONLY domain -> report ();

  // load in the time slices
  ui.resize(NSLICE * THREE);
  fi.resize(NSLICE * THREE);
  u0.resize(NSLICE * THREE);

  AuxField::couple(domain->u[1], domain->u[2], FORWARD);

  for(int field_i = 0; field_i < THREE/*NFIELD*/; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);

    *ui[field_i] = *domain->u[field_i];
  }

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0, session);

  AuxField::couple(domain->u[1], domain->u[2], INVERSE);

  domain->dump();
  delete domain;

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

  //if   (argc != 1) message (routine, "no session definition file", ERROR);
  //else             session = *argv;
  session = *argv;

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
