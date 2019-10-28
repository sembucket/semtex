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

#include "rpo_cart_utils.h"

#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

static char prog[] = "rpo";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

#define NFIELD 3
#define TESTING 1
//#define YMIN (0.0)
//#define YMAX (1.0)
//#define ZMAX (0.7391982714328925) // L = 2.pi/beta
//#define NELS_X 12
//#define NELS_Y 9

//#define XMIN (-1.0)
//#define XMAX (+1.0)

static PetscErrorCode RPOVecNormL2_Hookstep(void* ctx,Vec v,PetscScalar* norm) {
  Context* context = (Context*)ctx;
  PetscInt nDofs_l = Geometry::nZProc() * context->n_mesh_sum;
  PetscInt ind_i;
  double norm_sq, norm_l_sq;
  PetscScalar* vArray;
  Vec vl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl);

  VecScatterBegin(context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);

  norm_l_sq = 0.0;
  VecGetArray(vl, &vArray);
  for(ind_i=0; ind_i<nDofs_l; ind_i++) {
    norm_l_sq += vArray[ind_i]*vArray[ind_i];
  }
  VecRestoreArray(vl, &vArray);

  norm_sq = 0.0;
  MPI_Allreduce(&norm_l_sq, &norm_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *norm = sqrt(norm_sq);

  VecDestroy(&vl);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDot_Hookstep(void* ctx,Vec v1,Vec v2,PetscScalar* dot) {
  Context* context = (Context*)ctx;
  PetscInt nDofs_l = Geometry::nZProc() * context->n_mesh_sum;
  PetscInt ind_i;
  double dot_l;
  PetscScalar *v1Array, *v2Array;
  Vec vl1, vl2;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl1);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl2);

  VecScatterBegin(context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);

  dot_l = 0.0;
  VecGetArray(vl1, &v1Array);
  VecGetArray(vl2, &v2Array);
  for(ind_i=0; ind_i<nDofs_l; ind_i++) {
    dot_l += v1Array[ind_i]*v2Array[ind_i];
  }
  VecRestoreArray(vl1, &v1Array);
  VecRestoreArray(vl2, &v2Array);

  *dot = 0.0;
  MPI_Allreduce(&dot_l, dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecDestroy(&vl1);
  VecDestroy(&vl2);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDiff_Hookstep(void* ctx,Vec y,Vec F,PetscScalar h) {
  Context* context = (Context*)ctx;
  PetscInt nDofs_l = Geometry::nZProc() * context->n_mesh_sum;
  PetscInt ind_i;
  PetscScalar *yArray, *FArray;
  Vec yl, Fl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &yl);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &Fl);

  VecScatterBegin(context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(yl, &yArray);
  VecGetArray(Fl, &FArray);
  for(ind_i=0; ind_i<nDofs_l; ind_i++) {
    yArray[ind_i] = (yArray[ind_i] - FArray[ind_i])/h;
  }
  VecRestoreArray(yl, &yArray);
  VecRestoreArray(Fl, &FArray);

  VecScatterBegin(context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd  (context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&yl);
  VecDestroy(&Fl);

  PetscFunctionReturn(0);
}

void build_constraints(Context* context, Vec x_delta, double* f_phi, double* f_tau) {
  int          elOrd       = Geometry::nP() - 1;
  int          nNodesX     = Femlib::ivalue("NELS_X")*elOrd;
  int          index;
  int          nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int          nl          = 3 * nDofsCube_l;
  int          plane_j;
  int          pt_j, el_j, dof_j;
  double       p_y;
  double       k_z;
  double       f_phi_l, f_tau_l;
  double*      rz          = new double[nl];
  double*      rt          = new double[nl];
  double*      data_r      = new double[nNodesX * Femlib::ivalue("NELS_Y")*elOrd];
  double*      data_i      = new double[nNodesX * Femlib::ivalue("NELS_Y")*elOrd];
  PetscScalar* xArray;
  Vec          xl;
  int          nStep;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArray(xl, &xArray);

  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < 3; field_i++) {
      *context->domain->u[field_i] = *context->u0[field_i];
    }

    //delete context->analyst;
    //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);

    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);

    context->domain->time = 0.0;
    context->domain->step = 0;
    Femlib::value("t", 0.0);

    integrate(convective, context->domain, context->bman, context->analyst, context->ff);
    Femlib::ivalue("N_STEP", nStep);

    for(int field_i = 0; field_i < 3; field_i++) {
      *context->domain->u[field_i] -= *context->u0[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    }
  } 

  for(int dof_i = 0; dof_i < nl; dof_i++) { rz[dof_i] = rt[dof_i] = 0.0; }

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      elements_to_logical(context->u0[field_i]->plane(plane_i+0), data_r);
      elements_to_logical(context->u0[field_i]->plane(plane_i+1), data_i);
      dof_j = 0;
      for(int dof_i = 0; dof_i < Femlib::ivalue("NELS_Y") * elOrd * nNodesX; dof_i++) {
        // omit the boundaries
        if(dof_i % nNodesX == 0 || dof_i / nNodesX == 0) continue;

        k_z  = Femlib::value("BETA") * (plane_j / 2);

        index = field_i * nDofsCube_l + (plane_i+0) * context->nDofsPlane + dof_j;
        rz[index] = -k_z * data_i[dof_i];
        index = field_i * nDofsCube_l + (plane_i+1) * context->nDofsPlane + dof_j;
        rz[index] = +k_z * data_r[dof_i];

        dof_j++;
      }
    }
    if(!context->travelling_wave) {
      index = field_i * nDofsCube_l;

      for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
        elements_to_logical(context->domain->u[field_i]->plane(plane_i), data_r);
        for(int dof_i = 0; dof_i < Femlib::ivalue("NELS_Y") * elOrd * nNodesX; dof_i++) {
          // omit the boundaries
          if(dof_i % nNodesX == 0 || dof_i / nNodesX == 0) continue;

          rt[index++] = data_r[dof_i];
        }
      }
    }
  }

  f_phi_l = f_tau_l = 0.0;
  for(int dof_i = 0; dof_i < nl; dof_i++) {
    f_phi_l -= rz[dof_i] * xArray[dof_i];
    if(!context->travelling_wave) f_tau_l -= rt[dof_i] * xArray[dof_i];
  }
  MPI_Allreduce(&f_phi_l, f_phi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_tau_l, f_tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecRestoreArray(xl, &xArray);

  delete[] rz;
  delete[] rt;
  delete[] data_r;
  delete[] data_i;
  VecDestroy(&xl);
}

void _build_constraints(Context* context, Vec x_delta, double* f_phi, double* f_tau) {
  int          nl          = Geometry::nZProc() * context->n_mesh_sum;
  int          plane_j, dof_j, index;
  double       k_z;
  double       f_phi_l, f_tau_l;
  double*      rz          = new double[nl];
  double*      rt          = new double[nl];
  double*      data_r      = new double[context->n_mesh_max];
  double*      data_i      = new double[context->n_mesh_max];
  PetscScalar* xArray;
  Vec          xl;
  int          nStep;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArray(xl, &xArray);

  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < 3; field_i++) {
      *context->domain->u[field_i] = *context->u0[field_i];
    }

    //delete context->analyst;
    //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);

    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);
    context->domain->time = 0.0;
    context->domain->step = 0;
    Femlib::value("t", 0.0);

    integrate(convective, context->domain, context->bman, context->analyst, context->ff);
    Femlib::ivalue("N_STEP", nStep);

    for(int field_i = 0; field_i < 3; field_i++) {
      *context->domain->u[field_i] -= *context->u0[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    }
  } 

  for(int dof_i = 0; dof_i < nl; dof_i++) { rz[dof_i] = rt[dof_i] = 0.0; }

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;

      elements_to_vector(context, field_i, context->u0[field_i]->plane(plane_i+0), data_r, true);
      elements_to_vector(context, field_i, context->u0[field_i]->plane(plane_i+1), data_i, true);

      dof_j = 0;
      for(int dof_i = 0; dof_i < context->n_mesh[field_i]; dof_i++) {
        k_z  = Femlib::value("BETA") * (plane_j / 2);

        index = LocalIndex(context, field_i, plane_i+0, dof_i);
        rz[index] = -k_z * data_i[dof_i];
        index = LocalIndex(context, field_i, plane_i+1, dof_i);
        rz[index] = +k_z * data_r[dof_i];

        dof_j++;
      }
    }
    if(!context->travelling_wave) {
      for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
        elements_to_vector(context, field_i, context->u0[field_i]->plane(plane_i), data_r, true);
        for(int dof_i = 0; dof_i < context->n_mesh[field_i]; dof_i++) {
          index = LocalIndex(context, field_i, plane_i, dof_i);
          rt[index] = data_r[dof_i];
        }
      }
    }
  }

  f_phi_l = f_tau_l = 0.0;
  for(int dof_i = 0; dof_i < nl; dof_i++) {
    f_phi_l -= rz[dof_i] * xArray[dof_i];
    if(!context->travelling_wave) f_tau_l -= rt[dof_i] * xArray[dof_i];
  }
  MPI_Allreduce(&f_phi_l, f_phi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_tau_l, f_tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecRestoreArray(xl, &xArray);

  delete[] rz;
  delete[] rt;
  delete[] data_r;
  delete[] data_i;
  VecDestroy(&xl);
}

PetscErrorCode _snes_jacobian(SNES snes, Vec x, Mat J, Mat P, void* ctx) {
  Context* context = (Context*)ctx;

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  J, MAT_FINAL_ASSEMBLY);

  return 0;
}

PetscErrorCode _snes_function(SNES snes, Vec x, Vec f, void* ctx) {
  Context* context = (Context*)ctx;
  real_t dt;
  int field_i, ksp_i, snes_i, nSteps;
  real_t f_norm, x_norm;
  double dummy[2];
  double zero = 0.0;
  Vec x_snes;
  char filename[100];
  KSP ksp;
  PetscViewer viewer;
  KSPConvergedReason reason;

  SNESGetKSP(snes, &ksp);
  KSPGetConvergedReason(ksp, &reason); // if reason != 0 then gmres solve has completed
  if(!Geometry::procID()) cout << "\tksp converged reason: " << reason << endl;

  SNESGetSolution(snes, &x_snes);
  VecNorm(x_snes, NORM_2, &x_norm);
  if(!Geometry::procID()) cout << "\tx_snes: " << x_snes << ", |x_snes|: " << x_norm << endl;
  if(!reason) {
    if(!Geometry::procID()) cout << "\tunpacking constraints velocity from new_x";
    _UnpackX(context, context->u0, &dummy[0], &dummy[1], x_snes);
    MPI_Bcast(&dummy[1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else {
    // write the current state vector
    sprintf(filename, "x_curr_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(x, viewer);
    PetscViewerDestroy(&viewer);
  }

  if(!reason && context->build_dx) {
    // get the current estimate of dx
    SNESGetLastOrthoVec(snes, context->x_delta); if(!Geometry::procID()) cout << "\tusing last orthonormal vec as dx\n";
  }
  context->build_dx = true;

  _UnpackX(context, context->ui, &context->phi_i, &context->tau_i, x);

  // update the starting time for this slice
  context->domain->time = 0.0;
  context->domain->step = 0;
  Femlib::value("t", 0.0);

  // initialise the flow map fields with the solution fields
  for(field_i = 0; field_i < NFIELD; field_i++) {
    *context->domain->u[field_i] = *context->ui[field_i];
  }

  // solve the flow map for time tau_i
  MPI_Bcast(&context->tau_i, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&context->phi_i, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(!context->travelling_wave) {
    nSteps = context->tau_i / context->dt0;
    dt = context->tau_i / nSteps;
    Femlib::value("D_T", dt);
    Femlib::ivalue("N_STEP", nSteps);
  }

  if(!Geometry::procID()) {
    cout << "\trun time: " << Femlib::ivalue("N_STEP") * dt << "\tnstep: " << Femlib::ivalue("N_STEP") << "\tdt: " << Femlib::value("D_T") << endl;
    cout << scientific << "\ttau:   " << context->tau_i << "\tphi:   " << context->c_scale * context->phi_i * Femlib::value("BETA") << endl;
  }

  // don't want to call the dns analysis, use custom integrate routine instead
  //delete context->analyst;
  //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
#ifdef TESTING
  if(!reason) context->domain->dump();
#endif
  integrate(convective, context->domain, context->bman, context->analyst, context->ff);
#ifdef TESTING
  if(!reason) context->domain->dump();
#endif
  context->domain->update_bcs = true;

  _phase_shift_z(context, context->c_scale * context->phi_i * Femlib::value("BETA"), -1.0, context->domain->u);

  // set the residual vector
  for(field_i = 0; field_i < NFIELD; field_i++) {
    *context->fi[field_i]  = *context->domain->u[field_i];
    /*if(context->travelling_wave)*/ *context->fi[field_i] -= *context->ui[field_i];
    //*context->fi[field_i] -= *context->ui[field_i];
  }

  if(!reason) { // within  fgmres
    _build_constraints(context, context->x_delta, &context->f_phi, &context->f_tau);
    _RepackX(context, context->fi, context->f_phi, context->f_tau, f);
    if(!Geometry::procID()) cout << "\trepacking constraints as f_phi: "   << context->f_phi << ", f_tau: "   << context->f_tau << endl;
  } else {      // outside fgmres
    _RepackX(context, context->fi, zero, zero, f); if(!Geometry::procID()) cout << "\trepacking with zero constraints " << zero << endl;
  }

  RPOVecNormL2_Hookstep(context, x, &x_norm);
  RPOVecNormL2_Hookstep(context, f, &f_norm);
  context->iteration++;
  if(!Geometry::procID()) cout << "\t" << context->iteration << ":\tevaluating function, |x|: " << x_norm << "\t|f|: " << f_norm << endl;

  return 0;
}

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<AuxField*> ui, vector<AuxField*> fi, vector<AuxField*> u0, char* session)
{
  Context* context = new Context;
  Vec x, f;
  Mat P;
  KSP ksp;
  SNES snes;
  double norm;
  PetscBool is_fgmres;

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

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
  context->phi_i    = 0.0;
  context->tau_i    = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  context->iteration = Femlib::ivalue("RPO_LOAD_VEC");
  context->travelling_wave = Femlib::ivalue("TRAV_WAVE");
  context->dt0      = Femlib::value("D_T");
  context->build_dx = false;

  // add dofs for phi and tau for each time slice
  build_addToVector(context, context->domain->u);
  context->nDofsPlane = context->n_mesh[0];
  //context->nDofsSlice = Geometry::nZ() * context->n_mesh_sum + 1;
  //if(!Femlib::ivalue("TRAV_WAVE")) context->nDofsSlice++;
  context->localSize  = Geometry::nZProc() * context->n_mesh_sum;
  if(!Geometry::procID()) {
    context->localSize++;
    if(!Femlib::ivalue("TRAV_WAVE")) context->localSize++;
  }
  MPI_Allreduce(&context->localSize, &context->nDofsSlice, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  // assign the coordinate weights
  build_coordWeights(context);

  _assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &f);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &context->x_delta);
  VecZeroEntries(context->x_delta);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, context->localSize, context->localSize, context->nDofsSlice, context->nDofsSlice);
  MatMPIAIJSetPreallocation(P, 1, PETSC_NULL, 1, PETSC_NULL);
  MatSetOptionsPrefix(P, "P_");
  MatSetFromOptions(P);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESSetFunction(snes, f,    _snes_function, (void*)context);
  SNESSetJacobian(snes, P, P, _snes_jacobian, (void*)context);
  SNESGetKSP(snes, &ksp);
  KSPSetType(ksp, KSPGMRES);
  SNESSetType(snes, SNESNEWTONTR);
  SNESSetFromOptions(snes);

  PetscObjectTypeCompare((PetscObject)ksp, KSPFGMRES, &is_fgmres);
  if(!is_fgmres) {
    if(!Geometry::procID()) cout << "ERROR: KSP type must be FGMRES for hookstep\n";
    abort();
  }

  // set the custom operators to discount the constraint dofs
  KSPSetNorm_Hookstep(ksp,(void*)context,RPOVecNormL2_Hookstep);
  KSPSetDot_Hookstep(ksp,(void*)context,RPOVecDot_Hookstep);
  KSPSetDiff_Hookstep(ksp,(void*)context,RPOVecDiff_Hookstep);

  context->uBar = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'b');
  base_profile(context, context->domain->u[2], Femlib::value("BASE_PROFILE_SCALE"), context->uBar);
  *context->ui[2] -= *context->uBar;
  velocity_scales(context);
  *context->ui[2] += *context->uBar;

  _RepackX(context, context->ui, context->phi_i, context->tau_i, x);

  RPOVecNormL2_Hookstep(context, x, &norm);
  if(!Geometry::procID()) cout << "|x_0|: " << norm << endl;
  if(norm < 1.0e-6) {
    if(!Geometry::procID()) cout << "ERROR: initial state vector norm is SMALL! "
                                 << "Are you sure you loaded the initial condition correctly??\n";
    abort();
  }

  context->c_scale = context->tau_i / norm;
  if(!Geometry::procID()) printf("c scale: %g\n", context->c_scale);

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
  _UnpackX(context, context->ui, &context->phi_i, &context->tau_i, x);

  if(!Geometry::procID()) cout << "rpo solve complete.\n";
  if(!Geometry::procID()) cout << "\tshift phi:   " << context->c_scale * context->phi_i * Femlib::value("BETA") << endl;
  if(!Geometry::procID()) cout << "\tshift tau:   " << context->tau_i << endl;

  VecDestroy(&x);
  VecDestroy(&f);
  MatDestroy(&P);
  VecDestroy(&context->x_delta);
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
  vector<AuxField*>   u0;  // Additional fields for building the constraints
  char*            fname;
  BoundarySys*     bndry;

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  domain -> restart ();
  //ROOTONLY domain -> report ();
  
  // load in the time slices
  ui.resize(3);
  fi.resize(3);
  u0.resize(3);
  for(int field_i = 0; field_i < NFIELD; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);

    *ui[field_i] = *domain->u[field_i];
  }

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0, session);

  domain->dump();

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
