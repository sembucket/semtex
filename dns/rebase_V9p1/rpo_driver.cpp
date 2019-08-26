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
//#define RM_2FOLD_SYM

static PetscErrorCode RPOVecNormL2_Hookstep(void* ctx,Vec v,PetscScalar* norm) {
  Context* context = (Context*)ctx;
  PetscInt nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  PetscInt ind_i;
  double norm_orig, norm_sq, norm_l_sq = 0.0;
  PetscScalar* vArray;
  Vec vl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl);

  VecScatterBegin(context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);

#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==0) {
#endif

  //norm_l_sq = 0.0;
  VecGetArray(vl, &vArray);
  for(ind_i=0; ind_i<3*nDofsCube_l; ind_i++) {
    norm_l_sq += vArray[ind_i]*vArray[ind_i];
  }
  VecRestoreArray(vl, &vArray);

#ifdef RM_2FOLD_SYM
  }
#endif

  norm_sq = 0.0;
  MPI_Allreduce(&norm_l_sq, &norm_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *norm = sqrt(norm_sq);

  VecNorm(v,NORM_2,&norm_orig);
  if(!Geometry::procID())printf("\tSNES TR - modified norm - |x|: %g, |x_orig|: %g\n",*norm,norm_orig);

  VecDestroy(&vl);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDot_Hookstep(void* ctx,Vec v1,Vec v2,PetscScalar* dot) {
  Context* context = (Context*)ctx;
  PetscInt nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  PetscInt ind_i;
  double dot_l = 0.0;
  PetscScalar *v1Array, *v2Array;
  Vec vl1, vl2;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl1);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl2);

  VecScatterBegin(context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);

#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==0) {
#endif

  //dot_l = 0.0;
  VecGetArray(vl1, &v1Array);
  VecGetArray(vl2, &v2Array);
  for(ind_i=0; ind_i<3*nDofsCube_l; ind_i++) {
    dot_l += v1Array[ind_i]*v2Array[ind_i];
  }
  VecRestoreArray(vl1, &v1Array);
  VecRestoreArray(vl2, &v2Array);

#ifdef RM_2FOLD_SYM
  }
#endif

  *dot = 0.0;
  MPI_Allreduce(&dot_l, dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecDestroy(&vl1);
  VecDestroy(&vl2);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDiff_Hookstep(void* ctx,Vec y,Vec F,PetscScalar h) {
  Context* context = (Context*)ctx;
  PetscInt nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  PetscInt ind_i;
  PetscScalar *yArray, *FArray;
  Vec yl, Fl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &yl);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &Fl);

  VecScatterBegin(context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);

#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==0) {
#endif

  VecGetArray(yl, &yArray);
  VecGetArray(Fl, &FArray);
  for(ind_i=0; ind_i<3*nDofsCube_l; ind_i++) {
    yArray[ind_i] = (yArray[ind_i] - FArray[ind_i])/h;
  }
  VecRestoreArray(yl, &yArray);
  VecRestoreArray(Fl, &FArray);

#ifdef RM_2FOLD_SYM
  }
#endif

  VecScatterBegin(context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd  (context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&yl);
  VecDestroy(&Fl);

  PetscFunctionReturn(0);
}

void build_constraints(Context* context, Vec x_delta, double* f_theta, double* f_phi, double* f_tau) {
  int          elOrd       = Geometry::nP() - 1;
  int          nNodesX     = context->nElsX * elOrd;
  int          nModesX     = context->nModesX;
  int          index;
  int          nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int          nl          = context->nField * nDofsCube_l;
  int          plane_j, mode_j;
  int          pt_j, el_j;
  double       k_x, k_z;
  double       f_theta_l, f_phi_l, f_tau_l;
  double*      rx          = new double[nl];
  double*      rz          = new double[nl];
  double*      rt          = new double[nl];
  double*      data_r      = new double[context->nDofsPlane];
  double*      data_i      = new double[context->nDofsPlane];
  PetscScalar* xArray;
  Vec          xl;
  int          nStep;

  f_theta_l = f_phi_l = f_tau_l = 0.0;
  
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);

#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==0) {
#endif

  VecGetArray(xl, &xArray);

  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] = *context->u0[field_i];
    }

    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);
    if(!Geometry::procID()) cout << "\ttime step in constraints evaluation: " << scientific << Femlib::value("D_T") << endl;

    context->domain->time = 0.0;
    context->domain->step = 0;
    Femlib::value("t", 0.0);

    //delete context->analyst;
    //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
    integrate(convective, context->domain, context->bman, context->analyst, context->ff);
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] -= *context->u0[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    } 
    Femlib::ivalue("N_STEP", nStep);
  }

  for(int dof_i = 0; dof_i < nl; dof_i++) { rx[dof_i] = rz[dof_i] = rt[dof_i] = 0.0; }

  for(int field_i = 0; field_i < THREE; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      SEM_to_Fourier(plane_i, context, context->u0[field_i], data_r, data_i);
      for(int node_j = 0; node_j < context->nElsY * elOrd; node_j++) {
        for(int mode_i = 0; mode_i < nModesX; mode_i++) {
          mode_j = (mode_i <= nModesX/2) ? mode_i : mode_i - nModesX; // fftw ordering of complex data

          k_x = (2.0 * M_PI / (context->xmax /*-XMIN*/)) * (mode_j);

          index = LocalIndex(context, field_i, plane_i+0, mode_i, node_j);
          rx[index] = -k_x * data_i[node_j * nModesX + mode_i];

          index = LocalIndex(context, field_i, plane_i+1, mode_i, node_j);
          rx[index] = +k_x * data_r[node_j * nModesX + mode_i];
        }
      }
    }

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      SEM_to_Fourier(plane_i, context, context->u0[field_i], data_r, data_i);
      for(int dof_i = 0; dof_i < context->nElsY * elOrd * nModesX; dof_i++) {
        k_z = 1.0 * (plane_j / 2);

        index = LocalIndex(context, field_i, plane_i+0, dof_i%nModesX, dof_i/nModesX);
        rz[index] = -k_z * data_i[dof_i];
        index = LocalIndex(context, field_i, plane_i+1, dof_i%nModesX, dof_i/nModesX);
        rz[index] = +k_z * data_r[dof_i];
      }
    }

    if(!context->travelling_wave) {
      for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
        SEM_to_Fourier(plane_i, context, context->domain->u[field_i], data_r, data_i);

        for(int dof_i = 0; dof_i < context->nDofsPlane; dof_i++) {
          index = LocalIndex(context, field_i, plane_i+0, dof_i%nModesX, dof_i/nModesX);
          rt[index] = data_r[dof_i];

          index = LocalIndex(context, field_i, plane_i+1, dof_i%nModesX, dof_i/nModesX);
          rt[index] = data_i[dof_i];
        }
      }
    }
  }

  //f_theta_l = f_phi_l = f_tau_l = 0.0;
  for(int dof_i = 0; dof_i < nl; dof_i++) {
    f_theta_l -= rx[dof_i] * xArray[dof_i];
    f_phi_l   -= rz[dof_i] * xArray[dof_i];
    if(!context->travelling_wave) f_tau_l -= rt[dof_i] * xArray[dof_i];
  }

  VecRestoreArray(xl, &xArray);

#ifdef RM_2FOLD_SYM
  }
#endif

  MPI_Allreduce(&f_theta_l, f_theta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_phi_l,   f_phi,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_tau_l,   f_tau,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete[] rx;
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
  int field_i;
  double f_norm, x_norm, dx_norm, runTime, dummy[3], zero = 0.0;
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
    if(!Geometry::procID()) cout << "\tunpacking constraints velocity from new_x\n";
    UnpackX(context, context->u0, &dummy[0], &dummy[1], &dummy[2], x_snes);
    // get the current estimate of dx
    if(context->build_dx) SNESGetLastOrthoVec(snes, context->x_delta);
    //if(context->build_dx) KSPBuildSolution(ksp, context->x_delta, NULL); if(!Geometry::procID())cout<<"\tusing ksp solution as dx\n";
  } else {
    // write the current state vector
    sprintf(filename, "x_curr_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(x, viewer);
    PetscViewerDestroy(&viewer);
  }

  context->build_dx = true;
  RPOVecNormL2_Hookstep(context, context->x_delta, &dx_norm);
  RPOVecNormL2_Hookstep(context, x, &x_norm);
  if(!Geometry::procID()) cout << "\t|dx|: " << dx_norm << "\t|x|: " << x_norm << "\t|dx|/|x|: " << dx_norm/x_norm << endl;

  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  // update the starting time for this slice
  context->domain->time = 0.0;
  context->domain->step = 0;
  Femlib::value("t", 0.0);

  // initialise the flow map fields with the solution fields
  for(field_i = 0; field_i < context->nField; field_i++) {
    *context->domain->u[field_i] = *context->ui[field_i];
  }

  // solve the flow map for time tau_i
  MPI_Bcast(&context->tau_i[0]  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&context->theta_i[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&context->phi_i[0]  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  dt = context->tau_i[0] / Femlib::ivalue("N_STEP");
  Femlib::value("D_T", dt);
  runTime = Femlib::ivalue("N_STEP") * dt;

  if(!Geometry::procID()) {
    cout << "\trun time: " << runTime << "\tnstep: " << Femlib::ivalue("N_STEP") << "\tdt: " << Femlib::value("D_T") << endl;
    cout << scientific << "\ttau:   " << context->tau_i[0] 
                       << "\ttheta: " << context->c_scale * context->theta_i[0] * (2.0*M_PI/context->xmax)
                       << "\tphi:   " << context->c_scale * context->phi_i[0] << endl;
  }

#ifdef TESTING
  if(!reason) context->domain->dump();
#endif
  // don't want to call the dns analysis, use custom integrate routine instead
  //delete context->analyst;
  //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
  integrate(convective, context->domain, context->bman, context->analyst, context->ff);
#ifdef TESTING
  if(!reason) context->domain->dump();
#endif

  // phase shift in theta (axial direction)
  //phase_shift_x(context, context->theta_i[0] * (2.0 * M_PI / context->xmax), -1.0, context->domain->u);
  phase_shift_x(context, context->c_scale * context->theta_i[0] * (2.0 * M_PI / context->xmax), -1.0, context->domain->u);
  // phase shift in phi (azimuthal direction)
  //phase_shift_z(context, context->phi_i[0], -1.0, context->domain->u);
  phase_shift_z(context, context->c_scale * context->phi_i[0], -1.0, context->domain->u);

  // set the residual vector
  for(field_i = 0; field_i < context->nField; field_i++) {
    *context->fi[field_i]  = *context->domain->u[field_i];
    //*context->fi[field_i] -= *context->ui[field_i];
    //*context->fi[field_i] -= *context->u0[field_i];
  }

  if(!reason) { // within  fgmres
    build_constraints(context, context->x_delta, &context->f_theta, &context->f_phi, &context->f_tau);

    if(!Geometry::procID()) cout << "\trepacking constraints as f_theta: " << context->f_theta 
                                                          << ", f_phi: "   << context->f_phi 
                                                          << ", f_tau: "   << context->f_tau << endl;
    RepackX(context, context->fi, &context->f_theta, &context->f_phi, &context->f_tau, f);
  } else {      // outside fgmres
    RepackX(context, context->fi, &zero, &zero, &zero, f); if(!Geometry::procID())cout<<"\trepacking with zero constraints "<<zero<<endl;
  }

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
  FEML* file_i;
  char session_i[100];

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

  context->theta_i = new real_t[NSLICE];
  context->phi_i   = new real_t[NSLICE];
  context->tau_i   = new real_t[NSLICE];

  context->nModesX = nNodesX;
  context->xmax    = Femlib::value("XMAX");
  if(!Geometry::procID())cout<<"NELS_X: "<<context->nElsX<<", NELS_Y: "<<context->nElsY<<", XMAX: "<<context->xmax<<endl;
  if( context->travelling_wave && !Geometry::procID()) cout << "      travelling wave solution\n";
  if(!context->travelling_wave && !Geometry::procID()) cout << "NOT a travelling wave solution\n";
  context->iteration = 0;

  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    context->theta_i[slice_i] = 0.0;
    context->phi_i[slice_i]   = 0.0;
    context->tau_i[slice_i]   = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  }

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

  // set the velocity scales
  //context->u_scale[0] = 4.8721670664372931; context->u_scale[1] = 37.842234715773266; context->u_scale[2] = 42.182138079742955;
  //context->u_scale[0] = 0.30937660; context->u_scale[1] = 5.7519965; context->u_scale[2] = 3.7164474;
  context->u_scale[0] = Femlib::value("U_SCALE");
  context->u_scale[1] = Femlib::value("V_SCALE");
  context->u_scale[2] = Femlib::value("W_SCALE");
  //context->c_scale    = Femlib::value("C_SCALE");
  if(!Geometry::procID()) printf("u scales: %g, %g, %g\n", context->u_scale[0], context->u_scale[1], context->u_scale[2]);

  // add dofs for theta and tau for each time slice
#ifdef X_FOURIER
  context->nDofsPlane = context->nModesX*context->nElsY*elOrd;
#else
  context->nDofsPlane = context->nElsX*elOrd*context->nElsY*elOrd;
#endif

  context->nDofsSlice = context->nField * Geometry::nZ() * context->nDofsPlane + 1; // add azimuthal phase shift dof

#ifdef X_FOURIER
  context->nDofsSlice++; // add axial phase shift dof
#endif

  if(!context->travelling_wave) context->nDofsSlice++; // add temporal phase shift dof

  context->localSize  = NSLICE * NFIELD * Geometry::nZProc() * context->nDofsPlane;

  if(!Geometry::procID()) context->localSize += NSLICE; // azimuthal phase shift dof
#ifdef X_FOURIER
  if(!Geometry::procID()) context->localSize += NSLICE; // axial phase shift dof
#endif
  if(!context->travelling_wave && !Geometry::procID()) context->localSize += NSLICE; // temporal phase shift dof

#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==1) context->localSize = 0;
  MPI_Allreduce(&context->localSize, &context->nDofsSlice, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &f);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &context->x_delta);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, context->localSize, context->localSize, NSLICE * context->nDofsSlice, NSLICE * context->nDofsSlice);
  MatMPIAIJSetPreallocation(P, 1, PETSC_NULL, 1, PETSC_NULL);
  MatSetOptionsPrefix(P, "P_");
  MatSetFromOptions(P);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESGetKSP(snes, &ksp);
  SNESSetFunction(snes, f,    _snes_function, (void*)context);
  SNESSetJacobian(snes, P, P, _snes_jacobian, (void*)context);
  KSPSetType(ksp, KSPFGMRES);
  SNESSetType(snes, SNESNEWTONTR);
  //SNESSetNPCSide(snes, PC_LEFT);
  SNESSetFromOptions(snes);

  KSPSetNorm_Hookstep(ksp,(void*)context,RPOVecNormL2_Hookstep);
  KSPSetDot_Hookstep(ksp,(void*)context,RPOVecDot_Hookstep);
  KSPSetDiff_Hookstep(ksp,(void*)context,RPOVecDiff_Hookstep);

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  RepackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);
  VecNorm(x, NORM_2, &norm);
  if(!Geometry::procID()) cout << "|x_0|: " << norm << endl;
  if(norm < 1.0e-4) {
    if(!Geometry::procID()) cout << "ERROR: initial state vector norm is SMALL! "
                                 << "Are you sure you loaded the initial condition correctly??\n";
    abort();
  }
  VecZeroEntries(context->x_delta);
  // set the shift scale
#ifdef RM_2FOLD_SYM
  context->c_scale = context->tau_i[0] / norm;
  //context->c_scale *= 100.0;
  context->c_scale *= 55.37123711636364;
#else
  context->c_scale = 5179.3037;
#endif
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

    context->iteration = Femlib::ivalue("RPO_LOAD_VEC");
    if(!Geometry::procID()) cout << "loading vectors at iteration: " << context->iteration << endl;

    sprintf(filename, "x_curr_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    VecZeroEntries(x);
    VecLoad(x, viewer);
    PetscViewerDestroy(&viewer);
  }

  SNESSolve(snes, NULL, x);
  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  if(!Geometry::procID()) cout << "rpo solve complete.\n";
  if(!Geometry::procID()) cout << "\tshift theta: " << context->theta_i[0] * (2.0*M_PI/context->xmax) << endl;
  if(!Geometry::procID()) cout << "\tshift phi:   " << context->phi_i[0] << endl;
  if(!Geometry::procID()) cout << "\tshift tau:   " << context->tau_i[0] << endl;

  VecDestroy(&x);
  VecDestroy(&f);
  MatDestroy(&P);
  VecDestroy(&context->x_delta);
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
  char             session_i[100];
  FEML*            file_i;
  char*            fname;
  BoundarySys*     bndry;

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  //domain -> restart ();
  //ROOTONLY domain -> report ();

  // load in the time slices
  ui.resize(NSLICE * THREE);
  fi.resize(NSLICE * THREE);
  u0.resize(NSLICE * THREE);
  delete file;
  delete domain;
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();
  for(int field_i = 0; field_i < THREE/*NFIELD*/; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);

    *ui[field_i] = *domain->u[field_i];
  }
/*
  delete file_i;
  delete domain;

  // allocate the temporary fields
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();
*/

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0, session);
/*
  delete file_i;
  delete domain;

  // dump the output
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  for(int field_i = 0; field_i < THREE; field_i++) {
    *domain->u[field_i] = *ui[field_i];
  }
*/
  domain->dump();
  delete file_i;
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
