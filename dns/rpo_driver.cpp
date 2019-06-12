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

#include "rpo_utils.h"
#include "rpo_preconditioner.h"

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

#define TESTING
//#define ANALYTICALLY_FORCED
#define X_FOURIER
#define TRAVELLING_WAVE // constant phase speed, ct, so shifting for z and t leads to over determined system! 

#define NFIELD 3
//#define NFIELD 4
//#define NFIELD 2
#define XMIN 0.0
#define YMIN 0.0
#define YMAX 1.0
#define NELS_X 30
#define NELS_Y 8
//#define NELS_X 40
//#define NELS_Y 16
#define NSLICE 1
#define THREE 3
#define VEL_MAJOR 1

#ifdef ANALYTICALLY_FORCED
#define XMAX (2.0*M_PI)
#else
#define XMAX 4.053667940115862
#endif

double fieldnorm(Field* field) {
  double norm_l = 0.0;
  double norm_g = 0.0;

  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
    real_t* plane = field->plane(plane_i);
    for(int dof_i = 0; dof_i < Geometry::planeSize(); dof_i++) {
      norm_l += plane[dof_i]*plane[dof_i];
    }
  }
  MPI_Allreduce(&norm_l, &norm_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sqrt(norm_g);
}

void build_constraints(Context* context, Vec x_delta, double* f_theta, double* f_phi, double* f_tau) {
  int          elOrd       = Geometry::nP() - 1;
  int          nNodesX     = NELS_X*elOrd;
  int          nModesX     = context->nModesX;
  int          index;
  int          nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int          nl          = context->nField * nDofsCube_l;
  int          plane_j, field_j;
  int          pt_j, el_j;
  double       p_y;
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

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArray(xl, &xArray);

  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] = *context->u0[field_i];
    }
    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);
    if(!Geometry::procID()) cout << "time step in constraints evaluation: " << scientific << Femlib::value("D_T") << endl;
    //delete context->analyst;
    //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
    integrate(skewSymmetric, context->domain, context->bman, context->analyst, context->ff);
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] -= *context->u0[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    } 
    Femlib::ivalue("N_STEP", nStep);
  }

  for(int dof_i = 0; dof_i < nl; dof_i++) { rx[dof_i] = rz[dof_i] = rt[dof_i] = 0.0; }

  field_j = 0;
#ifdef VEL_MAJOR
  for(int field_i = 0; field_i < THREE; field_i++) {
    // omit the radial velocity
    if(context->nField < 3 && field_i == 1) continue;

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->u0[field_i], data_r);
      for(int node_j = 0; node_j < context->nElsY * elOrd; node_j++) {
        for(int mode_i = 0; mode_i < nModesX; mode_i++) {
          index = field_j * nDofsCube_l + plane_i * context->nDofsPlane + node_j * nModesX + mode_i;
          k_x = (2.0 * M_PI / (context->xmax /*- XMIN*/)) * (mode_i / 2);

          if(mode_i % 2 == 0) { // real part
            rx[index] = -k_x * data_r[node_j * nModesX + mode_i + 1];
          } else {              // imag part
            rx[index] = +k_x * data_r[node_j * nModesX + mode_i - 1];
          }
        }
      }
    }
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      SEM_to_Fourier(plane_i+0, context, context->u0[field_i], data_r);
      SEM_to_Fourier(plane_i+1, context, context->u0[field_i], data_i);
      for(int dof_i = 0; dof_i < context->nElsY * elOrd * nModesX; dof_i++) {
        el_j = dof_i / (nModesX * elOrd);
        pt_j = (dof_i / nModesX) % elOrd;
        p_y  = context->elmt[el_j]->_ymesh[pt_j*(elOrd+1)];
        //if(fabs(p_y) < 1.0e-6) p_y = 1.0;
        if(fabs(p_y) < 1.0e-6) {
          k_z = 0.0;
        } else {
          k_z  = (1.0 / p_y) * context->phi_i[0] * (plane_j / 2);
        }
        //k_z = 1.0 * (plane_j / 2);

        index = field_j * nDofsCube_l + (plane_i+0) * context->nDofsPlane + dof_i;
        rz[index] = -k_z * data_i[dof_i];
        index = field_j * nDofsCube_l + (plane_i+1) * context->nDofsPlane + dof_i;
        rz[index] = +k_z * data_r[dof_i];
      }
    }

    if(!context->travelling_wave) {
      for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
        SEM_to_Fourier(plane_i, context, context->domain->u[field_i], data_r);
        for(int dof_i = 0; dof_i < context->nDofsPlane; dof_i++) {
          index = field_j * nDofsCube_l + plane_i * context->nDofsPlane + dof_i;
          rt[index] = data_r[dof_i];
        }
      }
    }

    field_j++;
  }
#else
  double* data_u_r = new double[context->nDofsPlane];
  double* data_v_r = new double[context->nDofsPlane];
  double* data_w_r = new double[context->nDofsPlane];
  double* data_u_i = new double[context->nDofsPlane];
  double* data_v_i = new double[context->nDofsPlane];
  double* data_w_i = new double[context->nDofsPlane];

  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
    SEM_to_Fourier(plane_i, context, context->u0[0], data_u_r);
    SEM_to_Fourier(plane_i, context, context->u0[1], data_v_r);
    SEM_to_Fourier(plane_i, context, context->u0[2], data_w_r);

    for(int node_j = 0; node_j < context->nElsY * elOrd; node_j++) {
      for(int mode_i = 0; mode_i < nModesX; mode_i++) {
        index = THREE * ( plane_i * context->nDofsPlane + node_j * nModesX + mode_i );
        k_x = (2.0 * M_PI / (context->xmax /*- XMIN*/)) * (mode_i / 2);

        if(mode_i % 2 == 0) {
          rx[index+0] = -k_x * data_u_r[node_j * nModesX + mode_i + 1];
          rx[index+1] = -k_x * data_v_r[node_j * nModesX + mode_i + 1];
          rx[index+2] = -k_x * data_w_r[node_j * nModesX + mode_i + 1];
        } else {
          rx[index+0] = +k_x * data_u_r[node_j * nModesX + mode_i - 1];
          rx[index+1] = +k_x * data_v_r[node_j * nModesX + mode_i - 1];
          rx[index+2] = +k_x * data_w_r[node_j * nModesX + mode_i - 1];
        }
      }
    }
  }
  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
    plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;

    SEM_to_Fourier(plane_i+0, context, context->u0[0], data_u_r);
    SEM_to_Fourier(plane_i+1, context, context->u0[0], data_u_i);
    SEM_to_Fourier(plane_i+0, context, context->u0[1], data_v_r);
    SEM_to_Fourier(plane_i+1, context, context->u0[1], data_v_i);
    SEM_to_Fourier(plane_i+0, context, context->u0[2], data_w_r);
    SEM_to_Fourier(plane_i+1, context, context->u0[2], data_w_i);

    for(int dof_i = 0; dof_i < context->nElsY * elOrd * nModesX; dof_i++) {
      el_j = dof_i / (nModesX * elOrd);
      pt_j = (dof_i / nModesX) % elOrd;
      p_y  = context->elmt[el_j]->_ymesh[pt_j*(elOrd+1)];
      if(fabs(p_y) < 1.0e-6) {
        k_z = 0.0;
      } else {
        k_z  = (1.0 / p_y) * context->phi_i[0] * (plane_j / 2);
      }

      index = THREE * ( (plane_i+0) * context->nDofsPlane + dof_i ) + 0;
      rz[index] = +k_z * data_u_i[dof_i];
      index = THREE * ( (plane_i+1) * context->nDofsPlane + dof_i ) + 0;
      rz[index] = -k_z * data_u_r[dof_i];

      index = THREE * ( (plane_i+0) * context->nDofsPlane + dof_i ) + 1;
      rz[index] = +k_z * data_v_i[dof_i];
      index = THREE * ( (plane_i+1) * context->nDofsPlane + dof_i ) + 1;
      rz[index] = -k_z * data_v_r[dof_i];

      index = THREE * ( (plane_i+0) * context->nDofsPlane + dof_i ) + 2;
      rz[index] = +k_z * data_w_i[dof_i];
      index = THREE * ( (plane_i+1) * context->nDofsPlane + dof_i ) + 2;
      rz[index] = -k_z * data_w_r[dof_i];
    }
  }

  if(!context->travelling_wave) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->domain->u[0], data_u_r);
      SEM_to_Fourier(plane_i, context, context->domain->u[1], data_v_r);
      SEM_to_Fourier(plane_i, context, context->domain->u[2], data_w_r);

      for(int dof_i = 0; dof_i < context->nDofsPlane; dof_i++) {
        index = THREE * ( plane_i * context->nDofsPlane + dof_i );
        rt[index++] = data_u_r[dof_i];
        rt[index++] = data_v_r[dof_i];
        rt[index++] = data_w_r[dof_i];
      }
    }
  }

  delete[] data_u_r; delete[] data_u_i;
  delete[] data_v_r; delete[] data_v_i;
  delete[] data_w_r; delete[] data_w_i;
#endif

  f_theta_l = f_phi_l = f_tau_l = 0.0;
  for(int dof_i = 0; dof_i < nl; dof_i++) {
    f_theta_l -= rx[dof_i] * xArray[dof_i];
    f_phi_l   -= rz[dof_i] * xArray[dof_i];
    if(!context->travelling_wave) f_tau_l -= rt[dof_i] * xArray[dof_i];
  }
  MPI_Allreduce(&f_theta_l, f_theta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_phi_l,   f_phi,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_tau_l,   f_tau,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecRestoreArray(xl, &xArray);

  if(!Geometry::procID()) cout << scientific << "\tconstraints: theta: " << f_theta[0] << "\tphi: " << f_phi[0] << endl;

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

  if(context->build_PC) {
    build_preconditioner_ffs(context, P);
    //context->build_PC = false;
  }
  return 0;
}

PetscErrorCode _snes_function(SNES snes, Vec x, Vec f, void* ctx) {
  Context* context = (Context*)ctx;
  real_t dt;
  int field_i, mode_i, mode_j, dof_i, nStep;
  double f_norm, x_norm, dx_norm_2, dx_norm, u_norm;
  double f_theta[1], f_phi[1], f_tau[1], runTime, dummy[3];
  Vec dx_test;
  bool update = false;
  double shift_scale = (!context->iteration) ? 1.0 : context->shift_scale;
  char filename[100];
  PetscViewer viewer;

  // create the \delta x vector for use in the assembly of the constraints into the residual vector
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &dx_test);
  VecZeroEntries(dx_test);
  VecAXPY(dx_test, +1.0, x);
  VecAXPY(dx_test, -1.0, context->x_prev);
  VecNorm(dx_test, NORM_2, &dx_norm);
  VecNorm(x, NORM_2, &x_norm);

  if(context->iteration == 0) update = true;
  if(dx_norm/x_norm > 1.0e-6) update = true;

  if(!Geometry::procID()) cout << "\tupdate test; |dx|: " << dx_norm << "\t|x|: " << x_norm << "\t|dx|/|x|: " << dx_norm/x_norm << endl;
  VecNorm(context->x_delta, NORM_2, &dx_norm_2);
  if(!Geometry::procID()) cout << "\tactual |dx|: " << dx_norm_2 << endl;
  if(update) {
    if(!Geometry::procID()) cout << "\tunpacking solution vector" << endl;
    VecCopy(dx_test, context->x_delta);

    //UnpackX(context, context->u0, &dummy[0], &dummy[1], &dummy[2], context->x_prev, 1.0);
    UnpackX(context, context->u0, &dummy[0], &dummy[1], &dummy[2], x, 1.0);

    // write the current state vector
    sprintf(filename, "x_curr_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(x, viewer);
    PetscViewerDestroy(&viewer);
    // write the previous state vector
    sprintf(filename, "x_prev_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(context->x_prev, viewer);
    PetscViewerDestroy(&viewer);
  } else {
    if(!Geometry::procID()) cout << "\t existing solution vector" << endl;
  }
  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x, 1.0);

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
                       << "\ttheta: " << context->theta_i[0] * (2.0*M_PI/context->xmax) 
                       << "\tphi:   " << context->phi_i[0] << endl;
  }

  // don't want to call the dns analysis, use custom integrate routine instead
  //delete context->analyst;
  //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
  integrate(skewSymmetric, context->domain, context->bman, context->analyst, context->ff);
#ifdef TESTING
  for(field_i = 0; field_i < THREE; field_i++) *context->write_i->u[field_i] = *context->ui[field_i];
  *context->write_i->u[THREE] = *context->domain->u[THREE]; //dump the pressure as a sanity check
  context->write_i->dump();
#endif

  // phase shift in theta (axial direction)
  phase_shift_x(context, context->theta_i[0] * (2.0 * M_PI / context->xmax), +1.0, context->domain->u);
  // phase shift in phi (azimuthal direction)
  phase_shift_z(context, context->phi_i[0], +1.0, context->domain->u);

  // set the residual vector
  for(field_i = 0; field_i < context->nField; field_i++) {
    *context->fi[field_i]  = *context->domain->u[field_i];
    *context->fi[field_i] -= *context->ui[field_i];
  }

  build_constraints(context, context->x_delta, f_theta, f_phi, f_tau);

  RepackX(context, context->fi, f_theta, f_phi, f_tau, f, 1.0);
  //RepackX(context, context->fi, f_theta, f_phi, f_tau, f, context->shift_scale);

  // if not assembling the radial velocity into the rpo solver, then 
  // updating this from the DNS at the end of the simulation time
  //if(context->nField < 3) *context->ui[1] = *context->domain->u[1];

  VecNorm(x, NORM_2, &x_norm);
  VecNorm(f, NORM_2, &f_norm);
  context->iteration++;
  if(!Geometry::procID()) cout << "\t" << context->iteration << ":\tevaluating function, |x|: " << x_norm << "\t|f|: " << f_norm << endl;

  VecCopy(x, context->x_prev);
  VecDestroy(&dx_test);

  return 0;
}

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, vector<Field*> ui, vector<Field*> fi, vector<Field*> u0) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  real_t dx, er, es, ex, ey;
  double norm;
  const real_t* qx;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j;
  bool found;
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
  context->is_s     = NULL;
  context->is_u     = NULL;
  context->is_p     = NULL;
#ifdef X_FOURIER
  context->x_fourier = true;
#else
  context->x_fourier = false;
#endif
#ifdef TRAVELLING_WAVE
  context->travelling_wave = true;
#else
  context->travelling_wave = false;
#endif
  context->nElsX   = NELS_X;
  context->nElsY   = NELS_Y;

#ifdef TESTING
  sprintf(session_i, "tube8_tw_rpo_1");
  file_i = new FEML(session_i);
  context->write_i = new Domain(file_i, context->elmt, context->bman);
  delete file_i;
  sprintf(session_i, "tube8_tw_rpo_2");
  file_i = new FEML(session_i);
  context->write_f = new Domain(file_i, context->elmt, context->bman);
  delete file_i;
#endif

  context->theta_i = new real_t[NSLICE];
  context->phi_i   = new real_t[NSLICE];
  context->tau_i   = new real_t[NSLICE];

  //context->nModesX = nNodesX/2;
  context->nModesX = nNodesX;
  context->xmax    = XMAX;
  context->iteration = 0;
  context->dx_norm_prev = 0.0;

  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    context->theta_i[slice_i] = 0.0;//-0.1;
    context->phi_i[slice_i] = 0.0;
    context->tau_i[slice_i] = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  }

  // setup the fourier mapping data
  context->el = new int_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->r  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->s  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];

  dx = (XMAX - XMIN)/nNodesX;
  Femlib::quadrature(&qx, 0, 0, 0, elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < NELS_X*elOrd*NELS_Y*elOrd; pt_i++) {
    pt_x = pt_i%(NELS_X*elOrd);
    pt_y = pt_i/(NELS_X*elOrd);
    el_x = pt_x/elOrd;
    el_y = pt_y/elOrd;
    el_i = el_y*NELS_X + el_x;
    ex = XMIN + pt_x*dx;
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

  // add dofs for theta and tau for each time slice
#ifdef X_FOURIER
  context->nDofsPlane = context->nModesX*NELS_Y*elOrd;
#else
  context->nDofsPlane = NELS_X*elOrd*NELS_Y*elOrd;
#endif
  context->nDofsSlice = context->nField * Geometry::nZ() * context->nDofsPlane + 1; // add azimuthal phase shift dof
#ifdef X_FOURIER
  context->nDofsSlice++; // add axial phase shift dof
#endif
#ifndef TRAVELLING_WAVE
  context->nDofsSlice++; // add temporal phase shift dof
#endif
  context->localSize  = NSLICE * NFIELD * Geometry::nZProc() * context->nDofsPlane;
  if(!Geometry::procID()) context->localSize += NSLICE; // azimuthal phase shift dof
#ifdef X_FOURIER
  if(!Geometry::procID()) context->localSize += NSLICE; // axial phase shift dof
#endif
#ifndef TRAVELLING_WAVE
  if(!Geometry::procID()) context->localSize += NSLICE; // temporal phase shift dof
#endif

  assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &f);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &context->x_prev);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, NSLICE * context->nDofsSlice, &context->x_delta);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, context->localSize, context->localSize, NSLICE * context->nDofsSlice, NSLICE * context->nDofsSlice);
  //MatMPIAIJSetPreallocation(P, 2*context->nModesX, PETSC_NULL, 2*context->nModesX, PETSC_NULL);
  MatMPIAIJSetPreallocation(P, 1, PETSC_NULL, 1, PETSC_NULL);
  MatSetOptionsPrefix(P, "P_");
  MatSetFromOptions(P);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESSetFunction(snes, f,    _snes_function, (void*)context);
  SNESSetJacobian(snes, P, P, _snes_jacobian, (void*)context);
  SNESGetKSP(snes, &ksp);
  KSPSetType(ksp, KSPGMRES);
  SNESSetType(snes, SNESNEWTONTR);
  SNESSetNPCSide(snes, PC_LEFT);
  SNESSetFromOptions(snes);

  context->snes = snes;
  context->shift_scale = 1.0;//1.0e+4;
  RepackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x, 1.0);
  VecCopy(x, context->x_prev);
  VecNorm(x, NORM_2, &norm);
  if(!Geometry::procID()) cout << "|x_0|: " << norm << endl;
  if(norm < 1.0e-4) {
    if(!Geometry::procID()) cout << "ERROR: initial state vector norm is SMALL! "
                                 << "Are you sure you loaded the initial condition correctly??\n";
    MPI_Barrier(MPI_COMM_WORLD);
    abort();
  }
  context->shift_scale = 1.0;//1.0e+4;
  SNESSolve(snes, NULL, x);
  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x, 1.0);
  //UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x, context->shift_scale);

  if(!Geometry::procID()) cout << "rpo solve complete.\n";
  if(!Geometry::procID()) cout << "\tshift theta: " << context->theta_i[0] * (2.0*M_PI/context->xmax) << endl;
  if(!Geometry::procID()) cout << "\tshift phi:   " << context->phi_i[0] << endl;
  if(!Geometry::procID()) cout << "\tshift tau:   " << context->tau_i[0] << endl;

  // write out the phase shifted solution to the dns solve at time tau (phase shift already applied
  // at the end of the residual assembly routine
  sprintf(session_i, "tube8_tw_shift");
  file_i = new FEML(session_i);
  Domain* dom = new Domain(file_i, context->elmt, context->bman);
  for(int field_i = 0; field_i < THREE; field_i++) dom->u[field_i] = context->domain->u[field_i];
  dom->dump();
  delete file_i;
  delete dom;

  VecDestroy(&x);
  VecDestroy(&f);
  MatDestroy(&P);
  VecDestroy(&context->x_prev);
  VecDestroy(&context->x_delta);
  delete[] context->el;
  delete[] context->r;
  delete[] context->s;
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
  vector<Field*>   ui;  // Solution fields for velocities, pressure at the i time slices
  vector<Field*>   fi;  // Solution fields for flow maps at the i time slices
  vector<Field*>   u0;  // Initial guess for the solution velocities
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

#ifdef TRAVELLING_WAVE
  if(!Geometry::procID()) cout << "travelling wave solution\n";
#else
  if(!Geometry::procID()) cout << "NOT a travelling wave solution\n";
#endif

#ifdef VEL_MAJOR
  if(!Geometry::procID()) cout << "velocity component major indexing\n";
#else
  if(!Geometry::procID()) cout << "velocity component minor indexing\n";
#endif
  
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
    // this is the array of fields that contains the initial condition!!!
    ui[field_i] = domain->u[field_i];

    strcpy ((fname = new char [strlen (bman -> field()) + 1]), bman -> field());
    bndry = new BoundarySys(bman, elmt, fname[0]);
    fi[field_i] = new Field(bndry, new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, fname[0]);

    strcpy ((fname = new char [strlen (bman -> field()) + 1]), bman -> field());
    bndry = new BoundarySys(bman, elmt, fname[0]);
    u0[field_i] = new Field(bndry, new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, fname[0]);
  }
  delete file_i;
  delete domain;

  // allocate the temporary fields
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0);
  delete file_i;
  delete domain;

  // dump the output
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  for(int field_i = 0; field_i < THREE; field_i++) {
    *domain->u[field_i] = *ui[field_i];
  }
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
