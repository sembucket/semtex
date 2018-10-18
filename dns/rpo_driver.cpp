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
//#include "rpo_integrate.h"
#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

static char prog[] = "dns";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

static int myRank;

struct Context {
    int              nSlice;
    int              nDofsSlice;
    int              nDofsPlane;
    int              it;
    Mesh*            mesh;
    vector<Element*> elmt;
    Domain*          domain;
    BCmgr*           bman;
    DNSAnalyser*     analyst;
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
    // parallel vector scattering data
    int              localSize;
    int              localShift;
    IS               isl;
    IS               isg;
    VecScatter       ltog;
};

#define XMIN 0.0
#define XMAX (2.0*M_PI)
#define YMIN 0.0
#define YMAX 0.5
#define NELS_X 10
#define NELS_Y 12

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

void UnpackX(Context* context, vector<Field*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, slice_i, field_i, index;
  int nZ = Geometry::nZProc();
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  AuxField* field;
  PetscScalar *xArray;
  real_t* data = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->ltog, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->ltog, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArray(xl, &xArray);

  index = 0;
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    for(field_i = 0; field_i < context->domain->nField(); field_i++) {
      field = fields[slice_i * context->domain->nField() + field_i];

      for(kk = 0; kk < nZ; kk++) {
        SEM_to_Fourier(kk, context, field, data);

        // skip over redundant real dofs
        for(jj = 0; jj < NELS_Y*elOrd; jj++) {
          for(ii = 0; ii < nModesX; ii++) {
            data[jj*nNodesX + ii] = xArray[index++];
          }
        }
      }
    }
    // phase shift data lives on the 0th processors part of the vector
    if(!myRank) {
      theta[slice_i] = xArray[index++];
      phi[slice_i]   = xArray[index++];
      tau[slice_i]   = xArray[index++];
    }
  }
  VecRestoreArray(xl, &xArray);
  VecDestroy(&xl);

  delete[] data;
}

void RepackX(Context* context, vector<Field*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, slice_i, field_i, index;
  int nZ = Geometry::nZProc();
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  AuxField* field;
  PetscScalar *xArray;
  real_t* data = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecGetArray(xl, &xArray);

  index = 0;
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    for(field_i = 0; field_i < context->domain->nField(); field_i++) {
      field = fields[slice_i * context->domain->nField() + field_i];

      for(kk = 0; kk < nZ; kk++) {
        // skip over redundant real dofs
        for(jj = 0; jj < NELS_Y*elOrd; jj++) {
          for(ii = 0; ii < nModesX; ii++) {
            xArray[index++] = data[jj*nNodesX + ii];
          }
        }

        Fourier_to_SEM(kk, context, field, data);
      }
    }
    // phase shift data lives on the 0th processors part of the vector
    if(!myRank) {
      xArray[index++] = theta[slice_i];
      xArray[index++] = phi[slice_i];
      xArray[index++] = tau[slice_i];
    }
  }
  VecRestoreArray(xl, &xArray);
  VecScatterBegin(context->ltog, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->ltog, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecDestroy(&xl);

  delete[] data;
}

PetscErrorCode _snes_jacobian(SNES snes, Vec x, Mat J, Mat P, void* ctx) {
  Context* context = (Context*)ctx;
  int index = 0;
  int nZ = Geometry::nZProc();
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  int el_j, pt_j;
  int shift_proc;
  const real_t *qw, *DV, *DT;
  real_t waveNum, val, det;
  real_t *drdx, *dsdx, *drdy, *dsdy;

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

  if(context->it) return 0;

  MatZeroEntries(P);

  Femlib::quadrature(0, &qw, &DV, &DT, elOrd+1, GLJ, 0.0, 0.0);

  // assemble the schur complement preconditioner
  for(int slice_i = 0; slice_i < context->nSlice; slice_i++) {
    for(int kk = 0; kk < nZ; kk++) {
      index  = slice_i*context->nDofsSlice + kk;
      // global shift
      index += context->localShift;

      for(int field_i = 0; field_i < context->domain->nField(); field_i++) {
        for(int jj = 0; jj < NELS_Y*elOrd; jj++) {
          el_j = jj/elOrd;
          pt_j = (jj%elOrd)*(elOrd+1);
          context->elmt[el_j]->lTog(drdx, dsdx, drdy, dsdy);
          det = 1.0/(drdx[pt_j]*dsdy[pt_j] - drdy[pt_j]*dsdx[pt_j]);
          for(int ii = 0; ii < nModesX; ii++) {
            waveNum = (2.0*M_PI*ii)/(XMAX - XMIN);
            val  = -1.0*waveNum*waveNum;
            val *= dsdy[pt_j]*dsdy[pt_j]*DV[pt_j]*DT[pt_j]*qw[pt_j];
            val *= det;
            // assume contributions from both elements are the same for nodes on element boundaries
            if(pt_j == 0) val *= 2.0;
            MatSetValues(P, 1, &index, 1, &index, &val, ADD_VALUES);
            index++;
          }
        }
      }
    }
    // set preconditioner values for the theta, phi and tau phase shifts
    if(!myRank) {
      val = 1.0;
      MatSetValues(P, 1, &index, 1, &index, &val, ADD_VALUES);
      index++;
      MatSetValues(P, 1, &index, 1, &index, &val, ADD_VALUES);
      index++;
      MatSetValues(P, 1, &index, 1, &index, &val, ADD_VALUES);
      index++;
    }
  }
  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);

  context->it++;

  return 0;
}

PetscErrorCode _snes_function(SNES snes, Vec x, Vec f, void* ctx) {
  Context* context = (Context*)ctx;
  const real_t dt = Femlib::value ("D_T");
  int nField = context->domain->nField();
  int slice_i, field_i, mode_i, dof_i, nStep;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  real_t ckt, skt, rTmp, cTmp;
  register real_t* data_r;
  register real_t* data_c;
  real_t* data_f = new real_t[NELS_Y*elOrd*nModesX];

  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    // initialise the flow map fields with the solution fields
    for(field_i = 0; field_i < nField; field_i++) {
        *context->fi[slice_i * nField + field_i] = *context->ui[slice_i * nField + field_i];
        *context->domain->u[field_i] = *context->fi[slice_i * nField + field_i];
    }

    // solve the flow map for time tau_i
    nStep = (int)(context->tau_i[slice_i]/dt);
    // don't really want to call the dns analysis, use custom integrate routine instead?
    //integrate(skewSymmetric, context->domain, context->bman, context->ff, nStep);
    Femlib::ivalue("N_STEP", nStep);
    integrate(skewSymmetric, context->domain, context->bman, context->analyst, context->ff);

    // phase shift in theta (axial direction)
    for(mode_i = 1; mode_i < Geometry::nZ()/2; mode_i++) {
      for(field_i = 0; field_i < nField; field_i++) {
        SEM_to_Fourier(mode_i, context, context->domain->u[field_i], data_f);
        for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
          for(int mode_k = 1; mode_k < nModesX/2; mode_k++) {
            ckt = cos(mode_k*context->theta_i[slice_i]);
            skt = sin(mode_k*context->theta_i[slice_i]);
            rTmp = ckt*data_f[pt_y*nModesX+2*mode_k+0] - skt*data_f[pt_y*nModesX+2*mode_k+1];
            cTmp = skt*data_f[pt_y*nModesX+2*mode_k+0] + ckt*data_f[pt_y*nModesX+2*mode_k+1];
            data_f[pt_y*nModesX+2*mode_k+0] = rTmp;
            data_f[pt_y*nModesX+2*mode_k+1] = cTmp;
          }
        }
        Fourier_to_SEM(mode_i, context, context->domain->u[field_i], data_f);
      }
    }

    // phase shift in phi (azimuthal direction)
    for(mode_i = 1; mode_i < Geometry::nZ()/2; mode_i++) {
      ckt = cos(mode_i*context->phi_i[slice_i]);
      skt = sin(mode_i*context->phi_i[slice_i]);

      for(field_i = 0; field_i < nField; field_i++) {
        data_r = context->domain->u[field_i]->plane(2*mode_i+0);
        data_c = context->domain->u[field_i]->plane(2*mode_i+1);

        for(dof_i = 0; dof_i < NELS_Y*elOrd*nModesX; dof_i++) {
          rTmp = ckt*data_r[dof_i] - skt*data_c[dof_i];
          cTmp = skt*data_r[dof_i] + ckt*data_c[dof_i];
          data_r[dof_i] = rTmp;
          data_c[dof_i] = cTmp;
        }
      }
    }

    // set f
    for(field_i = 0; field_i < nField; field_i++) {
      *context->fi[slice_i * nField + field_i]  = *context->domain->u[field_i];
      *context->fi[slice_i * nField + field_i] *= -1.0;
      *context->fi[slice_i * nField + field_i] += *context->ui[slice_i * nField + field_i];
    }
  }

  RepackX(context, context->fi, context->theta_i, context->phi_i, context->tau_i, f);

  delete[] data_f;

  return 0;
}

void rpo_solve(int nSlice, Mesh* mesh, vector<Element*> elmt, BCmgr* bman, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<Field*> ui, vector<Field*> fi) 
{
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  //BoundarySys* bsys;
  //const NumberSys* nsys;
  real_t dx, er, es;
  real_t *xcoords, *ycoords;
  const real_t* qx;
  int_t pt_x, pt_y, el_x, el_y, el_i;
  const bool guess = true;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  Vec x, f, xl;
  Mat J, P;
  SNES snes;

  context->nSlice  = nSlice;
  context->mesh    = mesh;
  context->elmt    = elmt;
  context->domain  = domain;
  context->bman    = bman;
  context->analyst = analyst;
  context->ff      = FF;
  context->ui      = ui;
  context->fi      = fi;

  //bsys = ui[DOF]->bsys();
  //nsys = bsys->Nsys(0);
  //context->nDofsPlane = nsys->nGlobal() + Geometry::nInode();
  context->nDofsPlane = (((NELS_X*elOrd)/2 + 2)*NELS_Y*elOrd);
  // add dofs for theta and tau for each time slice
  context->nDofsSlice = context->domain->nField() * Geometry::nZ() * context->nDofsPlane + 3;
  context->it = 0;

  context->theta_i = new real_t[nSlice];
  context->phi_i   = new real_t[nSlice];
  context->tau_i   = new real_t[nSlice];
  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    context->theta_i[slice_i] = 0.0;
    context->phi_i[slice_i] = 0.0;
    context->tau_i[slice_i] = Femlib::ivalue("D_T");
  }

  // setup the fourier mapping data
  context->el = new int_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->x  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->y  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->r  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  context->s  = new real_t[NELS_X*elOrd*NELS_Y*elOrd];

  dx = (XMAX - XMIN)/(NELS_X*elOrd);
  Femlib::quadrature(&qx, 0, 0, 0  , elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < NELS_X*elOrd*NELS_Y*elOrd; pt_i++) {
    pt_x = pt_i%(NELS_X*elOrd);
    pt_y = pt_i/(NELS_X*elOrd);
    el_x = pt_i%NELS_X;
    el_y = pt_i/NELS_X;
    el_i = el_y*NELS_X + el_x;
    elmt[el_i]->gCoords(xcoords, ycoords);
    context->x[pt_i] = XMIN + pt_x*dx;
    // element y size increases with distance from the boundary
    context->y[pt_i] = ycoords[pt_y%elOrd];
    
    //for(el_i = 0; el_i < mesh->nEl(); el_i++) {
    // pass er and es by reference?
    if(elmt[el_i]->locate(context->x[pt_i], context->y[pt_i], er, es, &work[0], guess)) {
      context->el[pt_i] = el_i;
      context->r[pt_i] = er;
      context->s[pt_i] = es;
      //break;
    } else {
      cout << "ERROR! element does not contain point" << endl;
      abort();
    }
    //}
  }

  context->localSize = context->domain->nField() * Geometry::nZProc() * context->nDofsPlane;
  // store phase shifts on the 0th processor
  if(!myRank) context->localSize += 3;

  context->localShift = myRank * context->localSize;
  if(!myRank) context->localShift += 3;

  VecCreateMPI(MPI_COMM_WORLD, context->localSize, nSlice*context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, nSlice*context->nDofsSlice, &f);

  // create the local to global scatter object
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  ISCreateStride(MPI_COMM_SELF, context->localSize, 0, 1, &context->isl);
  ISCreateStride(MPI_COMM_WORLD, context->localSize, context->localShift, 1, &context->isg);
  VecScatterCreate(x, context->isg, xl, context->isl, &context->ltog);
  VecDestroy(&xl);

  MatCreate(MPI_COMM_WORLD, &J);
  MatSetType(J, MATMPIAIJ);
  MatSetSizes(J, context->localSize, context->localSize, nSlice*context->nDofsSlice, nSlice*context->nDofsSlice);
  MatMPIAIJSetPreallocation(J, 4*nSlice, PETSC_NULL, 4*nSlice, PETSC_NULL);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, context->localSize, context->localSize, nSlice*context->nDofsSlice, nSlice*context->nDofsSlice);
  MatMPIAIJSetPreallocation(P, 4*nSlice, PETSC_NULL, 4*nSlice, PETSC_NULL);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESSetFunction(snes, f, _snes_function, (void*)context);
  SNESSetJacobian(snes, J, P, _snes_jacobian, (void*)context);
  SNESSetFromOptions(snes);

  RepackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);
  SNESSolve(snes, NULL, x);
  RepackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  VecDestroy(&x);
  VecDestroy(&f);
  MatDestroy(&J);
  MatDestroy(&P);
  VecScatterDestroy(&context->ltog);
  ISDestroy(&context->isl);
  ISDestroy(&context->isg);
  delete[] context->el;
  delete[] context->x;
  delete[] context->y;
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
  int              nSlice = 8;
  vector<Field*>   ui; // Solution fields for velocities, pressure at the i time slices
  vector<Field*>   fi; // Solution fields for flow maps at the i time slices
  vector<Field*>   uTmp;
  char             session_i[100];

  PetscInitialize(&argc, &argv, (char*)0, help);

  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);

  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
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
  rpo_solve(nSlice, mesh, elmt, bman, domain, analyst, FF, ui, fi);

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
