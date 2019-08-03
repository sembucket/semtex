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
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

#include <fftw3.h>

#include "rpo_utils.h"
#include "rpo_preconditioner.h"

#define THREE 3
//#define RM_2FOLD_SYM

void velocity_scales(Context* context) {
  double Ku[3], Ku_bar;
  double fac[3] = {0.5, 0.2, 0.2};
  vector<AuxField*> vel;
  int elOrd = Geometry::nP() - 1;
  int nModesX = Femlib::ivalue("NELS_X")*elOrd;
  real_t* data_r = new real_t[context->nDofsPlane];
  real_t* data_i = new real_t[context->nDofsPlane];

  vel.resize(3);
  for(int field_i = 0; field_i < 3; field_i++) {
    vel[field_i] = context->u0[field_i];
  }

  // integrate the energy in each component separately
  for(int field_i = 0; field_i < 3; field_i++) {
    *vel[0] = 0.0;
    *vel[1] = 0.0;
    *vel[2] = 0.0;
    *vel[field_i] = *context->ui[field_i];
    vel[field_i]->transform(INVERSE);
    context->fi[0]->innerProduct(vel, vel);
    vel[field_i]->transform(FORWARD);
    context->fi[0]->transform(FORWARD);
    *context->fi[0] *= 0.5;
    Ku[field_i] = 2.0 * M_PI * context->fi[0]->integral(0);
    MPI_Bcast(&Ku[field_i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(!Geometry::procID()) cout << field_i << ": ke: " << Ku[field_i] << endl;
  }

  // now integrate the mean flow
  *context->u0[0] = *context->ui[0];
  *context->u0[1] = 0.0;
  *context->u0[2] = 0.0;
  SEM_to_Fourier(0, context, context->ui[0], data_r, data_i);
  for(int dof_i = 0; dof_i < context->nDofsPlane; dof_i++) {
    if(Geometry::procID() || dof_i%nModesX != 0) {
      data_r[dof_i]  = 0.0;
    } 
    data_i[dof_i] = 0.0;
  }
  Fourier_to_SEM(0, context, context->u0[0], data_r, data_i, 0);
  vel[0]->transform(INVERSE);
  context->fi[0]->innerProduct(vel, vel);
  vel[0]->transform(FORWARD);
  context->fi[0]->transform(FORWARD);
  *context->fi[0] *= 0.5;
  Ku_bar = 2.0 * M_PI * context->fi[0]->integral(0);
  MPI_Bcast(&Ku_bar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(!Geometry::procID()) cout << "Ku_bar:   " << Ku_bar << endl;

  Ku[0] -= Ku_bar;
  for(int field_i = 0; field_i < THREE; field_i++) {
    context->u_scale[field_i] = sqrt(fac[field_i] * Ku_bar / Ku[field_i]);
    if(!Geometry::procID()) cout << field_i << " velocity scale: " << context->u_scale[field_i] << endl;
  }

  delete[] data_r;
  delete[] data_i;
}

void data_transpose(real_t* data, int nx, int ny) {
  real_t* temp = new real_t[nx*ny];

  for(int iy = 0; iy < ny; iy++) {
    for(int ix = 0; ix < nx; ix++) {
      temp[ix*ny + iy] = data[iy*nx + ix];
    }
  }
  for(int ii = 0; ii < nx*ny; ii++) {
    data[ii] = temp[ii];
  }

  delete[] temp;
}

void elements_to_logical(int nex, int ney, real_t* data_els, real_t* data_log) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int pt_x, pt_y;
  int index = -1;

  for(int el_y = 0; el_y < ney; el_y++) {
    for(int el_x = 0; el_x < nex; el_x++) {
      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        index++;

        // skip over right and top edges for each element, as these are redundant
        if(pt_i%(elOrd+1) == elOrd || pt_i/(elOrd+1) == elOrd) continue;

        pt_x = el_x*elOrd + pt_i%(elOrd+1);
        pt_y = el_y*elOrd + pt_i/(elOrd+1);

        data_log[pt_y*nex*elOrd + pt_x] = data_els[index];
      }
    }
  }
}

void logical_to_elements(int nex, int ney, real_t* data_log, real_t* data_els, int plane_i, int field_i) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int shift_els, pt_r, pt_s, pt_x, pt_y;
  
  for(int el_y = 0; el_y < ney; el_y++) {
    for(int el_x = 0; el_x < nex; el_x++) {
      shift_els = (el_y*nex + el_x)*nodes_per_el;

      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        pt_r = pt_i%(elOrd+1);
        pt_s = pt_i/(elOrd+1);

        pt_x = el_x*elOrd + pt_r;
        pt_y = el_y*elOrd + pt_s;
        // assume periodic in x
        if(pt_x == nex*elOrd) pt_x = 0;
        // assume dirichlet bcs on the outer boundary
        if(pt_y == ney*elOrd) continue;
        // homogeneous dirichlet bcs along the axis
        //if(pt_y == 0 &&                 plane_i >= 4) continue;
        //if(pt_y == 0 && field_i == 0 && plane_i >= 2) continue;
        //if(pt_y == 0 && field_i >= 1 && plane_i <= 1) continue;

        data_els[shift_els+pt_i] = data_log[pt_y*nex*elOrd + pt_x];
      }
    }
  }
}

void SEM_to_Fourier(int plane_k, Context* context, AuxField* us, real_t* data_r, real_t* data_i) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nModes = context->nModesX;
  int pt_i;
  Element* elmt;

  for(int pt_y = 0; pt_y < ney*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < nModes; pt_x++) {
      pt_i = pt_y*nModes + pt_x;
      elmt = context->elmt[context->el[pt_i]];

      context->data_s[pt_x][0] = us->probe(elmt, context->r[pt_i], context->s[pt_i], plane_k+0);
      context->data_s[pt_x][1] = us->probe(elmt, context->r[pt_i], context->s[pt_i], plane_k+1);

      context->data_f[pt_x][0] = context->data_f[pt_x][1] = 0.0;
    }

    fftw_execute(context->trans_fwd);

    for(int pt_x = 0; pt_x < nModes; pt_x++) {
      data_r[pt_y*nModes+pt_x] = context->data_f[pt_x][0];
      data_i[pt_y*nModes+pt_x] = context->data_f[pt_x][1];
    }
  }
}

void Fourier_to_SEM(int plane_k, Context* context, AuxField* us, real_t* data_r, real_t* data_i, int field_i) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nModes = context->nModesX;
  int mode_l;
  const real_t *qx;
  double dx = (context->xmax)/nex; // assume constant
  double xr, theta;
  real_t* temp_r = new real_t[ney*elOrd*nModes];
  real_t* temp_i = new real_t[ney*elOrd*nModes];

  Femlib::quadrature(&qx, 0, 0, 0  , elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_y = 0; pt_y < ney*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < nModes; pt_x++) {
      // coordinate in real space
      xr = (pt_x/elOrd)*dx + 0.5*(1.0 + qx[pt_x%elOrd])*dx; // + XMIN
      theta = 2.0*M_PI*xr/(context->xmax /*-XMIN*/);

      temp_r[pt_y*nModes + pt_x] = data_r[pt_y*nModes];
      temp_i[pt_y*nModes + pt_x] = 0.0; // nyquist frequency
      for(int mode_k = 1; mode_k < nModes; mode_k++) {
        mode_l = (mode_k <= nModes/2) ? mode_k : mode_k - nModes; // fftw ordering of complex data

        temp_r[pt_y*nModes + pt_x] += (data_r[pt_y*nModes + mode_k]*cos(mode_l*theta) - data_i[pt_y*nModes + mode_k]*sin(mode_l*theta));
        temp_i[pt_y*nModes + pt_x] += (data_r[pt_y*nModes + mode_k]*sin(mode_l*theta) + data_i[pt_y*nModes + mode_k]*cos(mode_l*theta));
      }
      temp_r[pt_y*nModes + pt_x] /= nModes; // rescale!!
      temp_i[pt_y*nModes + pt_x] /= nModes; // rescale!!
    }
  }

  logical_to_elements(nex, ney, temp_r, us->plane(plane_k+0), Geometry::procID()*Geometry::nZProc()+plane_k+0, field_i);
  logical_to_elements(nex, ney, temp_i, us->plane(plane_k+1), Geometry::procID()*Geometry::nZProc()+plane_k+1, field_i);

  delete[] temp_r;
  delete[] temp_i;
}

int LocalIndex(Context* context, int field_i, int plane_i, int point_x, int point_y) {
  int shift_1 = (plane_i/2) * 3 * 2 * context->nDofsPlane;
  int shift_2 =  point_y    * 3 * 2 * context->nModesX;
  int shift_3 =  point_x    * 3 * 2;

  return shift_1 + shift_2 + shift_3 + 2*field_i + plane_i%2;

  //int nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  //return field_i * nDofsCube_l + (plane_i+0) * context->nDofsPlane + point_y * context->nModesX + point_x;
}

void UnpackX(Context* context, vector<AuxField*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int mode_i, mode_l, index;
  real_t* data_r;
  real_t* data_i;
  AuxField* field;
  const PetscScalar *xArray;
  Vec xl;
  double scale, rh, dr;
  int vecSize;

  data_r = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];
  data_i = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);

  VecScatterBegin(context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);

#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==0) {
#endif

  VecGetArrayRead(xl, &xArray);

  for(int field_i = 0; field_i < THREE; field_i++) {
    field = fields[field_i];

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      for(int point_y = 0; point_y < ney*elOrd; point_y++) {
        for(int point_x = 0; point_x < context->nModesX; point_x++) {
          mode_i = (Geometry::procID()*Geometry::nZProc() + plane_i)/2;
          mode_l = (point_x <= context->nModesX/2) ? point_x : point_x - context->nModesX; // fftw ordering of complex data
          scale = 4.0/(4.0 + (2.0*M_PI/context->xmax)*fabs(mode_l) + fabs(mode_i));
          rh = 0.5*(context->rad_coords[point_y+1] + context->rad_coords[point_y]);
          dr = context->rad_weights[point_y];
          scale *= sqrt(2.0 * M_PI * rh * dr);
          scale *= context->u_scale[field_i];

          index = LocalIndex(context, field_i, plane_i+0, point_x, point_y);
          data_r[point_y*context->nModesX+point_x] = xArray[index] / scale;

          index = LocalIndex(context, field_i, plane_i+1, point_x, point_y);
          data_i[point_y*context->nModesX+point_x] = xArray[index] / scale;
        }
      }

      Fourier_to_SEM(plane_i, context, field, data_r, data_i, field_i);
    }
  }

  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    index = THREE * nDofsCube_l;

    theta[0] = xArray[index++];
    phi[0]   = xArray[index++];
    if(!context->travelling_wave) {
      tau[0] = xArray[index++];
    }
  }

  VecRestoreArrayRead(xl, &xArray);

#ifdef RM_2FOLD_SYM
  }
#endif

  VecDestroy(&xl);
  delete[] data_r;
  delete[] data_i;
}

void RepackX(Context* context, vector<AuxField*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int mode_i, mode_l, index;
  real_t* data_r;
  real_t* data_i;
  AuxField* field;
  PetscScalar *xArray;
  Vec xl;
  double scale, dr, rh;

  data_r = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];
  data_i = (nNodesX > context->nModesX) ? new real_t[ney*elOrd*nNodesX] : new real_t[ney*elOrd*context->nModesX];

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);

#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==0) {
#endif

  VecGetArray(xl, &xArray);
  for(int field_i = 0; field_i < THREE; field_i++) {
    field = fields[field_i];

    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      SEM_to_Fourier(plane_i, context, field, data_r, data_i);

      for(int point_y = 0; point_y < ney*elOrd; point_y++) {
        for(int point_x = 0; point_x < context->nModesX; point_x++) {
          mode_i = (Geometry::procID()*Geometry::nZProc() + plane_i)/2;
          mode_l = (point_x <= context->nModesX/2) ? point_x : point_x - context->nModesX; // fftw ordering of complex data
          scale = 4.0/(4.0 + (2.0*M_PI/context->xmax)*fabs(mode_l) + fabs(mode_i));
          rh = 0.5*(context->rad_coords[point_y+1] + context->rad_coords[point_y]);
          dr = context->rad_weights[point_y];
          scale *= sqrt(2.0 * M_PI * rh * dr);
          scale *= context->u_scale[field_i];

          index = LocalIndex(context, field_i, plane_i+0, point_x, point_y);
          xArray[index] = data_r[point_y*context->nModesX+point_x] * scale;

          index = LocalIndex(context, field_i, plane_i+1, point_x, point_y);
          xArray[index] = data_i[point_y*context->nModesX+point_x] * scale;
        }
      }
    }
  }

  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    index = THREE * nDofsCube_l;

    xArray[index++]   = theta[0];
    xArray[index++]   = phi[0];
    if(!context->travelling_wave) {
      xArray[index++] = tau[0];
    }
  }

  VecRestoreArray(xl, &xArray);

#ifdef RM_2FOLD_SYM
  }
#endif

  VecScatterBegin(context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&xl);
  delete[] data_r;
  delete[] data_i;
}

// define a vec scatter object for mapping from parallel global vectors to semtex data
void assign_scatter_semtex(Context* context) {
  int   nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int   nDofsCube_g = Geometry::nZ()     * context->nDofsPlane;
  int   nShifts     = 1;
  int   ind_i       = 0;
  int*  inds;
  IS    isl, isg;
  Vec   vl, vg;

  if(context->x_fourier)        nShifts++;
  if(!context->travelling_wave) nShifts++;

#ifndef RM_2FOLD_SYM
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

    // assign the phase shifts from the 0th processor
    if(!Geometry::procID()) {
      for(int ind_j = 0; ind_j < nShifts; ind_j++) {
        inds[ind_i] = slice_i * context->nDofsSlice + context->nField * nDofsCube_g + ind_j;
        ind_i++;
      }
    }
  }
#endif

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nSlice * context->nDofsSlice, &vg);

  ISCreateStride(MPI_COMM_SELF, context->localSize, 0, 1, &isl);
#ifdef RM_2FOLD_SYM
  ISCreateStride(MPI_COMM_WORLD, context->localSize, 0, 1, &isg);
#else
  ISCreateGeneral(MPI_COMM_WORLD, context->localSize, inds, PETSC_COPY_VALUES, &isg);
#endif

  VecScatterCreate(vg, isg, vl, isl, &context->global_to_semtex);
  
  VecDestroy(&vl);
  VecDestroy(&vg);
  ISDestroy(&isl);
  ISDestroy(&isg);
#ifndef RM_2FOLD_SYM
  delete[] inds;
#endif
}

void phase_shift_x(Context* context, double theta, double sign, vector<Field*> fields) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nModesX = context->nModesX;
  int mode_i, field_i, pt_y, mode_k, mode_l;
  double ckt, skt, rTmp, cTmp;
  real_t *data_r, *data_i;

  data_r = new real_t[ney*elOrd*nModesX];
  data_i = new real_t[ney*elOrd*nModesX];

  for(mode_i = 0; mode_i < Geometry::nZProc(); mode_i += 2) {
    for(field_i = 0; field_i < THREE; field_i++) {
      SEM_to_Fourier(mode_i, context, fields[field_i], data_r, data_i);

      for(int pt_y = 0; pt_y < ney*elOrd; pt_y++) {
        for(int mode_k = 0; mode_k < nModesX; mode_k++) {
          mode_l = (mode_k <= nModesX/2) ? mode_k : mode_k - nModesX; // fftw ordering of complex data

          ckt = cos(sign * mode_l * theta);
          skt = sin(sign * mode_l * theta);

          rTmp = +ckt*data_r[pt_y*nModesX+mode_k] + skt*data_i[pt_y*nModesX+mode_k];
          cTmp = -skt*data_r[pt_y*nModesX+mode_k] + ckt*data_i[pt_y*nModesX+mode_k];
          data_r[pt_y*nModesX+mode_k] = rTmp;
          data_i[pt_y*nModesX+mode_k] = cTmp;
        }
      }
      Fourier_to_SEM(mode_i, context, fields[field_i], data_r, data_i, field_i);
    }
  }
  delete[] data_r;
  delete[] data_i;
}

void phase_shift_z(Context* context, double phi,   double sign, vector<Field*> fields) {
  int nex = context->nElsX;
  int ney = context->nElsY;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = nex*elOrd;
  int nModesX = context->nModesX;
  int mode_i, mode_j, field_i, dof_i;
  double ckt, skt, rTmp, cTmp;
  register real_t* data_r;
  register real_t* data_c;

  for(mode_i = 0; mode_i < Geometry::nZProc()/2; mode_i++) {
    mode_j = Geometry::procID() * (Geometry::nZProc()/2) + mode_i;
    if(!mode_j) continue;

    ckt = cos(sign * mode_j * phi);
    skt = sin(sign * mode_j * phi);

    for(field_i = 0; field_i < THREE; field_i++) {
      data_r = fields[field_i]->plane(2*mode_i+0);
      data_c = fields[field_i]->plane(2*mode_i+1);

      for(dof_i = 0; dof_i < nex*ney*(elOrd+1)*(elOrd+1); dof_i++) {
        rTmp = +ckt*data_r[dof_i] + skt*data_c[dof_i];
        cTmp = -skt*data_r[dof_i] + ckt*data_c[dof_i];
        data_r[dof_i] = rTmp;
        data_c[dof_i] = cTmp;
      }
    }
  }
}
