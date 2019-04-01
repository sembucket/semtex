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

#include "rpo_utils.h"
#include "rpo_preconditioner.h"

#define NFIELD 3
#define NELS_X 30
#define NELS_Y 8

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
        // don't do axis for now
        if(pt_y == NELS_Y*elOrd) continue;

        data_els[shift_els+pt_i] = data_log[pt_y*NELS_X*elOrd + pt_x];
      }
    }
  }
}

/*
data_f = data_f[num_nodes_y][num_modes_x],
where num_modes_x is equation to the number of imaginary + the number of real components
*/
void SEM_to_Fourier(int plane_k, Context* context, Field* us, real_t* data_f) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int pt_i;
  Element* elmt;
  real_t* data_r = new real_t[NELS_Y*elOrd*nNodesX];

  for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < nNodesX; pt_x++) {
      pt_i = pt_y*nNodesX + pt_x;
      elmt = context->elmt[context->el[pt_i]];
      data_r[pt_y*nNodesX+pt_x] = us->probe(elmt, context->r[pt_i], context->s[pt_i], plane_k);
    }
  }
  // semtex fft works on strided data, so transpose the plane before applying
  data_transpose(data_r, nNodesX, NELS_Y*elOrd);
  dDFTr(data_r, nNodesX, NELS_Y*elOrd, +1);
  data_transpose(data_r, NELS_Y*elOrd, nNodesX);

  for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < context->nModesX; pt_x++) {
      data_f[pt_y*context->nModesX + pt_x] = data_r[pt_y*nNodesX+pt_x];
    }
  }

  delete[] data_r;
}

/*
data_f = data_f[num_nodes_y][num_modes_x],
where num_modes_x is equation to the number of imaginary + the number of real components
*/
void Fourier_to_SEM(int plane_k, Context* context, Field* us, real_t* data_f) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int pt_r;
  double dx, xr, theta;
  real_t* temp = new real_t[context->nModesX];
  real_t* data_r = new real_t[NELS_Y*elOrd*nNodesX];
  const real_t *qx;

  Femlib::quadrature(&qx, 0, 0, 0  , elOrd+1, GLJ, 0.0, 0.0);

  dx = (context->xmax /*- XMIN*/)/NELS_X;

  for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < context->nModesX; pt_x++) {
      temp[pt_x] = data_f[pt_y*context->nModesX + pt_x];
    }
    // fourier interpolation to GLL grid
    for(int pt_x = 0; pt_x < nNodesX; pt_x++) {
      // coordinate in real space
      xr = /*XMIN +*/ (pt_x/elOrd)*dx + 0.5*(1.0 + qx[pt_x%elOrd])*dx;
      // coordinate in fourier space
      theta = 2.0*M_PI*xr/(context->xmax /*- XMIN*/);

      data_r[pt_y*nNodesX + pt_x] = temp[0];
      // ignore the nyquist frequency (entry [1])
      // all modes are scaled by 2.0, except the mean
      for(int mode_k = 1; mode_k < context->nModesX/2; mode_k++) {
        data_r[pt_y*nNodesX + pt_x] += 2.0*temp[2*mode_k+0]*cos(mode_k*theta);
        data_r[pt_y*nNodesX + pt_x] -= 2.0*temp[2*mode_k+1]*sin(mode_k*theta);
      }
    }
  }

  logical_to_elements(data_r, us->plane(plane_k));

  delete[] temp;
  delete[] data_r;
}

void UnpackX(Context* context, vector<Field*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, ll, slice_i, field_i, index;
  int nNodesX = NELS_X*elOrd;
  Field* field;
  const PetscScalar *xArray;
  real_t* data = (nNodesX > context->nModesX) ? new real_t[NELS_Y*elOrd*nNodesX] : new real_t[NELS_Y*elOrd*context->nModesX];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);

  index = 0;
  VecGetArrayRead(xl, &xArray);
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    for(field_i = 0; field_i < NFIELD; field_i++) {
      field = fields[slice_i * NFIELD + field_i];

      for(kk = 0; kk < Geometry::nZProc(); kk++) {
        for(jj = 0; jj < NELS_Y*elOrd; jj++) {
          if(context->x_fourier) {
            for(ii = 0; ii < context->nModesX; ii++) {
              data[jj*context->nModesX + ii] = xArray[index];
              index++;
            }
          } else {
            for(ii = 0; ii < nNodesX; ii++) {
              data[jj*nNodesX + ii] = xArray[index];
              index++;
            }
          }
        }

        if(context->x_fourier) {
          Fourier_to_SEM(kk, context, field, data);
        } else {
          logical_to_elements(data, field->plane(kk));
        }
      }
    }
    // phase shift data lives on the 0th processors part of the vector
    if(!Geometry::procID()) {
      if(context->x_fourier) {
        theta[slice_i] = xArray[index] / context->x_norm;
        index++;
      }
      phi[slice_i]     = xArray[index] / context->x_norm;
      index++;
      if(!context->travelling_wave) {
        tau[slice_i]   = xArray[index] / context->x_norm;
        index++;
      }
    }
  }
  VecRestoreArrayRead(xl, &xArray);

  VecDestroy(&xl);
  delete[] data;
}

void RepackX(Context* context, vector<Field*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, slice_i, field_i, index;
  int nNodesX = NELS_X*elOrd;
  Field* field;
  PetscScalar *xArray;
  real_t* data = (nNodesX > context->nModesX) ? new real_t[NELS_Y*elOrd*nNodesX] : new real_t[NELS_Y*elOrd*context->nModesX];
  Vec xl;
  double norm_l = 0.0, norm_g;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);
  VecGetArray(xl, &xArray);

  index = 0;
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    for(field_i = 0; field_i < NFIELD; field_i++) {
      field = fields[slice_i * NFIELD + field_i];

      for(kk = 0; kk < Geometry::nZProc(); kk++) {
        if(context->x_fourier) {
          SEM_to_Fourier(kk, context, field, data);
        } else {
          elements_to_logical(field->plane(kk), data);
        }

        // skip over redundant real dofs
        for(jj = 0; jj < NELS_Y*elOrd; jj++) {
          if(context->x_fourier) {
            for(ii = 0; ii < context->nModesX; ii++) {
              xArray[index] = data[jj*context->nModesX + ii];
              norm_l += xArray[index]*xArray[index];
              index++;
            }
          } else {
            for(ii = 0; ii < nNodesX; ii++) {
              xArray[index] = data[jj*nNodesX + ii];
              norm_l += xArray[index]*xArray[index];
              index++;
            }
          }
        }
      }
      MPI_Allreduce(&norm_l, &norm_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if(!Geometry::procID()) cout << field_i << " |x|: " << sqrt(norm_g) << endl;
    }
    // phase shift data lives on the 0th processors part of the vector
    if(!Geometry::procID()) {
      if(context->x_fourier) {
        xArray[index] = theta[slice_i] * context->x_norm;
        index++;
      }
      xArray[index]   = phi[slice_i]   * context->x_norm;
      index++;
      if(!context->travelling_wave) {
        xArray[index] = tau[slice_i] * context->x_norm;
        index++;
      }
    }
  }
  VecRestoreArray(xl, &xArray);

  VecScatterBegin(context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&xl);
  delete[] data;
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

void phase_shift_x(Context* context, double theta, double sign, vector<Field*> fields) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = context->nModesX;
  int mode_i, field_i, pt_y, mode_k;
  double ckt, skt, rTmp, cTmp;
  real_t* data_f = new real_t[NELS_Y*elOrd*nModesX];

  for(mode_i = 0; mode_i < Geometry::nZProc(); mode_i++) {
    for(field_i = 0; field_i < context->domain->nField(); field_i++) {
      SEM_to_Fourier(mode_i, context, fields[field_i], data_f);
      for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
        for(int mode_k = 1; mode_k < nModesX/2; mode_k++) {
          ckt = cos(sign * mode_k * theta);
          skt = sin(sign * mode_k * theta);
          rTmp = ckt*data_f[pt_y*nModesX+2*mode_k+0] - skt*data_f[pt_y*nModesX+2*mode_k+1];
          cTmp = skt*data_f[pt_y*nModesX+2*mode_k+0] + ckt*data_f[pt_y*nModesX+2*mode_k+1];
          data_f[pt_y*nModesX+2*mode_k+0] = rTmp;
          data_f[pt_y*nModesX+2*mode_k+1] = cTmp;
        }
      }
      Fourier_to_SEM(mode_i, context, fields[field_i], data_f);
    }
  }
  delete[] data_f;
}

void phase_shift_z(Context* context, double phi,   double sign, vector<Field*> fields) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
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

    for(field_i = 0; field_i < context->domain->nField(); field_i++) {
      data_r = fields[field_i]->plane(2*mode_i+0);
      data_c = fields[field_i]->plane(2*mode_i+1);

      //for(dof_i = 0; dof_i < NELS_Y*elOrd*nModesX; dof_i++) {
      for(dof_i = 0; dof_i < NELS_X *NELS_Y*(elOrd+1)*(elOrd+1); dof_i++) {
        rTmp = ckt*data_r[dof_i] - skt*data_c[dof_i];
        cTmp = skt*data_r[dof_i] + ckt*data_c[dof_i];
        data_r[dof_i] = rTmp;
        data_c[dof_i] = cTmp;
      }
    }
  }
}
