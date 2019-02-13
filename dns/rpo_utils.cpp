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

#define X_FOURIER
//#define NFIELD 4
#define NFIELD 3
#define XMIN 0.0
#define XMAX (2.0*M_PI)
#define YMIN 0.0
#define YMAX 1.0
#define NELS_X 30
#define NELS_Y 7
#define NSLICE 16
#define NSTEPS 40

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
void SEM_to_Fourier(int plane_k, Context* context, Field* us, real_t* data_f, int nModesX) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  //int nModesX = nNodesX/2 + 2;
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
    for(int pt_x = 0; pt_x < nModesX; pt_x++) {
      data_f[pt_y*nModesX + pt_x] = data_r[pt_y*nNodesX+pt_x];
    }
  }

  delete[] data_r;
}

/*
data_f = data_f[num_nodes_y][num_modes_x],
where num_modes_x is equation to the number of imaginary + the number of real components
*/
void Fourier_to_SEM(int plane_k, Context* context, Field* us, real_t* data_f, int nModesX) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  //int nModesX = nNodesX/2 + 2;
  int pt_r;
  double dx, xr, theta;
  real_t* temp = new real_t[nModesX];
  real_t* data_r = (nNodesX > nModesX) ? new real_t[NELS_Y*elOrd*nNodesX] : new real_t[NELS_Y*elOrd*nModesX];
  const real_t *qx;

  Femlib::quadrature(&qx, 0, 0, 0  , elOrd+1, GLJ, 0.0, 0.0);

  dx = (XMAX - XMIN)/NELS_X;

  for(int pt_y = 0; pt_y < NELS_Y*elOrd; pt_y++) {
    for(int pt_x = 0; pt_x < nModesX; pt_x++) {
      temp[pt_x] = data_f[pt_y*nModesX + pt_x];
    }
    // fourier interpolation to GLL grid
    for(int pt_x = 0; pt_x < nNodesX; pt_x++) {
      // coordinate in real space
      xr = XMIN + (pt_x/elOrd)*dx + 0.5*(1.0 + qx[pt_x%elOrd])*dx;
      // coordinate in fourier space
      theta = 2.0*M_PI*xr/(XMAX - XMIN);

      data_r[pt_y*nNodesX + pt_x] = temp[0];
      // ignore the nyquist frequency (entry [1])
      // all modes are scaled by 2.0, except the mean
      for(int mode_k = 1; mode_k < nModesX/2; mode_k++) {
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
  int ii, jj, kk, slice_i, field_i, index;
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2;// + 2;
  Field* field;
  const PetscScalar *xArray;
  real_t* data = (nNodesX > nModesX) ? new real_t[NELS_Y*elOrd*nNodesX] : new real_t[NELS_Y*elOrd*nModesX];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->ltog, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->ltog, x, xl, INSERT_VALUES, SCATTER_FORWARD);

  index = 0;
  VecGetArrayRead(xl, &xArray);
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    for(field_i = 0; field_i < NFIELD; field_i++) {
      field = fields[slice_i * NFIELD + field_i];

      for(kk = 0; kk < Geometry::nZProc(); kk++) {
        // skip over redundant real dofs
        for(jj = 0; jj < NELS_Y*elOrd; jj++) {
#ifdef X_FOURIER
          for(ii = 0; ii < nModesX; ii++) {
            data[jj*nModesX + ii] = xArray[index++];
#else
          for(ii = 0; ii < nNodesX; ii++) {
            data[jj*nNodesX + ii] = xArray[index++];
#endif
          }
        }

#ifdef X_FOURIER
        Fourier_to_SEM(kk, context, field, data, nModesX);
#else
        logical_to_elements(data, field->plane(kk));
#endif
      }
    }
    // phase shift data lives on the 0th processors part of the vector
    //if(!Geometry::procID()) {
    if(Geometry::procID() == slice_i) {
      theta[slice_i] = xArray[index++];
#ifdef X_FOURIER
      phi[slice_i]   = xArray[index++];
#endif
      tau[slice_i]   = xArray[index++];
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
  int nModesX = nNodesX/2;// + 2;
  Field* field;
  PetscScalar *xArray;
  real_t* data = (nNodesX > nModesX) ? new real_t[NELS_Y*elOrd*nNodesX] : new real_t[NELS_Y*elOrd*nModesX];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecGetArray(xl, &xArray);

  index = 0;
  for(slice_i = 0; slice_i < context->nSlice; slice_i++) {
    for(field_i = 0; field_i < NFIELD; field_i++) {
      field = fields[slice_i * NFIELD + field_i];

      for(kk = 0; kk < Geometry::nZProc(); kk++) {
#ifdef X_FOURIER
        SEM_to_Fourier(kk, context, field, data, nModesX);
#else
        elements_to_logical(field->plane(kk), data);
#endif

        // skip over redundant real dofs
        for(jj = 0; jj < NELS_Y*elOrd; jj++) {
#ifdef X_FOURIER
          for(ii = 0; ii < nModesX; ii++) {
            xArray[index] = data[jj*nModesX + ii];
#else
          for(ii = 0; ii < nNodesX; ii++) {
            xArray[index] = data[jj*nNodesX + ii];
#endif
            index++;
          }
        }
      }
    }
    // phase shift data lives on the 0th processors part of the vector
    //if(!Geometry::procID()) {
    if(Geometry::procID() == slice_i) {
      xArray[index++] = theta[slice_i];
#ifdef X_FOURIER
      xArray[index++] = phi[slice_i];
#endif
      xArray[index++] = tau[slice_i];
    }
  }
  VecRestoreArray(xl, &xArray);

  VecScatterBegin(context->ltog, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->ltog, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&xl);
  delete[] data;
}
