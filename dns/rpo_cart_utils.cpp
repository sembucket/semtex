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

#include "rpo_cart_utils.h"

#define NELS_X 12
#define NELS_Y 9

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

void UnpackX(Context* context, vector<Field*> fields, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, ll, field_i, index;
  int nNodesX = NELS_X*elOrd;
  Field* field;
  const PetscScalar *xArray;
  real_t* data = new real_t[NELS_Y*elOrd*nNodesX];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);

  index = 0;
  VecGetArrayRead(xl, &xArray);
  for(field_i = 0; field_i < 3; field_i++) {
    field = fields[field_i];

    for(kk = 0; kk < Geometry::nZProc(); kk++) {
      // skip over redundant real dofs
      for(jj = 0; jj < NELS_Y*elOrd; jj++) {
        for(ii = 0; ii < nNodesX; ii++) {
          data[jj*nNodesX + ii] = xArray[index++];
        }
      }

      logical_to_elements(data, field->plane(kk));
    }
  }
  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    *phi = xArray[index++] / context->x_norm;
    *tau = xArray[index++] / context->x_norm;
  }
  VecRestoreArrayRead(xl, &xArray);

  VecDestroy(&xl);
  delete[] data;
}

void RepackX(Context* context, vector<Field*> fields, real_t  phi, real_t  tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, field_i, index;
  int nNodesX = NELS_X*elOrd;
  Field* field;
  PetscScalar *xArray;
  real_t* data = new real_t[NELS_Y*elOrd*nNodesX];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);
  VecGetArray(xl, &xArray);

  index = 0;
  for(field_i = 0; field_i < 3; field_i++) {
    field = fields[field_i];

    for(kk = 0; kk < Geometry::nZProc(); kk++) {
      elements_to_logical(field->plane(kk), data);

      // skip over redundant real dofs
      for(jj = 0; jj < NELS_Y*elOrd; jj++) {
        for(ii = 0; ii < nNodesX; ii++) {
          xArray[index] = data[jj*nNodesX + ii];
          index++;
        }
      }
    }
  }
  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    xArray[index++] = phi * context->x_norm;
    xArray[index++] = tau * context->x_norm;
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
  int   nShifts     = 2;
  int   ind_i       = 0;
  int*  inds        = new int[context->localSize];
  IS    isl, isg;
  Vec   vl, vg;

  context->lShift = new int[3];

  for(int field_i = 0; field_i < 3; field_i++) {
    context->lShift[field_i] = field_i * nDofsCube_g + Geometry::procID() * nDofsCube_l;

    for(int ind_j = 0; ind_j < nDofsCube_l; ind_j++) {
      inds[ind_i] = context->lShift[field_i] + ind_j;
      ind_i++;
    }
  }

  // assign the phase shifts from the 0th processor
  if(!Geometry::procID()) {
    for(int ind_j = 0; ind_j < nShifts; ind_j++) {
      inds[ind_i] = 3 * nDofsCube_g + ind_j;
      ind_i++;
    }
  }

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &vg);

  ISCreateStride(MPI_COMM_SELF, context->localSize, 0, 1, &isl);
  ISCreateGeneral(MPI_COMM_WORLD, context->localSize, inds, PETSC_COPY_VALUES, &isg);

  VecScatterCreate(vg, isg, vl, isl, &context->global_to_semtex);
  
  VecDestroy(&vl);
  VecDestroy(&vg);
  ISDestroy(&isl);
  ISDestroy(&isg);
  delete[] inds;
}

void phase_shift_z(Context* context, double phi, double sign, vector<Field*> fields) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
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

      for(dof_i = 0; dof_i < NELS_X*NELS_Y*(elOrd+1)*(elOrd+1); dof_i++) {
        rTmp = ckt*data_r[dof_i] - skt*data_c[dof_i];
        cTmp = skt*data_r[dof_i] + ckt*data_c[dof_i];
        data_r[dof_i] = rTmp;
        data_c[dof_i] = cTmp;
      }
    }
  }
}
