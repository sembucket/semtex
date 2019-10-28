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

#define VEL_MAJOR

void elements_to_logical(real_t* data_els, real_t* data_log) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int pt_x, pt_y;
  int index = -1;

  for(int el_y = 0; el_y < Femlib::ivalue("NELS_Y"); el_y++) {
    for(int el_x = 0; el_x < Femlib::ivalue("NELS_X"); el_x++) {
      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        index++;

        // skip over right and top edges for each element
        if(pt_i%(elOrd+1) == elOrd || pt_i/(elOrd+1) == elOrd) continue;

        pt_x = el_x*elOrd + pt_i%(elOrd+1);
        pt_y = el_y*elOrd + pt_i/(elOrd+1);

        data_log[pt_y*Femlib::ivalue("NELS_X")*elOrd + pt_x] = data_els[index];
      }
    }
  }
}

void logical_to_elements(real_t* data_log, real_t* data_els) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int shift_els, pt_r, pt_s, pt_x, pt_y;
  
  for(int el_y = 0; el_y < Femlib::ivalue("NELS_Y"); el_y++) {
    for(int el_x = 0; el_x < Femlib::ivalue("NELS_X"); el_x++) {
      shift_els = (el_y*Femlib::ivalue("NELS_X") + el_x)*nodes_per_el;

      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        pt_r = pt_i%(elOrd+1);
        pt_s = pt_i/(elOrd+1);

        pt_x = el_x*elOrd + pt_r;
        pt_y = el_y*elOrd + pt_s;
        // omit boundaries (bottom left)
        if(pt_x == 0 || pt_y == 0) continue;

        data_els[shift_els+pt_i] = data_log[pt_y*Femlib::ivalue("NELS_X")*elOrd + pt_x];
      }
    }
  }
}

#ifdef VEL_MAJOR
void UnpackX(Context* context, vector<Field*> fields, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, ll, field_i, index;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  Field* field;
  const PetscScalar *xArray;
  real_t* data = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  Vec xl;
  double scale;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);

  index = 0;
  VecGetArrayRead(xl, &xArray);
  for(field_i = 0; field_i < 3; field_i++) {
    field = fields[field_i];

    for(kk = 0; kk < Geometry::nZProc(); kk++) {
      ll = ( Geometry::procID() * Geometry::nZProc() + kk ) / 2;
      for(jj = 0; jj < Femlib::ivalue("NELS_Y")*elOrd; jj++) {
        for(ii = 0; ii < nNodesX; ii++) {
          // omit boundaries
          if(ii == 0 || jj == 0) continue;

          scale = 2.0/(2.0 + fabs(ll)) * context->coord_weights[0][jj*nNodesX + ii];

          data[jj*nNodesX + ii] = xArray[index++] / scale;
        }
      }

      logical_to_elements(data, field->plane(kk));
    }
  }
  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    *phi = xArray[index++];
    if(!context->travelling_wave) *tau = xArray[index++];
  }
  VecRestoreArrayRead(xl, &xArray);

  VecDestroy(&xl);
  delete[] data;
}
#else
void UnpackX(Context* context, vector<Field*> fields, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, ll, index;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  const PetscScalar *xArray;
  real_t* data_u = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  real_t* data_v = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  real_t* data_w = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);

  index = 0;
  VecGetArrayRead(xl, &xArray);
  for(kk = 0; kk < Geometry::nZProc(); kk++) {
    ll = ( Geometry::procID() * Geometry::nZProc() + kk ) / 2;
    // skip over redundant real dofs
    for(jj = 0; jj < Femlib::ivalue("NELS_Y")*elOrd; jj++) {
      for(ii = 0; ii < nNodesX; ii++) {
        data_u[jj*nNodesX + ii] = xArray[index++];
        data_v[jj*nNodesX + ii] = xArray[index++];
        data_w[jj*nNodesX + ii] = xArray[index++];
      }
    }

    logical_to_elements(data_u, fields[0]->plane(kk));
    logical_to_elements(data_v, fields[1]->plane(kk));
    logical_to_elements(data_w, fields[2]->plane(kk));
  }
  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    *phi = xArray[index++];
    if(!context->travelling_wave) *tau = xArray[index++];
  }
  VecRestoreArrayRead(xl, &xArray);

  VecDestroy(&xl);
  delete[] data_u;
  delete[] data_v;
  delete[] data_w;
}
#endif

#ifdef VEL_MAJOR
void RepackX(Context* context, vector<Field*> fields, real_t  phi, real_t  tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, ll, field_i, index;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  Field* field;
  PetscScalar *xArray;
  real_t* data = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  Vec xl;
  double norm_l = 0.0, norm_g;
  double scale;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);
  VecGetArray(xl, &xArray);

  index = 0;
  for(field_i = 0; field_i < 3; field_i++) {
    field = fields[field_i];

    for(kk = 0; kk < Geometry::nZProc(); kk++) {
      ll = ( Geometry::procID() * Geometry::nZProc() + kk ) / 2;

      elements_to_logical(field->plane(kk), data);

      for(jj = 0; jj < Femlib::ivalue("NELS_Y")*elOrd; jj++) {
        for(ii = 0; ii < nNodesX; ii++) {
          // omit boundaries
          if(ii == 0 || jj == 0) continue;

          scale = 2.0/(2.0 + fabs(ll)) * context->coord_weights[0][jj*nNodesX + ii];

          xArray[index] = data[jj*nNodesX + ii] * scale;
          norm_l += xArray[index]*xArray[index];
          index++;
        }
      }
    }
    MPI_Allreduce(&norm_l, &norm_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(!Geometry::procID()) cout << field_i << " |x|: " << sqrt(norm_g) << endl;
  }
  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    xArray[index++] = phi;
    if(!context->travelling_wave) xArray[index++] = tau;
  }
  VecRestoreArray(xl, &xArray);

  VecScatterBegin(context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&xl);
  delete[] data;
}
#else
void RepackX(Context* context, vector<Field*> fields, real_t  phi, real_t  tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, ll, index;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  PetscScalar *xArray;
  real_t* data_u = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  real_t* data_v = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  real_t* data_w = new real_t[Femlib::ivalue("NELS_Y")*elOrd*nNodesX];
  Vec xl;
  double norm_l_u = 0.0, norm_g_u;
  double norm_l_v = 0.0, norm_g_v;
  double norm_l_w = 0.0, norm_g_w;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);
  VecGetArray(xl, &xArray);

  index = 0;

  for(kk = 0; kk < Geometry::nZProc(); kk++) {
    ll = ( Geometry::procID() * Geometry::nZProc() + kk ) / 2;

    elements_to_logical(fields[0]->plane(kk), data_u);
    elements_to_logical(fields[1]->plane(kk), data_v);
    elements_to_logical(fields[2]->plane(kk), data_w);

    // skip over redundant real dofs
    for(jj = 0; jj < Femlib::ivalue("NELS_Y")*elOrd; jj++) {
      for(ii = 0; ii < nNodesX; ii++) {
        xArray[index++] = data_u[jj*nNodesX + ii];
        xArray[index++] = data_v[jj*nNodesX + ii];
        xArray[index++] = data_w[jj*nNodesX + ii];
        norm_l_u += data_u[jj*nNodesX + ii] * data_u[jj*nNodesX + ii];
        norm_l_v += data_v[jj*nNodesX + ii] * data_v[jj*nNodesX + ii];
        norm_l_w += data_w[jj*nNodesX + ii] * data_w[jj*nNodesX + ii];
      }
    }
  }
  MPI_Allreduce(&norm_l_u, &norm_g_u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&norm_l_v, &norm_g_v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&norm_l_w, &norm_g_w, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(!Geometry::procID()) cout << " |u|: " << sqrt(norm_g_u) << endl;
  if(!Geometry::procID()) cout << " |v|: " << sqrt(norm_g_v) << endl;
  if(!Geometry::procID()) cout << " |w|: " << sqrt(norm_g_w) << endl;

  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    xArray[index++] = phi;
    if(!context->travelling_wave) xArray[index++] = tau;
  }
  VecRestoreArray(xl, &xArray);

  VecScatterBegin(context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&xl);
  delete[] data_u;
  delete[] data_v;
  delete[] data_w;
}
#endif

// define a vec scatter object for mapping from parallel global vectors to semtex data
#ifdef VEL_MAJOR
void assign_scatter_semtex(Context* context) {
  int   nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int   nDofsCube_g = Geometry::nZ()     * context->nDofsPlane;
  int   nShifts     = (!context->travelling_wave) ? 2 : 1;
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
#else
void assign_scatter_semtex(Context* context) {
  int   nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int   nDofsCube_g = Geometry::nZ()     * context->nDofsPlane;
  int   nShifts     = (!context->travelling_wave) ? 2 : 1;
  int   ind_i       = 0;
  int   start       = Geometry::procID() * nDofsCube_l * context->nField;
  int*  inds        = new int[context->localSize];
  IS    isl, isg;
  Vec   vl, vg;

  context->lShift = new int[3];

  for(int ind_j = 0; ind_j < context->nField * nDofsCube_l; ind_j++) {
    inds[ind_i++] = start + ind_j;
  }

  // assign the phase shifts from the 0th processor
  if(!Geometry::procID()) {
    for(int ind_j = 0; ind_j < nShifts; ind_j++) {
      inds[ind_i++] = 3 * nDofsCube_g + ind_j;
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
#endif

void phase_shift_z(Context* context, double phi, double sign, vector<Field*> fields) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  int mode_i, mode_j, field_i, dof_i;
  double ckt, skt, rTmp, cTmp;
  double* data_r = new double[nNodesX*Femlib::ivalue("NELS_Y")*elOrd];
  double* data_c = new double[nNodesX*Femlib::ivalue("NELS_Y")*elOrd];

  for(mode_i = 0; mode_i < Geometry::nZProc()/2; mode_i++) {
    mode_j = Geometry::procID() * (Geometry::nZProc()/2) + mode_i;
    if(!mode_j) continue;

    ckt = cos(sign * mode_j * phi);
    skt = sin(sign * mode_j * phi);

    for(field_i = 0; field_i < context->domain->nField(); field_i++) {
      elements_to_logical(fields[field_i]->plane(2*mode_i+0), data_r);
      elements_to_logical(fields[field_i]->plane(2*mode_i+1), data_c);

      for(dof_i = 0; dof_i < Femlib::ivalue("NELS_X")*Femlib::ivalue("NELS_Y")*elOrd*elOrd; dof_i++) {
        rTmp = +ckt*data_r[dof_i] + skt*data_c[dof_i];
        cTmp = -skt*data_r[dof_i] + ckt*data_c[dof_i];
        data_r[dof_i] = rTmp;
        data_c[dof_i] = cTmp;
      }

      logical_to_elements(data_r, fields[field_i]->plane(2*mode_i+0));
      logical_to_elements(data_c, fields[field_i]->plane(2*mode_i+1));
    }
    //logical_to_elements(data_r, fields[field_i]->plane(2*mode_i+0));
    //logical_to_elements(data_c, fields[field_i]->plane(2*mode_i+1));
  }
  delete[] data_r;
  delete[] data_c;
}

/***********************************************************************************/

int elBndryIndex(int pt) {
  int nP = Geometry::nP();

  if(pt == nP*nP-1)   return 4*(nP-1) - 1;                   // top right
  if(pt < nP)         return pt;                             // bottom
  if(pt%nP == 0)      return nP - 1 + 2*(pt/nP) - 1;         // left
  if(pt%nP == nP-1)   return nP - 1 + 2*(pt/nP);             // right
  if(pt >= nP*(nP-1)) return nP - 1 + 2*(pt/nP) - 1 + pt%nP; // top

  return -1;
}

void build_addToVector(Context* context, vector<Field*> fields) {
  const NumberSys* numSys;
  const int_t* btog;
  const int_t* bmask;
  int np2 = Geometry::nP() * Geometry::nP();
  int* inserted;
  int gid, el_bndry;
  int index;
  int n_bndry[3];

  context->n_mesh_sum = context->n_mesh_max = 0;

  context->addToVector = new int*[3];

  for(int fd_i = 0; fd_i < 3; fd_i++) {
    context->addToVector[fd_i] = new int[Geometry::nElmt() * np2];
    
    //numSys = fields[fd_i]->_bsys->Nsys(0); // 0th mode
    numSys = fields[fd_i]->_bsys->Nsys(Geometry::procID() * (Geometry::nZProc()/2) * Femlib::ivalue ("BETA"));
    btog = numSys->btog();                 // boundary to global index
    bmask = numSys->bmask();               // '1' if node is an essential bc

    for(int pt_i = 0; pt_i < Geometry::nElmt() * np2; pt_i++) context->addToVector[fd_i][pt_i] = -1;

    // number of boundary nodes for this field (including bcs)
    n_bndry[fd_i] = 0;
    for(int pt_i = 0; pt_i < Geometry::nElmt() * Geometry::nExtElmt(); pt_i++) 
      if(btog[pt_i] > n_bndry[fd_i] && !bmask[pt_i]) 
        n_bndry[fd_i] = btog[pt_i];

      n_bndry[fd_i]++;
if(!Geometry::procID())cout<<fd_i<<": num_bndry: "<<n_bndry[fd_i]<<endl;

    inserted = new int[n_bndry[fd_i]];
    for(int pt_i = 0; pt_i < n_bndry[fd_i]; pt_i++) inserted[pt_i] = -1;

    index = 0;
    for(int el_i = 0; el_i < Geometry::nElmt(); el_i++) {
      for(int pt_i = 0; pt_i < np2; pt_i++) {
        // element boundary node
        if((el_bndry = elBndryIndex(pt_i)) != -1) {
          // not an essential bc node
          if(!bmask[el_i*Geometry::nExtElmt() + el_bndry]) {
            gid = btog[el_i*Geometry::nExtElmt() + el_bndry];
            if(inserted[gid] == -1) {
              inserted[gid] = index;
              context->addToVector[fd_i][el_i*np2 + pt_i] = index;
              index++;
            } else {
              context->addToVector[fd_i][el_i*np2 + pt_i] = inserted[gid];
            }
          }
        // element internal node
        } else {
          context->addToVector[fd_i][el_i*np2 + pt_i] = index;
          index++;
        }
      }
    }
    // number of unique non bc dofs for this field in the spectral element mesh
    context->n_mesh[fd_i] = index;
    context->n_mesh_sum += index;
    if(index > context->n_mesh_max) context->n_mesh_max = index;
    cout << Geometry::procID() << ":\tnumber of non-bc mesh dofs for field: " << fd_i << ": " << context->n_mesh[fd_i] << endl;

    delete[] inserted;
  }

  context->n_mesh_sum_proc = new int[Geometry::nProc()];
  context->n_mesh_sum_proc[Geometry::procID()] = context->n_mesh_sum;
  for(int proc_i = 0; proc_i < Geometry::nProc(); proc_i++) {
    MPI_Bcast(&context->n_mesh_sum_proc[proc_i], 1, MPI_INT, proc_i, MPI_COMM_WORLD);
  }

  cout << Geometry::procID() << ":\tn_mesh_sum: " << context->n_mesh_sum << "\t, n_mesh_max: " << context->n_mesh_max << endl;
}

// NOTE: assumes addToVector has already been built
void build_coordWeights(Context* context) {
  int index;
  int np2 = Geometry::nP() * Geometry::nP();
  const real_t* q4;
  context->coord_weights = new real_t*[3];

  for(int fd_i = 0; fd_i < 3; fd_i++) {
    context->coord_weights[fd_i] = new real_t[context->n_mesh_max];

    for(int dof_i = 0; dof_i < context->n_mesh[fd_i]; dof_i++) context->coord_weights[fd_i][dof_i] = 0.0;

    for(int el_i = 0; el_i < Geometry::nElmt(); el_i++) {
      q4 = context->domain->elmt[el_i]->GetQ4();
      for(int pt_i = 0; pt_i < np2; pt_i++) {
        index = context->addToVector[fd_i][el_i*np2+pt_i];
        if(index != -1) {
          context->coord_weights[fd_i][index] += q4[pt_i];
        }
      }
    }

    for(int el_i = 0; el_i < Geometry::nElmt(); el_i++) {
      for(int pt_i = 0; pt_i < np2; pt_i++) {
        index = context->addToVector[fd_i][el_i*np2+pt_i];
        if(index != -1 && context->coord_weights[fd_i][pt_i] < 1.0e-8) {
          context->coord_weights[fd_i][index] = context->domain->elmt[el_i]->area() / np2;
cout<<Geometry::procID()<<": updating coord weight, field: "<<fd_i<<", index: "<<index<<",\tnew weight: "<<context->coord_weights[fd_i][index]<<endl;
        }
      }
    }
  }
}

void elements_to_vector(Context* context, int field_i, real_t* data_els, real_t* data_vec, bool fwd) {
  int index;
  int np2 = Geometry::nP() * Geometry::nP();

  for(int el_i = 0; el_i < Geometry::nElmt(); el_i++) {
    for(int pt_i = 0; pt_i < np2; pt_i++) {
      index = context->addToVector[field_i][el_i*np2+pt_i];
      if(index != -1) {
        if(fwd) data_vec[index] = data_els[el_i*np2+pt_i];
        else    data_els[el_i*np2+pt_i] = data_vec[index];
      }
    }
  }
}

int LocalIndex(Context* context, int field_i, int plane_i, int point_i) {
/*
  int shift_1 = (plane_i/2) * 3 * 2 * context->nDofsPlane;
  int shift_2 =  point_i    * 3 * 2;

  return shift_1 + shift_2 + 2*field_i + plane_i%2;
*/
  int index = 0;

  for(int shift_i = 0; shift_i < field_i; shift_i++) index += Geometry::nZProc() * context->n_mesh[shift_i];
  index += (plane_i * context->n_mesh[field_i] + point_i);

  return index;
}

void _UnpackX(Context* context, vector<AuxField*> fields, real_t* phi, real_t* tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, kk, ll, fd_i, index;
  const PetscScalar *xArray;
  double scale;
  real_t* data_u = new double[context->n_mesh_max];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x, xl, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArrayRead(xl, &xArray);
  for(kk = 0; kk < Geometry::nZProc(); kk++) {
    ll = ( Geometry::procID() * Geometry::nZProc() + kk ) / 2;

    for(fd_i = 0; fd_i < 3; fd_i++) {
      for(ii = 0; ii < context->n_mesh_max; ii++) data_u[ii] = 0.0;

      for(ii = 0; ii < context->n_mesh[fd_i]; ii++) {
        scale  = 2.0/(2.0 + fabs(ll)) * context->coord_weights[fd_i][ii];
        scale *= context->u_scale[fd_i];

        index = LocalIndex(context, fd_i, kk, ii);
        data_u[ii] = xArray[index] / scale;

        elements_to_vector(context, fd_i, fields[fd_i]->plane(kk), data_u, false);
      }
    }
  }

  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    index = Geometry::nZProc() * context->n_mesh_sum;

    *phi = xArray[index++];
    if(!context->travelling_wave) *tau = xArray[index++];
  }
  VecRestoreArrayRead(xl, &xArray);

  *fields[2] += *context->uBar;

  VecDestroy(&xl);
  delete[] data_u;
}

void _RepackX(Context* context, vector<AuxField*> fields, real_t phi, real_t tau, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, kk, ll, fd_i, index;
  PetscScalar *xArray;
  double scale;
  real_t* data_u = new real_t[context->n_mesh_max];
  Vec xl;

  *fields[2] -= *context->uBar;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecZeroEntries(xl);
  VecGetArray(xl, &xArray);

  for(kk = 0; kk < Geometry::nZProc(); kk++) {
    ll = ( Geometry::procID() * Geometry::nZProc() + kk ) / 2;

    for(fd_i = 0; fd_i < 3; fd_i++) {
      for(ii = 0; ii < context->n_mesh_max; ii++) data_u[ii] = 0.0;

      elements_to_vector(context, fd_i, fields[fd_i]->plane(kk), data_u, true);

      for(ii = 0; ii < context->n_mesh[fd_i]; ii++) {
        scale  = 2.0/(2.0 + fabs(ll)) * context->coord_weights[fd_i][ii];
        scale *= context->u_scale[fd_i];

        index = LocalIndex(context, fd_i, kk, ii);
        xArray[index] = data_u[ii] * scale;
      }
    }
  }

  // phase shift data lives on the 0th processors part of the vector
  if(!Geometry::procID()) {
    index = Geometry::nZProc() * context->n_mesh_sum;

    xArray[index++] = phi;
    if(!context->travelling_wave) xArray[index++] = tau;
  }
  VecRestoreArray(xl, &xArray);

  VecScatterBegin(context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  *fields[2] += *context->uBar;

  VecDestroy(&xl);
  delete[] data_u;
}

void _phase_shift_z(Context* context, double phi, double sign, vector<Field*> fields) {
  int elOrd = Geometry::nP() - 1;
  int mode_i, mode_j, field_i, dof_i;
  double ckt, skt, rTmp, cTmp;
  double* data_r = new double[context->n_mesh_max];
  double* data_c = new double[context->n_mesh_max];

  for(mode_i = 0; mode_i < Geometry::nZProc()/2; mode_i++) {
    mode_j = Geometry::procID() * (Geometry::nZProc()/2) + mode_i;
    if(!mode_j) continue;

    ckt = cos(sign * mode_j * phi);
    skt = sin(sign * mode_j * phi);

    for(field_i = 0; field_i < 3; field_i++) {
      for(dof_i = 0; dof_i < context->n_mesh_max; dof_i++) { data_r[dof_i] = data_c[dof_i] = 0.0; }

      elements_to_vector(context, field_i, fields[field_i]->plane(2*mode_i+0), data_r, true);
      elements_to_vector(context, field_i, fields[field_i]->plane(2*mode_i+1), data_c, true);

      for(dof_i = 0; dof_i < context->n_mesh[field_i]; dof_i++) {
        rTmp = +ckt*data_r[dof_i] + skt*data_c[dof_i];
        cTmp = -skt*data_r[dof_i] + ckt*data_c[dof_i];
        data_r[dof_i] = rTmp;
        data_c[dof_i] = cTmp;
      }

      elements_to_vector(context, field_i, fields[field_i]->plane(2*mode_i+0), data_r, false);
      elements_to_vector(context, field_i, fields[field_i]->plane(2*mode_i+1), data_c, false);
    }
  }
  delete[] data_r;
  delete[] data_c;
}

void velocity_scales(Context* context) {
  double Ku[3], Ku_bar;
  double fac[3] = {0.2, 0.2, 0.5};
  double Lz = Femlib::value("TWOPI / BETA");

  // integrate the energy in each component separately
  for(int field_i = 0; field_i < 3; field_i++) {
    *context->u0[0] = 0.0;
    *context->u0[1] = 0.0;
    *context->u0[2] = 0.0;
    *context->u0[field_i] = *context->ui[field_i];
    context->u0[field_i]->transform(INVERSE);
    context->fi[0]->innerProduct(context->u0, context->u0);
    context->u0[field_i]->transform(FORWARD);
    context->fi[0]->transform(FORWARD);
    *context->fi[0] *= 0.5;
    Ku[field_i] = Lz * context->fi[0]->integral(0);
    MPI_Bcast(&Ku[field_i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(!Geometry::procID()) cout << field_i << ": ke: " << Ku[field_i] << endl;
  }

  // now integrate the mean flow
  *context->u0[0] = 0.0;
  *context->u0[1] = 0.0;
  *context->u0[2] = 0.0;
  if(!Geometry::procID()) {
    for(int dof_i = 0; dof_i < context->n_mesh[2]; dof_i++) {
      context->u0[2]->plane(0)[dof_i] = context->ui[2]->plane(0)[dof_i];
    }
  }
  context->u0[2]->transform(INVERSE);
  context->fi[2]->innerProduct(context->u0, context->u0);
  context->u0[2]->transform(FORWARD);
  context->fi[2]->transform(FORWARD);
  *context->fi[2] *= 0.5;
  Ku_bar = Lz * context->fi[2]->integral(0);
  MPI_Bcast(&Ku_bar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(!Geometry::procID()) cout << "Ku_bar:   " << Ku_bar << endl;
  Ku[2] -= Ku_bar;
  if(!Geometry::procID()) cout << "Ku_prime: " << Ku[2]  << endl;

  if(Ku_bar < 1.0e-6) Ku_bar = 1.0;
  for(int field_i = 0; field_i < 3; field_i++) {
    context->u_scale[field_i] = sqrt(fac[field_i] * Ku_bar / Ku[field_i]);
    if(!Geometry::procID()) cout << field_i << " velocity scale: " << context->u_scale[field_i] << endl;
  }
}

void base_profile(Context* context, AuxField* uz, real_t scale, AuxField* uBar) {
  real_t energy;
  double Lz = Femlib::value("TWOPI / BETA");
  if(!Geometry::procID()) cout << "computing base profile for scale: " << scale << endl;

  *uBar = 0.0;

  if(scale < 1.0e-8) return;

  if(!Geometry::procID()) {
    for(int dof_i = 0; dof_i < context->n_mesh[2]; dof_i++) {
      uBar->plane(0)[dof_i] = scale * uz->plane(0)[dof_i];
    }
  }

  // check the energies
  *context->u0[0] = 0.0;
  *context->u0[1] = 0.0;
  *context->u0[2] = *uz;
  context->u0[2]->transform(INVERSE);
  context->fi[2]->innerProduct(context->u0, context->u0);
  context->fi[2]->transform(FORWARD);
  *context->fi[2] *= 0.5;
  if(!Geometry::procID()) cout << "Ku:       " << Lz * context->fi[2]->integral(0) << endl;

  *context->u0[2] = *uBar;
  context->u0[2]->transform(INVERSE);
  context->fi[2]->innerProduct(context->u0, context->u0);
  context->fi[2]->transform(FORWARD);
  *context->fi[2] *= 0.5;
  if(!Geometry::procID()) cout << "Ku_bar:   " << Lz * context->fi[2]->integral(0) << endl;

  *context->u0[2] = *uz;
  *context->u0[2] -= *uBar;
  context->u0[2]->transform(INVERSE);
  context->fi[2]->innerProduct(context->u0, context->u0);
  context->fi[2]->transform(FORWARD);
  *context->fi[2] *= 0.5;
  if(!Geometry::procID()) cout << "Ku_prime: " << Lz * context->fi[2]->integral(0) << endl;
  *context->u0[2] += *uBar;
}

void _assign_scatter_semtex(Context* context) {
  int   nDofs_l     = Geometry::nZProc() * context->n_mesh_sum;
  int   nDofs_g;//     = Geometry::nZ()     * context->n_mesh_sum;
  int   nShifts     = (!context->travelling_wave) ? 2 : 1;
  int   ind_i       = 0;
  int   start       = 0;//       = Geometry::procID() * nDofs_l;
  int*  inds        = new int[context->localSize];
  IS    isl, isg;
  Vec   vl, vg;

  if(Geometry::procID()) {
    for(int proc_i = 0; proc_i < Geometry::procID(); proc_i++) {
      start += Geometry::nZProc() * context->n_mesh_sum_proc[proc_i];
    }
  }
  MPI_Allreduce(&context->n_mesh_sum, &nDofs_g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  nDofs_g *= Geometry::nZProc();

  for(int ind_j = 0; ind_j < nDofs_l; ind_j++) {
    inds[ind_i++] = start + ind_j;
  }

  // assign the phase shifts from the 0th processor
  if(!Geometry::procID()) {
    for(int ind_j = 0; ind_j < nShifts; ind_j++) {
      inds[ind_i++] = nDofs_g + ind_j;
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
