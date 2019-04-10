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

#define NELS_X 30
#define NELS_Y 8
#define NFIELD 3

void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

double N0(double x, double y) { return 0.25*(1.0-x)*(1.0-y); }
double N1(double x, double y) { return 0.25*(1.0+x)*(1.0-y); }
double N2(double x, double y) { return 0.25*(1.0+x)*(1.0+y); }
double N3(double x, double y) { return 0.25*(1.0-x)*(1.0+y); }
double Ni(int i, double x, double y) {
  if     (i==0) return N0(x, y);
  else if(i==1) return N1(x, y);
  else if(i==2) return N2(x, y);
  else if(i==3) return N3(x, y);
  else abort();
}

double dN0dx(double x, double y) { return -0.25*(1.0-y); }
double dN1dx(double x, double y) { return +0.25*(1.0-y); }
double dN2dx(double x, double y) { return +0.25*(1.0+y); }
double dN3dx(double x, double y) { return -0.25*(1.0+y); }
double dNidx(int i, double x, double y) {
  if     (i==0) return dN0dx(x, y);
  else if(i==1) return dN1dx(x, y);
  else if(i==2) return dN2dx(x, y);
  else if(i==3) return dN3dx(x, y);
  return 0.0;
}

double dN0dy(double x, double y) { return -0.25*(1.0-x); }
double dN1dy(double x, double y) { return -0.25*(1.0+x); }
double dN2dy(double x, double y) { return +0.25*(1.0+x); }
double dN3dy(double x, double y) { return +0.25*(1.0-x); }
double dNidy(int i, double x, double y) {
  if     (i==0) return dN0dy(x, y);
  else if(i==1) return dN1dy(x, y);
  else if(i==2) return dN2dy(x, y);
  else if(i==3) return dN3dy(x, y);
  return 0.0;
}

double** f2r_transform(double lx, int nr, int nf) {
  double kj, xi;
  double** A = new double*[nr];

  for(int ii = 0; ii < nr; ii++) {
    xi = ii * (lx / nr);

    A[ii] = new double[nf];
    for(int jj = 0; jj < nf; jj++) {
      kj = (jj/2) * 2.0 * M_PI / lx;
      A[ii][jj] = (jj%2==0) ? cos(+kj * xi) : sin(+kj * xi);
    }
  }
  return A;
}

double** r2f_transform(double lx, int nr, int nf) {
  double ki, xj;
  double** A = new double*[nf];

  for(int ii = 0; ii < nf; ii++) {
    ki = (ii/2) * 2.0 * M_PI / lx;

    A[ii] = new double[nr];
    for(int jj = 0; jj < nr; jj++) {
      xj = jj * (lx / nr);
      A[ii][jj] = (ii%2==0) ? cos(-ki * xj) : sin(-ki * xj);
    }
  }
  return A;
}

double** mat_mat_mult(double** A, double** B, int ni, int nj, int nk) {
  double** C = new double*[ni];

  for(int ii = 0; ii < ni; ii++) {
    C[ii] = new double[nj];

    for(int jj = 0; jj < nj; jj++) {
      C[ii][jj] = 0.0;

      for(int kk = 0; kk < nk; kk++) {
        C[ii][jj] += A[ii][kk] * B[kk][jj];
      }
    }
  }
  return C;
}

double** real_space_pc(double lx, int nr, int nf, double* data_x, double* data_y, double* data_z) {
  int il, ir;
  double CX = 0.25;
  double dx = lx / nr;
  double** F2R = f2r_transform(lx, nr, nf);
  double** LAP;
  double** PRECON;

  // assemble the real space finite difference gradient operator
  LAP = new double*[nr];
  for(int ii = 0; ii < nr; ii++) {
    LAP[ii] = new double[nr];

    for(int jj = 0; jj < nr; jj++) LAP[ii][jj] = 0.0;

    il = (ii-1+nr)%nr;
    ir = (ii+1)%nr;
    LAP[ii][il] = 1.0/dx;
    LAP[ii][ii] = 2.0/dx;
    LAP[ii][ir] = 1.0/dx;

    // preconditioned convective terms
    if(data_x) LAP[ii][ii] += CX*(data_x[ir] - data_x[il])/(2.0*dx);
    if(data_y) LAP[ii][ii] += CX*(data_y[ir] - data_y[il])/(2.0*dx);
    if(data_z) LAP[ii][ii] += CX*(data_z[ir] - data_z[il])/(2.0*dx);
  }

  // assemble the preconditioner
  PRECON = mat_mat_mult(LAP, F2R, nr, nf, nr);

  // x modes along the axis - return the identity
  if(lx < 1.0e-6) {
    for(int ii = 0; ii < nf; ii++) { 
      for(int jj = 0; jj < nr; jj++) PRECON[ii][jj] = 0.0;
      PRECON[ii][ii] = 1.0; 
    }
  }

  for(int ii = 0; ii < nr; ii++) { delete[] F2R[ii]; } delete[] F2R;
  for(int ii = 0; ii < nr; ii++) { delete[] LAP[ii]; } delete[] LAP;

  return PRECON;
}

void schur_complement_constraints(Context* context, double* schur) {
  int     elOrd           = Geometry::nP() - 1;
  int     nNodesX         = NELS_X*elOrd;
  int     nModesX         = context->nModesX;
  double  schur_local[3];
  double  k_x, k_z;
  int     plane_j;
  int     pt_j, el_j;
  double  p_y;
  int     index;
  int     pRow, pCols[99999];
  int     nDofsCube_l     = Geometry::nZProc() * context->nDofsPlane;
  int     nl              = context->nField * nDofsCube_l;
  int     nStep;
  double* rxl             = new double[nl];
  double* rzl             = new double[nl];
  double* rtl             = new double[nl];
  double* cxl             = new double[nl];
  double* czl             = new double[nl];
  double* ctl             = new double[nl];
  double* data_r          = new double[NELS_Y*elOrd*nModesX];
  double* data_i          = new double[NELS_Y*elOrd*nModesX];

  // integrate the state and phase shifted state forwards
  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] = *context->ui[field_i];
    }
    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);
    //delete context->analyst;
    //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
    integrate(skewSymmetric, context->domain, context->bman, context->analyst, context->ff);
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] -= *context->ui[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    } 
    Femlib::ivalue("N_STEP", nStep);
  }

  // set the constraints as a schur complement (rows)
  for(int dof_i = 0; dof_i < nl; dof_i++) { 
    rxl[dof_i] = rzl[dof_i] = rtl[dof_i] = 0.0;
    cxl[dof_i] = czl[dof_i] = ctl[dof_i] = 0.0; 
  }

  for(int field_i = 0; field_i < context->nField; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->ui[field_i], data_r);
      for(int node_j = 0; node_j < NELS_Y * elOrd; node_j++) {
        for(int mode_i = 0; mode_i < nModesX; mode_i++) {
          index = field_i * nDofsCube_l + plane_i * context->nDofsPlane + node_j * nModesX + mode_i;
          //k_x = (context->xmax / 2.0 / M_PI) * context->theta_i[0] * (mode_i / 2);
          k_x = (2.0 * M_PI / context->xmax) * (mode_i / 2);

          if(mode_i % 2 == 0) {
            rxl[index] = -k_x * data_r[node_j * nModesX + mode_i + 1];
          } else {
            rxl[index] = +k_x * data_r[node_j * nModesX + mode_i - 1];
          }
        }
      }
    }
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i+=2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      SEM_to_Fourier(plane_i+0, context, context->ui[field_i], data_r);
      SEM_to_Fourier(plane_i+1, context, context->ui[field_i], data_i);
      for(int dof_i = 0; dof_i < NELS_Y * elOrd * nModesX; dof_i++) {
        el_j = dof_i / (nModesX * elOrd);
        pt_j = (dof_i / nModesX) % elOrd;
        p_y  = context->elmt[el_j]->_ymesh[pt_j*(elOrd+1)];
        if(fabs(p_y) < 1.0e-6) {
          k_z = 0.0;
        } else {
          // assume a radius of 1, so scale by (2*pi) / (2*pi*r)
          k_z  = (1.0 / p_y) * (plane_j / 2);
        }

        index = field_i * nDofsCube_l + (plane_i+0) * context->nDofsPlane + dof_i;
        rzl[index] = -k_z * data_i[dof_i];
        index = field_i * nDofsCube_l + (plane_i+1) * context->nDofsPlane + dof_i;
        rzl[index] = +k_z * data_r[dof_i];
      }
    }
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->domain->u[field_i], data_r);
      for(int dof_i = 0; dof_i < NELS_Y * elOrd * nModesX; dof_i++) {
        index = field_i * nDofsCube_l + plane_i * context->nDofsPlane + dof_i;
        rtl[index] = data_r[dof_i];
      }
    }
  }

  // phase shifted state vector
  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] = *context->ui[field_i];
      phase_shift_x(context, context->theta_i[0], -1.0, context->domain->u);
      phase_shift_z(context, context->phi_i[0],   -1.0, context->domain->u);
      *context->uj[field_i] = *context->domain->u[field_i];
    }
    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);
    //delete context->analyst;
    //context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
    integrate(skewSymmetric, context->domain, context->bman, context->analyst, context->ff);
    for(int field_i = 0; field_i < context->nField; field_i++) {
      *context->domain->u[field_i] -= *context->uj[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    } 
    Femlib::ivalue("N_STEP", nStep);
  }

  // set the constraints as a schur complement (columns)
  for(int field_i = 0; field_i < context->nField; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->uj[field_i], data_r);
      for(int node_j = 0; node_j < NELS_Y * elOrd; node_j++) {
        for(int mode_i = 0; mode_i < nModesX; mode_i++) {
          index = field_i * nDofsCube_l + plane_i * context->nDofsPlane + node_j * nModesX + mode_i;
          k_x = (context->xmax / 2.0 / M_PI) * context->theta_i[0] * (mode_i / 2);

          if(mode_i % 2 == 0) {
            cxl[index] = +k_x * data_r[node_j * nModesX + mode_i + 1];
          } else {
            cxl[index] = -k_x * data_r[node_j * nModesX + mode_i - 1];
          }
        }
      }
    }
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      SEM_to_Fourier(plane_i+0, context, context->uj[field_i], data_r);
      SEM_to_Fourier(plane_i+1, context, context->uj[field_i], data_i);
      for(int dof_i = 0; dof_i < NELS_Y * elOrd * nModesX; dof_i++) {
        el_j = dof_i / (nModesX * elOrd);
        pt_j = (dof_i / nModesX) % elOrd;
        p_y  = context->elmt[el_j]->_ymesh[pt_j*(elOrd+1)];
        if(fabs(p_y) < 1.0e-6) {
          k_z = 0.0;
        } else {
          // assume a radius of 1, so scale by (2*pi) / (2*pi*r)
          k_z  = (1.0 / p_y) * (plane_j / 2);
        }

        index = field_i * nDofsCube_l + (plane_i+0) * context->nDofsPlane + dof_i;
        czl[index] = +k_z * data_i[dof_i];
        index = field_i * nDofsCube_l + (plane_i+1) * context->nDofsPlane + dof_i;
        czl[index] = -k_z * data_r[dof_i];
      }
    }
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      SEM_to_Fourier(plane_i, context, context->domain->u[field_i], data_r);
      for(int dof_i = 0; dof_i < NELS_Y * elOrd * nModesX; dof_i++) {
        index = field_i * nDofsCube_l + plane_i * context->nDofsPlane + dof_i;
        ctl[index] = data_r[dof_i];
      }
    }
  }

  // add the constraints to the preconditioner
  /*for(int field_i = 0; field_i < context->nField; field_i++) {
    pRow = context->nField * Geometry::nZ() * context->nDofsPlane;
    for(int dof_i = 0; dof_i < nDofsCube_l; dof_i++) {
      pCols[dof_i] = context->lShift[0][field_i] + dof_i;
    }
    MatSetValues(P, 1, &pRow, nDofsCube_l, pCols, &rxl[field_i * nDofsCube_l], INSERT_VALUES);
    MatSetValues(P, nDofsCube_l, pCols, 1, &pRow, &cxl[field_i * nDofsCube_l], INSERT_VALUES);

    pRow++;
    MatSetValues(P, 1, &pRow, nDofsCube_l, pCols, &rzl[field_i * nDofsCube_l], INSERT_VALUES);
    MatSetValues(P, nDofsCube_l, pCols, 1, &pRow, &czl[field_i * nDofsCube_l], INSERT_VALUES);

    pRow++;
    MatSetValues(P, 1, &pRow, nDofsCube_l, pCols, &rtl[field_i * nDofsCube_l], INSERT_VALUES);
    MatSetValues(P, nDofsCube_l, pCols, 1, &pRow, &ctl[field_i * nDofsCube_l], INSERT_VALUES);
  }*/

  // add diagonal entries where required
  schur_local[0] = schur_local[1] = schur_local[2] = 0.0;
  for(int dof_i = 0; dof_i < nl; dof_i++) {
    //schur_local[0] += rxl[dof_i] * cxl[dof_i];
    schur_local[0] -= rxl[dof_i] * cxl[dof_i];
    //schur_local[1] += rzl[dof_i] * czl[dof_i];
    schur_local[1] -= rzl[dof_i] * czl[dof_i];
    if(!context->travelling_wave)
      //schur_local[2] += rtl[dof_i] * ctl[dof_i];
      schur_local[2] -= rtl[dof_i] * ctl[dof_i];
//cout << "rxl: " << rxl[dof_i] << 
//      "\trzl: " << rzl[dof_i] << 
//      "\trtl: " << rtl[dof_i] << 
//      "\tczl: " << cxl[dof_i] << 
//      "\tczl: " << czl[dof_i] << 
//      "\tctl: " << ctl[dof_i] << endl;
  }
  MPI_Allreduce(schur_local, schur, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete[] rxl;
  delete[] rzl;
  delete[] rtl;
  delete[] cxl;
  delete[] czl;
  delete[] ctl;
  delete[] data_r;
  delete[] data_i;
}

void assemble_K(Context* context, int plane_i, Mat K, Mat P) {
  int     elOrd         = Geometry::nP() - 1;
  int     nNodesX       = NELS_X*elOrd;
  int     nModesX       = context->nModesX;
  int     row_k[3], col_k[3], nCols;
  int     el_i;
  int     plane_j;
  int     row_dof, col_dof;
  const int* cols;
  const double* vals;
  int     pRow, pCols[99999];
  double  pVals[99999];
  double  k_x, k_z;
  double  yi[2], delta_y, det;
  double  qx[]          = {+1.0, +1.0};
  double  qy[]          = {-1.0, +1.0};
  double  kinvis        = Femlib::value("KINVIS");
  double  beta          = Femlib::value("BETA");
  double  alpha         = context->xmax * beta / (2.0*M_PI) / (context->nModesX / 2);
  double  dt            = Femlib::value("D_T");
  double  **PCX, **PCZ;
  double* data_f        = new double[NELS_Y*elOrd*nModesX];

  plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
  k_z = (plane_j / 2) / beta;
  if(plane_j < 2) k_z = 1.0;

  // same number of nodes as modes in x 
  SEM_to_Fourier(plane_i, context, context->ui[0], data_f);

  for(int elmt_y = 0; elmt_y < elOrd*NELS_Y; elmt_y++) {
    el_i = (elmt_y / elOrd) * NELS_X;

    yi[0] = context->elmt[el_i]->_ymesh[(elmt_y%elOrd+0)*(elOrd+1)];
    yi[1] = context->elmt[el_i]->_ymesh[(elmt_y%elOrd+1)*(elOrd+1)];
    
    delta_y = 0.5 * fabs(yi[1] - yi[0]);

    PCX = real_space_pc(2.0*M_PI, nModesX, nModesX, NULL, NULL, NULL);
    PCZ = real_space_pc(2.0*M_PI*yi[0], Geometry::nZ(), Geometry::nZ(), &data_f[elmt_y*nModesX], NULL, NULL);

    for(int mode_x = 0; mode_x < nModesX; mode_x++) {
      pVals[0] = pVals[1] = pVals[2] = pVals[3] = 0;
      for(int row = 0; row < 2; row++) {
        row_k[0] = (elmt_y+row)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary
        row_k[1] = row_k[0] + context->nDofsPlane;
        row_k[2] = row_k[1] + context->nDofsPlane;

        if(elmt_y+row == elOrd*NELS_Y) continue; // outer boundary (dirichlet bcs)

        for(int col = 0; col < 2; col++) {
          col_k[0] = (elmt_y+col)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary
          col_k[1] = col_k[0] + context->nDofsPlane;
          col_k[2] = col_k[1] + context->nDofsPlane;

          if(elmt_y+col == elOrd*NELS_Y) continue; // outer boundary (dirichlet bcs)

          // radial direction terms
          for(int pnt = 0; pnt < 2; pnt++) {
            //if(!elmt_y && !pnt) { // center axis
              //det = alpha * delta_y * beta * (2.0*M_PI/yi[1]);
            //} else {
              //det = alpha * delta_y * beta * (2.0*M_PI/yi[pnt]);
              det = delta_y * (2.0*M_PI*yi[pnt]);
            //}

            // mass matrix terms
            MatSetValue(K, row_k[0], row_k[0], det, ADD_VALUES);
            MatSetValue(K, row_k[1], row_k[1], det, ADD_VALUES);
            MatSetValue(K, row_k[2], row_k[2], det, ADD_VALUES);

            pVals[2*row+col] += -(dt * kinvis) * det * dNidy(row+1, qx[pnt], qy[pnt]) / delta_y * 
                                                       dNidy(col+1, qx[pnt], qy[pnt]) / delta_y;
          } //pnt
        } // col
      } // row

      //if(!elmt_y) { // center axis
        //det = alpha * delta_y * beta * (2.0*M_PI/yi[1]);
      //} else {
        //det = alpha * delta_y * beta * (2.0*M_PI/yi[0]);
        det = delta_y * (2.0*M_PI*yi[0]);
      //}

      // set the radial (finite element) laplacian for this (1d) element
      row_k[0] = row_k[1];
      row_k[1] = row_k[0] + nModesX;
      MatSetValues(K, 2, row_k, 2, row_k, pVals, ADD_VALUES);

      // set the axial terms
      pRow = elmt_y*nModesX + mode_x;
      for(int col = 0; col < nModesX; col++) {
        pCols[col] = elmt_y*nModesX + col;
        pVals[col] = (dt * kinvis) * det * PCX[mode_x][col];
      }
      MatSetValues(K, 1, &pRow, nModesX, pCols, pVals, ADD_VALUES);

      // set the azimuthal terms
      // TODO: do this outside of the plane_i loop so that this may be applied to all modes
      pRow = elmt_y*nModesX + mode_x + 2*context->nDofsPlane;
      pVals[0] = (dt * kinvis) * det * PCZ[plane_j][plane_j];
      MatSetValue(K, pRow, pRow, pVals[0], ADD_VALUES);
    } // mode_x

    for(int ii = 0; ii < nModesX;        ii++) { delete[] PCX[ii]; } delete[] PCX;
    for(int ii = 0; ii < Geometry::nZ(); ii++) { delete[] PCZ[ii]; } delete[] PCZ;
  } // elmt_y
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  K, MAT_FINAL_ASSEMBLY);

  // [u,u] block
  for(int row_i = 0; row_i < 3*context->nDofsPlane; row_i++) {
    row_dof = row_i/context->nDofsPlane;
    MatGetRow(K, row_i, &nCols, &cols, &vals);
    pRow = row_i%context->nDofsPlane + plane_i*context->nDofsPlane + context->lShift[0][row_dof];

    for(int col_i = 0; col_i < nCols; col_i++) {
      col_dof = cols[col_i]/context->nDofsPlane;
      pCols[col_i] = cols[col_i]%context->nDofsPlane + plane_i*context->nDofsPlane + context->lShift[0][col_dof];
    }
    MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
    MatRestoreRow(K, row_i, &nCols, &cols, &vals);
  //}

  // add in the -ve identity for the flow map to the next slice   
  //for(int row_i = 0; row_i < 3*context->nDofsPlane; row_i++) {
    //row_dof = row_i/context->nDofsPlane;
    //pRow = row_i%context->nDofsPlane + plane_i*context->nDofsPlane + context->lShift[0][row_dof];
    //MatSetValue(P, pRow, pRow, -1.0, ADD_VALUES);
  }

  delete[] data_f;
}

// as above, but for a fourier-fourier-sem discretisation
//
// fourier modes in the axial and azimuthal dimensions, and 
// spectral elements in the radial dimension
//
// only use basis functions N_1(1,y) and N_2(1,y) as defined above
//
// References:
//
//   Griffith (2009) "An accurate and efficient method for the incompressible Navier-Stokes
//   equations using the projection method as a preconditioner" J. Comp. Phys. 228, 7565-7595
//
//   Elman, Howle, Shadid, Shuttleworth, Tuminaro (2008) "A taxonomy and comparison of parallel 
//   block multi-level preconditioners  for the incompressible Navier-Stokes equations" J. Comp.
//   Phys. 227 1790-1808
//
//   formulated as
//
//   [ K      0    ][ I  K^{-1}G ] = [ K  G ]
//   [ D -DK^{-1}G ][ 0     I    ]   [ D  0 ]
//
void build_preconditioner_ffs(Context* context, Mat P) {
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  int nModesX = context->nModesX;
  int nCols, pRow;
  int nProws, nPcols;
  const int* cols;
  const double* vals;
  Mat K;
  double schur[3];

  MatCreateSeqAIJ(MPI_COMM_SELF, 3*context->nDofsPlane, 3*context->nDofsPlane, nModesX + 4, NULL, &K);
  MatSetOptionsPrefix(K, "K_");
  MatSetFromOptions(K);
  MatZeroEntries(K);
  MatZeroEntries(P);

  for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
    assemble_K(context, plane_i, K, P);
  } // plane_i

  schur_complement_constraints(context, schur);

  if(!Geometry::procID()) {
    pRow = context->nField * context->nDofsPlane * Geometry::nZ() + 0;
    if(fabs(schur[0]) < 1.0e-6) schur[0] = 1.0;
    MatSetValue(P, pRow, pRow, schur[0], INSERT_VALUES);

    pRow = context->nField * context->nDofsPlane * Geometry::nZ() + 1;
    if(fabs(schur[1]) < 1.0e-6) schur[1] = 1.0;
    MatSetValue(P, pRow, pRow, schur[1], INSERT_VALUES);

    if(!context->travelling_wave) {
      pRow = context->nField * context->nDofsPlane * Geometry::nZ() + 2;
      if(fabs(schur[2]) < 1.0e-6) schur[2] = 1.0;
      MatSetValue(P, pRow, pRow, schur[2], INSERT_VALUES);
    }
  }
  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  P, MAT_FINAL_ASSEMBLY);

  // test the matrix
  {
    int mi, mf;
    MatGetOwnershipRange(P, &mi, &mf);
    cout << "[" << Geometry::procID() << "]\t ownership range: "<< mi << "\t->\t" << mf << endl;
    for(int mm = mi; mm < mf; mm++) {
      MatGetRow(P, mm, &nCols, &cols, &vals);
      //if(nCols != 1) {
      //  cout << "[" << Geometry::procID() << "]\t incorrect no. cols in row:   " << mm << "\t, nCols: " << nCols   << endl;
      //}
      //if(cols[0] != mm) {
      //  cout << "[" << Geometry::procID() << "]\t row and column do not match: " << mm << "\t, col:   " << cols[0] << endl;
      //}
      //if(fabs(vals[0] - 1.0) > 1.0e-6) {
      //if(fabs(vals[0]) < 1.0e-6) {
      //  cout << "[" << Geometry::procID() << "]\t incorrect diagonal entry:    " << mm << "\t, val:   " << vals[0] << endl;
      //}
      MatRestoreRow(P, mm, &nCols, &cols, &vals);
    }
  }

  //rpo_set_fieldsplits(context);

  MatDestroy(&K);
}

void build_preconditioner_I(Context* context, Mat P) {
  Vec d;

  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &d);
  VecSet(d, 1.0);
  MatDiagonalSet(P, d, INSERT_VALUES);
  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  P, MAT_FINAL_ASSEMBLY);
  VecDestroy(&d);
}

// NOTE this implementation assumes that each slice fits exactly
// onto one processor!
void rpo_set_fieldsplits(Context* context) {
  KSP ksp;
  PC  pc;
  int uSize = context->nField * Geometry::nZ() * context->nDofsPlane;
  int pSize = 3;

  if(context->is_s) return;

  context->is_s = new IS[1];
  context->is_u = new IS[1];
  context->is_p = new IS[1];

  ISCreateStride(MPI_COMM_WORLD, uSize, Geometry::procID() * (uSize + pSize) + 0,     1, &context->is_u[0]);
  ISCreateStride(MPI_COMM_WORLD, pSize, Geometry::procID() * (uSize + pSize) + uSize, 1, &context->is_p[0]);
  //ISCreateStride(MPI_COMM_WORLD, uSize, 0,     1, &context->is_u[0]);
  //ISCreateStride(MPI_COMM_WORLD, pSize, uSize, 1, &context->is_p[0]);

  SNESGetKSP(context->snes, &ksp);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCFIELDSPLIT);
  PCFieldSplitSetIS(pc, "u", context->is_u[0]);
  PCFieldSplitSetIS(pc, "p", context->is_p[0]);
  PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
  PCSetUp(pc);
}

/*
void rpo_set_fieldsplits(Context* context) {
  PC pc, pc_i, pc_j;
  KSP ksp, *ksp_i, *ksp_j;
  int n_split, m_split, is, iu, ip, nDofs_u, nDofs_p, nSliceProc, startSlice;
  int *inds_s, *inds_u, *inds_p;

  if(context->is_s) return;

  context->is_s = new IS[context->nSlice];
  context->is_u = new IS[context->nSlice];
  context->is_p = new IS[context->nSlice];

  SNESGetKSP(context->snes, &ksp);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCFIELDSPLIT);
  PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);

  nDofs_u = 3 * Geometry::nZ() * context->nDofsPlane;
  nDofs_p = 3;
  if(NFIELD == 4) nDofs_p += Geometry::nZ() * context->nDofsPlane;
  nSliceProc = context->nSlice / Geometry::nProc();
  startSlice = nSliceProc * Geometry::procID();

  for(int slice_i = 0; slice_i < context->nSlice; slice_i++) {
    if(Geometry::procID() == slice_i) {
      ISCreateStride(MPI_COMM_SELF, context->nDofsSlice, slice_i * context->nDofsSlice, 1, &context->is_s[slice_i]);
    } else {
      ISCreateStride(MPI_COMM_SELF, 0, 0, 1, &context->is_s[slice_i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    PCFieldSplitSetIS(pc, NULL, context->is_s[slice_i]);
  }
  PCSetUp(pc);

  // setup a schur complement between the velocity dofs and the phase shift dofs
  PCFieldSplitGetSubKSP(pc, &n_split, &ksp_i);
  for(int slice_i = 0; slice_i < context->nSlice; slice_i++) {
    KSPGetPC(ksp_i[slice_i], &pc_i);
    PCSetType(pc_i, PCFIELDSPLIT);

    if(Geometry::procID() == slice_i) {
      ISCreateStride(MPI_COMM_SELF, nDofs_u, 0,       1, &context->is_u[slice_i]);
      ISCreateStride(MPI_COMM_SELF, nDofs_p, nDofs_u, 1, &context->is_p[slice_i]);
    } else {
      ISCreateStride(MPI_COMM_SELF, 0, 0, 1, &context->is_u[slice_i]);
      ISCreateStride(MPI_COMM_SELF, 0, 0, 1, &context->is_p[slice_i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    PCFieldSplitSetIS(pc_i, "u", context->is_u[slice_i]);
    PCFieldSplitSetIS(pc_i, "p", context->is_p[slice_i]);

    PCFieldSplitSetType(pc_i, PC_COMPOSITE_SCHUR);
    PCSetUp(pc_i);
  }

  for(int slice_i = 0; slice_i < context->nSlice; slice_i++) {
    KSPGetPC(ksp_i[slice_i], &pc_i);
    PCFieldSplitGetSubKSP(pc_i, &m_split, &ksp_j);
    KSPSetType(ksp_j[0], KSPGMRES);
    KSPSetType(ksp_j[1], KSPGMRES);
    KSPGetPC(ksp_j[0], &pc_j);
    PCSetType(pc_j, PCJACOBI);
    KSPGetPC(ksp_j[1], &pc_j);
    PCSetType(pc_j, PCJACOBI);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
*/
