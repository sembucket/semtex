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

#include "rpo_preconditioner.h"

#define XMIN 0.0
#define XMAX (2.0*M_PI)
#define YMIN 0.0
#define YMAX 1.0
#define NELS_X 30
#define NELS_Y 7
#define NFIELD 3

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

double** real_space_pc(double lx, int nr, int nf) {
  int jl, jr;
  double dx = lx / nr;
  double** F2R = f2r_transform(lx, nr, nf);
  double** R2F = r2f_transform(lx, nr, nf);
  double** LAP;
  double** R2F_LAP;
  double** PRECON;

  // assemble the real space finite difference gradient operator
  LAP = new double*[nr];
  for(int ii = 0; ii < nr; ii++) {
    LAP[ii] = new double[nr];

    for(int jj = 0; jj < nr; jj++) LAP[ii][jj] = 0.0;

    for(int jj = 0; jj < nr; jj++) {
      jl = (jj-1+nr)%nr;
      jr = (jj+1)%nr;
      LAP[ii][jl] -= 1.0/dx;
      LAP[ii][jj] += 2.0/dx;
      LAP[ii][jr] -= 1.0/dx;
    }
  }

  // assemble the preconditioner
  R2F_LAP = mat_mat_mult(R2F, LAP, nf, nr, nr);
  PRECON  = mat_mat_mult(R2F_LAP, F2R, nf, nf, nr);

  for(int ii = 0; ii < nf; ii++) { delete[] R2F[ii]; }     delete[] R2F;
  for(int ii = 0; ii < nr; ii++) { delete[] F2R[ii]; }     delete[] F2R;
  for(int ii = 0; ii < nr; ii++) { delete[] LAP[ii]; }     delete[] LAP;
  for(int ii = 0; ii < nf; ii++) { delete[] R2F_LAP[ii]; } delete[] R2F_LAP;

  return PRECON;
}

// Schur complement preconditioning for incompressible Navier-Stokes
// References:
//
//   Griffith (2009) "An accurate and efficient method for the incompressible Navier-Stokes
//   equations using the projection method as a preconditioner" J. Comp. Phys. 228, 7565-7595
//
//   Elman, Howle, Shadid, Shuttleworth, Tuminaro (2008) "A taxonomy and comparison of parallel 
//   block multi-level preconditioners  for the incompressible Navier-Stokes equations" J. Comp.
//   Phys. 227 1790-1808
//
void build_preconditioner(int nSlice, int nDofsSlice, int nDofsPlane, int localSize, int localShift, vector<Element*> elmt, Mat P) {
  int elOrd = Geometry::nP() - 1;
  int nZloc = Geometry::nZProc();
  int rank = Geometry::procID();
  int nx = elOrd*NELS_X;
  int ny = elOrd*NELS_Y;
  int row_i, row_j, col_i, col_j;
  int row_x, row_y, col_x, col_y;
  int row_k[3], col_k[3];
  int nCols;
  int offset_u, offset_v, offset_w, offset_p;
  const int* cols;
  const double* vals;
  int pRow, pCols[99];
  int nProws, nPcols;
  double pVals[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double xi[4], yi[4], det;
  double qx[] = {-1.0, +1.0, +1.0, -1.0};
  double qy[] = {-1.0, -1.0, +1.0, +1.0};
  double delta_x, delta_y;
  double kinvis = Femlib::value("KINVIS");
  double dt = Femlib::value("D_T");
  double beta = Femlib::value("BETA"); // L/(2\pi)
  double k_z;
  Mat G, K, S;
  Mat Kii_inv, D, LK, LKR, DKinv;

  MatCreateSeqAIJ(MPI_COMM_SELF, 3*nx*ny, 1*nx*ny, 16, NULL, &G);
  MatCreateSeqAIJ(MPI_COMM_SELF, 3*nx*ny, 3*nx*ny, 64, NULL, &K);
  MatCreateSeqAIJ(MPI_COMM_SELF, 3*nx*ny, 3*nx*ny,  1, NULL, &Kii_inv);

  MatZeroEntries(G);
  MatZeroEntries(K);
  MatZeroEntries(Kii_inv);

  MatZeroEntries(P);

  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    for(int plane_i = 0; plane_i < nZloc; plane_i++) {
      offset_u = localShift + slice_i*4*nZloc*nx*ny + 0*nZloc*nx*ny + plane_i*nx*ny;
      offset_v = localShift + slice_i*4*nZloc*nx*ny + 1*nZloc*nx*ny + plane_i*nx*ny;
      offset_w = localShift + slice_i*4*nZloc*nx*ny + 2*nZloc*nx*ny + plane_i*nx*ny;
      offset_p = localShift + slice_i*4*nZloc*nx*ny + 3*nZloc*nx*ny + plane_i*nx*ny;

      k_z = ((rank * nZloc + plane_i) / 2) / beta;

      // loop over the SPECTRAL elements
      for(int se = 0; se < NELS_X*NELS_Y; se++) {
        int se_x = se%NELS_X;
        int se_y = se/NELS_X;

        // loop over FINITE elements within the given SPECTRAL element
        for(int fe = 0; fe < elOrd*elOrd; fe++) {
          int fe_x = fe%elOrd;
          int fe_y = fe/elOrd;

          // global coordinates (counter clockwise)
          xi[0] = elmt[se]->_xmesh[(fe_y+0)*(elOrd+1)+fe_x+0];
          xi[1] = elmt[se]->_xmesh[(fe_y+0)*(elOrd+1)+fe_x+1];
          xi[2] = elmt[se]->_xmesh[(fe_y+1)*(elOrd+1)+fe_x+1];
          xi[3] = elmt[se]->_xmesh[(fe_y+1)*(elOrd+1)+fe_x+0];

          yi[0] = elmt[se]->_ymesh[(fe_y+0)*(elOrd+1)+fe_x+0];
          yi[1] = elmt[se]->_ymesh[(fe_y+0)*(elOrd+1)+fe_x+1];
          yi[2] = elmt[se]->_ymesh[(fe_y+1)*(elOrd+1)+fe_x+1];
          yi[3] = elmt[se]->_ymesh[(fe_y+1)*(elOrd+1)+fe_x+0];

          delta_x = 0.5 * fabs(xi[1] - xi[0]);
          delta_y = 0.5 * fabs(yi[3] - yi[0]);
          det  = delta_x * delta_y;
          det *= beta;

          // grad matrix (weak form)
          for(int row = 0; row < 4; row++) {
            row_x = row%2;
            row_y = row/2;

            row_i = se_x*elOrd + fe_x + row_x;
            row_j = se_y*elOrd + fe_y + row_y;

            // periodic in z
            if(row_i == elOrd*NELS_X) continue;

            // skip the first row of nodes (these are dirichlet bcs)
            if(!row_j) continue;
            row_j--;

            row_k[0] = row_j*nx + row_i;
            row_k[1] = row_k[0] + nx*ny;
            row_k[2] = row_k[1] + nx*ny;

            for(int col = 0; col < 4; col++) {
              col_x = col%2;
              col_y = col/2;

              col_i = se_x*elOrd + fe_x + col_x;
              col_j = se_y*elOrd + fe_y + col_y;

              // periodic in z
              if(col_i == elOrd*NELS_X) continue;

              // skip the first row of nodes (these are dirichlet bcs)
              if(!col_j) continue;
              col_j--;

              col_k[0] = col_j*nx + col_i;

              for(int pnt = 0; pnt < 4; pnt++) {
                pVals[0] = -dt * det / delta_x * dNidx(pnt, qx[pnt], qy[pnt]) * Ni(pnt, qx[pnt], qy[pnt]);
                pVals[1] = -dt * det / delta_y * dNidy(pnt, qx[pnt], qy[pnt]) * Ni(pnt, qx[pnt], qy[pnt]);
                pVals[2] =  dt * det * k_z;

                MatSetValues(G, 3, row_k, 1, col_k, pVals, ADD_VALUES);
              }
            }
          }

          // vector laplacian matrix
          for(int row = 0; row < 4; row++) {
            row_x = row%2;
            row_y = row/2;

            row_i = se_x*elOrd + fe_x + row_x;
            row_j = se_y*elOrd + fe_y + row_y;

            // periodic in z
            if(row_i == elOrd*NELS_X) continue;

            // skip the first row of nodes (these are dirichlet bcs)
            if(!row_j) continue;
            row_j--;

            row_k[0] = row_j*nx + row_i;
            row_k[1] = row_k[0] + nx*ny;
            row_k[2] = row_k[1] + nx*ny;

            for(int col = 0; col < 4; col++) {
              col_x = col%2;
              col_y = col/2;

              col_i = se_x*elOrd + fe_x + col_x;
              col_j = se_y*elOrd + fe_y + col_y;

              // periodic in z
              if(col_i == elOrd*NELS_X) continue;

              // skip the first row of nodes (these are dirichlet bcs)
              if(!col_j) continue;
              col_j--;

              col_k[0] = col_j*nx + col_i;
              col_k[1] = col_k[0] + nx*ny;
              col_k[2] = col_k[1] + nx*ny;

              for(int pnt = 0; pnt < 4; pnt++) {
                pVals[0] = -(dt / kinvis) * det / delta_x * dNidx(pnt, qx[pnt], qy[pnt]) / delta_x * dNidx(pnt, qx[pnt], qy[pnt]);
                pVals[1] = -(dt / kinvis) * det / delta_x * dNidx(pnt, qx[pnt], qy[pnt]) / delta_y * dNidy(pnt, qx[pnt], qy[pnt]);
                pVals[3] = -(dt / kinvis) * det / delta_y * dNidy(pnt, qx[pnt], qy[pnt]) / delta_x * dNidx(pnt, qx[pnt], qy[pnt]);
                pVals[4] = -(dt / kinvis) * det / delta_y * dNidy(pnt, qx[pnt], qy[pnt]) / delta_y * dNidy(pnt, qx[pnt], qy[pnt]);
                pVals[8] = -(dt / kinvis) * det * k_z * k_z;

                MatSetValues(K, 3, row_k, 3, col_k, pVals, ADD_VALUES);
                MatSetValue(K, row_k[0], row_k[0], det, ADD_VALUES);
                MatSetValue(K, row_k[1], row_k[1], det, ADD_VALUES);
                MatSetValue(K, row_k[2], row_k[2], det, ADD_VALUES);

                MatSetValue(Kii_inv, row_k[0], row_k[0], 1.0/(det+pVals[0]), ADD_VALUES);
                MatSetValue(Kii_inv, row_k[1], row_k[1], 1.0/(det+pVals[4]), ADD_VALUES);
                MatSetValue(Kii_inv, row_k[2], row_k[2], 1.0/(det+pVals[8]), ADD_VALUES);
              }
            }
          }
        }
      }
      MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(  G, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(  K, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(Kii_inv, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(  Kii_inv, MAT_FINAL_ASSEMBLY);

      // create the Schur complement operator
      MatTranspose(G, MAT_INITIAL_MATRIX, &D);
      MatMatMult(Kii_inv, K, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &LK);
      MatMatMult(LK, Kii_inv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &LKR);
      MatMatMult(D, LKR, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DKinv);
      MatMatMult(DKinv, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &S);
      // schur complement should be scaled by -1, however this is ommitted in 
      // order to account for the G^{T}G term in the azimutal direction (-k_z^2)
      //MatScale(S, -1.0);

      // [u,u] block
      for(int row_i = 0; row_i < 3*nx*ny; row_i++) {
        MatGetRow(K, row_i, &nCols, &cols, &vals);
        pRow = row_i + localShift + slice_i*4*nZloc*nx*ny;
        for(int col_i = 0; col_i < nCols; col_i++) {
          pCols[col_i] = cols[col_i] + localShift + slice_i*4*nZloc*nx*ny;
        }
        MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
        MatRestoreRow(K, row_i, &nCols, &cols, &vals);
      }

      // [u,p] block
      for(int row_i = 0; row_i < 3*nx*ny; row_i++) {
        MatGetRow(G, row_i, &nCols, &cols, &vals);

        pRow = row_i + localShift + slice_i*4*nZloc*nx*ny;
        for(int col_i = 0; col_i < nCols; col_i++) {
          pCols[col_i] = cols[col_i] + offset_p;
          pVals[col_i] = vals[col_i];
        }

        // multiplying by ik_z in the azimuthal direction implies that we have
        // to add contributions to real and imaginary planes separately here
        if(row_i >= 2*nZloc*nx*ny && plane_i%2==0) {
          for(int col_i = 0; col_i < nCols; col_i++) {
            pCols[col_i] += nx*ny; // imaginary plane for this mode
            pVals[col_i] *= -1.0; 
          }
        } else if(row_i >= 2*nZloc*nx*ny && plane_i%2==1) {
          for(int col_i = 0; col_i < nCols; col_i++) {
            pCols[col_i] -= nx*ny; // real plane for this mode
          }
        }

        MatSetValues(P, 1, &pRow, nCols, pCols, pVals, INSERT_VALUES);
        MatRestoreRow(G, row_i, &nCols, &cols, &vals);
      }

      // [p,p] block
      for(int row_i = 0; row_i < nx*ny; row_i++) {
        MatGetRow(S, row_i, &nCols, &cols, &vals);
        pRow = row_i + offset_p;
        for(int col_i = 0; col_i < nCols; col_i++) {
          pCols[col_i] = cols[col_i] + offset_p;
        }
        MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
        MatRestoreRow(S, row_i, &nCols, &cols, &vals);
      }
    }
  }
  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  P, MAT_FINAL_ASSEMBLY);

  MatDestroy(&G);
  MatDestroy(&D);
  MatDestroy(&K);
  MatDestroy(&Kii_inv);
  MatDestroy(&S);
  MatDestroy(&LK);
  MatDestroy(&LKR);
  MatDestroy(&DKinv);
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
void build_preconditioner_ffs(int nSlice, int nDofsSlice, int nDofsPlane, int localSize, int** lShift, int* els, vector<Element*> elmt, Mat P,
                              SNES snes, IS* is_s, IS* is_u, IS* is_p) {
  int elOrd = Geometry::nP() - 1;
  int nZloc = Geometry::nZProc();
  int rank = Geometry::procID();
  int nNodesX = NELS_X*elOrd;
  int nModesX = nNodesX/2 + 2;
  int row_j, col_j;
  int row_x, row_y, col_x, col_y;
  int row_k[3], col_k[3];
  int nCols;
  int el_i, plane_j;
  int ri, rf;
  const int* cols;
  const double* vals;
  int pRow, pCols[99];
  int nProws, nPcols;
  int row_dof, col_dof;
  double pVals[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double yi[2], delta_y, det;
  double qx[] = {+1.0, +1.0};
  double qy[] = {-1.0, +1.0};
  double kinvis = Femlib::value("KINVIS");
  double dt = Femlib::value("D_T");
  double beta = Femlib::value("BETA");
  double alpha = (XMAX - XMIN)/(2.0*M_PI);
  double k_x, k_z;
  double one = 1.0;
  Mat G, K, S;
  Mat Kii_inv, D, LK, LKR, DKinv;
  int nDofsCube = Geometry::nZ() * nDofsSlice;
  double** PCX = real_space_pc(2.0*M_PI, nNodesX, nModesX);

  MatCreateSeqAIJ(MPI_COMM_SELF, 3*nDofsPlane, 1*nDofsPlane, 16, NULL, &G);
  MatCreateSeqAIJ(MPI_COMM_SELF, 1*nDofsPlane, 3*nDofsPlane, 16, NULL, &D);
  MatCreateSeqAIJ(MPI_COMM_SELF, 3*nDofsPlane, 3*nDofsPlane, 16, NULL, &K);
  MatCreateSeqAIJ(MPI_COMM_SELF, 3*nDofsPlane, 3*nDofsPlane,  1, NULL, &Kii_inv);

  MatZeroEntries(G);
  MatZeroEntries(K);
  MatZeroEntries(Kii_inv);

  MatZeroEntries(P);

  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    for(int plane_i = 0; plane_i < nZloc; plane_i++) {
      plane_j = rank * nZloc + plane_i;
      k_z = (plane_j / 2) / beta;
      if(plane_j < 2) k_z = 1.0;

      for(int elmt_y = 0; elmt_y < elOrd*NELS_Y; elmt_y++) {
        el_i = (elmt_y / elOrd) * NELS_X;

        yi[0] = elmt[el_i]->_ymesh[(elmt_y%elOrd+0)*(elOrd+1)];
        yi[1] = elmt[el_i]->_ymesh[(elmt_y%elOrd+1)*(elOrd+1)];
    
        delta_y = 0.5 * fabs(yi[1] - yi[0]);
        det = alpha * delta_y * beta;

        for(int mode_x = 0; mode_x < nModesX; mode_x++) {
          k_x = (mode_x / 2) / alpha;
          if(mode_x < 2) k_x = 1.0;

          // grad and div matrices
          for(int row = 0; row < 2; row++) {
            row_k[0] = (elmt_y+row-1)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary
            row_k[1] = row_k[0] + nDofsPlane;
            row_k[2] = row_k[1] + nDofsPlane;

            if(row_k[0] < 0) continue;

            for(int col = 0; col < 2; col++) {
              col_k[0] = (elmt_y+col-1)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary

              if(col_k[0] < 0) continue;

              for(int pnt = 0; pnt < 2; pnt++) {
                // grad
                pVals[0] = +dt * det * k_x;
                pVals[1] = -dt * det / delta_y * dNidy(pnt+1, qx[pnt], qy[pnt]) * Ni(pnt+1, qx[pnt], qy[pnt]);
                pVals[2] = +dt * det * k_z;
                MatSetValues(G, 3, row_k, 1, col_k, pVals, ADD_VALUES);

                // div
                pVals[0] *= -1.0;
                pVals[2] *= -1.0;
                MatSetValues(D, 1, col_k, 3, row_k, pVals, ADD_VALUES);
              }
            }
          }

          // vector laplacian matrix
          for(int row = 0; row < 2; row++) {
            row_k[0] = (elmt_y+row-1)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary
            row_k[1] = row_k[0] + nDofsPlane;
            row_k[2] = row_k[1] + nDofsPlane;

            if(row_k[0] < 0) continue;

            for(int col = 0; col < 2; col++) {
              col_k[0] = (elmt_y+col-1)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary
              col_k[1] = col_k[0] + nDofsPlane;
              col_k[2] = col_k[1] + nDofsPlane;

              if(col_k[0] < 0) continue;

              for(int pnt = 0; pnt < 2; pnt++) {
                //pVals[0] = -(dt / kinvis) * det * k_x * k_x;
                //pVals[4] = -(dt / kinvis) * det / delta_y * dNidy(pnt+1, qx[pnt], qy[pnt]) / delta_y * dNidy(pnt+1, qx[pnt], qy[pnt]);
                //pVals[8] = -(dt / kinvis) * det * k_z * k_z;
                pVals[0] = (dt / kinvis) * PCX[mode_x][mode_x];
                pVals[4] = 0.0;
                pVals[8] = 0.0;

                MatSetValues(K, 3, row_k, 3, col_k, pVals, ADD_VALUES);
                MatSetValue(K, row_k[0], row_k[0], det, ADD_VALUES);
                MatSetValue(K, row_k[1], row_k[1], det, ADD_VALUES);
                MatSetValue(K, row_k[2], row_k[2], det, ADD_VALUES);

                MatSetValue(Kii_inv, row_k[0], row_k[0], 1.0/(det+pVals[0]), ADD_VALUES);
                MatSetValue(Kii_inv, row_k[1], row_k[1], 1.0/(det+pVals[4]), ADD_VALUES);
                MatSetValue(Kii_inv, row_k[2], row_k[2], 1.0/(det+pVals[8]), ADD_VALUES);
              }
            }
          }
        }
      }
      MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(  D, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(  G, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(  K, MAT_FINAL_ASSEMBLY);
      MatAssemblyBegin(Kii_inv, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(  Kii_inv, MAT_FINAL_ASSEMBLY);

      // create the Schur complement operator
      MatScale(D, -1.0);
      MatMatMult(Kii_inv, K, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &LK);
      MatMatMult(LK, Kii_inv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &LKR);
      //MatMatMult(D, LKR, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DKinv);
      MatMatMult(D, Kii_inv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DKinv);
      MatMatMult(DKinv, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &S);
      // schur complement should be scaled by -1, however this is ommitted in 
      // order to account for the G^{T}G term in the azimutal direction (-k_z^2)
      MatScale(S, -1.0);

      // [u,u] block
      for(int row_i = 0; row_i < 3*nDofsPlane; row_i++) {
        row_dof = row_i/nDofsPlane;
        MatGetRow(K, row_i, &nCols, &cols, &vals);
        pRow = row_i%nDofsPlane + plane_i*nDofsPlane + lShift[slice_i][row_dof];

        for(int col_i = 0; col_i < nCols; col_i++) {
          col_dof = cols[col_i]/nDofsPlane;
          pCols[col_i] = cols[col_i]%nDofsPlane + plane_i*nDofsPlane + lShift[slice_i][col_dof];
        }
        //MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
MatSetValues(P, 1, &pRow, 1, &pRow, &one, INSERT_VALUES);
        MatRestoreRow(K, row_i, &nCols, &cols, &vals);
      }

      if(NFIELD == 4) {
        // [p,u] block
        for(int row_i = 0; row_i < nDofsPlane; row_i++) {
          MatGetRow(D, row_i, &nCols, &cols, &vals);
          pRow = row_i + plane_i*nDofsPlane + lShift[slice_i][3];
          for(int col_i = 0; col_i < nCols; col_i++) {
            col_j = cols[col_i]%nDofsPlane;
            col_dof = cols[col_i]/nDofsPlane;
            pCols[col_i] = col_j + plane_i*nDofsPlane + lShift[slice_i][col_dof];
          }
          if(plane_i%2 == 0) {
            pRow += nDofsPlane; // imaginary plane for this mode
            for(int col_i = 0; col_i < nCols; col_i++) {
              pVals[col_i] *= -1.0;
            }
          } else {
            pRow -= nDofsPlane; // real plane for this mode
          }
          //MatSetValues(P, 1, &pRow, nCols, pCols, pVals, INSERT_VALUES);
pCols[0] = row_i + plane_i*nDofsPlane + lShift[slice_i][0];
pCols[1] = row_i + plane_i*nDofsPlane + lShift[slice_i][1];
pCols[2] = row_i + plane_i*nDofsPlane + lShift[slice_i][2];
pVals[0] = 1.0;
pVals[1] = 1.0;
pVals[2] = 1.0;
//MatSetValues(P, 1, &pRow, 3, pCols, pVals, INSERT_VALUES);
          MatRestoreRow(D, row_i, &nCols, &cols, &vals);
        }

        // [p,p] block
        for(int row_i = 0; row_i < nDofsPlane; row_i++) {
          MatGetRow(S, row_i, &nCols, &cols, &vals);
          pRow = row_i + plane_i*nDofsPlane + lShift[slice_i][3];
          for(int col_i = 0; col_i < nCols; col_i++) {
            pCols[col_i] = cols[col_i] + plane_i*nDofsPlane + lShift[slice_i][3];
          }
          //MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
//        MatSetValues(P, 1, &pRow, 1, &pRow, &one, INSERT_VALUES);
          MatRestoreRow(S, row_i, &nCols, &cols, &vals);
        }
      }
    }

    // add diagonal entries where required
    if(!Geometry::procID()) {
      for(int row_i = 0; row_i < 3; row_i++) {
        pRow = slice_i * nDofsSlice + NFIELD * nDofsPlane * Geometry::nZ() + row_i;
        MatSetValues(P, 1, &pRow, 1, &pRow, &one, INSERT_VALUES);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  P, MAT_FINAL_ASSEMBLY);

  rpo_set_fieldsplits(snes, is_s, is_u, is_p, nSlice, nDofsPlane, nDofsSlice, lShift);

  MatDestroy(&G);
  MatDestroy(&D);
  MatDestroy(&K);
  MatDestroy(&Kii_inv);
  MatDestroy(&S);
  MatDestroy(&LK);
  MatDestroy(&LKR);
  MatDestroy(&DKinv);
  for(int ii = 0; ii < nModesX; ii++) { delete[] PCX[ii]; }
  delete[] PCX;
}

void rpo_set_fieldsplits(SNES snes, IS* is_s, IS* is_u, IS* is_p, int nSlice, int nDofsPlane, int nDofsSlice, int** lShift) {
  PC pc, pc_i, pc_j;
  KSP ksp, *ksp_i, *ksp_j;
  int n_split, m_split, is, iu, ip, nDofs_u, nDofs_p, nSliceProc, startSlice;
  int *inds_s, *inds_u, *inds_p;

  if(is_s) return;

  is_s = new IS[nSlice];
  is_u = new IS[nSlice];
  is_p = new IS[nSlice];

  SNESGetKSP(snes, &ksp);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCFIELDSPLIT);
  PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);

  nDofs_u = 3 * Geometry::nZ() * nDofsPlane;
  nDofs_p = 3;
  if(NFIELD == 4) nDofs_p += Geometry::nZ() * nDofsPlane;
  nSliceProc = nSlice / Geometry::nProc();
  startSlice = nSliceProc * Geometry::procID();

  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    if(Geometry::procID() == slice_i) {
      ISCreateStride(MPI_COMM_SELF, nDofsSlice, slice_i * nDofsSlice, 1, &is_s[slice_i]);
    } else {
      ISCreateStride(MPI_COMM_SELF, 0, 0, 1, &is_s[slice_i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    PCFieldSplitSetIS(pc, NULL, is_s[slice_i]);
  }
  PCSetUp(pc);

  // setup a schur complement between the velocity dofs and the phase shift dofs
  PCFieldSplitGetSubKSP(pc, &n_split, &ksp_i);
  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    KSPGetPC(ksp_i[slice_i], &pc_i);
    PCSetType(pc_i, PCFIELDSPLIT);

    if(Geometry::procID() == slice_i) {
      ISCreateStride(MPI_COMM_SELF, nDofs_u, 0, 1, &is_u[slice_i]);
      ISCreateStride(MPI_COMM_SELF, nDofs_p, nDofs_u, 1, &is_p[slice_i]);
    } else {
      ISCreateStride(MPI_COMM_SELF, 0, 0, 1, &is_u[slice_i]);
      ISCreateStride(MPI_COMM_SELF, 0, 0, 1, &is_p[slice_i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    PCFieldSplitSetIS(pc_i, "u", is_u[slice_i]);
    PCFieldSplitSetIS(pc_i, "p", is_p[slice_i]);

    PCFieldSplitSetType(pc_i, PC_COMPOSITE_SCHUR);
    PCSetUp(pc_i);
  }

  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
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
