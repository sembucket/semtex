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

#define XMIN 0.0
#define XMAX (2.0*M_PI)
#define YMIN 0.0
#define YMAX 0.5
#define NELS_X 30
#define NELS_Y 7

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
}

// Schur complement preconditioning for incompressible Navier-Stokes
// Reference:
//   Griffith (2009) "An accurate and efficient method for the incompressible Navier-Stokes
//   equations using the projection method as a preconditioner" J. Comp. Phys. 228, 7565-7595
void build_preconditioner(int nSlice, int nDofsSlice, int nDofsPlane, int localSize, int localShift, vector<Element*> elmt, Mat P) {
  int elOrd = Geometry::nP() - 1;
  int nx = elOrd*NELS_X;
  int ny = elOrd*NELS_Y;
  int ii, jj, kk;
  int row_i, row_j, col_i, col_j;
  int row_x, row_y, col_x, col_y;
  int row_k[2], col_k[2];
  int nCols;
  int offset;
  const int* cols;
  const double* vals;
  int pRow, pCols[99];
  double xi[4], yi[4], det;
  double qx[] = {-1.0, +1.0, +1.0, -1.0};
  double qy[] = {-1.0, -1.0, +1.0, +1.0};
  double dx[4];
  double delta_x, delta_y;
  double kinvis = Femlib::value("KINVIS");
  double dt = Femlib::value("D_T");
  Mat G, K, S;
  Mat Kii_inv, D, LK, LKR, DKinv;

  MatCreateSeqAIJ(MPI_COMM_SELF, 2*nx*ny, 1*nx*ny, 16, NULL, &G);
  MatCreateSeqAIJ(MPI_COMM_SELF, 2*nx*ny, 2*nx*ny, 64, NULL, &K);
  MatCreateSeqAIJ(MPI_COMM_SELF, 2*nx*ny, 2*nx*ny,  1, NULL, &Kii_inv);

  MatZeroEntries(G);
  MatZeroEntries(K);
  MatZeroEntries(Kii_inv);

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
      det = delta_x * delta_y;

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

        row_k[0] = 2*(row_j*nx + row_i) + 0;
        row_k[1] = 2*(row_j*nx + row_i) + 1;

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
            dx[0] = -dt * delta_x * dNidx(pnt, qx[pnt], qy[pnt]) * Ni(pnt, qx[pnt], qy[pnt]);
            dx[1] = -dt * delta_y * dNidy(pnt, qx[pnt], qy[pnt]) * Ni(pnt, qx[pnt], qy[pnt]);

            MatSetValues(G, 2, row_k, 1, col_k, dx, ADD_VALUES);
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

        row_k[0] = 2*(row_j*nx + row_i) + 0;
        row_k[1] = 2*(row_j*nx + row_i) + 1;

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

          col_k[0] = 2*(col_j*nx + col_i) + 0;
          col_k[1] = 2*(col_j*nx + col_i) + 1;

          for(int pnt = 0; pnt < 4; pnt++) {
            dx[0]  = -(dt / kinvis) * det * delta_x * dNidx(pnt, qx[pnt], qy[pnt]) * delta_x * dNidx(pnt, qx[pnt], qy[pnt]);
            dx[1]  = -(dt / kinvis) * det * delta_x * dNidx(pnt, qx[pnt], qy[pnt]) * delta_y * dNidy(pnt, qx[pnt], qy[pnt]);
            dx[2]  = -(dt / kinvis) * det * delta_y * dNidy(pnt, qx[pnt], qy[pnt]) * delta_x * dNidx(pnt, qx[pnt], qy[pnt]);
            dx[3]  = -(dt / kinvis) * det * delta_y * dNidy(pnt, qx[pnt], qy[pnt]) * delta_y * dNidy(pnt, qx[pnt], qy[pnt]);

            MatSetValues(K, 2, row_k, 2, col_k, dx, ADD_VALUES);

            MatSetValue(K, row_k[0], row_k[0], det, ADD_VALUES);
            MatSetValue(K, col_k[1], col_k[1], det, ADD_VALUES);

            MatSetValue(Kii_inv, row_k[0], row_k[0], 1.0/(det+dx[0]), ADD_VALUES);
            MatSetValue(Kii_inv, col_k[1], col_k[1], 1.0/(det+dx[3]), ADD_VALUES);
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
  MatScale(S, -1.0);

  // assemble operators into the preconditioner
  MatZeroEntries(P);
  {
    int m, n;
    Vec v;
    MatGetSize(P, &m, &n);
    VecCreateMPI(MPI_COMM_WORLD, localSize, m, &v);
    VecSet(v, 1.0);
    MatDiagonalSet(P, v, INSERT_VALUES);
    VecDestroy(&v);
  }

  for(int slice_i = 0; slice_i < nSlice; slice_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
      offset = slice_i * nDofsSlice + plane_i * nDofsPlane + localShift;

      // [u,u] block
      for(int row_i = 0; row_i < 2*nx*ny; row_i++) {
        MatGetRow(K, row_i, &nCols, &cols, &vals);
        pRow = row_i + offset;
        for(int col_i = 0; col_i < nCols; col_i++) {
          pCols[col_i] = cols[col_i] + offset;
        }
        MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
        MatRestoreRow(K, row_i, &nCols, &cols, &vals);
      }

      // [u,p] block
      for(int row_i = 0; row_i < 2*nx*ny; row_i++) {
        MatGetRow(G, row_i, &nCols, &cols, &vals);
        pRow = row_i + offset;
        for(int col_i = 0; col_i < nCols; col_i++) {
          pCols[col_i] = cols[col_i] + 2*nx*ny + offset;
        }
        MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
        MatRestoreRow(G, row_i, &nCols, &cols, &vals);
      }

      // [p,p] block
      for(int row_i = 0; row_i < nx*ny; row_i++) {
        MatGetRow(S, row_i, &nCols, &cols, &vals);
        pRow = row_i + 2*nx*ny + offset;
        for(int col_i = 0; col_i < nCols; col_i++) {
          pCols[col_i] = cols[col_i] + 2*nx*ny + offset;
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
