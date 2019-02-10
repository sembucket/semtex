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
  int il, ir;
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

  for(int ii = 0; ii < nr; ii++) { delete[] F2R[ii]; }     delete[] F2R;
  for(int ii = 0; ii < nr; ii++) { delete[] LAP[ii]; }     delete[] LAP;

  return PRECON;
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
  double** PCX = real_space_pc(2.0*M_PI, nModesX, nModesX);
  double **PCZ1, **PCZ2;

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

        PCZ1 = real_space_pc(2.0*M_PI*yi[0], Geometry::nZ(), Geometry::nZ());
        PCZ2 = real_space_pc(2.0*M_PI*yi[1], Geometry::nZ(), Geometry::nZ());

        for(int mode_x = 0; mode_x < nModesX; mode_x++) {
          k_x = (mode_x / 2) / alpha;
          if(mode_x < 2) k_x = 1.0;

          // grad and div matrices
          for(int row = 0; row < 2; row++) {
            row_k[0] = (elmt_y+row)*nModesX + mode_x;
            row_k[1] = row_k[0] + nDofsPlane;
            row_k[2] = row_k[1] + nDofsPlane;

            if(row_k[0] == elOrd*NELS_Y) continue; // outer boundary (dirichlet bcs)

            for(int col = 0; col < 2; col++) {
              col_k[0] = (elmt_y+col)*nModesX + mode_x;

              if(col_k[0] == elOrd*NELS_Y) continue; // outer boundary (dirichlet bcs)

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
            row_k[0] = (elmt_y+row)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary
            row_k[1] = row_k[0] + nDofsPlane;
            row_k[2] = row_k[1] + nDofsPlane;

            if(elmt_y+row == elOrd*NELS_Y) continue; // outer boundary (dirichlet bcs)

            for(int col = 0; col < 2; col++) {
              col_k[0] = (elmt_y+col)*nModesX + mode_x; // elmt_y=0 along dirichlet boundary
              col_k[1] = col_k[0] + nDofsPlane;
              col_k[2] = col_k[1] + nDofsPlane;

              if(elmt_y+col == elOrd*NELS_Y) continue; // outer boundary (dirichlet bcs)

              for(int pnt = 0; pnt < 2; pnt++) {
                //pVals[0] = -(dt * kinvis) * det * k_x * k_x;
                //pVals[4] = -(dt * kinvis) * det / delta_y * dNidy(pnt+1, qx[pnt], qy[pnt]) / delta_y * dNidy(pnt+1, qx[pnt], qy[pnt]);
                //pVals[8] = -(dt * kinvis) * det * k_z * k_z;
                pVals[0] = (dt * kinvis) * PCX[mode_x][mode_x];
                pVals[4] = 0.0;
                if(row == 0 && col == 0)
                  pVals[8] = (dt * kinvis) * PCZ1[plane_j][plane_j];
                if(row == 1 && col == 1)
                  pVals[8] = (dt * kinvis) * PCZ2[plane_j][plane_j];
                //pVals[8] = 0.0;

                MatSetValues(K, 3, row_k, 3, col_k, pVals, ADD_VALUES);
                MatSetValue(K, row_k[0], row_k[0], det, ADD_VALUES);
                MatSetValue(K, row_k[1], row_k[1], det, ADD_VALUES);
                MatSetValue(K, row_k[2], row_k[2], det, ADD_VALUES);

                MatSetValue(Kii_inv, row_k[0], row_k[0], 1.0/(det+pVals[0]), ADD_VALUES);
                MatSetValue(Kii_inv, row_k[1], row_k[1], 1.0/(det+pVals[4]), ADD_VALUES);
                MatSetValue(Kii_inv, row_k[2], row_k[2], 1.0/(det+pVals[8]), ADD_VALUES);
              } // pnt
            } // col
          } // row
        } // mode_x

        for(int ii = 0; ii < Geometry::nZ(); ii++) {
          delete[] PCZ1[ii];
          delete[] PCZ2[ii];
        }
        delete[] PCZ1;
        delete[] PCZ2;
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
        MatSetValues(P, 1, &pRow, nCols, pCols, vals, INSERT_VALUES);
//MatSetValues(P, 1, &pRow, 1, &pRow, &one, INSERT_VALUES);
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
