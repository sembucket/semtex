#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// gcc -g test_eig.c -lm -I/usr/local/acml4.4.0/gfortran64/include/ -L/usr/local/acml4.4.0/gfortran64/lib/libacml.so
// g++ -g test_eig.cpp -lm -L/home/dlee0035/soft/lib/libopenblas.so -I/home/dlee0035/soft/include/

#include "cfemdef.h"
#include "blas.h"
#include "lapack.h"

//#include <cblas.h>
//#include <lapacke.h>
//#include <acml.h>

#define F77name(x) x ## _
//#define F77name(x) x ## 

#define NX 64
#define LX 1.0
#define _K 1.0
#define DT 0.1

double** LinearOp() {
  int ii, jj, ip1, im1;
  double dx = LX / NX;
  double fac = /*DT **/ _K / dx / dx;
  double** A = (double**)malloc(NX*sizeof(double*));

  for(ii = 0; ii < NX; ii++) {
    A[ii] = (double*)malloc(NX*sizeof(double));
    for(jj = 0; jj < NX; jj++) A[ii][jj] = 0.0;

    ip1 = (ii == NX-1) ? 0    : ii+1;
    im1 = (ii ==    0) ? NX-1 : ii-1;

    A[ii][im1] = -1.0*fac;
    A[ii][ii]  = +2.0*fac;
    A[ii][ip1] = -1.0*fac;
  }
  return A;
}

double* RHS() {
  int ii;
  double dx = LX / NX, xi;
  double ki = 2.0 * M_PI;
  double* b = (double*)malloc(NX*sizeof(double));

  for(ii = 0; ii < NX; ii++) {
    xi = ii*dx;
    b[ii] = -ki * ki * cos(ki * xi);
  }
  return b;
}

void residual(double** A, double* x, double* b, double* f) {
  int ii, jj;

  for(ii = 0; ii < NX; ii++) {
    f[ii] = -b[ii];
    for(jj = 0; jj < NX; jj++) {
      f[ii] += A[ii][jj]*x[jj];
    }
  }
}

double norm(int nl, double* v) {
  int ii;
  double norm_l = 0.0;

  for(ii = 0; ii < nl; ii++) {
    norm_l += v[ii]*v[ii];
  }
  return sqrt(norm_l);
}

double dot(int nl, double* u, double* v) {
  int ii;
  double dot_l = 0.0;

  for(ii = 0; ii < nl; ii++) {
    dot_l += u[ii]*v[ii];
  }
  return dot_l;
}

double normvec(int nl, double* v, double* vn) {
  int ii;
  double norm_q = norm(nl, v);

  for(ii = 0; ii < nl; ii++) {
    vn[ii] = v[ii] / norm_q;
  }
  return norm_q;
}

void mgs(int nk, int nl, double** AK, double** QK) {
  int ii, ki, kj;
  double dot_ij, dot_jj;

  normvec(nl, AK[0], QK[0]);

  for(ki = 1; ki < nk; ki++) {
    for(ii = 0; ii < nl; ii++) QK[ki][ii] = AK[ki][ii];

    for(kj = 0; kj < ki; kj++) {
      dot_ij = dot(nl, QK[ki], QK[kj]);
      dot_jj = dot(nl, QK[kj], QK[kj]); // should be 1!

      for(ii = 0; ii < nl; ii++) QK[ki][ii] -= (dot_ij / dot_jj) * QK[kj][ii];
    }

    normvec(nl, QK[ki], QK[ki]);
  }
}

// Generate an array of orthonormalised vectors, QK, and the Hessenberg, HK
void orthonorm(int nk, int nl, double** AK, double** QK, double** HK) {
  int ii, ki, kj;
  double* v = (double*)malloc(nl*sizeof(double));

  normvec(nl, AK[0], QK[0]);

  for(ki = 0; ki < nk; ki++) {
    for(ii = 0; ii < nl; ii++) v[ii] = AK[ki+1][ii];

    for(kj = 0; kj < ki+1; kj++) {
      HK[kj][ki] = dot(nl, QK[kj], v);

      for(ii = 0; ii < nl; ii++) v[ii] -= HK[kj][ki] * QK[kj][ii];
    }

    HK[ki+1][ki] = normvec(nl, v, QK[ki+1]);
  }

  free(v);
}

void rpo_solve() {
  int ier, lwork;
  int ki, kj, ji, nk = NX/*20*/;
  //int kdim, one = 1;
  double **QK, **HK, *AK, *x_o, *x_k, *f_o, *f_k, *wr, *wi, *ZK, *rwork, norm_q, eps = 1.0e-8;
  double** OP = LinearOp();
  double* rhs = RHS();
  char _N = 'N', _V = 'V';

  /* RUN THE ARNOLDI ITERATION */

  // 0: initialisation
  QK = (double**)malloc((nk+1)*sizeof(double*));
  HK = (double**)malloc((nk+1)*sizeof(double*));
  for(ki = 0; ki < nk+1; ki++) {
    QK[ki] = (double*)malloc(NX*sizeof(double));
    for(ji = 0; ji < NX; ji++) QK[ki][ji] = 0.0;
    HK[ki] = (double*)malloc(nk*sizeof(double));
    for(kj = 0; kj < nk; kj++) HK[ki][kj] = 0.0;
  }
  x_o = (double*)malloc(NX*sizeof(double));
  for(ji = 0; ji < NX; ji++) x_o[ji] = 0.0;
  x_k = (double*)malloc(NX*sizeof(double));
  f_o = (double*)malloc(NX*sizeof(double));
  f_k = (double*)malloc(NX*sizeof(double));

  // 1: first vector
  residual(OP, x_o, rhs, f_o);
  normvec(NX, f_o, QK[0]);

  // 2: iterate
  ki = 0;
  do {
    // x + h q_k, h << 1
    for(ji = 0; ji < NX; ji++) x_k[ji] = x_o[ji] + eps * QK[ki][ji];

    // next krylov vector: (f(x + hq) - f(x))/h
    residual(OP, x_k, rhs, f_k);
    for(ji = 0; ji < NX; ji++) f_k[ji] = (f_k[ji] - f_o[ji]) / eps;

    // orthonormalise (modified gram-schmidt)
    for(ji = 0; ji < NX; ji++) QK[ki+1][ji] = f_k[ji];

    for(kj = 0; kj < ki+1; kj++) {
      HK[kj][ki] = dot(NX, QK[kj], QK[ki+1]);

      for(ji = 0; ji < NX; ji++) QK[ki+1][ji] -= HK[kj][ki] * QK[kj][ji];
    }

    HK[ki+1][ki] = normvec(NX, QK[ki+1], QK[ki+1]);

    if(HK[ki+1][ki] < 1.0e-10) break;
    ki++;
  } while(ki < nk);
  const int one = 1;
  const int kdim = ki;

  // 3: eigenvalues
  AK = (double*)malloc(kdim*kdim*sizeof(double));
  for(ki = 0; ki < kdim; ki++) {
    for(kj = 0; kj < kdim; kj++) {
      AK[ki*kdim+kj] = HK[ki][kj];
    }
  }
  ZK = (double*)malloc(kdim*kdim*sizeof(double));
  wr = (double*)malloc(kdim*sizeof(double));
  wi = (double*)malloc(kdim*sizeof(double));
  lwork = 4*kdim;
  rwork = (double*)malloc(lwork*sizeof(double));
  
  //F77name(dgeev) (&_N,&_V,&kdim,AK,&kdim,wr,wi,0,&one,ZK,&kdim,rwork,&lwork,&ier);
  F77name(dgeev) (&_N,&_V,kdim,AK,kdim,wr,wi,0,one,ZK,kdim,rwork,lwork,ier);
  //dgeev(_N,_V,kdim,AK,kdim,wr,wi,0,one,ZK,kdim,rwork,lwork,ier);
  if(ier) printf("Hessenberg eigenvalue error!\n");

  printf("eigenvalues: \n");
  for(ki = 0; ki < kdim; ki++) {
    printf("\t%12.10e\t+\t%12.10ei\n", wr[ki],  wi[ki]);
  }

  // cleanup
  for(ki = 0; ki < nk+1; ki++) free(QK[ki]);
  free(QK);
  free(AK);
  free(ZK);
  free(x_o);
  free(x_k);
  free(f_o);
  free(f_k);
  free(wr);
  free(wi);
  free(rwork);
  for(ji = 0; ji < NX; ji++) free(OP[ji]);
  free(OP);
  free(rhs);
}

int main (int argc, char** argv) {
  rpo_solve();

  0;
}
