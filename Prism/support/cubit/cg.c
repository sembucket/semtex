
/* CG...see the header file for a description */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cg.h"

/* BLAS-type operations for arrays */

static double dot (int n, const double *x, const double *y) {
  double sum = 0.;
  int i;
  for(i = 0; i < n; i++)
    sum += x[i]*y[i];
  return sum;
}

static void copy (int n, const double *x, double *y) {
  int i;
  for(i = 0; i < n; i++)
    y[i] = x[i];
}

static void scal (int n, double d, double *x) {
  int i;
  for(i = 0; i < n; i++)
    x[i] *= d;
}

static void axpy (int n, double alpha, const double *x, double *y) {
  int i;
  for(i = 0; i < n; i++)
    y[i] += alpha*x[i];
}

/* ------------------------------------------------------------------------- */

#define NORM(n,x)      (*norm)(n,x)
#define M_SOLVE(n,x,y) (*M_solve)(n,x,y)
#define A_APPLY(n,x,y) (*A_apply)(n,x,y)

int CG (int n, double *x, const double *b, int *max_iter, double *tol,
	double (*norm)(), double (*dot)(), void (*M_solve)(), 
	void (*A_apply)())
{
  double resid;
  int i, iter;

  double alpha;
  double beta;
  double rho;
  double rho_1;

  double *p = (double*) malloc(n*sizeof(double));
  double *z = (double*) malloc(n*sizeof(double));
  double *q = (double*) malloc(n*sizeof(double));
  double *r = (double*) malloc(n*sizeof(double));

  double normb = NORM(n,b);
  if (normb == 0.0)
    normb = 1.0;

  
  A_APPLY(n, x, r);
  for (i = 0; i < n; i++)
    r[i] = b[i] - r[i];
  
  if ((resid = NORM(n,r)/normb) <= *tol) {
    *tol = resid;
    *max_iter = 0;
    goto cleanup;
  }

#ifndef NDEBUG
  printf ("Conjugate Gradient Iteration:\n");
  printf ("normb = %g, residual = %g, tolerance = %g\n", normb, resid, *tol);
  printf ("\titer    alpha       beta      resid\n");
#endif

  for (iter = 1; iter <= *max_iter; iter++) {
    M_SOLVE(n,r,z);
    rho = dot(n, r, z);

    if (iter == 1)
      copy (n, z, p);
    else {
      beta = rho / rho_1;
      for (i = 0; i < n; i++)
	p[i] = z[i] + beta * p[i];
    }

    A_APPLY(n,p,q);
    alpha = rho / dot(n, p, q);

    axpy (n, alpha, p, x);
    axpy (n,-alpha, q, r);

    if ((resid = NORM(n,r)/normb) <= *tol) {
      *tol = resid;
      *max_iter = iter;
    }

    if ((iter % 10) == 0)
      printf ("\t%4d   %#8.6g   %#8.6g   %#8.6g\n", iter, alpha, beta, resid);

    rho_1 = rho;
  }

  printf ("\t%4d   %#8.6g   %#8.6g   %#8.6g\n", iter, alpha, beta, resid);

cleanup:
  free (p);
  free (z);
  free (q);
  free (r);

  return *tol != resid;
}


#undef NORM
#undef M_SOLVE
#undef A_APPLY
