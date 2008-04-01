/*****************************************************************************
 * LU: L-U factorization & solution routines.  Numerical Recipes 2e.
 *****************************************************************************/

static char RCSid[] = "$Id$";


#include <stdio.h>
#include <math.h>
#include "panel.h"

    
void ludcmp(double **a, int n, int *indx, double *d)
/* ========================================================================= *
 * Given a matrix a[1..n][1..n] this routine replaces it by the LU decomp-   *
 * osition of a row-wise permutation of itself.  a & n are input.  a is      *
 * output, arranged (as shown in the book!); indx[1..n] is an output         *
 * vector which records the row permutation effected by partial pivoting;    *
 * d is output as +/-1 depending on whether the number of row interchanges   *
 * was even or odd, respectively.                                            *
 * ========================================================================= */
{
  int     i, imax, j, k;
  double  big, dum, sum, temp;
  double *vv;


  vv = dvector(1, n);
  *d = 1.0;
  for (i=1; i<=n; i++) {
    big = 0.0;
    for (j=1; j<=n; j++)
      if ((temp=fabs(a[i][j])) > big) big = temp;
    if (big == 0.0) message("ludcmp()", "Singular matrix", ERROR);
    vv[i] = 1.0 / big;
  }

  for (j=1; j<=n; j++) {
    for (i=1; i<j; i++) {
      sum = a[i][j];
      for (k=1; k<i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }

    big = 0.0;
    for (i=j; i<=n; i++) {
      sum = a[i][j];
      for (k=1; k<j; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big  = dum;
	imax = i;
      }
    }

    if (j != imax) {
      for (k=1; k<=n; k++) {
	dum = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }

    indx[j] = imax;
    if (a[j][j] == 0.0) a[j][j] = EPSm20;
    if (j != n) {
      dum = 1.0 / (a[j][j]);
      for (i=j+1; i<=n; i++) a[i][j] *= dum;
    }
  }

  free_dvector(vv, 1);
}





void lubksb(double **a, int n, int *indx, double *b)
/* ========================================================================= *
 * Solves the set of n linear equations Ax = b.  Here a[1..n][1..n] is       *
 * input, not as the matrix A, but as its row-permuted LU decomposition,     *
 * previously performed by ludcmp.  The permutation information is stored    *
 * in indx[1..n].  b[1..n] is input as the RHS vector and is returned with   *
 * the solution vector.                                                      *
 * ========================================================================= */
{
  int    i, ii=0, ip, j;
  double sum;
  
  
  for (i=1; i<=n; i++) {
    ip    = indx[i];
    sum   = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j=ii; j<=i-1; j++) sum -= a[i][j]*b[j];
    else if (sum) ii = i;
    b[i] = sum;
  }

  for (i=n; i>=1; i--) {
    sum = b[i];
    for (j=i+1; j<=n; j++) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i];
  }
}
