/*
 * Convert an ISO field file into a TECPLOT file
 *
 * usage: iso2tec input.fld > output.tec
 *
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "globals.h"

real magnitude (CVF, int, int, int, int);

int main (int argc, char *argv[])
{
  CVF U;
  CVF Q;
  CF          Work;
  header               U_info;
  int*              Dimension;
  complex*              Wtab;
  int                  N, Npts;
  FILE                 *fp;
  real                *u, *v, *w;
  register int         c, i, j, k;

  fp = efopen (argv[1], "r");
  if (fread(&U_info, 1, sizeof(header), fp) != sizeof(header))
    fprintf (stderr, "Can't read header from input file\n");
  
  /* Check for a '\n' in the title string */

  { char *p = strchr(U_info.Title,'\n'); if (p) *p = '\0'; }

  fprintf (stderr, "Title:  %s\n", U_info.Title);
  fprintf (stderr, "Re   :  %f\n", 1./U_info.K_Visc);
  fprintf (stderr, "N    :  %d\n", U_info.N_Grid);
  fprintf (stderr, "dt   :  %f\n", U_info.Delta_T);
  fprintf (stderr, "step :  %d\n", U_info.N_Step);
  fprintf (stderr, "time :  %f\n", U_info.N_Step * U_info.Delta_T);
  

  /* Set up the problem size */

  Dimension    =  ivect(1, 3);
  Dimension[1] = (N = U_info.N_Grid);
  Dimension[2] =  N;
  Dimension[3] =  N / 2;
  Npts         =  Dimension[1] * Dimension[2] * Dimension[3];

  /* Allocate storage for the solution */

  Wtab = cvect (0, Dimension[3]-1);
  cfield (Dimension, &U);
  cfield (Dimension, &Q);
  cbox   (0, Dimension[1]-1, 0, Dimension[2]-1, 0, Dimension[3]-1, &Work);

  read_components (fp, U, Npts);
  fclose          (fp);

  /* Compute vorticity */

  curl (U, Q, Work, Dimension);

  /* Get ready to transform to PHYSICAL space */

  preFFT (Dimension[3], Wtab);
  for (c = 1; c <= 3; c++)
    rc3DFT (U[c], Dimension, Wtab, INVERSE);
  for (c = 1; c <= 3; c++)
    rc3DFT (Q[c], Dimension, Wtab, INVERSE);

  u = (real*) & U[1][0][0][0].Re;   /* Set up fast pointers */
  v = (real*) & U[2][0][0][0].Re;
  w = (real*) & U[3][0][0][0].Re;

  /* TECPLOT header */

  printf ("TITLE = %s\n", U_info.Title);
  printf ("VARIABLES = x y z u v w q\n");
  printf ("ZONE T = \"BOX\", I=%d, J=%d, K=%d\n", N, N, N);

  /* Output the solution and grid */

  for (k = 0; k < N; k++) {
    const real z = 2.*M_PI * k / N;
    for (j = 0; j < N; j++) {
      const real y = 2.*M_PI * j / N;
      for (i = 0; i < N; i++) {
	const real x = 2.*M_PI * i / N;

	printf ("%g %g %g ", x, y, z);
	printf ("%g %g %g ",
		u[ k + N * (j + i * N) ],
		v[ k + N * (j + i * N) ],
		w[ k + N * (j + i * N) ]);
	printf ("%g\n", magnitude (Q, i, j, k, N));
      }
    }
  }

  return 0;
}


real magnitude (CVF Q, int i, int j, int k, int N)
{
  real *u = (real*) &Q[1][0][0][0].Re;
  real *v = (real*) &Q[2][0][0][0].Re;
  real *w = (real*) &Q[3][0][0][0].Re;

  const int p = k + N * (j + N * i);

  return u[p]*u[p] + v[p]*v[p] * w[p]*w[p];
}


