/*****************************************************************************
 * Convert an ISO field file into a TECPLOT file
 *
 * usage: iso2tec input.fld > output.tec
 *
 * : iso2tec.c,v 2.1 1995/11/08 00:42:04 hmb Exp hmb $
 *****************************************************************************/

#include "iso.h"

real magnitude (CVF, int, int, int, int);


int main (int argc, char *argv[])
{
  CVF            U;
  CVF            Q;
  CF             Work;
  Param*         Info = (Param*) calloc (1, sizeof (Param));
  int*           Dim;
  complex*       Wtab;
  int            N, Npts;
  FILE*          fp;
  real          *u, *v, *w;
  register int   c, i, j, k;

  fp = efopen (argv[1], "r");
  readParam   (fp, Info);
  printParam  (stderr, Info, ": iso2tec.c,v $", ": 2.1 $");

  /* -- Set up the problem size */

  Dim    =  ivector (1, 3);
  Dim[1] = (N = Info -> modes);
  Dim[2] =  N;
  Dim[3] =  N / 2;
  Npts   =  Dim[1] * Dim[2] * Dim[3];

  /* -- Allocate storage for the solution */

  Wtab = cvector (0, Dim[3]-1);
  cfield (Dim, &U);
  cfield (Dim, &Q);
  cbox   (0, Dim[1]-1, 0, Dim[2]-1, 0, Dim[3]-1, &Work);

  readCVF (fp, U, Dim);
  fclose  (fp);

  /* -- Compute vorticity. */

  curl (U, Q, Work, Dim);

  /* -- Transform to PHYSICAL space. */

  preFFT (Wtab, Dim[3]);
  for (c = 1; c <= 3; c++) {
    rc3DFT (U[c], Dim, Wtab, INVERSE);
    rc3DFT (Q[c], Dim, Wtab, INVERSE);
  }

  u = &U[1][0][0][0].Re;
  v = &U[2][0][0][0].Re;
  w = &U[3][0][0][0].Re;

  /* -- TECPLOT header. */

  printf ("TITLE = ISO FIELD FILE\n");
  printf ("VARIABLES = x y z u v w q\n");
  printf ("ZONE T = \"BOX\", I=%d, J=%d, K=%d\n", N, N, N);

  /* -- Output the solution and grid. */

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

  return EXIT_SUCCESS;
}


real magnitude (CVF Q, int i, int j, int k, int N)
{
  real *u = &Q[1][0][0][0].Re;
  real *v = &Q[2][0][0][0].Re;
  real *w = &Q[3][0][0][0].Re;

  const int p = k + N * (j + N * i);

  return u[p]*u[p] + v[p]*v[p] + w[p]*w[p];
}


