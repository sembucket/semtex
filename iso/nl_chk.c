/*****************************************************************************
 * nl_chk: generate velocity fields composed of products of cos & sin,
 * use nonlinear operator to produce the negative of the nonlinear
 * terms, and add on the analytical nonlinear terms.  Results shoud be
 * zero everywhere.
 *
 * Copyright (C) 1995-1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define SIZE 32

int N, K, FourKon3;

static void test        (void (*)(CVF), void (*)(CVF), CVF, CVF, CVF,
			 CF, CF, const complex*, const complex*);
static void test1_IC    (CVF);
static void test1_error (CVF);
static void test2_IC    (CVF);
static void test2_error (CVF);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine.
 * ------------------------------------------------------------------------- */
{
  CVF           U, G_temp, work1;
  CVF*          G = calloc (1, sizeof (CVF*));
  CF            work2, work3;
  complex       *Wtab, *Stab;

  N = SIZE; K  = N / 2; FourKon3 = (4 * K) / 3;

  allocate (&U, G, 1, &work1, &work2, &work3, &Wtab, &Stab);
  preFFT   (Wtab, K);
  preShift (Stab, N);

  printf ("-- Taylor--Green vortex IC:\n");
  test   (test1_IC, test1_error, U, G[0], work1, work2, work3, Wtab, Stab);

  printf ("-- Second IC:\n");
  test   (test2_IC, test2_error, U, G[0], work1, work2, work3, Wtab, Stab);

  return EXIT_SUCCESS;
}


static void test (void           (*generate) (CVF),
		  void           (*check)    (CVF),
		  CVF            U    ,
		  CVF            NL   ,
		  CVF            work1,
		  CF             work2,
		  CF             work3,
		  const complex* Wtab ,
		  const complex* Stab )
/* ------------------------------------------------------------------------- *
 * Generic test wrapper: pass in functions and storage.
 * ------------------------------------------------------------------------- */
{
  int c;

  /* -- Create analytical velocity field U. */

  generate (U);

  fprintf (stderr, "Maximum velocity components        : %g  %g  %g\n",
	   amaxF (U[1]), amaxF (U[2]), amaxF (U[3]));

  for (c = 1; c <= 3; c++) { rc3DFT (U[c], Wtab, FORWARD); scaleFT (U[c]); }

  fprintf (stderr, "Solution energy                    : %g\n", energyF (U));

  /* -- Compute nonlinear terms NL = d(UiUj)/dxj. */

  nonlinear (U, NL, work2, work3, work1, Wtab, Stab);

  fprintf (stderr, "Maximum nonlinear components       : %g  %g  %g\n",
	   amaxF (NL[1]), amaxF (NL[2]), amaxF (NL[3]));

  fprintf (stderr, "Nonlinear terms' energy            : %g\n", energyF (NL));

  /* -- Transform nonlinear terms to PHYSICAL space. */

  for (c = 1; c <= 3; c++) rc3DFT (NL[c], Wtab, INVERSE);

  /* -- Subtract analytical solution. */

  check (NL);

  fprintf (stderr, "Maximum nonlinear error components : %g  %g  %g\n",
	   amaxF (NL[1]), amaxF (NL[2]), amaxF (NL[3]));

  /* -- Back to FOURIER space, just to compute error energy. */

  for (c = 1; c <= 3; c++) { rc3DFT (NL[c], Wtab, FORWARD); scaleFT (NL[c]); }

  fprintf (stderr, "Error energy                       : %g\n", energyF (NL));
}


static void test1_IC (CVF IC)
/* ------------------------------------------------------------------------- *
 * Generate initial conditions of the 3D Taylor--Green vortex, in
 * PHYSICAL space.
 *
 *   u =  sin x cos y cos z
 *   v = -cos x sin y cos z
 *   w =  0
 * ------------------------------------------------------------------------- */
{
  register real* u = &IC[1][0][0][0].Re;
  register real* v = &IC[2][0][0][0].Re;
  register real* w = &IC[3][0][0][0].Re;
  register int   i, j, k;

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] =  sin(x) * cos(y) * cos(z);
	v[k + N * (j + i * N)] = -cos(x) * sin(y) * cos(z);
	w[k + N * (j + i * N)] =  0.0;
      }
    }
  }
}


static void test1_error (CVF IC)
/* ------------------------------------------------------------------------- *
 * Generate the nonlinear terms in the Navier--Stokes equations for
 * the initial conditions of the 3D Taylor--Green vortex, in PHYSICAL
 * space.  Add them to input, which should be the negative of the
 * nonlinear terms (being generated for the RHS of the N-S equations).
 *
 *   u =  sin x cos y cos z
 *   v = -cos x sin y cos z
 *   w =  0
 *
 *   d(uu)/dx = 2 cos^2(y) sin(x) cos(x)            cos^2(z)
 *   d(uv)/dy = (sin^2(y) - cos^2(y)) cos(x) sin(x) cos^2(z)
 *   d(vu)/dx = (sin^2(x) - cos^2(x)) cos(y) sin(y) cos^2(z)
 *   d(vv)/dy = 2 cos^2(x) cos(y) sin(y)            cos^2(z)
 * ------------------------------------------------------------------------- */
{
  register real* u = &IC[1][0][0][0].Re;
  register real* v = &IC[2][0][0][0].Re;
  register real* w = &IC[3][0][0][0].Re;
  register int   i, j, k;
  real           UUx, UVy, VUx, VVy, cz2;

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	UUx = 2.0 * cos(x) * sin(x) * SQR(cos(y))           * SQR(cos(z));
	UVy = (SQR(sin(y)) - SQR(cos(y))) * cos(x) * sin(x) * SQR(cos(z));
	VUx = (SQR(sin(x)) - SQR(cos(x))) * cos(y) * sin(y) * SQR(cos(z));
	VVy = 2.0 * cos(y) * sin(y) * SQR(cos(x))           * SQR(cos(z));

	u[k + N * (j + i * N)] += (UUx + UVy);
	v[k + N * (j + i * N)] += (VUx + VVy);
      }
    }
  }
}


static void test2_IC (CVF IC)
/* ------------------------------------------------------------------------- *
 *   u = sin x cos y cos z
 *   v = cos x sin y cos z
 *   w = cos x cos y sin z
 * ------------------------------------------------------------------------- */
{
  register real* u = &IC[1][0][0][0].Re;
  register real* v = &IC[2][0][0][0].Re;
  register real* w = &IC[3][0][0][0].Re;
  register int   i, j, k;

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] = sin(x) * cos(y) * cos(z);
	v[k + N * (j + i * N)] = cos(x) * sin(y) * cos(z);
	w[k + N * (j + i * N)] = cos(x) * cos(y) * sin(z);
      }
    }
  }
}


static void test2_error (CVF IC)
/* ------------------------------------------------------------------------- *
 * For:
 *   u = sin x cos y cos z
 *   v = cos x sin y cos z
 *   w = cos x cos y sin z
 * IC[1]:
 *   d(uu)/dx = 2 cos(x) sin(x) cos^2(y)            cos^2(z)
 *   d(uv)/dy = cos(x) sin(x) (cos^2(y) - sin^2(y)) cos^2(z)
 *   d(uw)/dz = cos(x) sin(x) cos^2(y) (cos^2(z) - sin^2(z))
 * IC[2]:
 *   d(vu)/dx = (cos^2(x) - sin^2(x)) cos(y) sin(y) cos^2(z)
 *   d(vv)/dy = 2 cos^2(x) cos(y) sin(y)            cos^2(z)
 *   d(vw)/dz = cos^2(x) cos(y) sin(y) (cos^2(z) - sin^2(z))
 * IC[3]:
 *   d(wu)/dx = (cos^2(x) - sin^2(x)) cos^2(y) cos(z) sin(z)
 *   d(wv)/dy = cos^2(x) (cos^2(y) - sin^2(y)) cos(z) sin(z)
 *   d(ww)/dz = cos^2(x) cos^2(y) (cos^2(z) - sin^2(z))
 * ------------------------------------------------------------------------- */
{
  register real* u = &IC[1][0][0][0].Re;
  register real* v = &IC[2][0][0][0].Re;
  register real* w = &IC[3][0][0][0].Re;
  register int   i, j, k;
  real           UUx, UVy, UWz, VUx, VVy, VWz, WUx, WVy, WWz;

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	UUx = 2.0 * cos(x) * sin(x) * SQR(cos(y))           * SQR(cos(z));
	UVy = cos(x) * sin(x) * (SQR(cos(y)) - SQR(sin(y))) * SQR(cos(z));
	UWz = cos(x) * sin(x) * SQR(cos(y)) * (SQR(cos(z)) - SQR(sin(z)));

	VUx = (SQR(cos(x)) - SQR(sin(x))) * cos(y) * sin(y) * SQR(cos(z));
	VVy = 2.0 * SQR(cos(x)) * cos(y) * sin(y)           * SQR(cos(z));
	VWz = SQR(cos(x)) * cos(y) * sin(y) * (SQR(cos(z)) - SQR(sin(z)));
	
	WUx = (SQR(cos(x)) - SQR(sin(x))) * SQR(cos(y)) * cos(z) * sin(z);
	WVy = SQR(cos(x)) * (SQR(cos(y)) - SQR(sin(y))) * cos(z) * sin(z);
	WWz = 2.0 * SQR(cos(x)) * SQR(cos(y)) * cos(z) * sin(z);

	u[k + N * (j + i * N)] += (UUx + UVy + UWz);
	v[k + N * (j + i * N)] += (VUx + VVy + VWz);
	w[k + N * (j + i * N)] += (WUx + WVy + WWz);
      }
    }
  }
}
