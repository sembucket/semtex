/*****************************************************************************
 * DFT_chk.c: exercise 3D real--complex FFT routines.
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 * 
 * $Id$
 *****************************************************************************/

#include "iso.h"

int N, K, FourKon3;

#define SIZE 8

static char prog[] = "DFT_chk";


int main()
/* ------------------------------------------------------------------------- *
 * Test real DFT of SIZE*SIZE*SIZE data.
 * ------------------------------------------------------------------------- */
{
  int       c, i, j, k, Npts;
  CF        U, V;
  real*     u;
  real*     v;
  int*      Dim;
  complex*  Wtab;
  char      s[STR_MAX];
  int       seed = 1;

  /* -- Allocation. */

  N    = SIZE;
  K    = SIZE / 2;
  Npts = N * N * K;
  u    = cbox (0, N-1, 0, N-1, 0, K-1, &U);
  v    = cbox (0, N-1, 0, N-1, 0, K-1, &V);

  Wtab = cvector (0, K-1);
  preFFT  (Wtab, K);
  sprintf (s, "Dimensions: %d, %d, %d", N, N, K);
  message (prog, s, REMARK);

  message (prog, "Checking constant (1) ----------------------------", REMARK);

  zeroF (U);
  U[ 0][ 0][ 0].Re = 1.0;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= 1.0;
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x) ----------------------------------", REMARK);

  zeroF (U);
  U[  1][ 0][ 0].Re = 0.5;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(x);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking sin(x) ----------------------------------", REMARK);

  zeroF (U);
  U[  1][ 0][ 0].Im = -0.5;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= sin(x);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(y) ----------------------------------", REMARK);

  zeroF (U);
  U[ 0][  1][ 0].Re = 0.5;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(y);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking sin(y) ----------------------------------", REMARK);

  zeroF (U);
  U[ 0][ 1][ 0].Im =  -0.5;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= sin(y);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);
 
  message (prog, "Checking cos(z) ----------------------------------", REMARK);

  zeroF (U);
  U[ 0][ 0][ 1].Re = 0.5;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(z);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);
 
  message (prog, "Checking sin(z) ----------------------------------", REMARK);

  zeroF (U);
  U[ 0][ 0][ 1].Im = -0.5;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
 
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= sin(z);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)cos(y) ----------------------------", REMARK);

  zeroF (U);
  U[  1][  1][  0].Re = 0.25;
  U[N-1][  1][  0].Re = 0.25;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
 
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(x) * cos(y);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)sin(y) ----------------------------", REMARK);

  zeroF (U);
  U[  1][  1][  0].Im = -0.25;
  U[N-1][  1][  0].Im = -0.25;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
 
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(x) * sin(y);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)cos(z) ----------------------------", REMARK);

  zeroF (U);
  U[  1][  0][  1].Re = 0.25;
  U[N-1][  0][  1].Re = 0.25;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
 
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(x) * cos(z);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(y)cos(z) ----------------------------", REMARK);

  zeroF (U);
  U[  0][  1][  1].Re = 0.25;
  U[  0][N-1][  1].Re = 0.25;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
 
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(y) * cos(z);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)cos(y)cos(z) ----------------------", REMARK);

  zeroF (U);
  U[  1][  1][  1].Re = 0.125;
  U[N-1][  1][  1].Re = 0.125;
  U[  1][N-1][  1].Re = 0.125;
  U[N-1][N-1][  1].Re = 0.125;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
          
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(x) * cos(y) * cos(z);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)sin(y)cos(z) ----------------------", REMARK);

  zeroF (U);
  U[  1][  1][  1].Im = -0.125;
  U[N-1][  1][  1].Im = -0.125;
  U[  1][N-1][  1].Im =  0.125;
  U[N-1][N-1][  1].Im =  0.125;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
          
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(x) * sin(y) * cos(z);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking cos(2x)cos(2y)cos(2z) -------------------", REMARK);

  zeroF (U);
  U[  2][  2][  2].Re = 0.125;
  U[N-2][  2][  2].Re = 0.125;
  U[  2][N-2][  2].Re = 0.125;
  U[N-2][N-2][  2].Re = 0.125;

  rc3DFT  (U, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);
          
  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] -= cos(2.0*x) * cos(2.0*y) * cos(2.0*z);
      }
    }
  }

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  message (prog, "Checking full DFT with random numbers ------------", REMARK);

  for (i = 0; i < 2 * Npts; i++) v[i] = ran2PI (&seed);
  copyF (U, V);

  sprintf (s, "maximum physical component: %g", amaxF (U));
  message (prog, s, REMARK);  
  
  rc3DFT  (U, Wtab, FORWARD);
  scaleFT (U);
  rc3DFT  (U, Wtab, INVERSE);

  for (i = 0; i < 2 * Npts; i++) u[i] -= v[i];

  sprintf (s, "maximum physical error:     %g", amaxF (U));
  message (prog, s, REMARK);

  return (EXIT_SUCCESS);
}


  
