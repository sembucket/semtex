/*****************************************************************************
 * DFT_chk.c: exercise 3D real--complex FFT routines.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


#define SIZE 8

static char prog[] = "DFT_chk";


int main()
/* ------------------------------------------------------------------------- *
 * Test real DFT of SIZE*SIZE*SIZE data.
 * ------------------------------------------------------------------------- */
{
  int       c, i, j, k, Npts, N, K;
  CF        U;
  real*     u;
  int*      Dim;
  complex*  Wtab;
  char      s[STR_MAX];

  /* -- Allocation. */

  Dim  = ivector (1, 3);
  N    = Dim[1] = (Dim[2] = SIZE);
  K    = Dim[3] = SIZE / 2;
  Npts = Dim[1] * Dim[2] * Dim[3];
  u    = cbox (0, N-1, 0, N-1, 0, K-1, &U);

  Wtab = cvector (0, K-1);
  preFFT  (Wtab, K);
  sprintf (s, "Dimensions: %d, %d, %d", Dim[1], Dim[2], Dim[3]);
  message (prog, s, REMARK);

  message (prog, "Checking constant (1) ----------------------------", REMARK);

  zeroF (U, Dim);
  U[ 0][ 0][ 0].Re = 2.0;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x) ----------------------------------", REMARK);

  zeroF (U, Dim);
  U[  1][ 0][ 0].Re = 1.0;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking sin(x) ----------------------------------", REMARK);

  zeroF (U, Dim);
  U[  1][ 0][ 0].Im =  -1.0;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(y) ----------------------------------", REMARK);

  zeroF (U, Dim);
  U[ 0][  1][ 0].Re = 1.0;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking sin(y) ----------------------------------", REMARK);

  zeroF (U, Dim);
  U[ 0][ 1][ 0].Im =  -1.0;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);
 
  message (prog, "Checking cos(z) ----------------------------------", REMARK);

  zeroF (U, Dim);
  U[ 0][ 0][ 1].Re = 1.0;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);
 
  message (prog, "Checking sin(z) ----------------------------------", REMARK);

  zeroF (U, Dim);
  U[ 0][ 0][ 1].Im = -1.0;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)cos(y) ----------------------------", REMARK);

  zeroF (U, Dim);
  U[  1][  1][  0].Re = 0.5;
  U[N-1][  1][  0].Re = 0.5;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)sin(y) ----------------------------", REMARK);

  zeroF (U, Dim);
  U[  1][  1][  0].Im = -0.5;
  U[N-1][  1][  0].Im = -0.5;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)cos(z) ----------------------------", REMARK);

  zeroF (U, Dim);
  U[  1][  0][  1].Re = 0.5;
  U[N-1][  0][  1].Re = 0.5;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(y)cos(z) ----------------------------", REMARK);

  zeroF (U, Dim);
  U[  0][  1][  1].Re = 0.5;
  U[  0][N-1][  1].Re = 0.5;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)cos(y)cos(z) ----------------------", REMARK);

  zeroF (U, Dim);
  U[  1][  1][  1].Re = 0.25;
  U[N-1][  1][  1].Re = 0.25;
  U[  1][N-1][  1].Re = 0.25;
  U[N-1][N-1][  1].Re = 0.25;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(x)sin(y)cos(z) ----------------------", REMARK);

  zeroF (U, Dim);
  U[  1][  1][  1].Im = -0.25;
  U[N-1][  1][  1].Im = -0.25;
  U[  1][N-1][  1].Im = 0.25;
  U[N-1][N-1][  1].Im = 0.25;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  message (prog, "Checking cos(2x)cos(2y)cos(2z) -------------------", REMARK);

  zeroF (U, Dim);
  U[  2][  2][  2].Re = 0.25;
  U[N-2][  2][  2].Re = 0.25;
  U[  2][N-2][  2].Re = 0.25;
  U[N-2][N-2][  2].Re = 0.25;

  rc3DFT  (U, Dim, Wtab, INVERSE);
  sprintf (s, "maximum physical component: %g", amaxf (U, Dim));
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

  sprintf (s, "maximum physical error:     %g", amaxf (U, Dim));
  message (prog, s, REMARK);

  return (EXIT_SUCCESS);
}


  
