#ifndef ISO_H
#define ISO_H
/*****************************************************************************
 * iso.h: typedefs, prototypes for isotropic code & associated routines.
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 * 
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

extern int N, K, FourKon3;	/* Global sizes: N = 2*K = grid size.   */
				/* These are the only global variables. */

#define STR_MAX 256
#define TRUE    1
#define FALSE   0
#define FORWARD 1
#define INVERSE 0

#define SQR(a)       ((a) * (a))
#define MIN(a,b)     ((a) < (b) ? (a) : (b))
#define MAX(a,b)     ((a) > (b) ? (a) : (b))
#define CLAMP(t,a,b) (MAX(MIN((t),(b)),(a)))
#define RM(i,j,k)    ((k) + K * ((j) + N * (i)))
#define INSIDE       <= FourKon3
#define OUTSIDE      >  FourKon3
#define MAG(Z)       ((Z).Re*(Z).Re + (Z).Im*(Z).Im)
#define SHIFT(Z, W)  tempRe = (Z).Re; \
                     (Z).Re = tempRe*(W).Re - (Z).Im*(W).Im; \
                     (Z).Im = (Z).Im*(W).Re + tempRe*(W).Im

enum err_lev   { WARNING, ERROR, REMARK };
enum sgs_model { NONE,    SMAG,  S3     };

typedef double     real;	/* Global precision control. */
typedef struct     {real Re, Im;} complex;
typedef complex*** CF;		/* 3D array of complex.      */
typedef CF*        CVF;		/* 3-component vector of CF. */

typedef struct {
  char* session;		/* Name of simulation                  */
  FILE* fld_dmp;		/* Field file for writes               */
  FILE* his_dmp;		/* History file for energy spectra     */
  int   io_fld;			/* Number of timesteps between dumps   */
  int   io_his;			/* Number of timesteps between dumps   */
  int   chkpnt;			/* Flag checkpointing on/off           */
  int   ngrid;			/* Power of 2: grid is ngrid^3         */
  int   norder;			/* Timestepping order                  */
  int   nstep;			/* Maximum number of timesteps to take */
  int   step;			/* Number of steps taken so far        */
  real  dt;			/* Time step                           */
  real  time;			/* Evolution time                      */
  real  kinvis;			/* Kinematic viscosity                 */
} Param;

/* -- io.c */

void  message    (const char*, const char*, int);
FILE* efopen     (const char*, const char*);
void  readCVF    (FILE*,       CVF);
void  writeCVF   (FILE*, const CVF);
void  startup    (Param*);
void  printParam (FILE*, const Param*);
void  readParam  (FILE*,       Param*);
void  writeParam (FILE*, const Param*);
void  restart    (CVF, Param*);
void  analyze    (CVF, Param*, const complex*);
void  dump       (const CVF, Param*);
void  cleanup    (Param*);
void  format     (char*);

/* -- allocate.c */

real**   cfield   (CVF*);
int*     ivector  (int, int);
complex* cvector  (int, int);
real*    rvector  (int, int);
real*    cbox     (int, int, int, int, int, int, CF*);
void     allocate (CVF*, CVF*, const int, CVF*, CF*,CF*, complex**,complex**);

/* -- FFT.c */

void preFFT   (complex*, const int);
void preShift (complex*, const int);
void rc3DFT   (CF, const complex*, const int);
void scaleFT  (CF);
int  ispow2   (int);

/* -- nonlinear.c */

void nonlinear (CVF, CVF, CVF, CF, CF, const complex*, const complex*);

/* -- integrate.c */

void integrate (CVF, const CVF*, const Param*);

/* -- energy.c */

real energyP     (CVF V, const complex*);
real energyF     (const CVF);
real enstrophyF  (const CVF, CVF, CF);
real microF      (const CVF, CF);
real rmsEns      (const CVF);
real L2norm      (const CF);
real amaxF       (const CF);
void normaliseVF (CVF);
void energySpec  (const CVF, real*);

/* -- operator.c */

void deriv    (const CVF, const int, CF, const int);
void curl     (const CVF, CVF, CF);
void shift    (CF, const complex*, const int);
void convolve (const CF, const CF, const CF, const CF, CF, CF,
	       const complex*, const complex*);

/* -- pressure.c */

void pressure (CVF, CVF, CF, CVF, CF, const complex*, const complex*);

/* -- taylor.c */

void Taylor2D       (CVF, const int);
void Taylor2D_error (CVF, const Param*, const int);
void TaylorGreen    (CVF);
real Brachet        (const real);

/* -- misc.c */

void zeroF      (CF);
void zeroVF     (CVF);
void copyF      (CF,  const CF);
void copyVF     (CVF, const CVF);
void scaleF     (CF,  const real);
void scaleVF    (CVF, const real);
void setF       (CF,  const CF);
void setVF      (CVF, const CVF);
void addF       (CF,  const CF);
void addVF      (CVF, const CVF);
void subF       (CF,  const CF);
void subVF      (CVF, const CVF);
void truncateF  (CF);
void truncateVF (CVF);

void filterF    (CF, real*, real*);
void projectVF  (CVF, CF);

/* -- random.c */

void randomise (CVF, int*, CF);
real ran2PI    (int*);

/* -- filter.c */

void bvdFilter (const int, const int, const int, const real, real*);
void ispectrum (CVF, const real, const real,
		real (*)(const real, const real, const real));
real unknown01 (const real, const real, const real);
real vonKarman (const real, const real, const real);

#endif
