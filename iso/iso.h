/*****************************************************************************
 * iso.h: typedefs, prototypes for isotropic code & associated routines.
 * 
 * : iso.h,v 2.3 1995/11/24 19:28:02 hmb Exp $
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define STR_MAX     80
#define TRUE        1
#define FALSE       0
#define FORWARD     1
#define INVERSE     0

#define SQR(a)      ((a) * (a))
#define MIN(a, b)   ((a) < (b) ?     (a) : (b))
#define MAX(a, b)   ((a) > (b) ?     (a) : (b))

enum    err_lev     {WARNING, ERROR, REMARK};

typedef float       real;	/* Global precision control. */
typedef struct      {real Re, Im;} complex;
typedef complex***  CF;		/* 3D array of complex.      */
typedef CF*         CVF;	/* 3-component vector of CF. */

typedef struct {
  char    name[STR_MAX];  /*  Name of simulation                    */
  FILE*   output;         /*  Field file for writes                 */
  int     modes;          /*  Power of 2: grid is (2*Modes)^3       */
  int     stepSave;       /*  Number of timesteps between restarts  */
  int     stepMax;        /*  Maximum number of timesteps to take   */
  int     step;           /*  Number of steps taken so far          */
  real    Re;             /*  Re = 1 / Kinematic viscosity          */
  real    dt;             /*  Time step                             */
  real    time;           /*  Evolution time                        */
} Param;


/* -- io.c */

void  message     (const char*, const char*, int);
FILE* efopen      (const char*, const char*);
void  readCVF     (FILE*,       CVF, const int*);
void  writeCVF    (FILE*, const CVF, const int*);
void  startup     (FILE*,        Param*, const char*, const int);
void  printParam  (FILE*, const  Param*, const char *, const char *);
void  readParam   (FILE*,        Param*);
void  writeParam  (FILE*, const  Param*);
void  initialize  (       CVF,   Param*, const int*);
void  analyze     (       CVF,   Param*, const complex*, const int*);
void  dump        (const  CVF,   Param*, const int, const int*);
void  cleanup     (Param*, const int);

/* -- allocate.c */

real**    cfield     (int*, CVF*);
int*      ivector    (int, int);
complex*  cvector    (int, int);
real*     cbox       (int, int, int, int, int, int, CF*);
real**    cfield     (int*, CVF*);
void      allocate   (CVF*, CVF*, CVF*, CVF*, CF*, CF*,
		      complex**, complex**, int*);

/* -- FFT.c */

void preFFT   (complex*, const int);
void preShift (complex*, const int);
void rc3DFT   (CF, const int*, const complex*, const int);
void scaleFT  (CF, const int*);
int  ispow2   (int);

/* -- nonlinear.c */

void nonlinear (CVF, CVF, CF, CF, CVF,
		const complex*, const complex*, const int*);

/* -- integrate.c */

void integrate (CVF, const CVF, const CVF, const Param*, const int*);

/* -- energy.c */

real  energyP   (CVF   V,   const complex*, const int*);
real  energyF   (const CVF, const int*);
real  rmsEns    (const CVF, const int*);
real  L2norm    (const CF,  const int*);
real  amaxf     (const CF,  const int*);
void  normalize (      CVF, const int*);

/* -- truncation.c */

void  truncation (CF, const int*);

/* -- derivative.c */

void  deriv (const CVF, const int, CF, const int, const int*);
void  curl  (const CVF, CVF, CF, const int*);

/* -- pressure.c */

void  pressure (CVF, CVF, CF, CVF, CF, complex*, complex*, int*);

/* -- taylor.c */

void  Taylor2D       (CVF, const int*, const int);
void  Taylor2D_error (CVF, const int*, const Param*, const int);
void  TaylorGreen    (CVF, const int*);
real  Brachet        (const real);

/* -- misc.c */

void  zeroVF  (CVF, const int*);
void  zeroF   (CF,  const int*);
void  copyF   (CF,  const CF, const int*);
void  scaleF  (CF,  const real, const int*);
void  setF    (CF,  const CF, const int*);
void  addF    (CF,  const CF, const int*);
void  subF    (CF,  const CF, const int*);
void  project (CVF, CF, const int*);

/* -- random.c */

real  ran2PI (int*);

/* -- tophat.c */

void  tophat (int*, CVF, int);
