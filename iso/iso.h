/*****************************************************************************
 * globals.h: typedefs, prototypes for isotropic code & associated routines.
 * 
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define STRMAX      80
#define TRUE        1
#define FALSE       0
#define FORWARD     1
#define INVERSE     0
#define MAGIC       0x80e9f3ef

#define SQR(a)      ((a) * (a))
#define SGN(a)      ((a) < 0.0 ?    -1.0 : 1.0)
#define MIN(a, b)   ((a) < (b) ?     (a) : (b))
#define MAX(a, b)   ((a) > (b) ?     (a) : (b))

enum    err_lev     {WARNING, ERROR, REMARK};

typedef float       real;
typedef char        string[STRMAX];
typedef struct      {real Re, Im;} complex;
typedef complex***  CF;
typedef CF*         CVF;

typedef struct {
  int     Magic;    /* 0x80e9f3ef = 0x" iso" & 0x80808080    */
  string  Title;    /* Name of simulation                    */
  string  IC_File;  /* Used as root for restart file names   */
  int     N_Grid;   /* Power of 2: grid is N_Grid^3          */
  real    Delta_T;  /* Time step                             */
  int     N_Save;   /* Number of timesteps between restarts  */
  int     Max_Step; /* Maximum number of timesteps to take   */
  int     N_Step;   /* Number of steps taken so far          */
  real    K_Visc;   /* Kinematic viscosity                   */
} header;


/* -- io.c */

void  message         (const char*, const char*, int);
FILE* efopen          (string, string);
void  read_start_file (FILE*, header*);
void  read_field      (FILE*, CVF, const int);
void  write_field     (FILE*, CVF, const int);
void  make_file_name  (const header*, string);
void  write_restart   (FILE*, const header*, const CVF, const CVF, const int);
void  print_header    (FILE*, const header*);
void  read_header     (FILE*, header*);
void  write_header    (FILE*, const header*);



real**  cfield    (int*, CVF*);
void    tophat    (int*, CVF, int);
real    normalize (int*, complex*, real**, CVF);

int*      ivector (int, int);
complex*  cvector (int, int);
real*     cbox    (int, int, int, int, int, int, CF*);

real**  cfield           (int*, CVF*);
void    zero             (CVF, const int*);
void    allocate_storage (CVF*, CVF*, CVF*, CVF*,
			  CF*,  complex**, complex**, int*);

/* -- FFT.c */

void preFFT   (int, complex*);
void preShift (int, complex*);
void rc3DFT   (CF, const int*, const complex*, const int);
void scaleFT  (CF, const int*);
int  ispow2   (int);

/* -- nonlin.c */

void  nonlin (const CVF, CVF, CF, CVF,
	      const complex*, const complex*, const int*);

/* -- integrate.c */

void  integrate (CVF, const CVF, const CVF,
		 const header*, const int*, const int);

/* -- energy.c */

real  energyP      (const int*, CVF V, const complex*);
real  energyF      (const int*, const CVF);
real  genEnstrophy (const int*, const CVF);

/* -- truncation.c */

void  truncation (CF, const int*);

/* -- derivative.c */

void  deriv (const CVF, const int, CF, const int, const int*);
void  curl  (const CVF, CVF, CF, const int*);

/* -- pressure.c */

void  pressure (CVF, CVF, CF, CVF, CF, complex*, complex*, int*);

/* -- taylor.c */

void  Taylor2D       (const int*, CVF, const int);
void  Taylor2D_error (const int*, CVF, const header*, const int);
void  TaylorGreen    (const int*, CVF);
real  Brachet        (const real);

