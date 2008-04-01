/*****************************************************************************
 *                             P A N E L . H
 *****************************************************************************/

/* $Id$ */

#define TWOPI   6.28318530717958647692
#define PI_180  0.017453292519943295
#define EPSm20  1.0e-20
#define STR_MAX 2048

#define MAX(a, b)  ((a) > (b) ?     (a) : (b))
#define MIN(a, b)  ((a) < (b) ?     (a) : (b))
#define SQR(a)     ((a) * (a))

char    buf[STR_MAX];
enum    err_lev {WARNING, ERROR, REMARK};
void    message (const char *routine, const char *txt, int level);
FILE*   efopen  (const char *file,    const char *mode);

typedef struct point {
  double x, z;
} Point;

/* ------------------------------------------------------------------------- *
 * Routines in panel2d.c:
 * ------------------------------------------------------------------------- */

double*  dvector      (int, int);
int*     ivector      (int, int);
double** dmatrix      (int, int, int, int);
void     free_dvector (double *, int);
void     dzero        (int, double*, int);
void     dmxv         (double*, int, double*, int, double*);


/* ------------------------------------------------------------------------- *
 * Routines in lu.c:
 * ------------------------------------------------------------------------- */

void ludcmp (double **, int, int *, double *);
void lubksb (double **, int, int *, double *);


/* ------------------------------------------------------------------------- *
 * Routines in singular.c:
 * ------------------------------------------------------------------------- */

void vor2Dl(Point, Point, Point, Point, Point *, Point *, int);


