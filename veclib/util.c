/*****************************************************************************
 *                            FILE & I-O UTILITIES
 *
 * $Id$
 *****************************************************************************/


#include <sys/types.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "alplib.h"


int     _alpIreg[NVREG];	/* For FORTRAN linkage. */
char    _alpCreg[NVREG];
float   _alpSreg[NVREG];
double  _alpDreg[NVREG];

char     buf[STR_MAX];		/* A string for general use. */




void message (const char *routine, const char *text, int level)
/* ========================================================================= *
 * A simple error handler.
 * ========================================================================= */
{
  switch (level) {
  case WARNING:
    fprintf (stderr, "WARNING: %s: %s\n", routine, text); 
    break;
  case ERROR:
    fprintf (stderr, "ERROR: %s: %s\n", routine, text); 
    break;
  case REMARK:
    fprintf (stdout, "%s: %s\n", routine, text);
    break;
  default:
    fprintf (stderr, "bad error level in message: %d\n", level);
    exit (EXIT_FAILURE);
    break;
  }

  if (level == ERROR) exit (EXIT_FAILURE);
}





FILE *efopen(const char *file, const char *mode)
/* ========================================================================= *
 * fopen file, die if can't.
 * ========================================================================= */
{
  FILE *fp;


  if ( fp = fopen(file, mode) ) return fp;

  sprintf(buf, "can't open %s mode %s", file, mode);
  message("efopen()", buf, ERROR);
}





#if !defined(i860) && !defined(dclock)

double dclock (void)
/* ========================================================================= *
 * Double-precision timing routine.
 * ========================================================================= */
{
  static double tps = 1.0 / CLOCKS_PER_SEC;
  return (double) clock() * tps;
}





float sclock (void)
/* ========================================================================= *
 * Single-precision timing routine.
 * ========================================================================= */
{
  static float tps = 1.0F / CLOCKS_PER_SEC;
  return (float) clock() * tps;
}

#endif





void printDvector(FILE  *fp     ,
		  int    width  ,
		  int    prec   ,
		  int    ntot   ,
		  int    inc    ,
		  int    nfield , ...)
/* ========================================================================= *
 * Print up a variable number of dvectors on fp, in columns.
 * ========================================================================= */
{
  char      routine[] = "printDvector()";
  int       i, j, k;
  double  **u;
  va_list   ap;


  u = (double **) calloc (nfield, sizeof (double*));
  va_start (ap, nfield);
  for (i = 0; i < nfield; i++) u[i] = va_arg (ap, double*); 
  va_end (ap);

  k = (inc < 0) ? (-ntot + 1)*inc : 0;

  for (i = 0; i < ntot; i++) {
    for (j = 0; j < nfield; j++)
      if (fprintf(fp, "%*.*f", width, prec, u[j][k]) < 0)
	message (routine, "unable to write to file", ERROR);
    k += inc;
    fprintf(fp, "\n");
  }

  free(u);
}





void printIvector(FILE  *fp     ,
		  int    width  ,
		  int    ntot   ,
		  int    inc    ,
		  int    nfield , ...)
/* ========================================================================= *
 * Print up a variable number of ivectors on fp, in columns.
 * ========================================================================= */
{
  char       routine[] = "printIvector()";
  int        i, j, k;
  int      **u;
  va_list    ap;


  u = (int **) calloc (nfield, sizeof (int*));
  va_start (ap, nfield);
  for (i = 0; i < nfield; i++) u[i] = va_arg (ap, int*); 
  va_end (ap);

  k = (inc < 0) ? (-ntot + 1)*inc : 0;

  for (i = 0; i < ntot; i++) {
    for (j = 0; j < nfield; j++)
      if (fprintf (fp, "%*d", width, u[j][k]) < 0)
	message (routine, "couldn't write to file", ERROR);
    k += inc;
    fprintf (fp, "\n");
  }

  free (u);
}





void printSvector(FILE  *fp     ,
		  int    width  ,
		  int    prec   ,
		  int    ntot   ,
		  int    inc    ,
		  int    nfield , ...)
/* ========================================================================= *
 * Write (ASCII) a variable number of svectors on fp, in columns.
 * ========================================================================= */
{
  char      routine[] = "printSvector()";
  int       i, j, k;
  float   **u;
  va_list   ap;


  u = (float **) calloc (nfield, sizeof (float*));
  va_start (ap, nfield);
  for (i = 0; i < nfield; i++) u[i] = va_arg (ap, float*); 
  va_end (ap);

  k = (inc < 0) ? (-ntot + 1)*inc : 0;

  for (i = 0; i < ntot; i++) {
    for (j = 0; j < nfield; j++)
      if (fprintf (fp, "%*.*f", width, prec, u[j][k]) < 0)
	message (routine, "unable to write to file", ERROR);
    k += inc;
    fprintf (fp, "\n");
  }

  free (u);
}





