/*****************************************************************************
 *                            FILE & I-O UTILITIES                           *
 *****************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/
static char
  RCS_util[] = "$Id$";


#include <sys/types.h>
#include <malloc.h>
#include <stdio.h>
#include <time.h>
#include <alplib.h>


int     _lapIreg[NVREG];	/* For FORTRAN linkage. */
char    _lapCreg[NVREG];
float   _lapSreg[NVREG];
double  _lapDreg[NVREG];

char     buf[STR_MAX];		/* A string for general use. */




void message(char *routine, char *text, err_lev level)
/* ========================================================================= *
 * A general error handler.                                                  *
 * ========================================================================= */
{
  switch (level) {
  case WARNING:
    fprintf(stderr, "WARNING");
    break;
  case ERROR:
    fprintf(stderr, "ERROR  ");
    break;
  case REMARK:
    fprintf(stderr, "REMARK ");
    break;
  default:
    fprintf(stderr, "bad error level in message: %d\n", level);
    exit(1);
    break;
  }

  fprintf(stderr, " in %s: %s\n", routine, text);

  if (level == ERROR) exit(1);
}





FILE *efopen(char *file, char *mode)
/* ========================================================================= *
 * fopen file, die if can't.                                                 *
 * ========================================================================= */
{
  FILE *fp;


  if ( fp = fopen(file, mode) ) return fp;

  sprintf(buf, "can't open %s mode %s", file, mode);
  message("efopen()", buf, ERROR);
}





#if !defined(i860) && !defined(dclock)

double dclock(void)
/* ========================================================================= *
 * Double-precision timing routine.                                          *
 * ========================================================================= */
{
  static double tps = 1.0 / CLOCKS_PER_SEC;
  return (double) clock() * tps;
}


float sclock(void)
/* ========================================================================= *
 * Single-precision timing routine.                                          *
 * ========================================================================= */
{
  static float tps = 1.0F / CLOCKS_PER_SEC;
  return (float) clock() * tps;
}

#endif





void putDvector(FILE      *fp     ,
		int        width  ,
		int        prec   ,
		int        ntot   ,
		int        nfield , ...)
/* ========================================================================= *
 * Write (ASCII) a variable number of dvectors on fp, in columns.            *
 * ========================================================================= */
{
  int        i, j;
  double   **u;
  va_list    ap;


  u = (double **) calloc (nfield, sizeof(double*));
  va_start(ap, nfield);
  for (i=0; i<nfield; i++) u[i] = va_arg(ap, double*); 
  va_end(ap);

  for (i=0; i<ntot; i++) {
    for (j=0; j<nfield; j++) {
      fprintf(fp, "%*.*f", width, prec, u[j][i]);
    }
    fprintf(fp, "\n");
  }

  free(u);

}





void putIvector(FILE      *fp     ,
		int        width  ,
		int        ntot   ,
		int        nfield , ...)
/* ========================================================================= *
 * Write (ASCII) a variable number of ivectors on fp, in columns.            *
 * ========================================================================= */
{
  int        i, j;
  int      **u;
  va_list    ap;


  u = (int **) calloc (nfield, sizeof(int*));
  va_start(ap, nfield);
  for (i=0; i<nfield; i++) u[i] = va_arg(ap, int*); 
  va_end(ap);

  for (i=0; i<ntot; i++) {
    for (j=0; j<nfield; j++) {
      fprintf(fp, "%*d", width, u[j][i]);
    }
    fprintf(fp, "\n");
  }

  free(u);

}





void putSvector(FILE      *fp     ,
		int        width  ,
		int        prec   ,
		int        ntot   ,
		int        nfield , ...)
/* ========================================================================= *
 * Write (ASCII) a variable number of svectors on fp, in columns.            *
 * ========================================================================= */
{
  int       i, j;
  float   **u;
  va_list   ap;


  u = (float **) calloc (nfield, sizeof(float*));
  va_start(ap, nfield);
  for (i=0; i<nfield; i++) u[i] = va_arg(ap, float*); 
  va_end(ap);

  for (i=0; i<ntot; i++) {
    for (j=0; j<nfield; j++) {
      fprintf(fp, "%*.*f", width, prec, u[j][i]);
    }
    fprintf(fp, "\n");
  }

  free(u);

}





