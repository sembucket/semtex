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

char    buf[STR_MAX];           /* General I-O buffer.  */

int     _lap_ireg[NVREG];	/* For FORTRAN linkage. */
char    _lap_creg[NVREG];
float   _lap_sreg[NVREG];
double  _lap_dreg[NVREG];






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
