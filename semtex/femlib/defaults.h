/*****************************************************************************
 * DEFAULTS.H:  default parameter initializations for initial.y.             *
 *****************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/
/* $Id$ */

#include <femdef.h>


static struct {			/* Math constants */
  char  *name;
  double cval;
} consts[] = {
  "PHI"   ,   1.61803398874989484820 ,
  "E"     ,   2.71828182845904523536 ,
  "DEG"   ,  57.29577951308232087721 ,
  "TWOPI" ,   6.28318530717958647688 ,
  "PI"    ,   3.14159265358979323844 ,
  NULL,       0.0
};

static struct {			/* Default options (global flags) */
  char *name;
  int   oval;
} option[] = {
  "binary",  1 ,		/* Save field-files in binary form */
  "core"  ,  1 ,		/* Do matrix solves in memory      */
  NULL    ,  0
};

static struct {			/* Default integer parameters */
  char *name;
  int   ival;
} iparam[] = {
  "ndim"          ,  2            ,
  "nz"            ,  1            ,
  "nsteps"        ,  1            ,
  "io-step"       ,  0            ,
  "verbose"       ,  0            ,
  "problem-type"  ,  SNSE         ,
  "problem-form"  ,  PENALTY      ,
  "element-type"  ,  QUAD         ,
  "element-order" ,  2            ,
  NULL            ,  0
};

static struct {			/* Default double parameters */
  char  *name;
  double dval;
} dparam[] = {
  "delta-t"   ,    0.05,
  "tolabs"    ,    1.0e-8,
  "tolrel"    ,    1.0e-4,
  "kinvisc"   ,    1.0,
  "rho"       ,    1.0,
  NULL        ,    0.0
};
