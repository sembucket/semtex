/*****************************************************************************
 * DEFAULTS.H:  default parameter initializations for initial.y.
 *****************************************************************************/

/* $Id$ */

#include <femtype.h>
#include <femlib.h>


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
} option_init[] = {
  "BINARY"  ,  0         ,	/* Save field-files in binary form */
  "CORE"    ,  1         ,	/* Do matrix solves in memory      */
  "BASIS"   ,  GLL       ,	/* GLL basis functions             */
  "RULE"    ,  GL        ,	/* Gauss-Legendre quadrature rule  */
  "PROBLEM" ,  POTENTIAL ,	/* Potential (flow) problem        */
  NULL    ,  0
};

static struct {			/* Default integer parameters */
  char *name;
  int   ival;
} iparam_init[] = {
  "NDIM"     ,  2 ,
  "NZ"       ,  1 ,
  "NSTEPS"   ,  1 ,
  "IOSTEP"   ,  1 ,
  "VERBOSE"  ,  0 ,
  NULL       ,  0
};

static struct {			/* Default double parameters */
  char  *name;
  double dval;
} dparam_init[] = {
  "BETA"      ,    0.0,
  "DELTAT"    ,    0.05,
  "TOLABS"    ,    1.0e-8,
  "TOLREL"    ,    1.0e-4,
  "KINVIS"    ,    1.0,
  "RHO"       ,    1.0,
  NULL        ,    0.0
};
