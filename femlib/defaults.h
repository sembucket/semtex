/*****************************************************************************
 * DEFAULTS.H:  default parameter initializations for initial.y.
 *****************************************************************************/

/* $Id$ */

#include <femdef.h>


static struct {			/* Mathematical constants. */
  char  *name;
  double cval;
} consts[] = {
  "PHI"   ,   1.61803398874989484820 ,
  "E"     ,   2.71828182845904523536 ,
  "DEG"   ,  57.29577951308232087721 ,
  "TWOPI" ,   6.28318530717958647688 ,
  "PI"    ,   3.14159265358979323844 ,
  0       ,   0.0
};

static struct {			/* Default options (global flags). */
  char *name;
  int   oval;
} option_init[] = {
  "BINARY"  ,  0         ,
  "GEOM"    ,  CART2D    ,
  "BASIS"   ,  GLL       ,
  "RULE"    ,  GL        ,
  "PROBLEM" ,  HELMHOLTZ ,
  "VERBOSE" ,  0         ,
  0         ,  0
};

static struct {			/* Default integer parameters. */
  char *name;
  int   ival;
} iparam_init[] = {
  0 , 0
};

static struct {			/* Default floating point parameters. */
  char  *name;
  double dval;
} dparam_init[] = {
  "BETA"   , 1.0 ,
  "LAMBDA2", 0.0 ,
  "KINVIS" , 1.0 ,
  "RHO"    , 1.0 ,
  0        , 0.0
};
