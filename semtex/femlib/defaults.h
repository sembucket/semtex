/*****************************************************************************
 * DEFAULTS.H:  default parameter initializations for initial.y.
 * All variables are now in a single table, everything is double precision.
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

static struct {
  char*  name;
  double cval;
} consts[] = {

  /* -- Mathematical constants. */

  "PHI"      ,   1.61803398874989484820 ,
  "E"        ,   2.71828182845904523536 ,
  "DEG"      ,  57.29577951308232087721 ,
  "TWOPI"    ,   6.28318530717958647688 ,
  "PI"       ,   3.14159265358979323844 ,
  
  /* -- Default named parameters. */

  "D_T"      ,   0.01  ,
  "TOL_REL"  ,   1.0e-6,
  "TOL_ABS"  ,   1.0e-6,
  "BETA"     ,   1.0   ,	/* -- TWOPI / Lz. */
  "LAMBDA2"  ,   0.0   ,	/* -- Helmholtz constant. */
  "KINVIS"   ,   1.0   ,
  "RHO"      ,   1.0   ,

  /* -- Option switches. */

  "ITERATIVE",   0        ,
  "PROBLEM"  ,   HELMHOLTZ,
  "BINARY"   ,   0        ,
  "GEOMETRY" ,   CART2D   ,
  "BASIS"    ,   GLL      ,
  "RULE"     ,   LL       ,
  "VERBOSE"  ,   0        ,
  "CHKPOINT" ,   0        ,
  "OPTIMIZE" ,   1        ,
  
  /* -- Default integer values. */

  "IO_FLD"   ,   1000,
  "IO_HIS"   ,   10  ,
  "IO_CFL"   ,   0   ,
  "STEP_MAX" ,   500 ,
  "N_POLY"   ,   5   ,
  "N_TIME"   ,   2   ,
  "N_STEP"   ,   1   ,
  "N_Z"      ,   1   ,

  0          ,   0.0
};




