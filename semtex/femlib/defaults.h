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

  "PHI"         ,   1.61803398874989484820 ,
  "E"           ,   2.71828182845904523536 ,
  "DEG"         ,  57.29577951308232087721 ,
  "PI"          ,   3.14159265358979323844 ,
  "TWOPI"       ,   6.28318530717958647688 ,
  
  /* -- Default named parameters. */

  "D_T"         ,   0.01   ,	/* -- Time step.                     */
  "t"           ,   0.0    ,	/* -- Time.                          */
  "z"           ,   0.0    ,	/* -- z-plane location.              */
  "TOL_REL"     ,   1.0e-8 ,	/* -- Relative tolerance.            */
  "TOL_ABS"     ,   1.0e-8 ,	/* -- Absolute tolerance.            */
  "BETA"        ,   1.0    ,	/* -- TWOPI / Lz (Fourier constant). */
  "LAMBDA2"     ,   0.0    ,	/* -- Helmholtz constant.            */

  "KINVIS"      ,   1.0    ,	/* -- Kinematic viscosity.           */
  "RHO"         ,   1.0    ,	/* -- Density.                       */
  "GRAVITY"     ,   9.81   ,	/* -- Gravitational acceleration.    */
  "T_REF"       ,   288.15 ,	/* -- Reference temperature (15C).   */
  "PRANDTL"     ,   0.72   ,	/* -- Prandtl number for air at STP. */

  /* -- Option switches. */

  "ITERATIVE"   ,   0   ,
  "BINARY"      ,   0   ,
  "CYLINDRICAL" ,   0   ,
  "BASIS"       ,   GLL ,
  "RULE"        ,   LL  ,
  "VERBOSE"     ,   0   ,
  "CHKPOINT"    ,   0   ,
  
  /* -- Default integer values. */

  "IO_FLD"      ,   500 ,
  "IO_HIS"      ,   10  ,
  "IO_CFL"      ,   0   ,
  "STEP_MAX"    ,   500 ,
  "N_POLY"      ,   5   ,
  "N_TIME"      ,   2   ,
  "N_STEP"      ,   1   ,
  "N_Z"         ,   1   ,

  0             ,   0.0
};




