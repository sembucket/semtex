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

  "E"           ,   2.71828182845904523536 ,
  "DEG"         ,  57.29577951308232087721 ,
  "PI"          ,   3.14159265358979323844 ,
  "TWOPI"       ,   6.28318530717958647688 ,
  "EULER"       ,   0.57721566490153286061 ,
  "GOLDEN"      ,   1.61803398874989484820 ,
  
  /* -- Default named parameters. */

  "t"           ,   0.0    ,	/* -- Time.                               */
  "D_T"         ,   0.01   ,	/* -- Time step.                          */

  "TOL_REL"     ,   1.0e-8 ,	/* -- Relative tolerance (PCG)            */
  "TOL_ABS"     ,   1.0e-8 ,	/* -- Absolute tolerance.                 */
  "TOL_POS"     ,   1.0e-8 ,    /* -- Positional tolerance.               */

  "z"           ,   0.0    ,	/* -- z-plane location.                   */
  "BETA"        ,   1.0    ,	/* -- TWOPI / Lz (Fourier constant).      */
  "LAMBDA2"     ,   0.0    ,	/* -- Helmholtz constant.                 */

  "KINVIS"      ,   1.0    ,	/* -- Kinematic viscosity.                */
  "REFVIS"      ,   1.0    ,	/* -- Reference kinematic viscosity.      */
  "RHO"         ,   1.0    ,	/* -- Density.                            */
  "GRAVITY"     ,   9.81   ,	/* -- Gravitational acceleration.         */
  "T_REF"       ,   288.15 ,	/* -- Reference temperature (15C).        */
  "PRANDTL"     ,   0.72   ,	/* -- Prandtl number for air at STP.      */

  "C_SMAG"      ,   0.1114 ,	/* -- Smagorinsky's constant (RNG value). */
  "RNG_C"       ,   75.0   ,    /* -- a / gamma^4 for 3D RNG.             */
  "RNG_BIG"     ,   500.0  ,    /* -- Switch to plain Smag in RNG.        */

  "FFX"         ,   0.0    ,	/* -- Body force per unit mass (x).       */
  "FFY"         ,   0.0    ,	/* -- y component.                        */
  "FFZ"         ,   0.0    ,	/* -- z component.                        */

  /* -- Option switches. */

  "BASIS"       ,   GLL ,
  "RULE"        ,   LL  ,
  "ITERATIVE"   ,   0   ,
  "CYLINDRICAL" ,   0   ,
  "VERBOSE"     ,   0   ,
  "CHKPOINT"    ,   1   ,
  "AVERAGE"     ,   0   ,
  "SPAWN"       ,   0   ,
  "RNG"         ,   0   ,
  
  /* -- Default integer values. */

  "STEP_MAX"    ,   500 ,
  "IO_FLD"      ,   500 ,
  "IO_HIS"      ,   10  ,
  "IO_CFL"      ,   50  ,
  "N_POLY"      ,   5   ,
  "N_TIME"      ,   2   ,
  "N_STEP"      ,   1   ,
  "N_Z"         ,   1   ,
  "N_PROC"      ,   1   ,
  "I_PROC"      ,   0   ,

  0             ,   0.0
};



