#ifndef PRISM_CONSTANT_H
#define PRISM_CONSTANT_H

/* Prism's constants and compiled limits
 *
 * $Id$
 * ------------------------------------------------------------------------- */
 
#define _MAX_HP       32   /* Maximum number of history points          */
#define _MAX_NZ      128   /* Maximum number of z-planes                */
#define _MAX_TORDER    3   /* Maximum integration order (time)          */

typedef enum {                    /* ......... ACTION Flags .......... */
  Rotational,                     /* N(U) = U x curl U                 */
  SkewSymmetric,                  /* N(U) =[U . grad U + grad (U.U)]/2 */
  Stokes,                         /* N(U) = 0  [drive force only]      */
  Pressure,                       /* div (U'/dt)                       */
  Viscous,                        /* (U' - dt grad P) Re / dt          */
  Prep,                           /* Run the PREP phase of MakeF()     */
  Post,                           /* Run the POST phase of MakeF()     */
  Fourier,                        /* Transform Physical -> Fourier     */
  Physical,                       /* Transform Fourier  -> Physical    */
  Full,                           /* --------------------------------- */
  Depart,                         /* Parallel Transpose flags          */
  Arrive,                         /* --------------------------------- */
  eq_Helmholtz    = 0,            /* Helmholtz solver flags            */
  eq_Poisson      = 1,            /*                                   */
  eq_Laplace      = 2,            /* --------------------------------- */
  StifflyStable   = 0,            /* Integration methods               */
  AdamsBashforth  = 1,            /* --------------------------------- */
  Convective      = 9,            /* Special flag for convective form  */
  EndOfActionList = 999           /* Always keep as the last enum!     */
} ACTION;

#endif
