/*****************************************************************************
 * integration.C:  supply coefficients for discrete time integration schemes.
 *
 * Maximum time order supported is 3.
 *****************************************************************************/

static char 
RCSid[] = "$Id$";

#include <Fem.h>


const int Integration::OrderMax = 3;


void Integration::AdamsBashforth  (const int n, real* coeff)
// ---------------------------------------------------------------------------
// Adams-Bashforth predictor coefficients of order n.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::AdamsBashforth";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    break;
  case 2:
    coeff[0] =  1.5;
    coeff[1] = -0.5;
    break;
  case 3:
    coeff[0] = 23.0 / 12.0;
    coeff[1] = -4.0 /  3.0;
    coeff[2] =  5.0 / 12.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}


void Integration::StifflyStable (const int n, real* coeff)
// ---------------------------------------------------------------------------
// "Stiffly-stable" backwards differentiation coefficients of order n.
// NB: vector coeff must be of length n + 1.  First coefficient in each
// case applies to the new time level.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::StifflyStable";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    coeff[1] = -1.0;
    break;
  case 2:
    coeff[0] =  1.5;
    coeff[1] = -2.0;
    coeff[2] =  0.5;
    break;
  case 3:
    coeff[0] = 11.0 / 6.0;
    coeff[1] = -3.0;
    coeff[2] =  1.5;
    coeff[3] = -1.0 / 3.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}


void Integration::Extrapolation  (const int n, real* coeff)
// ---------------------------------------------------------------------------
// Coefficients of order n for explicit extrapolation to end of timestep.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::Extrapolation";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    break;
  case 2:
    coeff[0] =  2.0;
    coeff[1] = -1.0;
    break;
  case 3:
    coeff[0] =  3.0;
    coeff[1] = -3.0;
    coeff[2] =  1.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}
