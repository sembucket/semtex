///////////////////////////////////////////////////////////////////////////////
// integration.C:  supply coefficients for discrete time integration schemes.
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// Maximum time order supported is 3 (4 for implicit Adams--Moulton methods).
// Coefficients for all schemes can be found in Gear's book, "Numerical
// Initial Value Problems in Ordinary Differential Equations", 1971.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <Sem.h>


const integer Integration::OrderMax = 4;


void Integration::AdamsBashforth  (const integer n    ,
				   real*         coeff)
// ---------------------------------------------------------------------------
// Adams--Bashforth (predictor) coefficients of order n.  Gear, Table 7.3.
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
    coeff[0] =  23.0 / 12.0;
    coeff[1] = -16.0 / 12.0;
    coeff[2] =   5.0 / 12.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}


void Integration::AdamsMoulton (const integer n    ,
				real*         coeff)
// ---------------------------------------------------------------------------
// Adams--Moulton (corrector) coefficients of order n.  Gear, Table 7.5.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::AdamsMoulton";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    break;
  case 2:
    coeff[0] =  0.5;
    coeff[1] =  0.5;
    break;
  case 3:
    coeff[0] =  5.0 / 12.0;
    coeff[1] =  8.0 / 12.0;
    coeff[2] = -1.0 / 12.0;
    break;
  case 4:
    coeff[0] =  9.0 / 24.0;
    coeff[1] = 19.0 / 24.0;
    coeff[2] = -5.0 / 24.0;
    coeff[3] =  1.0 / 24.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}


void Integration::StifflyStable (const integer n    ,
				 real*         coeff)
// ---------------------------------------------------------------------------
// "Stiffly-stable" backwards differentiation coefficients of order n.
// NB: vector coeff must be of length n + 1.  First coefficient in each
// case applies to the new time level.  Gear, Table 11.1.
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


void Integration::Extrapolation  (const integer n    ,
				  real*         coeff)
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
