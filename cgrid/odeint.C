/*****************************************************************************
 * ODE integration driver.  Numerical Recipes, 2e.
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <Utility.h>
    
#define MAXSTP 10000

    


    
void odeint (double   *ystart                                                ,
	     int       nvar                                                  ,
	     double    x1                                                    ,
	     double    x2                                                    ,
	     double    eps                                                   ,
	     double    h1                                                    ,
	     double    hmin                                                  ,
	     int      *nok                                                   ,
	     int      *nbad                                                  ,
	     int       kmax                                                  ,
	     int      *kount                                                 ,
	     double  **state                                                 ,
	     double    dxsav                                                 ,
	     void    (*derivs) (double, double*, double*, void*)             ,
	     void     *params                                                ,
	     void    (*rkqc)   (double*, double*, int,     double*,  double  ,
 			        double,  double*, double*, double*,
			        void (*) (double, double*, double*, void*),
			        void*)                                       )
/* ========================================================================= *
 * Runge-Kutta driver with adaptive stepsize control.  Integrate starting
 * values ystart[1..nvar] from x1 to x2 with accuracy eps, storing intermed-
 * iate results in state[1..kmax][1..nvar+1].  h1 should be set as a guessed
 * first stepsize, hmin as the minimum allowed stepsize (can be zero).
 *
 * On output, nok & nbad are the number of good and bad (but retried and 
 * fixed) steps taken, and ystart is replaced by values at the end of the
 * integration interval.  derivs is the user-supplied routine for calculating
 * the RHS derivative, while rkqs is the name of the stepper routine to use.
 *
 * All arrays are base-1 indexed.
 * ========================================================================= */
{
  int      nstp, i;
  double   xsav, x, hnext, hdid, h;
  double  *yscal, *y, *dydx;

    
  yscal = dvector (1, nvar);
  y     = dvector (1, nvar);
  dydx  = dvector (1, nvar);

  x = x1;
  h = (x2 > x1) ? fabs (h1) : -fabs (h1);
  *nok = (*nbad) = *kount = 0;

  for (i = 1; i <= nvar; i++) y[i] = ystart[i];
  if (state) xsav = x - dxsav * 2.0;
  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    (*derivs) (x, y, dydx, params);
    for (i = 1; i <= nvar; i++)
      yscal[i] = 1.0;

    /* -- Alternative: yscal[i] = fabs(y[i]) + fabs (dydx[i]) + EPSm20; */

    if (state) {
      if (fabs(x - xsav) > fabs (dxsav)) {
	if (*kount < kmax - 1) {
	  ++*kount;
	  state[*kount][1] = x;
	  for (i = 1; i <= nvar; i++) state[*kount][i+1] = y[i];
	  xsav = x;
	}
      }
    }

    if ((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x;
    (*rkqc) (y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs, params);
    if (hdid == h) ++*nok; else ++*nbad;

    if ((x - x2) * (x2 - x1) >= 0.0) {
      for (i = 1; i <= nvar; i++) ystart[i] = y[i];
      if (state) {
	++*kount;
	state[*kount][1] = x;
	for (i = 1; i <= nvar; i++) state[*kount][i+1] = y[i];
      }
      freeDvector (dydx,  1);
      freeDvector (y,     1);
      freeDvector (yscal, 1);
      return;
    }

    if (fabs (hnext) <= hmin)
      message ("odeint", "Step size too small", ERROR);
    h = hnext;
  }
  message ("odeint", "Too many steps taken", ERROR);
}
    
#undef MAXSTP



