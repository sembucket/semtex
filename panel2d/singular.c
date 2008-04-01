/*****************************************************************************
 * SINGULAR.C:  Routines for singularity panels.
 *
 * In general, notation conforms to that used by Katz & Plotkin 1991.
 *****************************************************************************/

static char RCSid[] = "$Id$";

#include <math.h>
#include <stdio.h>
#include "panel.h"


void vor2Dl(Point  CP ,	/* Collocation point                         */
	    Point  LE ,	/* Leading edge coordinates                  */
	    Point  TE ,	/* Trailing edge coordinates                 */
	    Point  n  ,	/* Components of unit outward normal         */
	    Point *Va ,	/* Velocity components in global system      */
	    Point *Vb ,	/*                                           */
            int    IN )
/* ========================================================================= *
 * Influence coefficient calculation for velocity components measured in     *
 * panel coordinate system on the assumption of linear vorticity             *
 * distribution over the element.  As discussed in Katz & Plotkin \S 11.4.2, *
 * velocity components are split up according to the panel-end vorticity     *
 * strengths: this facilitates matrix solution for the unknowns (which are   *
 * just these panel-end strengths).                                          *
 *                                                                           *
 * NB: K&P eq (11.100) is in error: the term (\theta_{j+1} - \theta_{j})     *
 * should be (\theta_{j} - \theta_{j+1}).                                    *
 *                                                                           *
 * Obtain panel coordinates by subtracting LE value & rotating TE & CP.      * 
 * ========================================================================= */
{
  double  tmp,     den;
  double  thetaLE, thetaTE;
  double  rLE,     rTE;


  TE.x -= LE.x;  TE.z -= LE.z;
  CP.x -= LE.x;  CP.z -= LE.z;
  LE.x  = 0.0;   LE.z  = 0.0;

  TE.x = n.z * TE.x - n.x * TE.z;
  TE.z = 0.0;

  tmp  = n.z * CP.x - n.x * CP.z;
  CP.z = n.x * CP.x + n.z * CP.z;
  CP.x = tmp;

  /* Compute velocity components in panel coordinates. */

  if (IN) { /* Collocation point internal to panel. */

    Va->x = -0.5 * (CP.x - TE.x) / TE.x;  Va->z = -1.0 / TWOPI;
    Vb->x =  0.5 *  CP.x         / TE.x;  Vb->z =  1.0 / TWOPI;

  } else {  /* Collocation point outside panel.     */

    thetaLE = atan2(CP.z, CP.x       );
    thetaTE = atan2(CP.z, CP.x - TE.x);
    rLE     = hypot(CP.x,        CP.z);
    rTE     = hypot(CP.x - TE.x, CP.z);
    den     = 1.0 / (TWOPI * TE.x);

    Va->x  = (TE.x - CP.x) * (thetaTE - thetaLE);
    Va->x -= CP.z * log( rTE / rLE );
    Va->x *= den;

    Va->z  = (CP.x - TE.x) * log( rLE / rTE );
    Va->z -= TE.x - CP.z * (thetaTE - thetaLE);
    Va->z *= den;
    
    Vb->x  = CP.z * log( rTE / rLE );
    Vb->x += CP.x * (thetaTE - thetaLE);
    Vb->x *= den;
    
    Vb->z  = TE.x - CP.z * (thetaTE - thetaLE);
    Vb->z -= CP.x * log( rLE / rTE );
    Vb->z *= den;
  }

  /* Rotate velocity components back into global coordinates. */

  tmp   = n.z * Va->x + n.x * Va->z;
  Va->z = n.z * Va->z - n.x * Va->x;
  Va->x = tmp;

  tmp   = n.z * Vb->x + n.x * Vb->z;
  Vb->z = n.z * Vb->z - n.x * Vb->x;
  Vb->x = tmp;
}
