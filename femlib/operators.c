/*****************************************************************************
 * OPERATORS.C:  maintain operators for mesh points, derivatives.            *
 *****************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/
static char
  RCSid[] = "$Id$";

#include <stdio.h>
#include <linalp.h>
#include "femlib.h"





typedef struct meshopr {	/* ---- operator information structure  ---- */
  int              basis    ;	/* 'f' ==>  finite, 's' ==> spectral element */
  int              nk       ;	/* number of interpolant knot points         */
  int              nq       ;	/* number of quadrature points               */
  double          *knot     ;	/* Lagrange knots    on [-1, 1]              */
  double          *quad     ;	/* Quadrature points on [-1, 1]              */
  double          *weight   ;	/* Quadrature weights                        */
  double         **interp   ;	/* Interpolate: knots --> quadrature points  */
  double         **interpT  ;	/* Transpose of the above                    */
  double         **deriv    ;	/* Derivative:  knots --> quadrature points  */
  double         **derivT   ;	/* Transpose of the above                    */
  struct meshopr  *next     ;	/* link to next one                          */
} Meshopr;			/* ----------------------------------------- */





void  get_ops(char      basis ,	/* input: element basis: 'f' or 's'          */
	      int       nk ,	/* input: number of knot points              */
	      int       nq ,	/* input: number of quadrature points        */
	      double  **kp ,	/* pointer to knot point storage             */
	      double  **qp ,	/* pointer to quadrature point storage       */
	      double  **qw ,	/* pointer to quadrature weight storage      */
	      double ***in ,	/* pointer to interpolation matrix           */
	      double ***it ,	/* pointer to transposed interp matrix       */
	      double ***dr ,	/* pointer to derivative matrix              */
	      double ***dt )	/* pointer to transposed derivative matrix   */
/* ========================================================================= *
 * Maintain/return operators for finite elements with high- or low-order     *
 * basis functions.   All operators are defined for the interval [-1, 1].    *
 *                                                                           *
 * If the required operators are "in stock", they are returned from a list,  *
 * otherwise they are created and added to the list as a first step.         *
 * NULL pointers are not assigned a value.                                   *
 * ========================================================================= */
{
  int              found;
  static Meshopr  *head = NULL;
         Meshopr  *p;


  for (p=head; p; p=p->next) {
    found = p->basis == basis && p->nk == nk && p->nq == nq;
    if (found) break;
  }

  if (!found) {		/* Make more storage and operators. */

    p = (Meshopr *) calloc(1, sizeof(Meshopr));
    if (head) p -> next = head -> next;
    head = p;

    p -> basis = basis;
    p -> nk    = nk;
    p -> nq    = nq;

    switch (basis) {

    case STANDARD:		/* Standard finite-element basis functions. */

      p -> knot    = dvector(0, nk-1);
      p -> quad    = dvector(0, nq-1);
      p -> weight  = dvector(0, nq-1);
      p -> interp  = dmatrix(0, nq-1, 0, nk-1);
      p -> interpT = dmatrix(0, nk-1, 0, nq-1);
      p -> deriv   = dmatrix(0, nq-1, 0, nk-1);
      p -> derivT  = dmatrix(0, nk-1, 0, nq-1);

      uniknot  (nk, p->knot);
      zwgl     (p->quad, p->weight, nq);
      intmat_g (nk, p->knot, nq, p->quad, p->interp, p->interpT);
      dermat_g (nk, p->knot, nq, p->quad, p->deriv,  p->derivT );
      break;

    case GLL:			/* Spectral element basis functions. */

      if (nk != nq)
	message("GetOps()", "spectral element: need nk & nq identical", ERROR);

      p -> knot    = dvector(0, nk-1);
      p -> quad    = p -> knot;
      p -> weight  = dvector(0, nk-1);
      p -> interp  = NULL;	/* No interpolation is needed. */
      p -> interpT = NULL;
      p -> deriv   = dmatrix(0, nk-1, 0, nk-1);
      p -> derivT  = dmatrix(0, nk-1, 0, nk-1);

      zwgll (p->knot, p->weight, nk);
      dgll  (nk, p->knot, p->deriv, p->derivT);
      break;
      
    default:
      message("GetOps()", "unknown basis kind", ERROR);
      break;
    }
  }

  /* p now points to valid storage: return requested operators. */

  if (kp) *kp = p -> knot;
  if (qp) *qp = p -> quad;
  if (qw) *qw = p -> weight;
  if (in) *in = p -> interp;
  if (it) *it = p -> interpT;
  if (dr) *dr = p -> deriv;
  if (dt) *dt = p -> derivT;

}
