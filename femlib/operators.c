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





typedef struct quadopr {	/* ---- quadrature operator information  --- */
  int             basis   ;	/* STANDARD or GLL                           */
  int             np      ;	/* number of interpolant knot points         */
  int             nq      ;	/* number of quadrature points               */
  double         *knot    ;	/* Lagrange knots    on [-1, 1]              */
  double         *quad    ;	/* Quadrature points on [-1, 1]              */
  double         *weight  ;	/* Quadrature weights                        */
  double        **interp  ;	/* Interpolant operator: knots-->quadrature  */
  double        **interpT ;	/* Transpose of the above                    */
  double        **deriv   ;	/* Derivative operator:  knots-->quadrature  */
  double        **derivT  ;	/* Transpose of the above                    */
  struct quadopr *next    ;	/* link to next one                          */
} QuadOpr;			/* ----------------------------------------- */


typedef struct meshopr {	/* ------- mesh operator information ------- */
  int             oldbasis;	/* STANDARD or GLL                           */
  int             newbasis;	/* STANDARD or GLL                           */
  int             np      ;	/* Number of points on original mesh         */
  int             ni      ;	/* Number of points in interpolated mesh     */
  double         *mesh    ;	/* Location of points on interpolated mesh   */
  double        **interp  ;	/* Interpolant operator: original-->new      */
  double        **interpT ;	/* Transpose of the above                    */
  double        **deriv   ;	/* Derivative operator: original-->new       */
  double        **derivT  ;	/* Transpose of the above                    */
  struct meshopr *next    ;	/* link to next one                          */
} MeshOpr;			/* ----------------------------------------- */





void  quadOps(int basis ,	/* input: element basis: STANDARD or GLL     */
	      int np    ,	/* input: number of knot points              */
	      int nq    ,	/* input: number of quadrature points        */
	      double  **kp ,	/* pointer to knot point storage             */
	      double  **qp ,	/* pointer to quadrature point storage       */
	      double  **qw ,	/* pointer to quadrature weight storage      */
	      double ***in ,	/* pointer to interpolation matrix           */
	      double ***it ,	/* pointer to transposed interpolant matrix  */
	      double ***dr ,	/* pointer to derivative matrix              */
	      double ***dt )	/* pointer to transposed derivative matrix   */
/* ========================================================================= *
 * Maintain/return QUADRATURE operators for finite elements with high- or    *
 * low-order basis functions. Operators are defined on the interval [-1, 1]. *
 *                                                                           *
 * If the required operators are "in stock", they are returned from a list,  *
 * otherwise they are created and added to the list as a first step.         *
 * NULL input pointers are not assigned a value.                             *
 *                                                                           *
 * ========================================================================= */
{
  char            *routine = "quadOps()";
  int              found   = 0;
  static QuadOpr  *head    = NULL;
         QuadOpr  *p;


  for (p=head; p; p=p->next) {
    found = p->basis == basis && p->np == np && p->nq == nq;
    if (found) break;
  }


  if (!found) {		/* Make more storage and operators. */

    p = (QuadOpr *) calloc(1, sizeof(QuadOpr));
    if (head) p -> next = head -> next;
    head = p;

    p -> basis = basis;
    p -> np    = np;
    p -> nq    = (basis == STANDARD) ? nq : np;

    switch (basis) {

    case GLL:

      p -> knot    = dvector(0, np-1);
      p -> quad    = p -> knot;
      p -> weight  = dvector(0, np-1);
      p -> interp  = (double **) NULL;	/* No interpolation needed. */
      p -> interpT = (double **) NULL;
      p -> deriv   = dmatrix(0, np-1, 0, np-1);
      p -> derivT  = dmatrix(0, np-1, 0, np-1);

      zwgll (p->knot, p->weight, np);
      dgll  (np, p->knot, p->deriv, p->derivT);
      break;

    case STANDARD: default:

      p -> knot    = dvector(0, np-1);
      p -> quad    = dvector(0, nq-1);
      p -> weight  = dvector(0, nq-1);
      p -> interp  = dmatrix(0, nq-1, 0, np-1);
      p -> interpT = dmatrix(0, np-1, 0, nq-1);
      p -> deriv   = dmatrix(0, nq-1, 0, np-1);
      p -> derivT  = dmatrix(0, np-1, 0, nq-1);
      
      uniknot  (np, p->knot);
      zwgl     (p->quad, p->weight, nq);
      intmat_g (np, p->knot, nq, p->quad, p->interp, p->interpT);
      dermat_g (np, p->knot, nq, p->quad, p->deriv,  p->derivT );
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





void  meshOps(int oldbasis   ,	/* input: element basis: STANDARD or GLL     */
	      int newbasis   ,	/* input: desired basis: STANDARD or GLL     */
	      int np         ,	/* input: number of knot points              */
	      int ni         ,	/* input: number of interpolant points       */
	      double  **mesh ,	/* pointer to interpolated mesh storage      */
	      double ***in   ,	/* pointer to interpolation matrix           */
	      double ***it   ,	/* pointer to transposed interpolant matrix  */
	      double ***dr   ,	/* pointer to derivative matrix              */
	      double ***dt   )	/* pointer to transposed derivative matrix   */
/* ========================================================================= *
 * Maintain/return INTERPOLATION operators for finite elements with high- or *
 * low-order basis functions. Operators are defined on the interval [-1, 1]. *
 * It is legal to ask for interpolation onto the same mesh as input, in      *
 * which case the mesh storage and interpolation matrices are NULL.          *
 *                                                                           *
 * If the required operators are "in stock", they are returned from a list,  *
 * otherwise they are created and added to the list as a first step.         *
 * NULL input pointers are not assigned a value.                             *
 *                                                                           *
 * ========================================================================= */
{
  char            *routine = "meshOps()";
  int              found   = 0;
  static MeshOpr  *head    = NULL;
         MeshOpr  *p;
  double          *oldmesh;


  for (p=head; p; p=p->next) {
    found = p->oldbasis == oldbasis
      &&    p->newbasis == newbasis 
      &&    p->np       == np
      &&    p->ni       == ni;
    if (found) break;
  }


  if (!found) {		/* Make more storage and operators. */

    p = (MeshOpr *) calloc(1, sizeof(MeshOpr));
    if (head) p -> next = head -> next;
    head = p;

    p -> oldbasis = oldbasis;
    p -> newbasis = newbasis;
    p -> np       = np;
    p -> ni       = ni;

    oldmesh       = dvector(0, np-1);
    p -> mesh     = dvector(0, ni-1);
    p -> deriv    = dmatrix(0, np-1, 0, ni-1);
    p -> derivT   = dmatrix(0, ni-1, 0, np-1);
    
    if (oldbasis == newbasis && np == ni) {	/* Meshes are the same; */
      p -> interp  = (double **) NULL;          /* no interpolation.    */
      p -> interpT = (double **) NULL;
    } else {
      p -> interp  = dmatrix(0, np-1, 0, ni-1);
      p -> interpT = dmatrix(0, ni-1, 0, np-1);
    }   

    switch (oldbasis) {
    case GLL:
      switch (newbasis) {
      case GLL:
	if (np==ni) {		/* Meshes are the same.          */
	  jacgl(np-1, 0.0, 0.0, p->mesh);
	  dgll (np, p->mesh, p->deriv, p->derivT);
	} else {		/* Same basis, different orders. */
	  jacgl    (np-1, 0.0, 0.0, oldmesh);
	  jacgl    (np-1, 0.0, 0.0, p->mesh);
	  intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	  dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
	}
	break;
      case STANDARD: default:
	jacgl    (np-1, 0.0, 0.0, oldmesh);
	uniknot  (ni,             p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
	break;
      }
      break;
    case STANDARD: default:
      switch (newbasis) {
      case STANDARD:
	if (np==ni) {		/* Meshes are the same.          */
	  uniknot  (np, p->mesh);
	  dermat_k (np, p->mesh, p->deriv, p->derivT);
 	} else {		/* Same basis, different orders. */
	  uniknot  (np, oldmesh);
	  uniknot  (ni, p->mesh);
	  intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	  dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
	}
	break;
      case GLL: default:
	uniknot  (np,             oldmesh);
	jacgl    (ni-1, 0.0, 0.0, p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
	break;
      }	
      break;
    }
  }

  freeDvector(oldmesh, 0);

  /* p now points to valid storage: return requested operators. */

  if (mesh) *mesh = p -> mesh;
  if (in)   *in   = p -> interp;
  if (it)   *it   = p -> interpT;
  if (dr)   *dr   = p -> deriv;
  if (dt)   *dt   = p -> derivT;
}
