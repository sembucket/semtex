/*****************************************************************************
 * OPERATORS.C:  operators for mesh points, derivatives, quadratures.        *
 *****************************************************************************/

static char
RCSid02[] = "$Id$";

#include <stdio.h>
#include <malloc.h>
#include <alplib.h>
#include <femdef.h>
#include <femlib.h>


typedef struct quadopr {	/* ---- quadrature operator information  --- */
  int             rule    ;	/* quadrature rule: GL or LL                 */
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

static QuadOpr *Qhead = 0;


typedef struct meshopr {	/* ------- mesh operator information ------- */
  int             oldbasis;	/* STD or GLL                                */
  int             newbasis;	/* STD or GLL                                */
  int             np      ;	/* Number of points on original mesh         */
  int             ni      ;	/* Number of points in interpolated mesh     */
  double         *mesh    ;	/* Location of points on interpolated mesh   */
  double        **interp  ;	/* Interpolant operator: original-->new      */
  double        **interpT ;	/* Transpose of the above                    */
  double        **deriv   ;	/* Derivative operator: original-->new       */
  double        **derivT  ;	/* Transpose of the above                    */
  struct meshopr *next    ;	/* link to next one                          */
} MeshOpr;			/* ----------------------------------------- */

static MeshOpr *Mhead = 0;



void  quadOps(int rule     ,	/* input: quadrature rule: GL or LL          */
	      int np       ,	/* input: number of knot points              */
	      int nq       ,	/* input: number of quadrature points        */
	      double  **kp ,	/* pointer to knot point storage             */
	      double  **qp ,	/* pointer to quadrature point storage       */
	      double  **qw ,	/* pointer to quadrature weight storage      */
	      double ***in ,	/* pointer to interpolation matrix           */
	      double ***it ,	/* pointer to transposed interpolant matrix  */
	      double ***dr ,	/* pointer to derivative matrix              */
	      double ***dt )	/* pointer to transposed derivative matrix   */
/* ========================================================================= *
 * Maintain/return QUADRATURE operators for finite elements with GLL BASIS   *
 * FUNCTIONS.  Quadrature rule may be at nodes (LL, nq = np enforced) or at  *
 * Gauss(-Legendre) points (GL).                                             *
 *                                                                           *
 * Operators are defined on the interval [-1, 1].                            *
 *                                                                           *
 * If the required operators are "in stock", they are returned from a list,  *
 * otherwise they are created and added to the list as a first step.         *
 * NULL input pointers are not assigned a value.                             *
 *                                                                           *
 * ========================================================================= */
{
  char      routine[] = "quadOps";
  int       found   = 0;
  QuadOpr  *p;


  for (p=Qhead; p; p=p->next) {
    found = p->rule == rule && p->np == np && p->nq == nq;
    if (found) break;
  }

  if (!found) {

    if (rule != LL && rule != GL)
      message(routine, "unrecognized quadrature rule", ERROR);

    p = (QuadOpr *) calloc(1, sizeof(QuadOpr));
    if (Qhead) p -> next = Qhead;
    Qhead = p;

    p -> rule = rule;
    p -> np   = np;
    p -> nq   = (rule == GL) ? nq : np;

    if (rule == LL && np != nq)
      message(routine, "np != nq in LL rule...enforcing", WARNING);

    if (rule == LL) {

      p -> knot    = dvector(0, np-1);
      p -> quad    = p -> knot;
      p -> weight  = dvector(0, np-1);
      p -> interp  = (double **) 0;	/* No interpolation needed. */
      p -> interpT = (double **) 0;
      p -> deriv   = dmatrix(0, np-1, 0, np-1);
      p -> derivT  = dmatrix(0, np-1, 0, np-1);

      zwgll (p->knot, p->weight, np);
      dgll  (np, p->knot, p->deriv, p->derivT);

    } else {

      p -> knot    = dvector(0, np-1);
      p -> quad    = dvector(0, nq-1);
      p -> weight  = dvector(0, nq-1);
      p -> interp  = dmatrix(0, nq-1, 0, np-1);
      p -> interpT = dmatrix(0, np-1, 0, nq-1);
      p -> deriv   = dmatrix(0, nq-1, 0, np-1);
      p -> derivT  = dmatrix(0, np-1, 0, nq-1);

      jacgl    (np-1, 0.0, 0.0, p->knot);     /* Knots at Lobatto points.    */
      zwgl     (p->quad, p->weight, nq);      /* Quadrature at Gauss points. */
      intmat_g (np, p->knot, nq, p->quad, p->interp, p->interpT);
      dermat_g (np, p->knot, nq, p->quad, p->deriv,  p->derivT );

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





void  meshOps(int oldbasis   ,	/* input: element basis: STD or GLL          */
	      int newbasis   ,	/* input: desired basis: STD or GLL          */
	      int np         ,	/* input: number of knot points              */
	      int ni         ,	/* input: number of interpolant points       */
	      double  **mesh ,	/* pointer to interpolated mesh storage      */
	      double ***in   ,	/* pointer to interpolation matrix           */
	      double ***it   ,	/* pointer to transposed interpolant matrix  */
	      double ***dr   ,	/* pointer to derivative matrix              */
	      double ***dt   )	/* pointer to transposed derivative matrix   */
/* ========================================================================= *
 * Maintain/return INTERPOLATION operators for STD and GLL meshes.           *
 * Operators are defined on the interval [-1, 1].                            *
 *                                                                           *
 * It is legal to ask for interpolation onto the same mesh as input, in      *
 * which case the mesh storage and interpolation matrices are NULL.          *
 *                                                                           *
 * If the required operators are "in stock", they are returned from a list,  *
 * otherwise they are created and added to the list as a first step.         *
 * NULL input pointers are not assigned a value.                             *
 *                                                                           *
 * ========================================================================= */
{
  char      routine[] = "meshOps()";
  int       found = 0;
  MeshOpr  *p;
  double   *oldmesh;


  for (p=Mhead; p; p=p->next) {
    found = p->oldbasis == oldbasis
      &&    p->newbasis == newbasis 
      &&    p->np       == np
      &&    p->ni       == ni;
    if (found) break;
  }


  if (!found) {		/* Make more storage and operators. */

    p = (MeshOpr *) calloc(1, sizeof(MeshOpr));
    if (Mhead) p -> next = Mhead;
    Mhead = p;

    p -> oldbasis = oldbasis;
    p -> newbasis = newbasis;
    p -> np       = np;
    p -> ni       = ni;

    oldmesh       = dvector(0, np-1);
    p -> mesh     = dvector(0, ni-1);
    p -> deriv    = dmatrix(0, ni-1, 0, np-1);
    p -> derivT   = dmatrix(0, np-1, 0, ni-1);
    
    if (oldbasis == newbasis && np == ni) {	/* Meshes are the same; */
      p -> interp  = (double **) 0;             /* no interpolation.    */
      p -> interpT = (double **) 0;
    } else {
      p -> interp  = dmatrix(0, ni-1, 0, np-1);
      p -> interpT = dmatrix(0, np-1, 0, ni-1);
    }   

    if (oldbasis == GLL && newbasis == GLL) {

      if (np==ni) {  /* Meshes are the same. */
	jacgl (np-1, 0.0, 0.0, p->mesh);
	dgll  (np, p->mesh, p->deriv, p->derivT);
      } else {
	jacgl    (np-1, 0.0, 0.0, oldmesh);
	jacgl    (ni-1, 0.0, 0.0, p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
      }

    } else if (oldbasis == GLL && newbasis == STD) {

      jacgl    (np-1, 0.0, 0.0, oldmesh);
      uniknot  (ni,             p->mesh);
      intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
      dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );

    } else if (oldbasis == STD && newbasis == STD) {

      if (np==ni) {  /* Meshes are the same. */
	uniknot  (np, p->mesh);
	dermat_k (np, p->mesh, p->deriv, p->derivT);
      } else {
	uniknot  (np, oldmesh);
	uniknot  (ni, p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
      }

    } else if (oldbasis == STD && newbasis == GLL) {

      uniknot  (np,             oldmesh);
      jacgl    (ni-1, 0.0, 0.0, p->mesh);
      intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
      dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );

    } else

      message(routine, "basis function unrecognized as STD or GLL", ERROR);

    freeDvector(oldmesh, 0);
  }

  /* p now points to valid storage: return requested operators. */

  if (mesh) *mesh = p -> mesh;
  if (in)   *in   = p -> interp;
  if (it)   *it   = p -> interpT;
  if (dr)   *dr   = p -> deriv;
  if (dt)   *dt   = p -> derivT;
}
