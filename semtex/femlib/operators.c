/*****************************************************************************
 * operators.c:  operators for mesh points, derivatives, quadratures.
 *
 * All 2D matrices have row-major ordering.
 * Routines for single and double precision are supplied.
 *****************************************************************************/

static char
RCSid_ops[] = "$Id$";

#include <stdio.h>
#include <malloc.h>
#include <femdef.h>
#include <alplib.h>
#include <femlib.h>


typedef struct dquadopr {	/* ---- quadrature operator information  --- */
  integer          rule    ;	/* quadrature rule: GL or LL                 */
  integer          np      ;	/* number of interpolant knot points         */
  integer          nq      ;	/* number of quadrature points               */
  double*          knot    ;	/* Lagrange knots    on [-1, 1]              */
  double*          quad    ;	/* Quadrature points on [-1, 1]              */
  double*          weight  ;	/* Quadrature weights                        */
  double**         interp  ;	/* Interpolant operator: knots-->quadrature  */
  double**         interpT ;	/* Transpose of the above                    */
  double**         deriv   ;	/* Derivative operator:  knots-->quadrature  */
  double**         derivT  ;	/* Transpose of the above                    */
  struct dquadopr* next    ;	/* link to next one                          */
} dQuadOpr;			/* ----------------------------------------- */

static dQuadOpr* dQhead = 0;

typedef struct squadopr {	/* ---- quadrature operator information  --- */
  integer          rule    ;	/* quadrature rule: GL or LL                 */
  integer          np      ;	/* number of interpolant knot points         */
  integer          nq      ;	/* number of quadrature points               */
  float*           knot    ;	/* Lagrange knots    on [-1, 1]              */
  float*           quad    ;	/* Quadrature points on [-1, 1]              */
  float*           weight  ;	/* Quadrature weights                        */
  float**          interp  ;	/* Interpolant operator: knots-->quadrature  */
  float**          interpT ;	/* Transpose of the above                    */
  float**          deriv   ;	/* Derivative operator:  knots-->quadrature  */
  float**          derivT  ;	/* Transpose of the above                    */
  struct squadopr* next    ;	/* link to next one                          */
} sQuadOpr;			/* ----------------------------------------- */

static sQuadOpr *sQhead = 0;

typedef struct dmeshopr {	/* ------- mesh operator information ------- */
  integer          oldbasis;	/* STD or GLL                                */
  integer          newbasis;	/* STD or GLL                                */
  integer          np      ;	/* Number of points on original mesh         */
  integer          ni      ;	/* Number of points in interpolated mesh     */
  double*          mesh    ;	/* Location of points on interpolated mesh   */
  double**         interp  ;	/* Interpolant operator: original-->new      */
  double**         interpT ;	/* Transpose of the above                    */
  double**         deriv   ;	/* Derivative operator: original-->new       */
  double**         derivT  ;	/* Transpose of the above                    */
  struct dmeshopr* next    ;	/* link to next one                          */
} dMeshOpr;			/* ----------------------------------------- */

static dMeshOpr* dMhead = 0;

typedef struct smeshopr {	/* ------- mesh operator information ------- */
  integer          oldbasis;	/* STD or GLL                                */
  integer          newbasis;	/* STD or GLL                                */
  integer          np      ;	/* Number of points on original mesh         */
  integer          ni      ;	/* Number of points in interpolated mesh     */
  float*           mesh    ;	/* Location of points on interpolated mesh   */
  float**          interp  ;	/* Interpolant operator: original-->new      */
  float**          interpT ;	/* Transpose of the above                    */
  float**          deriv   ;	/* Derivative operator: original-->new       */
  float**          derivT  ;	/* Transpose of the above                    */
  struct smeshopr* next    ;	/* link to next one                          */
} sMeshOpr;			/* ----------------------------------------- */

static sMeshOpr* sMhead = 0;



void dQuadOps (const integer   rule, /* input: quadrature rule: GL or LL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   nq  , /* input: number of quadrature points   */
	       const double**  kp  , /* pointer to knot point storage        */
	       const double**  qp  , /* pointer to quadrature point storage  */
	       const double**  qw  , /* pointer to quadrature weight storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  ) /* pointer to transposed deriv matrix   */
/* ------------------------------------------------------------------------- *
 * Maintain/return QUADRATURE operators for finite elements with GLL BASIS
 * FUNCTIONS.  Quadrature rule may be at nodes (LL, nq = np enforced) or at
 * Gauss(-Legendre) points (GL).
 *
 * Operators are defined on the interval [-1, 1].
 *
 * If the required operators are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 * ------------------------------------------------------------------------- */
{
  char               routine[] = "dQuadOps";
  register integer   found = 0;
  register dQuadOpr* p;

  for (p = dQhead; p; p=p->next) {
    found = p->rule == rule && p->np == np && p->nq == nq;
    if (found) break;
  }

  if (!found) {

    if (rule != LL && rule != GL)
      message (routine, "unrecognized quadrature rule", ERROR);

    p = (dQuadOpr *) calloc (1, sizeof (dQuadOpr));
    if (dQhead) p -> next = dQhead;
    dQhead = p;

    p -> rule = rule;
    p -> np   = np;
    p -> nq   = (rule == GL) ? nq : np;

    if (rule == LL && np != nq)
      message (routine, "np != nq in LL rule...enforcing equality", WARNING);

    if (rule == LL) {

      p -> knot    = dvector (0, np-1);
      p -> quad    = p -> knot;
      p -> weight  = dvector (0, np-1);
      p -> interp  = (double**) 0;	/* No interpolation needed. */
      p -> interpT = (double**) 0;
      p -> deriv   = dmatrix (0, np-1, 0, np-1);
      p -> derivT  = dmatrix (0, np-1, 0, np-1);

      zwgll (p->knot, p->weight, np);
      dgll  (np, p->knot, p->deriv, p->derivT);

    } else {

      p -> knot    = dvector (0, np-1);
      p -> quad    = dvector (0, nq-1);
      p -> weight  = dvector (0, nq-1);
      p -> interp  = dmatrix (0, nq-1, 0, np-1);
      p -> interpT = dmatrix (0, np-1, 0, nq-1);
      p -> deriv   = dmatrix (0, nq-1, 0, np-1);
      p -> derivT  = dmatrix (0, np-1, 0, nq-1);

      jacgl    (np-1, 0.0, 0.0, p->knot);     /* Knots at Lobatto points.    */
      zwgl     (p->quad, p->weight, nq);      /* Quadrature at Gauss points. */
      intmat_g (np, p->knot, nq, p->quad, p->interp, p->interpT);
      dermat_g (np, p->knot, nq, p->quad, p->deriv,  p->derivT );

    }
  }

  /* -- p now points to valid storage: return requested operators. */

  if (kp) *kp = (const double* ) p -> knot;
  if (qp) *qp = (const double* ) p -> quad;
  if (qw) *qw = (const double* ) p -> weight;
  if (in) *in = (const double**) p -> interp;
  if (it) *it = (const double**) p -> interpT;
  if (dr) *dr = (const double**) p -> deriv;
  if (dt) *dt = (const double**) p -> derivT;
}


void dMeshOps (const integer   old , /* input: element basis: STD or GLL     */
	       const integer   new , /* input: desired basis: STD or GLL     */
	       const integer   np  , /* input: number of knot points         */
	       const integer   ni  , /* input: number of interpolant points  */
	       const double**  mesh, /* pointer to interpolated mesh storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  ) /* pointer to transposed deriv matrix   */
/* ------------------------------------------------------------------------- *
 * Maintain/return INTERPOLATION operators for STD and GLL meshes.
 * Operators are defined on the interval [-1, 1].
 *
 * It is legal to ask for interpolation onto the same mesh as input, in
 * which case the mesh storage and interpolation matrices are NULL.
 *
 * If the required operators are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 * ------------------------------------------------------------------------- */
{
  char               routine[] = "dMeshOps";
  register integer   found = 0;
  register dMeshOpr* p;
  double*            oldmesh;

  for (p = dMhead; p; p = p->next) {
    found = p->oldbasis == old
      &&    p->newbasis == new 
      &&    p->np       == np
      &&    p->ni       == ni;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */

    p = (dMeshOpr *) calloc (1, sizeof (dMeshOpr));
    if (dMhead) p -> next = dMhead;
    dMhead = p;

    p -> oldbasis = old;
    p -> newbasis = new;
    p -> np       = np;
    p -> ni       = ni;

    oldmesh       = dvector (0, np-1);
    p -> mesh     = dvector (0, ni-1);
    p -> deriv    = dmatrix (0, ni-1, 0, np-1);
    p -> derivT   = dmatrix (0, np-1, 0, ni-1);
    
    if (old == new && np == ni) {	/* Meshes are the same; */
      p -> interp  = (double**) 0;      /* no interpolation.    */
      p -> interpT = (double**) 0;
    } else {
      p -> interp  = dmatrix (0, ni-1, 0, np-1);
      p -> interpT = dmatrix (0, np-1, 0, ni-1);
    }   

    if (old == GLL && new == GLL) {

      if (np == ni) {  /* Meshes are the same. */
	jacgl (np-1, 0.0, 0.0, p->mesh);
	dgll  (np, p->mesh, p->deriv, p->derivT);
      } else {
	jacgl    (np-1, 0.0, 0.0, oldmesh);
	jacgl    (ni-1, 0.0, 0.0, p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
      }

    } else if (old == GLL && new == STD) {

      jacgl    (np-1, 0.0, 0.0,  oldmesh);
      uniknot  (ni,              p->mesh);
      intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
      dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );

    } else if (old == STD && new == STD) {

      if (np == ni) {  /* Meshes are the same. */
	uniknot  (np, p->mesh);
	dermat_k (np, p->mesh, p->deriv, p->derivT);
      } else {
	uniknot  (np, oldmesh);
	uniknot  (ni, p->mesh);
	intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
	dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );
      }

    } else if (old == STD && new == GLL) {

      uniknot  (np,              oldmesh);
      jacgl    (ni-1, 0.0, 0.0,  p->mesh);
      intmat_g (np, oldmesh, ni, p->mesh, p->interp, p->interpT);
      dermat_g (np, oldmesh, ni, p->mesh, p->deriv,  p->derivT );

    } else

      message (routine, "basis function unrecognized as STD or GLL", ERROR);

    freeDvector (oldmesh, 0);
  }

  /* p now points to valid storage: return requested operators. */

  if (mesh) *mesh = (const double* ) p -> mesh;
  if (in)   *in   = (const double**) p -> interp;
  if (it)   *it   = (const double**) p -> interpT;
  if (dr)   *dr   = (const double**) p -> deriv;
  if (dt)   *dt   = (const double**) p -> derivT;
}


void sQuadOps (const integer   rule ,
	       const integer   np   ,
	       const integer   nq   ,
	       const float**   kp   ,
	       const float**   qp   ,
	       const float**   qw   ,
	       const float***  in   ,
	       const float***  it   ,
	       const float***  dr   ,
	       const float***  dt   )
/* ------------------------------------------------------------------------- *
 * As for dQuadOps, but here single-precision operators are made.
 * ------------------------------------------------------------------------- */
{
  char               routine[] = "sQuadOps";
  register integer   found = 0;
  register sQuadOpr* p;

  for (p = sQhead; p; p = p->next) {
    found = p->rule == rule && p->np == np && p->nq == nq;
    if (found) break;
  }

  if (!found) {

    if (rule != LL && rule != GL)
      message (routine, "unrecognized quadrature rule", ERROR);

    p = (sQuadOpr *) calloc (1, sizeof (sQuadOpr));
    if (sQhead) p -> next = sQhead;
    sQhead = p;

    p -> rule = rule;
    p -> np   = np;
    p -> nq   = (rule == GL) ? nq : np;

    if (rule == LL && np != nq)
      message (routine, "np != nq in LL rule...enforcing equality", WARNING);

    if (rule == LL) {
      
      p -> knot    = svector (0, np-1);
      p -> quad    = p -> knot;
      p -> weight  = svector (0, np-1);
      p -> interp  = (float**) 0;	/* No interpolation needed. */
      p -> interpT = (float**) 0;
      p -> deriv   = smatrix (0, np-1, 0, np-1);
      p -> derivT  = smatrix (0, np-1, 0, np-1);

      {
	double*  dk = dvector (0, np-1);
	double*  dw = dvector (0, np-1);
	double** dd = dmatrix (0, np-1, 0, np-1);
	double** dt = dmatrix (0, np-1, 0, np-1);

	zwgll (dk, dw, np);
	dgll  (np, dk, dd, dt);

	vsngl (np,     dk, 1,  p -> knot,   1);
	vsngl (np,     dw, 1,  p -> weight, 1);
	vsngl (np*np, *dd, 1, *p -> deriv,  1);
	vsngl (np*np, *dt, 1, *p -> derivT, 1);
	
	freeDvector (dk, 0);
	freeDvector (dw, 0);
	freeDmatrix (dd, 0, 0);
	freeDmatrix (dt, 0, 0);
      }
    } else {

      p -> knot    = svector (0, np-1);
      p -> quad    = svector (0, nq-1);
      p -> weight  = svector (0, nq-1);
      p -> interp  = smatrix (0, nq-1, 0, np-1);
      p -> interpT = smatrix (0, np-1, 0, nq-1);
      p -> deriv   = smatrix (0, nq-1, 0, np-1);
      p -> derivT  = smatrix (0, np-1, 0, nq-1);

      {
	double*  dk  = dvector (0, np-1);
	double*  dq  = dvector (0, nq-1);
	double*  dw  = dvector (0, nq-1);
	double** din = dmatrix (0, nq-1, 0, np-1);
	double** dit = dmatrix (0, np-1, 0, nq-1);
	double** ddv = dmatrix (0, nq-1, 0, np-1);
	double** ddt = dmatrix (0, np-1, 0, nq-1);

	jacgl    (np-1, 0.0, 0.0, dk);
	zwgl     (dq, dw, nq);
	intmat_g (np, dk, nq, dq, din, dit);
	dermat_g (np, dk, nq, dq, ddv, ddt);

	vsngl (np,     dk,  1,  p -> knot,    1);
	vsngl (nq,     dq,  1,  p -> quad,    1);
	vsngl (nq,     dw,  1,  p -> weight,  1);
	vsngl (nq*np, *din, 1, *p -> interp,  1);
	vsngl (np*nq, *dit, 1, *p -> interpT, 1);
	vsngl (nq*np, *ddv, 1, *p -> deriv,   1);
	vsngl (np*nq, *ddt, 1, *p -> derivT,  1);

	freeDvector (dk,  0);
	freeDvector (dq,  0);
	freeDvector (dw,  0);
	freeDmatrix (din, 0, 0);
	freeDmatrix (dit, 0, 0);
	freeDmatrix (ddv, 0, 0);
	freeDmatrix (ddt, 0, 0);
      }
    }
  }

  if (kp) *kp = (const float* ) p -> knot;
  if (qp) *qp = (const float* ) p -> quad;
  if (qw) *qw = (const float* ) p -> weight;
  if (in) *in = (const float**) p -> interp;
  if (it) *it = (const float**) p -> interpT;
  if (dr) *dr = (const float**) p -> deriv;
  if (dt) *dt = (const float**) p -> derivT;

}


void sMeshOps (const integer   old ,
	       const integer   new ,
	       const integer   np  ,
	       const integer   ni  ,
	       const float**   mesh,
	       const float***  in  ,
	       const float***  it  ,
	       const float***  dr  ,
	       const float***  dt  )
/* ------------------------------------------------------------------------- *
 * As for dMeshOps, but make single-precision operators.
 * ------------------------------------------------------------------------- */
{
  char               routine[] = "sMeshOps";
  register integer   found = 0;
  register sMeshOpr* p;

  for (p = sMhead; p; p = p->next) {
    found = p->oldbasis == old
      &&    p->newbasis == new 
      &&    p->np       == np
      &&    p->ni       == ni;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */

    p = (sMeshOpr *) calloc (1, sizeof (sMeshOpr));
    if (sMhead) p -> next = sMhead;
    sMhead = p;

    p -> oldbasis = old;
    p -> newbasis = new;
    p -> np       = np;
    p -> ni       = ni;

    p -> mesh     = svector (0, ni-1);
    p -> deriv    = smatrix (0, ni-1, 0, np-1);
    p -> derivT   = smatrix (0, np-1, 0, ni-1);
    
    if (old == new && np == ni) {	/* Meshes are the same; */
      p -> interp  = (float**) 0;     /* no interpolation.    */
      p -> interpT = (float**) 0;
    } else {
      p -> interp  = smatrix (0, ni-1, 0, np-1);
      p -> interpT = smatrix (0, np-1, 0, ni-1);
    }   

    if (old == GLL && new == GLL) {

      if (np == ni) {  /* Meshes are the same; return derivative matrices. */

	double*  dm  = dvector (0, np-1);
	double** ddv = dmatrix (0, ni-1, 0, np-1);
	double** ddt = dmatrix (0, np-1, 0, ni-1);

	jacgl (np-1, 0.0, 0.0, dm);
	dgll  (np, dm, ddv, ddt);

	vsngl (np,     dm,  1,  p -> mesh,   1);
	vsngl (ni*np, *ddv, 1, *p -> deriv,  1);
	vsngl (np*ni, *ddt, 1, *p -> derivT, 1);

	freeDvector (dm,  0);
	freeDmatrix (ddv, 0, 0);
	freeDmatrix (ddt, 0, 0);

      } else {

	double*  om  = dvector (0, np-1);
	double*  dm  = dvector (0, ni-1);
	double** ddv = dmatrix (0, ni-1, 0, np-1);
	double** ddt = dmatrix (0, np-1, 0, ni-1);
	double** din = dmatrix (0, ni-1, 0, np-1);
	double** dit = dmatrix (0, np-1, 0, ni-1);
	
	jacgl    (np-1, 0.0, 0.0, om);
	jacgl    (ni-1, 0.0, 0.0, dm);
	intmat_g (np, om, ni, dm, din, dit);
	dermat_g (np, om, ni, dm, ddv, dit);

	vsngl (np,     dm,  1,  p -> mesh,    1);
	vsngl (ni*np, *ddv, 1, *p -> deriv,   1);
	vsngl (np*ni, *ddt, 1, *p -> derivT,  1);
	vsngl (ni*np, *din, 1, *p -> interp,  1);
	vsngl (np*ni, *dit, 1, *p -> interpT, 1);

	freeDvector (om,  0);
	freeDvector (dm,  0);
	freeDmatrix (ddv, 0, 0);
	freeDmatrix (ddt, 0, 0);
	freeDmatrix (din, 0, 0);
	freeDmatrix (dit, 0, 0);
      }

    } else if (old == GLL && new == STD) {

      double*  om  = dvector (0, np-1);
      double*  dm  = dvector (0, ni-1);
      double** ddv = dmatrix (0, ni-1, 0, np-1);
      double** ddt = dmatrix (0, np-1, 0, ni-1);
      double** din = dmatrix (0, ni-1, 0, np-1);
      double** dit = dmatrix (0, np-1, 0, ni-1);
      
      jacgl    (np-1, 0.0, 0.0,  om);
      uniknot  (ni,              dm);
      intmat_g (np, om, ni, dm, din, dit);
      dermat_g (np, om, ni, dm, ddv, ddt);

      vsngl (np,     dm,  1,  p -> mesh,    1);
      vsngl (ni*np, *ddv, 1, *p -> deriv,   1);
      vsngl (np*ni, *ddt, 1, *p -> derivT,  1);
      vsngl (ni*np, *din, 1, *p -> interp,  1);
      vsngl (np*ni, *dit, 1, *p -> interpT, 1);

      freeDvector (om,  0);
      freeDvector (dm,  0);
      freeDmatrix (ddv, 0, 0);
      freeDmatrix (ddt, 0, 0);
      freeDmatrix (din, 0, 0);
      freeDmatrix (dit, 0, 0);

    } else if (old == STD && new == STD) {

      if (np == ni) {  /* Meshes are the same. */

	double*  dm  = dvector (0, np-1);
	double** ddv = dmatrix (0, ni-1, 0, np-1);
	double** ddt = dmatrix (0, np-1, 0, ni-1);

	uniknot  (np, dm);
	dermat_k (np, dm, ddv, ddt);

	vsngl (np,     dm,  1,  p -> mesh,   1);
	vsngl (ni*np, *ddv, 1, *p -> deriv,  1);
	vsngl (np*ni, *ddt, 1, *p -> derivT, 1);

	freeDvector (dm,  0);
	freeDmatrix (ddv, 0, 0);
	freeDmatrix (ddt, 0, 0);

      } else {

	double*  om  = dvector (0, np-1);
	double*  dm  = dvector (0, ni-1);
	double** ddv = dmatrix (0, ni-1, 0, np-1);
	double** ddt = dmatrix (0, np-1, 0, ni-1);
	double** din = dmatrix (0, ni-1, 0, np-1);
	double** dit = dmatrix (0, np-1, 0, ni-1);

	uniknot  (np, om);
	uniknot  (ni, dm);
	intmat_g (np, om, ni, dm, din, dit);
	dermat_g (np, om, ni, dm, ddv, ddt);

	vsngl (np,     dm,  1,  p -> mesh,    1);
	vsngl (ni*np, *ddv, 1, *p -> deriv,   1);
	vsngl (np*ni, *ddt, 1, *p -> derivT,  1);
	vsngl (ni*np, *din, 1, *p -> interp,  1);
	vsngl (np*ni, *dit, 1, *p -> interpT, 1);

	freeDvector (om,  0);
	freeDvector (dm,  0);
	freeDmatrix (ddv, 0, 0);
	freeDmatrix (ddt, 0, 0);
	freeDmatrix (din, 0, 0);
	freeDmatrix (dit, 0, 0);
      }

    } else if (old == STD && new == GLL) {

      double*  om  = dvector (0, np-1);
      double*  dm  = dvector (0, ni-1);
      double** ddv = dmatrix (0, ni-1, 0, np-1);
      double** ddt = dmatrix (0, np-1, 0, ni-1);
      double** din = dmatrix (0, ni-1, 0, np-1);
      double** dit = dmatrix (0, np-1, 0, ni-1);

      uniknot  (np,              om);
      jacgl    (ni-1, 0.0, 0.0,  dm);
      intmat_g (np, om, ni, dm, din, dit);
      dermat_g (np, om, ni, dm, ddv, ddt);

      vsngl (np,     dm,  1,  p -> mesh,    1);
      vsngl (ni*np, *ddv, 1, *p -> deriv,   1);
      vsngl (np*ni, *ddt, 1, *p -> derivT,  1);
      vsngl (ni*np, *din, 1, *p -> interp,  1);
      vsngl (np*ni, *dit, 1, *p -> interpT, 1);

      freeDvector (om,  0);
      freeDvector (dm,  0);
      freeDmatrix (ddv, 0, 0);
      freeDmatrix (ddt, 0, 0);
      freeDmatrix (din, 0, 0);
      freeDmatrix (dit, 0, 0);

    } else

      message (routine, "basis function unrecognized as STD or GLL", ERROR);
  }

  /* -- p now points to valid storage: return requested operators. */

  if (mesh) *mesh = (const float* ) p -> mesh;
  if (in)   *in   = (const float**) p -> interp;
  if (it)   *it   = (const float**) p -> interpT;
  if (dr)   *dr   = (const float**) p -> deriv;
  if (dt)   *dt   = (const float**) p -> derivT;
}


void dIntpOps (const integer basis,  /* element basis: STD or GLL            */
	       const integer np   ,  /* number of knot points                */
	       const double  r    ,  /* location of r in [-1, 1]             */
	       const double  s    ,  /* location of s in [-1, 1]             */
	       double*       inr  ,  /* 1D shape function at r               */
	       double*       ins  ,  /* 1D shape function at s               */
	       double*       dvr  ,  /* 1D shape function derivative at r    */
	       double*       dvs  )  /* 1D shape function derivative at s    */
/* ------------------------------------------------------------------------- *
 * Return the interpolation/derivative vectors for tensor-product shape
 * functions evaluated at canonical coordinates (r, s), each in [-1. 1].
 *
 * Take no action for empty vectors.
 * ------------------------------------------------------------------------- */
{
  const double* kp;
  double        x[1], **dv, **dt;

  dv = dmatrix (0, 0, 0, np - 1);
  dt = dmatrix (0, np - 1, 0, 0);

  dQuadOps (basis, np, np, &kp, 0, 0, 0, 0, 0, 0);

  x[0] = r;
  if (inr) {
    intmat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, inr, 1);
  }
  if (dvr) {
    dermat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, dvr, 1);
  }

  x[0] = s;
  if (ins) {
    intmat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, ins, 1);
  }
  if (dvs) {
    dermat_g (np, kp, 1, x, dv, dt);
    dcopy    (np, *dv, 1, dvs, 1);
  }

  freeDmatrix (dv, 0, 0);
  freeDmatrix (dt, 0, 0);
}


void sIntpOps (const integer basis,
	       const integer np   ,
	       const float   r    ,
	       const float   s    ,
	       float*        inr  ,
	       float*        ins  ,
	       float*        dvr  ,
	       float*        dvs  )
/* ------------------------------------------------------------------------- *
 * Single-precision version of dIntpOps.
 * ------------------------------------------------------------------------- */
{
  const double* kp;
  double        x[1], **dv, **dt;

  dv = dmatrix (0, 0, 0, np - 1);
  dt = dmatrix (0, np - 1, 0, 0);

  dQuadOps (basis, np, np, &kp, 0, 0, 0, 0, 0, 0);

  x[0] = r;
  if (inr) {
    intmat_g (np, kp, 1, x, dv, dt);
    vsngl    (np, *dv, 1, inr, 1);
  }
  if (dvr) {
    dermat_g (np, kp, 1, x, dv, dt);
    vsngl    (np, *dv, 1, dvr, 1);
  }

  x[0] = s;
  if (ins) {
    intmat_g (np, kp, 1, x, dv, dt);
    vsngl    (np, *dv, 1, ins, 1);
  }
  if (dvs) {
    dermat_g (np, kp, 1, x, dv, dt);
    vsngl    (np, *dv, 1, dvs, 1);
  }

  freeDmatrix (dv, 0, 0);
  freeDmatrix (dt, 0, 0);
}
