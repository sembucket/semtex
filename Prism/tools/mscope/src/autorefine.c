/*
 * AutoRefine() 
 *
 * Compute the solution to a Helmholtz problem by adaptive refinement of
 * the computational domain.  Refinement is based on the L2 norm of the
 * solution gradient.
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "veclib/veclib.h"
#include "cubit/cubit.h"

#include "mscope.h"


/* ------------------------------------------------------------------------- *
 * AutoRefine()  -- Driver for the adaptive solution                         *
 * ------------------------------------------------------------------------- */
 
void AutoRefine 
(Domain *problem, double tol, int maxDepth, int numPasses, estimate_t type)
{
#if 0
  Mesh   *mesh  = problem->mesh;
  Matrix *A     = problem->A;
  Field  *U     = problem->U;
  Field  *F     = problem->F;

  error_t *eps  = Error_alloc (mesh, type);
  double   sum  = 0.;
  int      nref = 0;

  do {
    Matrix_update (A, mesh);
    Field_set     (F, problem->force->info[0]);
    Solve         (U, F, mesh, A);

    sum  = Error_compute (eps, U);
    Error_info (eps, tol, stdout);
    nref = Error_adapt (eps, tol, maxDepth);

    DoErase();    /* Redraw the mesh */
    DoGrid ();

  } while (nref > 0 && --numPasses);

  Error_free (eps);
#endif
}

