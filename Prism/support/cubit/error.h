#ifndef ERROR_H
#define ERROR_H

/* Error Estimates
 *
 * $Revision$
 *
 * Author:     R. D. Henderson
 *
 * This is a generic class of "error estimates" for scalar fields.  The 
 * error estimate can be used to drive adaptive mesh generation.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "cubit.h"
#include "mesh.h"

/* Macros */

#define ERROR_MAX(e) (e->max)
#define ERROR_MIN(e) (e->min)
#define ERROR_SUM(d) (d->sum)

typedef enum {       /* Various error estimates */
  E_SPECTRUM,        /* Refine based on the local polynomial spectrum */
  E_GRADIENT,        /* Refine based on the local solution gradient */
  E_REGRESS,         /* Refine based on spectral regression */
  E_DEFAULT          /* Dummy type */
} estimate_t;

typedef struct {     
  int      depth;        /* depth in the tree                      */
  double   size;         /* element size                           */
  double   error;        /* error estimate                         */

  struct {               /* local solution norm                    */
    double inf;
    double l2;
    double h1;
  } norm;

} eps_t;

typedef struct {         /* -------------------------------------- */
  double min;            /* minimum local error                    */
  double max;            /* maximum  "   "   "                     */
  double sum;            /* sum of all local errors (global error) */
  double scale;          /* scale for error estimate               */

  struct {               /* norm of the function we're refining on */
    double inf;
    double l2;
    double h1;
  } norm;

  eps_t  eps[_MAX_NEL];  /* array of local error estimates         */

  estimate_t  type;      /* type of error estimate                 */
  Mesh       *mesh;      /* mesh that defines the solution space   */
} error_t;

/* ------------------------------------------------------------------------- */


/* Prototypes */

error_t *Error_alloc   (Mesh *mesh, estimate_t type);
void     Error_free    (error_t *error);
void     Error_info    (const error_t *error, double tol, FILE *fp);

/* Compute the error estimate based on a given Field; returns global error */

double Error_compute  (error_t *error, const Field *u);

/* Driver for mesh adaption.  Return the number of mesh elements that were   *
 * refined during one pass of the adaptor.  Each element with a local error  *
 * greater than tol and a tree depth less than maxDepth is refined.          */

int Error_adapt   (error_t *error, double tol, int maxDepth);

#endif
