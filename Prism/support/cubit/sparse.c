/*
 * Sparse Matrix Solver
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

#include "vdb/vdb.h"
#include "veclib/veclib.h"

#include "cubit.h"
#include "cg.h"

static void show_vector (int, double*);

#define _MAX_NP (_MAX_NORDER * _MAX_NORDER)

/* External */

static Mesh   *A_mesh;
static Matrix *A_matrix;

/* Private Functions */

int elmt_get_data (Element *elmt, double *x, double *xloc);
int elmt_sum_data (Element *elmt, double *x, double *xloc);
static Helmholtz  (Element *elmt, Matrix *A, double *x, double *Ax);

/* Preconditioners */


static void M_diagonal (int n, const double *x, double *y)
{
  int i;
  for(i = 0; i < n; i++)
    y[i] = A_matrix->diag[i] * x[i];
}

static void M_identity (int n, const double *x, double *y)
{
  int i;
  for(i = 0; i < n; i++)
    y[i] = x[i];
}


/* Compute the 2-norm, |u| = sqrt(dot(u,u)) */

double norm (int n, double *u) {
  double inner_product();
  return sqrt(inner_product(n, u, u));
}

/* ------------------------------------------------------------------------- *
 * The following function solves the system:                                 *
 *                                                                           *
 *             D[u,{x,2}] - LAMBDA u + f = 0,                                *
 *                                                                           *
 * where LAMBDA > 0 is defined as a parameter and f is a given field.        *
 * ------------------------------------------------------------------------- */

void Solve (Field *U, Field *F, Mesh *mesh, Matrix *A)
{
  Element *elmt;

  const int    n   = A->bpts + A->ipts;
  const double tau = dparam("TOLABS");

  int max = n;
  double *g = (double*) calloc (n, sizeof(double));
  double *r = (double*) calloc (n, sizeof(double));
  double *v = (double*) calloc (n, sizeof(double));

  double tol   = 0.;
  double normb = 1.;
  double normf = Field_L2(F);
  if (normf < FLT_EPSILON)
    normf = 1.0;

  A_mesh   = mesh;
  A_matrix = A;

  Form_G   (U, mesh, A, g);
  Form_RHS (F, mesh, A, Laplace_IP, g, r);

  if (option("singular")) dsadd(n, -dsum(n,r,1)/n, r, 1, r, 1);

  /* NB: We really want || A(u) - F || < tau * || F ||, but the input to the *
  *  conjugate gradient solver already has the initial guess in it.  There-  *
  *  fore, the error tolerance is adjusted so that we don't try to converge  *
  *  to anything lower than the tolerance "tau" given above.                 */

  normb = norm(n,r);
  tol   = tau * normf/normb;

  CG    (n, v, r, &max, &tol, norm, inner_product, M_diagonal, Laplace_IP);
  dvadd (n, v, 1, g, 1, g, 1);

  for (elmt = mesh->head; elmt; elmt = elmt->next)
    elmt_get_data (elmt, g, FIELD_DATA(U,elmt->id));


  free (v);
  free (g);
  free (r);

  return;
}


/* ------------------------------------------------------------------------- *
 * Laplace_IP() -- Discrete Laplace Inner Product                            *
 *                                                                           *
 * Compute the effect of the discrete Laplacian A on a vector x.  The input  *
 * is x in flat-field storage, the output is Ax in flat-field storage.  The  *
 * number "n" is the number of degrees-of-freedom in the mesh.               *
 * ------------------------------------------------------------------------- */

void Laplace_IP (int n, double *x, double *Ax)
{
  Element *elmt;
  double   xloc[_MAX_NP], Axloc[_MAX_NP], *Avert, *Aedge;
  int i;

  const int np   = A_mesh->nr - 2;
  const int bpts = A_matrix->bpts;

  /* Initialization */

  Avert = Ax;
  Aedge = Ax + A_matrix->nvert;
  dzero (n, Ax, 1);

  /* ------------------------------------------------------------ */

  for (elmt = A_mesh->head; elmt; elmt = elmt->next) {
    elmt_get_data (elmt,  x,  xloc);
    Helmholtz     (elmt, A_matrix, xloc, Axloc);
    elmt_sum_data (elmt, Ax, Axloc);
  }

  /* ------------------------------------------------------------ */

  /* Combine vertices and edges */
  
  vdbcombine 
    (A_mesh->vertx, Avert, vdbdsum, sizeof(double),  1,    sizeof(double));
  vdbcombine 
    (A_mesh->edges, Aedge, vdbdsum, sizeof(double), np, np*sizeof(double));

  /* Mask boundary conditions */

  Mask_RHS (A_mesh, A_matrix, Ax);

  return;
}

/* ------------------------------------------------------------------------- *
 * inner_product()                                                           *
 *                                                                           *
 * Compute the distributed inner-product of two vectors.                     *
 * ------------------------------------------------------------------------- */

double inner_product (int n, double *u, double *v)
{
  const int np    = A_mesh->head->nr-2;
  const int nvert = A_matrix->nvert;
  const int nedge = A_matrix->nedge;
  const int *secr = A_matrix->secr;

  double sum = 0.;
  int i, j;

  for (i = 0; i < nvert; i++)              /* ...... vertices ..... */
    if (secr[i]) sum = sum + u[i] * v[i];

  u    += nvert;
  v    += nvert;
  secr += nvert;
  n    -= nvert;

  for (i = 0; i < nedge; i++) {            /* ....... edges ....... */
    if (secr[i]) {
      const int start = np * i;
      const int stop  = np + start;
      for (j = start; j < stop; j++)
	sum = sum + u[j] * v[j];
    }
  }

  u += np*nedge;
  v += np*nedge;
  n -= np*nedge;
  for (i = 0; i < n; i++)                  /* ..... interiors ..... */
    sum = sum + u[i] * v[i];

#ifdef PARALLEL
  { int local = sum; comm_dsum(1, &local, &sum); }
#endif

  return sum;
}

/* ------------------------------------------------------------------------- *
 * edge_xxx_data() -- Data Movement                                          *
 *                                                                           *
 *   (not to be confused with the movement of XXX data!)                     *
 *                                                                           *
 * The following two routines move data from flat-field storage to the local *
 * tensor product storage.  The "get" version just copies global data to the *
 * local buffer, while the "sum" version accumulates data into a shared      *
 * global memory.                                                            *
 *                                                                           *
 * In the nonconforming version, a projection matrix Z sits between the      *
 * global data and the local edge values.  It's stored in sparse format      *
 * with the length of a given row (the number of locally-coupled boundary    *
 * points) stored in a vector n.                                             *
 * ------------------------------------------------------------------------- */

#define INIT(id,nb,ni,local,global,n,bmap,imap,Z,xelmt) \
   const int  id  = elmt->id;                      \
   const int  nb  = elmt->nr * 4 - 4;              \
   const int  ni  = SQR(elmt->nr-2);               \
   int    **index = A_matrix->map[id].index;       \
   int    *npts   = A_matrix->map[id].npts;        \
   double **Z     = A_matrix->map[id].Z;           \
   int    *bmap   = elmt->emap;                    \
   int    *imap   = elmt->emap + nb;               

int elmt_get_data (Element *elmt, double *x, double *xloc)
{
  double xbuf[_MAX_NB], *xb, *xi, *xelmt;
  int i, p;

  INIT(id,nb,ni,local,global,n,bmap,imap,Z,xelmt);

  xb    = x;
  xi    = x  + A_matrix->bpts;
  xelmt = xi + ni * elmt->id;

  /* ..... multiply by Z */

  for (p = 0; p < nb; p++) {
    if (npts[p] == 1)
      xbuf[p] = xb[ index[p][0] ];
    else {
      xbuf[p] = 0.;
      for (i = 0; i < npts[p]; i++)
	xbuf[p] += Z[p][i] * xb[ index[p][i] ];
    }
  }
  
  for (p = 0; p < nb; p++)
    xloc[ bmap[p] ] = xbuf [p];
  for (p = 0; p < ni; p++)
    xloc[ imap[p] ] = xelmt[p];

  return 0;
}


int elmt_sum_data (Element *elmt, double *x, double *xloc)
{
  double   xbuf[_MAX_NB], *xb, *xi, *xelmt;
  int i, p;

  INIT(id,nb,ni,local,global,n,bmap,imap,Z,xelmt);

  xb    = x;
  xi    = x  + A_matrix->bpts;
  xelmt = xi + ni * elmt->id;

  for (p = 0; p < ni; p++)
    xelmt[ p ] += xloc[ imap[p] ];
  for (p = 0; p < nb; p++)
    xbuf [ p ]  = xloc[ bmap[p] ];
  
  /* ..... multiply by Z-transpose */

  for (p = 0; p < nb; p++) {
    if (npts[p] == 1)
      xb[ index[p][0] ] += xbuf[p];
    else {
      for (i = 0; i < npts[p]; i++)
	xb [ index[p][i] ] += Z[p][i] * xbuf[p];
    }
  }

  return 0;
}

#undef INIT

/* ------------------------------------------------------------------------- *
 * Spectral Element Operators                                                *
 *                                                                           *
 * The following routines are among the most important in the library, as    *
 * they define the spectral element operators: the matrices which act on the *
 * elemental data to generate the discrete Helmholtz equations.  There is a  *
 * function for the conjugate gradient solver and a function used to form    *
 * stiffness matrices for the static-condensation solver.  You can use a     *
 * different operator (variable coefficient, polar coordinates, etc.) by     *
 * providing a new operator function to the solvers.                         *
 *                                                                           *
 * Arguments:                                                                *
 *                                                                           *
 * Operator   ( Field   *U,    Template element array                        *
 *              BSystem *B,    Boundary system structure                     *
 *              double  *in,   Input array, default to U if (NULL)           *
 *              double  *out ) Output array, possible the same as "in"       *
 *                                                                           *
 * OperatorSC ( Element *U,    Template element structure (single!)          *
 *              BSystem *B,    Boundary system structure                     *
 *              double  **A,   Boundary-boundary coupling matrix             *
 *              double  **B,   Boundary-interior coupling matrix             *
 *              double  **C )  Interior-interior coupling matrix             *
 * ------------------------------------------------------------------------- */

static Helmholtz (Element *elmt, Matrix *A, double *x, double *Ax)
{
  const int nr   = elmt->nr;
  const int ns   = elmt->ns;
  const int nrns = nr * ns;
  int i, k;

  double lambda = A->constant;
  Matrix_CG  *G = A->CG + elmt->id;
  
  double **dr, **ds, *Pr, *Ps;

  tempVector (tmp, nrns + nrns + nrns);
  
  Pr = tmp + nrns;
  Ps = tmp + nrns + nrns;

  getops (nr, 0, 0, &dr, 0);             /* Get derivative matrices */
  getops (ns, 0, 0, &ds, 0);   

  /* Compute (r,s)-derivatives */

  dgemm ('T','N', nr, ns, nr, 1., *dr, nr, x, nr, 0., Pr, nr); 
  dgemm ('N','N', nr, ns, ns, 1., x, nr, *ds, ns, 0., Ps, nr); 

  if (G->g3) {
    for (i = 0; i < nrns; i++) {
      tmp[i]  = (*G->g1) [i] * Pr[i] + (*G->g3)[i] * Ps[i];
      Ps [i]  = (*G->g2) [i] * Ps[i] + (*G->g3)[i] * Pr[i];
      Ax [i]  = (*elmt->mass)[i] * x[i];
      Pr [i]  = tmp [i];
    }
  } else {
    for (i = 0; i < nrns; i++) {
      Pr [i] *= (*G->g1)[i];
      Ps [i] *= (*G->g2)[i];
      Ax [i]  = (*elmt->mass)[i] * x[i];
    }
  }
  
  /* Compute final derivatives */

  dgemm ('N', 'N', nr, ns, nr, 1., *dr, nr, Pr, nr, lambda, Ax, nr);
  dgemm ('N', 'T', nr, ns, ns, 1., Ps, nr, *ds, ns, 1., Ax, nr);

  freeVector (tmp); 
  return 0;
}


static void show_vector (int n, double *u)
{
  int i;

  for (i = 0; i < n; i++)
    printf ("%4d %g\n", i, u[i]);
}
