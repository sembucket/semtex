/*
 * error estimation and mesh refinement 
 *
 * $Revision$
 *
 * Author:    R. D. Henderson
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "cubit.h"
#include "element.h"
#include "mesh.h"
#include "field.h"
#include "error.h"

#include "veclib/veclib.h"

/* ------------------------------------------------------------------------- */

static void gradient (error_t *error, const Field *u) 
{
  Field   *tmp = Field_alloc(error->mesh);
  Element *elmt;

  double min = FLT_MAX;
  double max = -min;
  double sum = 0.;

  int i;

  Field_dx (u, tmp);
  for (elmt = error->mesh->head; elmt; elmt = elmt->next) {
    const int     id   = ELEMENT_ID(elmt);
    const int     npts = ELEMENT_NR(elmt)*ELEMENT_NS(elmt);
    const double *dx   = FIELD_DATA(tmp,id);
    const double *mass = *elmt->mass;

    error->eps[id].error = 0.;
    for (i = 0; i < npts; i++)
      error->eps[id].error += mass[i] * SQR(dx[i]);
  }

  Field_dy (u, tmp);
  for (elmt = error->mesh->head; elmt; elmt = elmt->next) {
    const int     id   = ELEMENT_ID(elmt);
    const int     npts = ELEMENT_NR(elmt)*ELEMENT_NS(elmt);
    const double *dy   = FIELD_DATA(tmp,id);
    const double *mass = *elmt->mass;
    
    for (i = 0; i < npts; i++)
      error->eps[id].error += mass[i] * SQR(dy[i]);
    error->eps[id].error = sqrt(error->eps[id].error);

    if (error->eps[id].error < min) min = error->eps[id].error;
    if (error->eps[id].error > max) max = error->eps[id].error;

    sum += SQR(error->eps[id].error);
  }

  error->min   = min;
  error->max   = max;
  error->sum   = sqrt(sum);
  error->scale = error->norm.h1;

  Field_free (tmp);
}

/* ------------------------------------------------------------------------- *
 * Compute the Legendre polynomial spectrum                                  *
 * ------------------------------------------------------------------------- */

static void spectrum_0 (error_t *error, const Field *u) 
{
  Element *elmt = MESH_HEAD(error->mesh);

  const int nr = MESH_NR(error->mesh);
  const int ns = MESH_NS(error->mesh);

  double min  = FLT_MAX;
  double max  = -min;
  double sum  = 0.;
  double hmin = dparam("HMIN");
  double **a  = dmatrix (0, ns-1, 0, nr-1);

  while (elmt) {
    const int id = ELEMENT_ID(elmt);
    int i, j;

    double trace = 0.;

    legcoef (nr, ns, FIELD_DATA(u,elmt->id), *a);

    j = nr-1; for (i = 0; i < ns; i++)
      trace += fabs(a[i][j]);
    i = ns-1; for (j = 0; j < nr; j++)
      trace += fabs(a[i][j]);
    error->eps[id].error = 
      trace - fabs(a[ns-1][nr-1]);  /* this was added twice */

    error->eps[id].error = trace * error->eps[id].size;

    if (error->eps[id].error < min) min = error->eps[id].error;
    if (error->eps[id].error > max) max = error->eps[id].error;

    sum += SQR(error->eps[id].error);
    elmt = elmt->next;
  }

  error->min   = min;
  error->max   = max;
  error->sum   = sqrt(sum);
  error->scale = error->norm.l2;
  
  free_dmatrix (a, 0, 0);
}

/* Same as above but we also multiply by the Jacobian */

static void spectrum_1 (error_t *error, const Field *u) 
{
  Element *elmt  = MESH_HEAD(error->mesh);
  const int nr   = MESH_NR(error->mesh);
  const int ns   = MESH_NS(error->mesh);
  const int npts = nr*ns;

  double min  = FLT_MAX;
  double max  = -min;
  double sum  = 0.;
  double hmin = dparam("HMIN");
  double **a  = dmatrix (0, ns-1, 0, nr-1);

  while (elmt) {
    const int id = ELEMENT_ID(elmt);
    int i, j;

    double trace = 0.;
    double buf[_MAX_NORDER*_MAX_NORDER];

    memcpy (buf, FIELD_DATA(u,elmt->id), npts*sizeof(double));

    for (i = 0; i < npts; i++)
      buf[i] *= (*elmt->jac)[i];

    legcoef (nr, ns, buf, *a);

    j = nr-1; for (i = 0; i < ns; i++)
      trace += fabs(a[i][j]);
    i = ns-1; for (j = 0; j < nr; j++)
      trace += fabs(a[i][j]);
    error->eps[id].error = 
      trace - fabs(a[ns-1][nr-1]);  /* this was added twice */

    error->eps[id].error = trace;

    if (error->eps[id].error < min) min = error->eps[id].error;
    if (error->eps[id].error > max) max = error->eps[id].error;

    sum += SQR(error->eps[id].error);
    elmt = elmt->next;
  }

  error->min   = min;
  error->max   = max;
  error->sum   = sqrt(sum);
  error->scale = error->norm.l2;
  
  free_dmatrix (a, 0, 0);
}

static void spectrum (error_t *error, const Field *u) {
  spectrum_1 (error,u);
}

/* ------------------------------------------------------------------------- *
 * Regression                                                                *
 * ------------------------------------------------------------------------- */

/* Fit the spectrum, const * exp(-sigma n), to values a[n1] to a[n2] */

static void compute_regress 
(int n1, int n2, double *a, double *trunc, double *sigma)
{
  double S, Sx, Sxx, Sxy, Syy, Sy, del, A, B;
  int i;

  S  = (double) n2 - n1 + 1;
  Sx = Sxx = Sxy = Syy = Sy = 0.;

  for (i = n1; i <= n2; i++) {
    double xval = i;
    double yval = log(a[i]);

    Sx  += xval;
    Sxx += xval * xval;
    Sxy += xval * yval;
    Syy += yval * yval;
    Sy  += yval;
  }
    
  del = ( S   * Sxx - Sx * Sx );
  A   = ( Sxx * Sy  - Sx * Sxy ) / del;
  B   = ( Sxy * S   - Sx * Sy  ) / del;
  
  *trunc = exp(A);
  *sigma = -B;
}

/* Extrapolate the error using the fit */

static double extrapolate (int n, double sigma)
{
  double term  = exp(-2.*sigma*n);
  double error = term;

  while (term > DBL_EPSILON) {
    n++;
    term   = exp(-2.*sigma*n);
    error += term;
  }

  return error;
}

static void regress (error_t *error, const Field *u)
{
  Element *elmt  = MESH_HEAD(error->mesh);
  const int np   = MESH_NR(error->mesh);
  const int npts = np*np;

  double min  = FLT_MAX;
  double max  = -min;
  double sum  = 0.;
  double **a  = dmatrix (0, np-1, 0, np-1);

  while (elmt) {
    const int id = ELEMENT_ID(elmt);
    int i, p;

    double trunc;
    double sigma;
    double buf[_MAX_NORDER*_MAX_NORDER];
    double b  [_MAX_NORDER];

    memcpy (buf, FIELD_DATA(u,elmt->id), npts*sizeof(double));

    for (i = 0; i < npts; i++)
      buf[i] *= (*elmt->jac)[i];

    legcoef (np, np, buf, *a);

    /* Collapse to a 1D spectrum */

    for (p = 0; p < np; p++) {
      b[p] = fabs(a[p][p]);
      for (i = 0; i < p-1; i++)
	b[p] += fabs(a[i][p]) + fabs(a[p][i]);
    }

    /* Compute the decay rate and truncation level */

    compute_regress (MAX(np-5,0), MAX(np-1,0), b, &trunc, &sigma);
    
    /* Extrapolate the results */

    if (sigma > 0.)
      error->eps[id].error = 
	sqrt(SQR(b[p-1]) + SQR(trunc) * extrapolate(np,sigma));
    else
      error->eps[id].error = 1.;

    if (error->eps[id].error < min) min = error->eps[id].error;
    if (error->eps[id].error > max) max = error->eps[id].error;

    sum += SQR(error->eps[id].error);
    elmt = elmt->next;
  }

  error->min   = min;
  error->max   = max;
  error->sum   = sqrt(sum);
  error->scale = error->norm.l2;
  
  free_dmatrix (a, 0, 0);
}


/* ------------------------------------------------------------------------- */

static void compute_norms (error_t *error, const Field *u)
{
  Field *tmp = Field_alloc(u->mesh);
  Element *elmt;
  int i;
  
  /* Initialize */

  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    const int k         = ELEMENT_ID   (elmt);
    error->eps[k].depth = ELEMENT_DEPTH(elmt);
    error->eps[k].size  = ELEMENT_SIZE (elmt);
    error->eps[k].error = 0.;
  }

  /* Compute LOCAL norms */

  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    const int npts = ELEMENT_NR(elmt)*ELEMENT_NS(elmt);
    const int k    = ELEMENT_ID(elmt);
    double   *val  = ELEMENT_DATA(elmt,u);
    double   *mass = ELEMENT_MASS(elmt);

    double inf = -FLT_MAX;
    double sum =  0.;

    for (i = 0; i < npts; i++) {
      inf  = MAX(inf, fabs(val[i]));
      sum += mass[i] * SQR(val[i]);
    }

    error->eps[k].norm.inf = inf;
    error->eps[k].norm.l2  = sqrt(sum);
    error->eps[k].norm.h1  = sum;
  }

  Field_dx (u, tmp);
  for (elmt = FIELD_HEAD(tmp); elmt; elmt = elmt->next) {
    const int npts = ELEMENT_NR(elmt)*ELEMENT_NS(elmt);
    const int k    = ELEMENT_ID(elmt);
    double   *val  = ELEMENT_DATA(elmt,tmp);
    double   *mass = ELEMENT_MASS(elmt);

    double sum = 0.;

    for (i = 0; i < npts; i++)
      sum += mass[i] * SQR(val[i]);

    error->eps[k].norm.h1 += sum;
  }

  Field_dy (u, tmp);
  for (elmt = FIELD_HEAD(tmp); elmt; elmt = elmt->next) {
    const int npts = ELEMENT_NR(elmt)*ELEMENT_NS(elmt);
    const int k    = ELEMENT_ID  (elmt);
    double   *val  = ELEMENT_DATA(elmt,tmp);
    double   *mass = ELEMENT_MASS(elmt);

    double sum = 0.;

    for (i = 0; i < npts; i++)
      sum += mass[i] * SQR(val[i]);

    error->eps[k].norm.h1 += sum;
    error->eps[k].norm.h1  = sqrt(error->eps[k].norm.h1);
  }
  
  /* Compute GLOBAL norms */

  error->norm.inf = 0.;
  error->norm.l2  = 0.;
  error->norm.h1  = 0.;

  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    const int k = ELEMENT_ID(elmt);

    error->norm.inf = MAX(error->norm.inf, error->eps[k].norm.inf);
    error->norm.l2 += SQR(error->eps[k].norm.l2);
    error->norm.h1 += SQR(error->eps[k].norm.h1);
  }

  error->norm.l2 = sqrt(error->norm.l2);
  error->norm.h1 = sqrt(error->norm.h1);

  dparam_set("norm:inf",   error->norm.inf);
  dparam_set("norm:ltwo",  error->norm.l2);
  dparam_set("norm:hone",  error->norm.h1);

  Field_free(tmp);
}

  
/* ------------------------------------------------------------------------- */

error_t *Error_alloc (Mesh *mesh, estimate_t type)
{
  error_t *error = (error_t*) calloc(1,sizeof(error_t));

  error->type = type;
  error->mesh = mesh;

  return error;
}

void Error_free (error_t *error) {
  free (error);
}

/* ------------------------------------------------------------------------- */

double Error_compute (error_t *error, const Field *u)
{
  /* Compute local and global solution norms */

  compute_norms(error,u);

  switch (error->type) {
  case E_SPECTRUM:
    spectrum (error, u);
    break;
  case E_GRADIENT:
    gradient (error, u);
    break;
  case E_REGRESS:
    regress  (error, u);
    break;
  default:
    break;
  }

  return error->sum;
}

/* ------------------------------------------------------------------------- */

int Error_adapt (error_t *error, double tol, int maxdepth)
{
  const double threshold = tol * error->scale;
  int nref = 0;
  Element *elmt;

  Mesh_resetFlags(error->mesh);

  for (elmt = MESH_HEAD(error->mesh); elmt; elmt = elmt->next) {
    const int k = ELEMENT_ID(elmt);

    if (error->eps[k].error > threshold && maxdepth > error->eps[k].depth) {
      ELEMENT_REFINE(elmt) = 1;
      nref++;
    }
  }

  Mesh_refine (error->mesh);
  return nref;
}

/* ------------------------------------------------------------------------- */

void Error_info (const error_t *error, double tol, FILE *fp)
{
  const double threshold = tol * error->scale;
  Element *elmt;

  fputs ("error estimate [", fp);
  switch (error->type) {
  case E_REGRESS:
    fputs ("regression", fp);
    break;
  case E_SPECTRUM:
    fputs ("spectrum", fp);
    break;
  case E_GRADIENT:
    fputs ("gradient", fp);
    break;
  default:
    fputs ("unknown", fp);
    break;
  }
  fputs ("]:\n", fp);

  fprintf (fp, "\tnorm  = [%g %g %g]\n", 
	   error->norm.inf, error->norm.l2, error->norm.h1);
  fprintf (fp, "\tscale = %g\n", error->scale);
  fprintf (fp, "\tmin   = %g\n", error->min / error->scale);
  fprintf (fp, "\tmax   = %g\n", error->max / error->scale);
  fprintf (fp, "\tsum   = %g\n", error->sum / error->scale);

  for (elmt = MESH_HEAD(error->mesh); elmt; elmt = elmt->next) {
    const int k = ELEMENT_ID(elmt);

    if (error->eps[k].error > threshold) {
      fprintf (fp, "[%3d] norm = [%g %g %g] err = %g\n", k, 
	       error->eps[k].norm.inf, 
	       error->eps[k].norm.l2,
	       error->eps[k].norm.h1,
	       error->eps[k].error);
    }
  }

  /* Store values as parameters */

  dparam_set ("eps:min", error->min / error->scale);
  dparam_set ("eps:max", error->max / error->scale);
  dparam_set ("eps:sum", error->sum / error->scale);
}

