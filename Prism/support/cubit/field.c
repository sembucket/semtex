/*
 * Functions that operate on a Field
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "cubit.h"
#include "field.h"

#include "veclib/veclib.h"

/* -------------------------------------------------------------------- * 
 * Field_alloc(), Field_realloc()                                       *
 *                                                                      *
 * Allocate memory for a Field                                          *
 * -------------------------------------------------------------------- */

Field *Field_alloc (Mesh *mesh)
{
  const int nr = MESH_NR(mesh);
  const int ns = MESH_NR(mesh);
  const int nz = MESH_NZ(mesh);
  Element *elmt;

  Field *f = (Field*) calloc (1, sizeof(Field));
  assert(f);

  f->data  = (double**) calloc(_MAX_NEL, sizeof(Element*));
  f->mesh  = mesh;

  for (elmt = MESH_HEAD(mesh); elmt; elmt = elmt->next)
    f->data[ELEMENT_ID(elmt)] = (double*) calloc (nr*ns*nz, sizeof(double));

  return f;
}

Field *Field_realloc (Field *f)
{
  Mesh* mesh = f->mesh;
  Element *elmt;

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    const int k = ELEMENT_ID(elmt);
    if (f->data[k] == NULL)
        f->data[k] = (double*) calloc(ELEMENT_NPTS(elmt), sizeof(double));
  }

  return f;
}

Field *Field_dup (const Field *u) {
  Field *v = Field_alloc(FIELD_MESH(u));
  Field_copy(u,v);
  return v;
}
  
void Field_free (Field *f) {
  int  i;
  for (i = 0; i < _MAX_NEL; i++)
    if (f->data[i])
      free (f->data[i]);
  free (f);
}


/* -------------------------------------------------------------------- * 
 * Field_count() -- Count Active Elements                               *
 * -------------------------------------------------------------------- */

int Field_count (const Field *u)
{
  Element *elmt = FIELD_HEAD(u);
  int n = 0;

  while (elmt) { n++; elmt = elmt->next; }

  return n;
}

int Field_npts (const Field *u) {
  return FIELD_NPTS(u);
}

/* -------------------------------------------------------------------- * 
 * Field_set()                                                          *
 *                                                                      *
 * Use the parser to define a Field's value                             *
 * -------------------------------------------------------------------- */

void Field_set (Field *u, char *expr)
{
  Mesh     *mesh = u->mesh;
  Element  *elmt = MESH_HEAD(mesh);
  const int nrns = MESH_NR(mesh) * MESH_NS(mesh);
  const int nz   = MESH_NZ(mesh);

  vector_def("x y", expr);

  if (nz == 1) {
    while (elmt) {
      vector_set(nrns, *elmt->xmesh, *elmt->ymesh, u->data[elmt->id]);
      elmt = elmt->next;
    }
  } else {
    const double dz = dparam("LZ")/nz;
    double z;
    int k;

    while (elmt) {
      double *uptr = u->data[elmt->id];
      for (k = 0, z = 0.; k < nz; k++, z += dz) {
	scalar_set("z", z);
	vector_set (nrns, *elmt->xmesh, *elmt->ymesh, uptr + k*nrns);
      }
      elmt = elmt->next;
    }
  }
}

/* Compute the integral L2-norm of a Field */

double Field_L2 (const Field *u)
{
  const int nrns = MESH_NR(u->mesh)*MESH_NS(u->mesh);
  const int nz   = MESH_NZ(u->mesh);
  double sum = 0.0;

  int i, k;
  Element *elmt;

  if (nz==1) {
    for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
      const double *mass = *elmt->mass;
      const double *uptr = FIELD_DATA(u,elmt->id);
      for (i = 0; i < nrns; i++, uptr++)
	sum += mass[i] * (*uptr) * (*uptr);
    }
  } else {
    for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
      const double *mass = *elmt->mass;
      const double *uptr = FIELD_DATA(u,elmt->id);
      for (i = 0; i < nrns; i++)
	for (k = 0; k < nz; k++, uptr++)
	  sum += mass[i] * (*uptr) * (*uptr);
    }
    sum *= dparam("LZ")/nz;
  }

  return sqrt(sum);
}

/* Compute the integral H1-norm of a Field */

double Field_H1 (const Field *u)
{
  const int nrns = MESH_NR(u->mesh)*MESH_NS(u->mesh);
  const int nz   = MESH_NZ(u->mesh);
  double sum = 0.0;
  double val;

  Field *tmp = Field_alloc(u->mesh);

  Element *elmt;
  int i;

  val = Field_L2(u);
  sum = val*val;

  Field_dx (u, tmp);
  val  = Field_L2(tmp);
  sum += val*val;
  
  Field_dy (u, tmp);
  val  = Field_L2(tmp);
  sum += val*val;

  if (nz > 1) {
    Field_dz(u, tmp);
    val  = Field_L2(tmp);
    sum += val*val;
  }

  
  Field_free (tmp);
  return sqrt(sum);
}


double Field_nrm2 (const Field *u)
{
  double sum = 0.0;
  Element *elmt = FIELD_HEAD(u);

  while (elmt) {
    double s = dnrm2(ELEMENT_NPTS(elmt), u->data[elmt->id], 1);
    sum += s*s;
    elmt = elmt->next;
  }

  return sqrt(sum);
}

void Field_scal (double alpha, Field *u)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);

  while (elmt) {
    dscal (npts, alpha, FIELD_DATA(u,elmt->id), 1);
    elmt = elmt->next;
  }
}

void Field_shift (double alpha, Field *u)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);
  int i;

  while (elmt) {
    double *uptr = FIELD_DATA(u,elmt->id);

    for (i = 0; i < npts; i++)
      uptr[i] += alpha;

    elmt = elmt->next;
  }
}

void Field_abs (const Field *u, Field *v)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);
  int i;

  while (elmt) {
    double *uptr = FIELD_DATA(u,elmt->id);
    double *vptr = FIELD_DATA(v,elmt->id);
    for (i = 0; i < npts; i++)
      vptr[i] = fabs(uptr[i]);
    elmt = elmt->next;
  }
}

void Field_axpy (double alpha, const Field *u, Field *v)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);

  while (elmt) {
    const int id = elmt->id;
    daxpy(npts, alpha, FIELD_DATA(u,id), 1, FIELD_DATA(v,id), 1);
    elmt = elmt->next;
  }
}

void Field_copy (const Field *u, Field *v)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);

  while (elmt) {
    const int id = elmt->id;
    memcpy (FIELD_DATA(v,id), FIELD_DATA(u,id), npts*sizeof(double));
    elmt = elmt->next;
  }
}

/* binary operations */

void Field_add (const Field *u, const Field *v, Field *w)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);
  int i;

  while (elmt) {
    const int id = elmt->id;
    double *uptr = FIELD_DATA(u,id);
    double *vptr = FIELD_DATA(v,id);
    double *wptr = FIELD_DATA(w,id);

    for (i = 0; i < npts; i++)
      wptr[i] = uptr[i] + vptr[i];

    elmt = elmt->next;
  }
}

void Field_sub (const Field *u, const Field *v, Field *w)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);
  int i;

  while (elmt) {
    const int id = elmt->id;
    double *uptr = FIELD_DATA(u,id);
    double *vptr = FIELD_DATA(v,id);
    double *wptr = FIELD_DATA(w,id);

    for (i = 0; i < npts; i++)
      wptr[i] = uptr[i] - vptr[i];

    elmt = elmt->next;
  }
}

void Field_mult (const Field *u, const Field *v, Field *w)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);
  int i;

  while (elmt) {
    const int id = elmt->id;
    double *uptr = FIELD_DATA(u,id);
    double *vptr = FIELD_DATA(v,id);
    double *wptr = FIELD_DATA(w,id);

    for (i = 0; i < npts; i++)
      wptr[i] = uptr[i] * vptr[i];

    elmt = elmt->next;
  }
}

void Field_div (const Field *u, const Field *v, Field *w)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);
  int i;

  while (elmt) {
    const int id = elmt->id;
    double *uptr = FIELD_DATA(u,id);
    double *vptr = FIELD_DATA(v,id);
    double *wptr = FIELD_DATA(w,id);

    for (i = 0; i < npts; i++)
      wptr[i] = uptr[i] / vptr[i];

    elmt = elmt->next;
  }
}

void Field_negative (const Field *u, Field *v)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);
  int i;

  while (elmt) {
    const int id = elmt->id;
    double *uptr = FIELD_DATA(u,id);
    double *vptr = FIELD_DATA(v,id);

    for (i = 0; i < npts; i++)
      vptr[i] = -uptr[i];

    elmt = elmt->next;
  }
}

double Field_min (const Field *u)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);

  double min = FLT_MAX;
  int i;

  while (elmt) {
    const double *val = FIELD_DATA(u,elmt->id);
    for (i = 0; i < npts; i++)
      min = min < val[i] ? min : val[i];
    elmt = elmt->next;
  }

  return min;
}

double Field_max (const Field *u)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);

  double max = -FLT_MAX;
  int i;

  while (elmt) {
    const double *val = FIELD_DATA(u,elmt->id);
    for (i = 0; i < npts; i++)
      max = max > val[i] ? max : val[i];
    elmt = elmt->next;
  }

  return max;
}

double Field_amax (const Field *u)
{
  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u);
  Element  *elmt = FIELD_HEAD(u);

  double amax = 0.;
  int i;

  while (elmt) {
    const double *val = FIELD_DATA(u,elmt->id);
    for (i = 0; i < npts; i++)
      amax = amax > fabs(val[i]) ? amax : fabs(val[i]);
    elmt = elmt->next;
  }

  return amax;
}

/* Integral of a field over the domain */

double Field_integral (const Field *u)
{
  double sum = 0.;

  const int nrns = FIELD_NR(u)*FIELD_NS(u);
  const int nz   = FIELD_NZ(u);

  Element *elmt;
  int i, k;

  if (nz==1) {
    for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
      double *uptr = FIELD_DATA(u,elmt->id);
      double *mass = *elmt->mass;
      for (i = 0; i < nrns; i++, uptr++)
	sum += mass[i] * (*uptr);
    }
  } else {
    for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
      double *uptr = FIELD_DATA(u,elmt->id);
      double *mass = *elmt->mass;
      for (i = 0; i < nrns; i++)
	for (k = 0; k < nz; k++, uptr++)
	  sum += mass[i] * (*uptr);
    } 
    sum *= dparam("LZ")/nz;
  }
  
  return sum;
}

/* ------------------------------------------------------------------------ *
 * dx, dy, dr, dz -- Compute spatial derivatives                            *
 * ------------------------------------------------------------------------ */

void Field_dx (const Field *u, Field *d)
{
  Element *elmt  = FIELD_HEAD(u);
  const int nrns = FIELD_NR(u)*FIELD_NS(u);
  const int nz   = FIELD_NZ(u);

  if (nz==1) {
    while (elmt) {
      Element_dx (elmt, FIELD_DATA(u,elmt->id), FIELD_DATA(d,elmt->id));
      elmt = elmt->next;
    }
  } else {
    while (elmt) {
      double *uptr = FIELD_DATA(u,elmt->id);
      double *dptr = FIELD_DATA(d,elmt->id);
      int k;
      for (k = 0; k < nz; k++)
	Element_dx (elmt, uptr + k*nrns, dptr + k*nrns);
      elmt = elmt->next;
    }
  }
}

void Field_dy (const Field *u, Field *d)
{
  Element *elmt  = FIELD_HEAD(u);
  const int nrns = FIELD_NR(u)*FIELD_NS(u);
  const int nz   = FIELD_NZ(u);

  if (nz==1) {
    while (elmt) {
      Element_dy (elmt, FIELD_DATA(u,elmt->id), FIELD_DATA(d,elmt->id));
      elmt = elmt->next;
    }
  } else {
    while (elmt) {
      double *uptr = FIELD_DATA(u,elmt->id);
      double *dptr = FIELD_DATA(d,elmt->id);
      int k;
      for (k = 0; k < nz; k++)
	Element_dy (elmt, uptr + k*nrns, dptr + k*nrns);
      elmt = elmt->next;
    }
  }
}

void Field_dz (const Field *u, Field *du)
{
  const int m = FIELD_NR(u)*FIELD_NS(u);
  const int n = FIELD_NZ(u);

  const double beta = dparam("BETA");
  
  Element *elmt;

  Field *tmp = Field_dup(u);
  Field_FFT(u,tmp,-1);

  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    double *uptr = FIELD_DATA(tmp,elmt->id);
    int i;

    /* First two planes are set to zero */

    dzero (2*m, uptr, 1);

    /* For remaining planes we have:
     *
     *    du/dz = I beta_k (ur + I ui) = beta_k (-ui + I ur)
     */

    for (i = 2; i < n; i += 2) {
      double scal = beta * (i >> 1);
      double *ur  = uptr + m*i;
      double *ui  = ur   + m;

      dswap (m, ur, 1, ui, 1);
      dscal (m, -scal, ur, 1);
      dscal (m,  scal, ui, 1);
    }
  }

  Field_FFT (tmp,du,-1);
  Field_free(tmp);
}

void Field_gradient (const Field *u, Field *du, int dir)
{
  switch (dir) {
  case 0:
    Field_dx(u,du);
    break;
  case 1:
    Field_dy(u,du);
    break;
  case 2:
    Field_dz(u,du);
    break;
  default:
    break;
  }
}

/* ------------------------------------------------------------------------ *
 * Field_avg() -- Weighted Average                                          *
 *                                                                          *
 * Use Field_avg() to solve the weighted-residual equation                  *
 *                                                                          *
 *                           (u,w) = (g,w)                                  *
 *                                                                          *
 * where "g" is input and "u" is output.  The end result is to compute a    *
 * weighted average of nodal values using the global mass matrix, giving    *
 * a solution "u" that is continuous across element boundaries.             *
 * ------------------------------------------------------------------------ */

void Field_avg (Field *u, Matrix *A)
{
  Element *elmt  = FIELD_HEAD(u);
  const int bpts = A->bpts;
  const int nb   = (elmt->nr + elmt->ns - 2) << 1;
  int p, i;

  tempVector (tmp, bpts);
  dzero(bpts, tmp, 1);

  while (elmt) {
    const int k = ELEMENT_ID(elmt);
    int *emap   = elmt->emap;
    int *npts   = A->map[k].npts;
    int **index = A->map[k].index;
    double **Z  = A->map[k].Z;

    for (p = 0; p < nb; p++)
      if (npts[p] == 1)
	tmp[index[p][0]] += (*elmt->mass)[emap[p]] * FIELD_DATA(u,k)[emap[p]];
      else
	for (i = 0; i < npts[p]; i++)
	  tmp[index[p][i]] += Z[p][i] * 
	    (*elmt->mass)[emap[p]] * FIELD_DATA(u,k)[emap[p]];

    elmt = elmt->next;
  }

  dvmul (bpts, A->massinv, 1, tmp, 1, tmp, 1);

  for (elmt = FIELD_HEAD(u); elmt ; elmt = elmt->next) {
    const int k = elmt->id;
    int *emap   = elmt->emap;
    int *npts   = A->map[k].npts;
    int **index = A->map[k].index;
    double **Z  = A->map[k].Z;

    for (p = 0; p < nb; p++)
      if (npts[p] == 1)
	FIELD_DATA(u,k)[emap[p]] = tmp[index[p][0]];
      else {
	FIELD_DATA(u,k)[emap[p]] = 0.;
	for (i = 0; i < npts[p]; i++)
	  FIELD_DATA(u,k)[emap[p]] += Z[p][i] * tmp[index[p][i]];
      }
  }

  freeVector (tmp);
}

/* ------------------------------------------------------------------------- *
 * Fourier transform of a field:  dir = -1 (Fourier), dir = +1 (Physical)    *
 * ------------------------------------------------------------------------- */

void Field_FFT (const Field *u, Field *v, int dir)
{
  const int m = FIELD_NR(u)*FIELD_NS(u);
  const int n = FIELD_NZ(u);

  Element *elmt;

  int i, j, ierr;

  const int irev  = 1;
  const int ireal = 1;

  static int     ncmplx = 0;
  static int    *bitrev = NULL;
  static double *factor = NULL;
  
  double **tmp = dmatrix(0, m-1, 0, n+1);

  /* Initialize the FFT */

  if (ncmplx != n/2) {
    if (factor) free (factor);
    if (bitrev) free (bitrev);

    ncmplx = n/2;
    factor = (double*) malloc(6*ncmplx*sizeof(double));
    bitrev = (int*) malloc(6*ncmplx*sizeof(int));
    fftdf(*tmp, ncmplx, 1, 1, 1, 0, factor, irev, bitrev, ierr, ireal);
  }

  for (i = 0; i < m; i++)
    for (j = n; j < n+2; j++)
      tmp[i][j] = 0.;

  /* Begin transform */

  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    double *uptr = FIELD_DATA(u,elmt->id);
    double *vptr = FIELD_DATA(v,elmt->id);

    /* Transpose the data */

    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	tmp[j][i] = uptr[i*m+j];
    
    fftdf(*tmp,ncmplx,1,m,ncmplx+1,dir,factor,irev,bitrev,ierr,ireal);

    if (dir == -1)
      dscal (m*(n+2), 1./(4.*ncmplx), *tmp, 1);
    
    /* Transpose back */
    
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	vptr[i*m+j] = tmp[j][i];
  }
    
  free_dmatrix(tmp, 0, 0);
}

int Field_show (const Field *u) {
  return Field_print(u,stdout);
}

int Field_print (const Field *u, FILE *fp)
{
  const int nr  = FIELD_NR(u);
  const int ns  = FIELD_NS(u);
  const int nz  = FIELD_NZ(u);
  const int nel = FIELD_NELMT(u);

  Element *elmt;
  int i, j, k;

  fprintf (fp, "Field %c:\n", FIELD_TYPE(u));

  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    const int     id   = ELEMENT_ID(elmt);
    const double *uptr = FIELD_DATA(u,id);

    fprintf (fp, "Element %d:\n", id);
    for (k = 0; k < nz; k++) {
      for (i = 0; i < ns; i++) {
	for (j = 0; j < nr; j++) {
	  fprintf (fp, " %#10.6g", uptr[j + ns*(i + nr*k)]);
	}
	fprintf (fp, "\n");
      }
      fprintf (fp, "--\n");
    }
  }

  return 0;
}
