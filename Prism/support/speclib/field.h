#ifndef FIELD_H
#define FIELD_H

/* Field
 *
 * $Id$
 *
 * Copyright 1994-1998 R. D. Henderson and Caltech
 *
 * This is the main "workhorse" class.  A Field represents a scalar function
 * of position, u[x,y].  In the case of "3D" Field, there is a collection
 * of M such objects that share the same 2D geometry: u_m[x,y], m <= M.  A
 * Field can be treated as either a continuous function (for differentiation,
 * integration, etc.) or a finite-dimensional vector of discrete values.
 *
 * A Field is allocated by reading its geometry definition from a file and
 * specifying its resolution parameters and type.  Only one Field needs to
 * be allocated for most applications, and additional Fields can be created
 * with the same geometry using Field_dup().
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "speclib/element.h"

typedef Element Field;

/* macros to access elements of a Field */

#define FIELD_NR(u) (u->nr)
#define FIELD_NS(u) (u->ns)
#define FIELD_NZ(u) (u->nz)
#define FIELD_NELMT(u) (Field_count(u))
#define FIELD_TYPE(u) (u->type)

/* Access a Field as a flat array.  The index "i" runs from 0 to n-1, where  *
 * n = Field_npts() is the total number of grid points in the mesh.          *
 *                                                                           *
 * Example (pointwise addition of two Fields):                               *
 *                                                                           *
 * const int n = Field_npts(u);                                              *
 * int i;                                                                    *
 * for (i = 0; i < n; i++)                                                   *
 *    FIELD_FLAT(u,i) += FIELD_FLAT(v,i);                                    *
 *                                                                           *
 * This loop is equivalent to calling Field_axpy(1., v, u).                  */

#define FIELD_FLAT(u,i) ((*u->base)[i])

/* Access a Field as a set of 2D frames.  The index "i" runs from 0 to n-1,  *
 * where n = NR*NS*NEL is the frame dimension; the index "j" runs from 0 to  *
 * m-1, where m = NZ is the number of frames.  The value of "n" can be set   *
 * by calling Field_frameSize(); the value of "m" can be set by calling      *
 * Field_frameCount().                                                       *
 *                                                                           *
 * Example (same as above):                                                  *
 *                                                                           *
 * const int n = Field_frameSize (u);                                        *
 * const int m = Field_frameCount(u);                                        *
 * int i, j;                                                                 *
 * for (j = 0; j < m; j++)                                                   *
 *   for (i = 0; i < n; i++)                                                 *
 *      FIELD_FRAME(u,n,i,j) += FIELD_FRAME(v,n,i,j);                        *
 *                                                                           *
 * Note that the index "i" runs fastest in memory.                           */

#define FIELD_FRAME(u,n,i,j) ((*u->base)[(j)*(n)+(i)])

/* Field allocation: */

Field*   Field_alloc     (char type, int nr, int ns, int nz, int nel);
void     Field_free      (Field *u);
Field*   Field_realloc   (Field *u, int nr, int ns, int nz, char type);
Field*   Field_build     (FILE *fp, int nr, int ns, int nz, char type);

int      Field_count     (const Field *u);
int      Field_npts      (const Field *u);
int      Field_frameSize (const Field *u);
int      Field_frameCount(const Field *u);
double*  Field_base      (const Field *u);
Element* Field_head      (const Field *u);
int      Field_size      (const Field *u);

Field*   Field_aux       (const Field *u, int nr, int ns, int nz, char type);
Field*   Field_dup       (const Field *u);
int      Field_project   (const Field *u, Field *v);
void     Field_set       (Field *u, const char *expr);

int      Field_setFrame      (Field *u, int m);
int      Field_setFrameMulti (int m, int nfields, ...);

/* Derivative operators */

void     Field_grad (const Field *u, Field *dx, Field *dy);
void     Field_dx   (const Field *u, Field *dx);
void     Field_dy   (const Field *u, Field *dy);
void     Field_dz   (const Field *u, Field *dz);
void     Field_dzhat(const Field *u, Field *du);

/* Reduction operations */

double   Field_nrm2 (const Field *u);
double   Field_asum (const Field *u);
double   Field_amax (const Field *u);
double   Field_min  (const Field *u);
double   Field_max  (const Field *u);
double   Field_L2   (const Field *u);
double   Field_H1   (const Field *u);
double   Field_dot  (const Field *u, const Field *v);

/* Integral over the domain */

double   Field_integral (Field *u);

/* Unary operations */

void     Field_scal (double d, Field *u);

/* Binary operations */

void     Field_axpy (double d, const Field *u, Field *v);
void     Field_copy (const Field *u, Field *v);
void     Field_swap (Field *u, Field *v);

/* Triad operations:
 *    
 *   mul         w  = u * v
 *   vvtp        w += u * v
 */

void     Field_mul  (const Field *u, const Field *v, Field *w);
void     Field_vvtp (const Field *u, const Field *v, Field *w);

/* This one is a bit shady.  It computes the Fourier transform of a 3D Field *
 * where the NZ data planes are treated as discrete samples of a periodic    *
 * function, either in time or space.  The parameter "dir" specifies the     *
 * direction of the transform: dir = -1 (forward to Fourier coefficients) or *
 * dir = +1 (reverse to Physical space coefficients).                        */

void Field_FFT (Field *u, int dir);

#endif


