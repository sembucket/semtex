#ifndef FIELD_H
#define FIELD_H

#include <stdio.h>
#include "mesh.h"

#define FIELD_TYPE(u)           (u->type)
#define FIELD_HEAD(u)           (u->mesh->head)
#define FIELD_DATA(u,k)         (u->data[k])
#define FIELD_NELMT(u)          (MESH_NELMT(u->mesh))
#define FIELD_MESH(u)           (u->mesh)
#define FIELD_NR(u)             (MESH_NR(u->mesh))
#define FIELD_NS(u)             (MESH_NS(u->mesh))
#define FIELD_NZ(u)             (MESH_NZ(u->mesh))
#define FIELD_NPTS(u)           (FIELD_NR(u)*FIELD_NS(u)*FIELD_NZ(u))

typedef struct Field {    /* ........ FIELD Definition ......... */
  char           type  ;   /* Variable type stored in this Field  */
  int            frame ;   /* Active frame number (3D)            */
  double**       data  ;   /* Data values for each element        */
  Mesh*          mesh  ;   /* Underlying mesh                     */
} Field;

Field*   Field_alloc     (Mesh *mesh);
Field*   Field_realloc   (Field *u);
Field*   Field_dup       (const Field *u);
void     Field_free      (Field *u);

void     Field_set       (Field *u, char *expr);
int      Field_count     (const Field *u);
int      Field_npts      (const Field *u);
int      Field_print     (const Field *u, FILE *fp);


/* unary */

void     Field_negative  (const Field *u, Field *v);

/* binary */

void     Field_add       (const Field *u, const Field *v, Field *w);
void     Field_sub       (const Field *u, const Field *v, Field *w);
void     Field_mult      (const Field *u, const Field *v, Field *w);
void     Field_div       (const Field *u, const Field *v, Field *w);

/* ternary */

/* scalar reduction */

double   Field_L2        (const Field *u);
double   Field_H1        (const Field *u);
double   Field_nrm2      (const Field *u);
double   Field_min       (const Field *u);
double   Field_max       (const Field *u);
double   Field_amax      (const Field *u);
double   Field_integral  (const Field *u);

/* misc */

void     Field_scal      (double alpha, Field *u);
void     Field_shift     (double alpha, Field *u);
void     Field_abs       (const Field *u, Field *v);
void     Field_axpy      (double alpha, const Field *u, Field *v);
void     Field_copy      (const Field *u, Field *v);

void     Field_dx        (const Field *u, Field *du);
void     Field_dy        (const Field *u, Field *du);
void     Field_dz        (const Field *u, Field *du);
void     Field_gradient  (const Field *u, Field *du, int dir);
void     Field_FFT       (const Field *u, Field *v, int dir);

#endif

