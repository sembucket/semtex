#ifndef SPECLIB_H
#define SPECLIB_H

/* S P E C L I B 
 *
 * $Id$
 *
 * Copyright (c) 1994-1998 R. D. Henderson and Caltech
 *
 * Prototype declarations for speclib.  This is the top-level include file,
 * and the only one that should be necessary in most source code.  All of
 * the other files for data type declarations are included by "speclib.h".
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdarg.h>

/* The following are limits that may or may not be compiled into the code. */
/* Most are just "suggestions" for decent operating parameters.            */

#define _MAX_NEL        1024   /* Elements in a single mesh                */
#define _MAX_NORDER       32   /* Maximum Polynomial order in space        */
#define _MAX_NB          128   /* Points along the edge of an element      */
#define _MAX_NC         1024   /* Points describing a curved side          */
#define _MAX_FIELDS       24   /* Maximum fields in a solution file        */
#define  UNSET      0xffffff   /* Flag for an unset parameter              */

typedef enum {                 /* --------- List of flags --------- */
  OFF            = 0,          /* General-purpose flags for on/off  */
  ON             = 1,          /*                                   */
  ERROR          = 2,          /* Non-recoverable error             */
  WARNING        = 3           /* Not too serious                   */
} FLAG;

/* Other headers */

#include "speclib/config.h"
#include "speclib/bc.h"
#include "speclib/edge.h"
#include "speclib/element.h"
#include "speclib/family.h"
#include "speclib/field.h"
#include "speclib/fieldfile.h"
#include "speclib/frame.h"
#include "speclib/operators.h"
#include "speclib/matrix.h"
#include "speclib/probe.h"
#include "speclib/splib.h"

/* Macros */

#define MAX(a,b) ( (b) < (a) ? (a) : (b) )
#define MIN(a,b) ( (b) > (a) ? (a) : (b) )
#define CLAMP(t,a,b)  (MAX(MIN(t,b),a))

/*.........................................................................*/

/* The initialization routine MUST be called before any library functions! */

int speclib_init(void);
int speclib_exit(void);

int speclib_warning (const char *fmt, ...);
int speclib_error   (const char *fmt, ...);
int speclib_options (FILE *fp);
int speclib_params  (FILE *fp);

/* Functions for dealing with isoparametric mesh generation: */
void     genmap           (Element *elmt);
void     normals          (Element *elmt);
double*  zmesh            (int m);

void     Field_local     (Field *U, BSystem *B, double *local, double *global);
void     Field_dsum      (Field *U, BSystem *B, double *local, double *global);
void     Field_davg      (Field *U, BSystem *B);

/* symbol table manager and parser: */
void     manager_init     (void);
int      option           (char *name);
int      option_set       (char *name, int value);
int      iparam           (char *name);
int      iparam_set       (char *name, int value);
double   dparam           (char *name);
double   dparam_set       (char *name, double value);
double   scalar           (char *expr);
double   scalar_set       (char *name, double value);
void     vector_def       (const char *vlist, const char *expr);
void     vector_set       (int   vsize, ... /* v1, v2, ..., f(v) */ );
void     show_symbols     (FILE *fp);
void     show_params      (FILE *fp);
void     show_options     (FILE *fp);

void param_set       (char *name, char *expr);
void param_exception (char *name, int type);

/* Functions for processing the input file: */
void     ReadParams       (FILE *rea);
void     ReadScalars      (FILE *rea);
void     ReadLogics       (FILE *rea);
void     ReadVPs          (FILE *rea);
Bedge*   ReadBCs          (FILE *rea, int group, Field *U);
Field*   ReadMesh         (FILE *rea);

char    *findSection      (char *name, char *buf, FILE *fp);


/* Utilities: */
#define  Element_get(U,k)  ((U) + k)
#define  Element_edge(U,l) ((U)->edges + l)
int      idcomp (const void *p1, const void *p2);
int      ecopy  (int n, double *x, int incx, double *y, int incy);

#ifdef DEBUG
void     show_field       (Field *U);
void     show_matrix      (double **m, int ra, int ca, int rb, int cb);
void     show_vector      (double *v, int imin, int imax);
void     show_bcs         (Bedge *Ubc);
void     show_normals     (Bedge *Ubc);
#endif

/* Elliptic solvers (Laplace, Poisson, and Helmholtz): */
void     Solve            (Field *U, Field *F, Bedge *Ubc, BSystem *B);
void     Solve_CG         (Field *U, Field *F, Bedge *Ubc, BSystem *B);
void     Solve_A          (Field *U, Field *F, Bedge *Ubc, BSystem *B, 
			                       Operator A);
void     Solve_CG_A       (Field *U, Field *F, Bedge *Ubc, BSystem *B, 
			                       Operator A, double  *M);
double   Form_RHS         (Field *U, Field *F, Bedge *Ubc, BSystem *B, 
                                               Operator A, double *r);
void     Mask_RHS         (Field *U, BSystem *B, double *r);

void     BC_set           (Bedge *list, Field *u, BSystem *B);


#endif      /* END OF SPECLIB.H DECLARATIONS */
