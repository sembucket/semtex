#ifndef CUBIT_H
#define CUBIT_H

#ifdef __cpluscplus
extern "C" {
#endif

/*
 * C U B I T
 *
 * Prototype declarations for Cubit.  This is the top-level include file,
 * and the only one that should be necessary in most source code.  All of
 * the other files for data type declarations are included by "cubit.h"
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ----------------------------------------------------------- */

#define  UNSET      0xffffff   /* Flag for an unset parameter              */

/* Other headers */

#include "cubit/cubitypes.h"
#include "cubit/splib.h"
#include "cubit/vertex.h"
#include "cubit/edge.h"
#include "cubit/alias.h"
#include "cubit/family.h"
#include "cubit/element.h"
#include "cubit/curve.h"
#include "cubit/bc.h"
#include "cubit/mesh.h"
#include "cubit/field.h"
#include "cubit/fieldfile.h"
#include "cubit/error.h"
#include "cubit/param.h"
#include "cubit/keyword.h"
#include "cubit/manager.h"
#include "cubit/probe.h"

extern char *cubit_err_msg;
extern int   cubit_err_num;

/* User call-back functions */

extern int (*user_refine)();      /* Refine an element      */
extern int (*user_bc)();          /* Compute BC information */
extern int (*user_prune)();       /* Prune an element       */
extern int (*user_perm)();        /* Permute two elements   */

/* Macros */

#define MAX(a,b) ( (b) < (a) ? (a) : (b) )
#define MIN(a,b) ( (b) > (a) ? (a) : (b) )
#define SQR(a)   ( (a) * (a) )
#define CLAMP(t,a,b)  (MAX(MIN(t,b),a))

/*.........................................................................*/

/*  Always call this first! */

int cubit_init (void);
int cubit_exit (void);

/*
 *  Functions from FAMILIES:
 *  -------------------------  */

Family*  Family_get       (Element *U);
Family*  Family_create    (Element *U);
int      Family_count     (Element *U);
int      Family_set       (Family  *F);
int      Family_add       (Family  *F, Element *U);
void     Family_reset     (void);
void     Family_disable   (void);

/*
 *  Functions from ISOMESH:
 *  ----------------------  */

void     genmap           (Element *U);
void     normals          (Element *U);
double*  zmesh            (int M);

/*
 *  Functions from SOLVE:
 *  --------------------- */

void   Solve        (Field *, Field *, Mesh *, Matrix *);
void   Laplace_IP   (int n, double *x, double *Ax);
double inner_product(int n, double *u, double *v);
double norm         (int n, double *u);

/*
 *  Functions from OPERATORS:
 *  -------------------------  */

void coef   (int N);
void getops (int N, double* *zp, double* *wp, double* **dp, double* **dpt);
void legcoef (int n, int m, const double *u, double *a);

/*
 *  Functions from SESSION:
 *  -----------------------  */

void     ReadParams       (FILE *rea);
Element* ReadMesh         (FILE *rea);
BC*      ReadBCs          (FILE *rea, int group, Element *list);
void     ReadKeys         (FILE *rea, Element *list);
char    *findSection      (const char *name, char *buf, FILE *fp);

/*
 *  Functions from RESIDUAL:
 *  ------------------------  */

void     Form_G           (Field *U, Mesh *mesh, Matrix *A, double *g);
double   Form_RHS         (Field *F, Mesh *mesh, Matrix *M, void (*A)(),
			   double *g, double *r);
void     Mask_RHS         (Mesh *mesh, Matrix *A, double *r);

void     cubit_err        (char *msg);


#ifdef __cpluscplus
}
#endif
#endif      /* END OF CUBIT.H DECLARATIONS */




