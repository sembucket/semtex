#ifndef MATRIX_H
#define MATRIX_H

/* Matrix
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 *
 * This is the data structure for global matrix systems.  There are really
 * two parts to this: (1) bookeeping for the global numbering system, and 
 * (2) storage of the matrix operator that defines a linear system.
 * ------------------------------------------------------------------------- */

#include "speclib/bc.h"
#include "speclib/edge.h"
#include "speclib/element.h"
#include "speclib/field.h"

typedef struct {             /* ----- Inner Product Matrices --- */
  double  **g1;              /*                                   */
  double  **g2;              /*                                   */
  double  **g3;              /*                                   */
} Matrix_IP;

typedef struct {             /* -- Static Condensation Matrices - */
  int      nrows;            /* Number of coupled boundary points */
  int     *rowmap;           /* Map to those boundary points      */
  double  *A_11;             /* Condensed boundary matrix         */
  double  *A_12;             /* Condensed coupling matrix         */
  double  *A_22;             /* Factored interior matrix          */
} Matrix_SC;

typedef struct bsys {        /* ------- BOUNDARY SYSTEM --------- */
  int      elements    ;     /* Number of elements in the mesh    */
  int      families    ;     /*    "      families                */
  int      bpts, ipts,       /* Number of global nodes            */
           bdof,             /*    "      global boundary dof     */
           ndof        ;     /*    "      global dof              */
  int      bandwidth   ;     /* Bandwidth of the boundary matrix  */
  int**    bmap        ;     /* Global numbering map for B-nodes  */

  /* NB:       bdof  <=   bpts   <=   ndof   <=  bpts + ipts,     //
  //                                                              //
  //      where bpts + ipts gives the total number of independent //
  //      nodes in mesh.   The degrees of freedom are determined  //
  //      by the boundary conditions (essential vs. natural).     */

  /* Mass Matrix */

  double*  mass        ;     /* Global Mass Matrix (diagonal)     */
  double*  massinv     ;     /* Global Mass Matrix (inverted)     */

  /* Stiffnes Matrix */

  int      singular    ;     /* Singularity flag                  */
  double   constant    ;     /* Helmholtz constant                */
  double   condition   ;     /* Matrix condition number           */

  Matrix_IP     *CG    ;     /* Inner Product Factors             */
  Matrix_SC     *SC    ;     /* Elemental Stiffness Matrices      */

  void*    Hp          ;     /* Packed coefficient matrix         */
  void*    other       ;     /* Anything else ?                   */
} BSystem;

typedef void (*Operator)
     (Element *U, BSystem *B, double *u, double *v);
typedef void (*OperatorSC) 
     (Element *U, BSystem *B, double **Abb, double **Abi, double **Aii);

/* Prototypes */

BSystem* Matrix_alloc (Field *U, Bedge *Ubc, char *session);
void     Matrix_free  (BSystem *B);

void Matrix_build     (Field *U, BSystem *B);
void Matrix_build_A   (Field *U, BSystem *B, OperatorSC A);
void Matrix_assemble  (Field *U, BSystem *B, OperatorSC A);

void Matrix_local (const BSystem *B, Field *u, double *local, double *global);
void Matrix_dsum  (const BSystem *B, Field *u, double *local, double *global);
void Matrix_davg  (const BSystem *B, Field *u);

#endif
