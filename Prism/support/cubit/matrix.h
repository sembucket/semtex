#ifndef MATRIX_H
#define MATRIX_H

#include "mesh.h"

/* Macros */

#define MATRIX_LAMBDA(m)   (m->constant)
#define MATRIX_KAPPA(m)    (m->condition)
#define MATRIX_SINGULAR(m) (m->singular)

typedef struct {             /* ------ Conjugate Gradients ------ */
  int     family, key;       /* element ID                        */
  double  **g1;              /* (dx/ds + dy/ds) factors           */
  double  **g2;              /* (dx/dr + dy/dr) factors           */
  double  **g3;              /* (dr/ds + ds/dr) factors           */
} Matrix_CG;

typedef struct {             /* ------ Static Condensation ------ */
  int     family, key;       /* element ID                        */
  double  *A_11;             /* Condensed boundary matrix         */
  double  *A_12;             /* Condensed coupling matrix         */
  double  *A_22;             /* Factored interior matrix          */
} Matrix_SC;

typedef struct Matrix {      /* -------- Global Matrix ---------- */
  int      elements    ;     /* Number of elements in the mesh    */
  int      bpts, ipts  ;     /* Number of global nodes            */
  int      nvert,            /* Number of global vertices         */
           nedge       ;     /*    "      global edges            */
  int      bandwidth   ;     /* Bandwidth of the boundary matrix  */
  int*     solve       ;     /* Global solve mask                 */
  int*     secr        ;     /* Global secretary node mask        */

  /* NB:       bdof  <=   bpts   <=   ndof   <=  bpts + ipts,     //
  //                                                              //
  //      where bpts + ipts gives the total number of independent //
  //      nodes in mesh.   The degrees of freedom are determined  //
  //      by the boundary conditions (essential vs. natural).     */

  /* Mapping between local and global values */

  struct map {
    int*        npts;        /* Number of points coupled to p     */
    int**       index;       /* Global indices of the points      */
    double**    Z;           /* Projection matrix                 */
  } *map;


  /* Mass Matrix */

  double*  mass        ;     /* Global Mass Matrix (diagonal)     */
  double*  massinv     ;     /* Global Mass Matrix (inverted)     */

  /* Stiffnes Matrix */

  int      singular    ;     /* Singularity flag                  */
  double   constant    ;     /* Helmholtz constant                */
  double   condition   ;     /* Matrix condition number           */

  Matrix_CG     *CG    ;     /* Conjugate Gradient Matrices       */
  Matrix_SC     *SC    ;     /* Static Condensation Matrices      */
  double        *diag  ;     /* Diagonal preconditioner           */

} Matrix;

/* ------------------------------------------------------------------------- */

Matrix* Matrix_alloc  (void);
void    Matrix_free   (Matrix *A);
void    Matrix_update (Matrix *A, Mesh *mesh);

#endif
