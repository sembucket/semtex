/*
 * Type definitions for mortar-based patching
 */

#ifndef  MORTAR_H
#define  MORTAR_H

#include "speclib/element.h"

#define M_MAP     1              /* Flag for Q.map    */
#define M_SOLVE   2              /* Flag for Q.solve  */
#define M_PROJ    4              /* Flag for Q.proj   */

/* Variable Types used for the Mortar Patching Code */

typedef struct branch  {         /* ...........  BRANCH  ............. */
  int            id       ;      /* branch number                      */
  double         s0       ;      /* mortar offset                      */
  double         gamma_s  ;      /* integration strip length           */
  double         **q      ;      /* projection operator                */
  struct mortar  *m       ;      /* pointer to a mortar                */
  struct branch  *next    ;      /* next branch along this segment     */
} Branch;

typedef struct segment {         /* ...........  SEGMENT  ............ */
  int            id       ;      /* segment number (from input file)   */
  int            type     ;      /* segment type (Master or Slave)     */
  int            n        ;      /* polynomial order along the segment */
  int            vleft    ;      /* flag for virtual vertices on the   */
  int            vright   ;      /* left and right of the segment.     */
  int            branches ;      /* number of branches (multiplicity)  */
  double         gamma_k  ;      /* segment length                     */
  double         *mesh    ;      /* nodal coordinates                  */
  double         *u       ;      /* solution                           */
  struct patch   *parent  ;      /* parent patch                       */
  struct edge    *edge    ;      /* parent edge                        */
  struct element *elmt    ;      /* parent element                     */
  struct branch  *blist   ;      /* branches to the mortars            */
  struct segment *next    ;      /* pointer to the next segment        */
} Segment;

typedef struct mortar  {         /* ...........  MORTAR  ............. */
  int            id       ;      /* mortar identification number       */
  int            n        ;      /* polynomial order along the mortar  */
  int            *solve   ;      /* solve mask                         */
  double         gamma_p  ;      /* mortar length in (x,y)-coords.     */
  double         *mesh    ;      /* nodal coordinates                  */
  double         *phi     ;      /* solution on the mortar             */
  struct element *elmt    ;      /* parent element                     */
  struct edge    *edge    ;      /* parent edge                        */
  struct mortar  *next    ;      /* pointer to the next mortar         */
} Mortar;

typedef struct patch   {         /* ............  PATCH  ............. */
  int            id       ;      /* patch ID number (from input file)  */
  int            nodes    ;      /* number of mortar nodes             */
  double         x0, y0   ;      /* patch origin in (x,y)-coords.      */
  struct segment *slaves  ;      /* slave segments (interp. at vertex) */
  struct segment *masters ;      /* master segments (exact continuity) */
  struct mortar  *mortars ;      /* the list of mortar segments        */
  struct patch   *next    ;      /* pointer to the next patch          */
} Patch;

typedef struct Qmap {            /* ......... Projection Map ......... */
  int      nq     ;              /* number of linked points            */
  int*     map    ;              /* global boundary map                */
  int*     solve  ;              /* solve mask for each point          */
  double** proj   ;              /* global->local projection matrix    */
} Qmap;


/*
 * Function prototypes
 */

Patch*     build_patches (Field *U, FILE *fp);

void       Project_s_m   (Field *U, Patch *plist, double *u);
void       Project_m_s   (Field *U, Patch *plist, double *u);

int*       Patch_list    (Field *U, Patch *P);
Qmap       Qmap_alloc    (Field *U, BSystem *M, int flags);
void       Qmap_free     (Qmap   q);

void       mult_QAQ      (Element *U, Qmap q, double **A, double **QAQ);
void       mult_QA       (Element *U, Qmap q, double **A, double **QA );
void       mult_AQ       (Element *U, Qmap q, double **A, double ** AQ);

#ifdef DEBUG
void       show_patches  (Patch *p);
#endif

/* Special library functions */

BSystem   *BndrySetup_q  (Field *U, Bedge *Ubc, const char *name);
void       MatrixSC_q    (Field *U, BSystem *M);
void       MatrixSC_A_q  (Field *U, BSystem *M, OperatorSC A);
void       Assemble_q    (Field *U, BSystem *B, OperatorSC A);
void       Solve_CG_q    (Field *U, Field *F, Bedge *Ubc, BSystem *B);
void       Solve_q       (Field *U, Field *F, Bedge *Ubc, BSystem *B);
void       Solve_A_q     (Field *U, Field *F, Bedge *Ubc, BSystem *B, 
			                      Operator A);
void       Solve_CG_A_q  (Field *U, Field *F, Bedge *Ubc, BSystem *B, 
                                              Operator A, double  *M);
double     Form_RHS_q    (Field *U, Field *F, Bedge *Ubc, BSystem *B, 
			                      Operator A, double *r);
void       Mask_RHS_q    (Field *U, BSystem *B, double *resid);

/* Translation Macros 
 *       Conforming >>>> Nonconforming   */

/* Old Library */

#define  BndrySetup      BndrySetup_q
#define  MatrixSC        MatrixSC_q
#define  Assemble        Assemble_q
#define  Solve           Solve_q
#define  Solve_CG        Solve_CG_q
#define  dssum           qssum
#define  dssum3d         qssum3d

/* New Library */

#define  Matrix_alloc    BndrySetup_q
#define  Matrix_build    MatrixSC_q
#define  Matrix_assemble Assemble_q
#define  Field_davg      qssum3d
#define  Field_davg_3D   qssum3d

#endif
