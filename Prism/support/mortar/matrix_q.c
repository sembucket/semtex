/*
 * Spectral Element Stiffness Matrix >> Nonconforming
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
#include <string.h>
#include <math.h>
#include <float.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "mortar/mortar.h"

#ifdef _CRAY
#undef  FLT_EPSILON
#define FLT_EPSILON 1.e-6  /* Lower error tolerance for the Cray */
#endif

/* Private functions */

static void     Post         (double a, BSystem *B, int row, int col);

static void     Memory       (Element *U, BSystem *B),
                Condense     (Element *U, BSystem *B),
                Factor       (Element *U, BSystem *B);

static Matrix_SC*
                alloc_SC     (Element *U, BSystem *B);
static Matrix_IP*
                build_geof   (Element *U, BSystem *B);
static double*  build_mass   (Element *U, BSystem *B);
static int    **build_bmap   (Element *U, BSystem *B, 
                                          int [][_MAX_NB], int [], int []),
                build_solve  (Element *U, Bedge *Ubc, int,
                                          int [][_MAX_NB], int [], int []),
                build_map    (Element *U, FILE *fp, int [][_MAX_NB]);

static int      bandwidth    (Element *U, BSystem *B);

/* ... Added for the MORTAR code ... */

static int      
  build_links  (Element *U, Bedge *Ubc, int, int [][_MAX_NB], int[]);


/* Contact mason to build the connectivity file */

static FILE *mason (const char *name)
{
  FILE *fp = (FILE*) NULL;
  char *p;
  char session[FILENAME_MAX];
  char fname  [FILENAME_MAX];

  /* Did somebody pass us <name>.mor? */

  if ((p=strstr(strcpy(session,name),".mor")))
    *p = '\0';

  /* Did somebody pass us <name>.rea? */

  else if ((p=strstr(session,".rea")))
    *p = '\0';

  /* Now we have <session> = <name> stripped of any .mor or .rea suffix.  *
   * Open the .mor file if it exists, otherwise open a pipe to mason.     */

  sprintf (fname,"%s.mor",session);

  if (!(fp=fopen(fname,"r"))) {
    char buf[BUFSIZ];
    sprintf(buf,"mason -n %d %s", iparam("NORDER"), session);
    fp = popen(buf,"r");
  }

  return fp;
}
       
/* ------------------------------------------------------------------------- *
 * BndrySetup() - Set up the boundary system                                 *
 *                                                                           *
 * This function does the initial setup for the boundary system of a spec-   *
 * tral element mesh, and returns a pointer to a BSystem structure.  The     *
 * function computes:                                                        *
 *                                                                           *
 *        - boundary system parameters (size, bandwidth, ...)                *
 *                                                                           *
 *        - numbering system for the boundary nodes                          *
 *                                                                           *
 *        - indexing arrays for packing/unpacking the element matrices       *
 *                                                                           *
 *        - patching connections within the mesh                             *
 *                                                                           *
 * BndrySetup should only be called once to create a new BSystem, which can  *
 * then be given to MatrixSC() to build the stiffness matrix.  Boundary      *
 * conditions are used only to decide which points must be included in the   *
 * boundary stiffness matrix, not to assign actual values.                   *
 *                                                                           *
 * To create a family of matrices, i.e., the same essential boundary condi-  *
 * tions but a different Helmholtz constant, you would use a loop like the   *
 * following:                                                                *
 *                                                                           *
 *       N = 10;                                                             *
 *       B = BndrySetup (U, Ubc, "example");                                 *
 *       for (n = 1; n < N; n++)  {                                          *
 *           B->constant = 1. + n*n;    // Set the constant = 1 + n**2       *
 *           MatrixSC (U, B);           // Compute the Helmholtz matrix      *
 *           SaveFrame(U, B, n);        // Save as frame number n            *
 *       }                                                                   *
 *                                                                           *
 *       LoadFrame (U, B, 5);           // Load frame U.5 and solve the      * 
 *       Solve     (U, F, Ubc, B);      // Helmholtz problem on it           *
 *                                                                           *
 *                                                                           *
 * NB: SaveFrame and LoadFrame must be used to store and initialize a given  *
 *     BSystem before a call to the Helmholtz solver.                        *
 * ------------------------------------------------------------------------- */

BSystem *BndrySetup_q (Element *U, Bedge *Ubc, const char *name)
{
  int solve  [_MAX_NEL * _MAX_NB];
  int scatter[_MAX_NEL] [_MAX_NB], bpts, nsolve;
  int links  [_MAX_NEL * _MAX_NB], nlinks;
  Patch* patch;

  const int nel = Field_count(U);
  BSystem*  sys = (BSystem*) calloc (1, sizeof(BSystem));

  /* Open the connectivity file */

  FILE* mor = mason(name);
  if (mor == (FILE *) NULL)
    speclib_error("couldn't open/create the connectivity (session.mor) file");
  
  /* Build the indexing arrays & initialize Patching */

  bpts   = build_map     (U, mor, scatter);
  patch  = build_patches (U, mor);
  nlinks = build_links   (U, Ubc, bpts, scatter, links);
  nsolve = build_solve   (U, Ubc, bpts, scatter, solve, links);
  
  fclose(mor);

  /* Set the members of the new boundary system structure */

  sys -> other     = patch;
  sys -> elements  = nel;
  sys -> families  = Family_count(U);
  sys -> bpts      = bpts;
  sys -> ipts      = (U->nr - 2) * (U->ns - 2) * nel;
  sys -> bdof      = nsolve;
  sys -> ndof      = nsolve + sys->ipts;

  sys -> bmap      = build_bmap (U, sys, scatter, solve, links);
  sys -> CG        = build_geof (U, sys);
  sys -> mass      = build_mass (U, sys);
  sys -> massinv   = dvector    (0, bpts-1);
  sys -> bandwidth = bandwidth  (U, sys);

  dvrecp (bpts, sys->mass, 1, sys->massinv, 1);

  sys -> singular  = UNSET;  /* All taken care of somewhere else */
  sys -> constant  = 0.;
  sys -> condition = 0.;

  sys -> SC        = NULL;
  sys -> Hp        = NULL;

  /* If the paramater MWORDS is set, then keep the size  *
   * of the boundary matrix (bdof x (bw + 1)) under this *
   * limit.                                              */

  if (iparam("MWORDS")) {
    const int bw = CLAMP (sys->bandwidth, 1, iparam("MWORDS")/nsolve - 1);

    if (sys->bandwidth > bw && option("verbose") > 1)
      fprintf (stderr, "\nReducing %c matrix size to %d diagonals\n", 
	       U->type, bw);
    sys->bandwidth = bw;
  }

  return  sys;
}


/* ------------------------------------------------------------------------- *
 *  Matrix Static Condensation                                               *
 *                                                                           * 
 *  This function computes the Helmholtz matrices for a given geometry       *
 *  and Helmholtz constant (lambda), packs the boundary matrix into a con-   *
 *  densed form, and then factors it into LU-form.                           *
 * ------------------------------------------------------------------------- */

void MatrixSC_q (Element *U, BSystem *M) 
{ 
  MatrixSC_A_q (U, M, HelmholtzSC);
  return;
}

void MatrixSC_A_q (Element *U, BSystem *B, OperatorSC A)
{
  const int    nel = B->elements;
  register int k;

  Memory        (U, B);       /* Allocate system memory */
  Assemble_q    (U, B, A);    /* Construct the stiffness matrix */
  Condense      (U, B);       /* Perform the static condensation */
  Factor        (U, B);       /* Factor into LU-form for the solvers */


  /* Check size parameter to make sure everything is within limits */

  for (k = 0; k < nel; k++) {
    if (B->SC[k].nrows > _MAX_NB)
      speclib_error("MatrixSC_A_q: "
		    "number of linked bnodes exceeds _MAX_NB\n");

#if defined(_CRAY)
    if (B->SC[k].nrows > 64)
      speclib_error("MatrixSC_A_q [CRAY]: "
		    "more than 64 bnodes -> definite shortloop problems\n");
#endif
  }
}


/* ------------------------------------------------------------------------- *
 * Assemble_q() - Global Matrix Assembly                                     *
 *                                                                           *
 * Construct the global stiffness matrix from the elemental matrices.        *
 * ------------------------------------------------------------------------  */

void Assemble_q (Element *U, BSystem *M, OperatorSC A)
{
  Patch  *P     = M->other;
  int    *emask = Patch_list (U, P);

  const int nrns  = U->nr * U->ns;
  const int nb    =(U->nr + U->ns - 2) << 1;
  const int ni    = nrns - nb;
  const int nel   = M->elements;

  int i, k, *rmap, nrows;
  double  **A_bb, **A_bi, **A_ii, *A_bb_dof, *A_bi_dof;

  /*  Allocate memory for the elemental matrices */

  A_bb  = dmatrix (0, nb-1, 0, nb-1);
  A_bi  = dmatrix (0, nb-1, 0, ni-1);
  A_ii  = dmatrix (0, ni-1, 0, ni-1);

  Family_reset();
  for (k = 0; k < nel; k++) {
    Family *fam = Family_get (&U[k]);

    /* Get the matrices */

    A (&U[k], M, A_bb, A_bi, fam->set ? NULL : A_ii);

    /* Save the boundary and coupling matrices */

    nrows    = M-> SC[k].nrows;
    rmap     = M-> SC[k].rowmap; 
    A_bb_dof = M-> SC[k].A_11;
    A_bi_dof = M-> SC[k].A_12;

    if (!emask[k]) {
      for (i = 0; i < nrows; i++, A_bi_dof += ni, A_bb_dof += nrows) {
	dcopy (ni,    A_bi[rmap[i]], 1, A_bi_dof, 1);
	dgathr(nrows, A_bb[rmap[i]], rmap, A_bb_dof);
      }
    } else {

      Qmap  q  = Qmap_alloc (& U[k], M, M_MAP | M_SOLVE | M_PROJ);
      const int nq = q.nq;

      double **A_bb_tmp = dmatrix (0, nq-1, 0, nq-1),
             **A_bi_tmp = dmatrix (0, nq-1, 0, ni-1);

      mult_QAQ (& U[k], q, A_bb, A_bb_tmp);
      mult_QA  (& U[k], q, A_bi, A_bi_tmp);
      
      for (i = 0; i < nrows; i++, A_bi_dof += ni, A_bb_dof += nrows) {
	dcopy (ni,    A_bi_tmp[rmap[i]], 1, A_bi_dof, 1);
	dgathr(nrows, A_bb_tmp[rmap[i]], rmap, A_bb_dof);
      }
    
      free_dmatrix (A_bb_tmp, 0, 0);
      free_dmatrix (A_bi_tmp, 0, 0);
      Qmap_free    (q);
    }

    /* Factor the interior matrix */

    if (!fam->set) {
      double* C = (M->SC)[k].A_22;
      int     p = 0;
      int     j, info;

      for (i = 0; i < ni; i++)
	for (j = 0; j <= i; j++)
	  C[p++] = A_ii[i][j];

      dpptrf ('U', ni, C, info);
      Family_set(fam);
    }
  }

  free (emask);
  free_dmatrix (A_bb, 0, 0);
  free_dmatrix (A_bi, 0, 0);
  free_dmatrix (A_ii, 0, 0);
  return;
}


/* -----------------  P R I V A T E    F U N C T I O N S  ----------------- */


/* ------------------------------------------------------------------------ *
 * Condense() - Static Condensation of Global Matrices                      *
 *                                                                          *
 * This function performs the re-arrangement of the global system matrices  *
 * into a form which allows the solution of the boundary and interior sys-  *
 * tems separately.  This technique, known as "static condensation", may be *
 * applied to any algebraic system which may be written in the form:        *
 *                                                                          *
 *                A  ub  +   B  ui = Fb                     (1a)            *
 *                                                                          *
 *                B' ub  +   C  ui = Fi                     (1b)            *
 *                                                                          *
 * which is then transformed into...                                        *
 *                                                                          *
 *                (A - B C^ B') ub = Fb - B C^ Fi           (2a)            *
 *                                                                          *
 *                            C ui = Fi - B' ub             (2b)            *
 *                                                                          *
 * In these equations, (ub|ui) is the solution at (boundary|interior) nodes *
 * of the spectral element mesh.  The "A" encompasses the global coupling   *
 * of the boundary nodes, the "B" matrix the coupling of boundary|interior, *
 * and "C" the local Helmholtz operator acting only on the interior nodes.  *
 * ------------------------------------------------------------------------ */

static void Condense (Element *U, BSystem *M)
{
  const int  nel       = M->elements;
  const int  ni        = M->ipts / nel;
  int**      bmap      = M->bmap;
  Patch*     P         = M->other;
  int*       emask     = Patch_list (U, P);
  Matrix_SC* SC        = M->SC;
  register int i, j, k;

  /* ------------------------------------------------------------------ *
   * Multiply B * inv(C) and store for future use, then post-multiply   *
   * by trans(B) and subtract from the boundary matrix.  Probably lousy *
   * memory performance because of the indirect addressing.             *
   * ------------------------------------------------------------------ */

  for (k = 0; k < nel; k++) {
    int  nrows = SC[k].nrows;
    int* rmap  = SC[k].rowmap;
    int  info;

    tempVector (B, nrows*ni);

    /* Change the "rmap" from local nodes to global nodes */

    if (emask[k]) {
      Qmap q = Qmap_alloc (& U[k], M, M_MAP);
      for (i = 0; i < nrows; i++) rmap[i] = q.map[rmap[i]];
      Qmap_free(q);
    } else {
      for (i = 0; i < nrows; i++) rmap[i] = bmap[k][rmap[i]];
    }

    /* Solve [ A_22 ][ X ] = [ A_12 ] */

    dcopy  (nrows * ni,     SC[k].A_12, 1, B, 1);
    dpptrs ('U', ni, nrows, SC[k].A_22,    B, ni, info);
    
    /* Compute A_11 - A_21 [ X ] */
    
    dgemm  ('T', 'N', nrows, nrows, ni, -1., SC[k].A_12, ni, B, ni,
	                                 1., SC[k].A_11, MAX(nrows,1));

    for (i = 0; i < nrows; i++) {                     /* ----- Post ---- */
      const int row = rmap[i];
      for (j = 0; j < nrows; j++) {
	const int col = rmap[j];
	if (col >= row) 
	  Post (SC[k].A_11 [i*nrows + j], M, row, col);
      }
    }

    /* Save [ A_12 ] [ A_22 ]^-1 */

    dcopy  (nrows * ni, B, 1, SC[k].A_12, 1);
    
    freeVector (B);
  }
    
  return;
}

/*
 * Allocate memory for the boundary stiffness matrix 
 */

static void Memory (Element *U, BSystem *M)
{
  int       bw      = M->bandwidth;
  const int nsolve  = MAX (M->bdof,1);
  size_t    unit_H  = sizeof(double);

  int size_Hp = (bw > 0) ? nsolve * (bw + 1) : nsolve * (nsolve + 1) / 2;

  /* Allocate storage for the matrices.  The following   *
   * line allocates the majority of storage used by the  *
   * code when computing the stiffness matrices.         */

  M -> SC  =  alloc_SC (U, M);
  M -> Hp  = (void*) calloc (size_Hp, unit_H);

  return;
}

static void Factor (Element *U, BSystem *M)
{
  int     nsolve  =  M -> bdof,
          bw      =  M -> bandwidth, info;
  size_t  unit_H  =  sizeof (double);
  void *  f_work  =  malloc (nsolve * 3 * unit_H);
  int     i_work [_MAX_NEL*_MAX_NB];

  if (bw) { 
    dpbtrf ('U', nsolve, bw, M->Hp, bw+1, info); 
    dpbcon ('U', nsolve, bw, M->Hp, bw+1, 1.,
	    M->condition, f_work, i_work, info);
  } else { 
    dpptrf ('U', nsolve,     M->Hp,       info); 
    dppcon ('U', nsolve,     M->Hp,       1.,
	    M->condition, f_work, i_work, info);
  }
    
  free (f_work);
  return;
}


/* Build the local to global mapping from the connectivity file */

static int build_map (Element *U, FILE *fp, int scatter[][_MAX_NB])
{
  int nel    = Field_count(U);
  int global = 0;
  char  buf [BUFSIZ];
  Edge *edg;
  int   element, face, offset, *p;

  register int k, n, np;

  rewind (fp);

  for (k = 0; k < nel; k++) {
    for (edg = U[k].edges; edg; edg = edg->next) {
      p  = &scatter[k][(*edg).bindex];
      np = edg->np;
      for (n = 0; n < np; n++, p++) {

	while (fgets (buf, 128, fp) && *buf == '#');  /* comments */

	if (sscanf(buf, "%d%d%d%d", &element, &face, &offset, p) != 4 || 
	    element != k+1 || offset != n+1)
	  speclib_error("connectivity file (session.mor) is invalid");
	
	global = MAX(global, *p); (*p)--;
      }
    }    
  }
  
  return global;
}

/* Create the global solve mask */

static int build_solve
  (Element *U, Bedge *Ubc, int bpts, int scatter[][_MAX_NB], int solve[], 
                                                             int links[])
{
  Bedge *bc;
  int n;

  ifill (bpts, 1, solve, 1);      /* All points are initially "unknown" */

  /* Loop through the list of boundary points and turn off all of those *
   * correponding to essential (Dirichlet) boundaries                   */

  for (bc = Ubc; bc; bc = bc->next) {
    if (BC_getType(bc->type)==DIRICHLET) {
      int np    = bc->edge->np;
      int *node = bc->edge->bindex + scatter[bc->elmt->id];
      while (np--) solve [ *(node++) ] = 0;
    }
  }

  /* Now perform an additional "and" with the link mask */

  for (n = 0; n < bpts; n++)
    solve [n] &= links[n];

  return icount (bpts, solve, 1);
}


/* Build the global numbering map */

static int **build_bmap
  (Element *U, BSystem *B, int scatter[][_MAX_NB], int solve[], int link[])
{
  const int nel  =  B->elements,
        bpts =  B->bpts,
        nb   = (U->nr + U->ns - 2) << 1;
  int **bmap =  imatrix (0, nel-1, 0, nb),
        bnodes [_MAX_NEL * _MAX_NB];

  register int k, n, known, unknown, unlinked;

  unknown  = 0;
  unlinked = B->bpts;
  known    = B->bdof;

  for  (n = 0; n < bpts; n++) 
    if (link[n])
      bnodes[n] = (solve[n]) ? unknown++ : known++;
    else
      bnodes[n] = (--unlinked);

  /* Assign to the global bmap */

  for (k = 0; k < nel; k++)
    for (n = 0; n <= nb; n++)
      bmap[k][n] = bnodes[scatter[k][n]];

  return bmap;
}

/*
 * Build the Inner Product Geometric Factors
 */

static Matrix_IP *build_geof (Element *U, BSystem *B)
{
  const int nel = B->elements;
  Matrix_IP *G = (Matrix_IP*) malloc (nel*sizeof(Matrix_IP));
  Family *fam;
  int k;

  Family_reset();
  for (k = 0; k < nel; k++) {
    fam = Family_get(&U[k]);

    if (fam->set) {
      memcpy (& G[k], & G[fam->parent->id], sizeof(Matrix_IP));
    } else {
      geofac (& U[k], & G[k]);
      Family_set(fam);
    }
  }

  return G;
}


/*
 * Form the global mass matrix 
 */

static double *build_mass (Element *U, BSystem *M)
{
  Patch *P      = M->other;
  int   *emask  = Patch_list (U, P);

  int    bpts   = M->bpts,
         nel    = M->elements,
         nb     =(U->nr + U->ns - 2) << 1,
       **bmap   = M->bmap;

  double *mass, tmp [_MAX_NB];
  register int i, k;

  if (!(mass = (double*) calloc (bpts, sizeof(double))))
    speclib_error("build_mass: out of memory");

  for (k = 0; k < nel; k++) {
    if (!emask[k]) {                           /* Conforming */
      dgathr (nb, *U[k].mass, U[k].emap, tmp);
      for (i = 0; i < nb; i++)
	mass [bmap[k][i]] += tmp[i];

    } else {                                   /* Nonconforming */

      Qmap  q  = Qmap_alloc (& U[k], M, M_MAP | M_SOLVE | M_PROJ );
      const int nq = q.nq;

      double ** B  = dmatrix (0, nb-1, 0, nb-1),
             **QBQ = dmatrix (0, nq-1, 0, nq-1);
      
      dzero (nb*nb, *B, 1);
      dgathr(nb, *U[k].mass, U[k].emap, tmp);

      for (i = 0; i < nb; i++)
	B [i][i] = tmp[i];

      mult_QAQ (& U[k], q, B, QBQ);

      /* Add the lumped matrix to the global array */

      for (i = 0; i < nq; i++)
	mass [ q.map[i] ] += dsum (nq, QBQ[i], 1);

      free_dmatrix (QBQ, 0, 0);
      free_dmatrix ( B , 0, 0);
      Qmap_free    (q);
    }
  }

  while (mass[--bpts] < FLT_EPSILON)   /* Eliminate the "unlinked" points */
    mass[bpts] = -1.;

  free (emask);
  return mass;
}

/*
 * Memory allocation for the elemental systems
 */

static Matrix_SC *alloc_SC (Element *U, BSystem *M)
{
  const int nel   = M->elements;
  const int nfam  = M->families;
  const int bdof  = M->bdof;
  const int nb    = (U->nr + U->ns - 2) << 1;
  const int ni    = M->ipts / nel;

  int**     bmap  = M->bmap;
  int*      rmap  = NULL;
  double*   F     = NULL;
  double*   A     = NULL;
  double*   C     = NULL;

  Patch*    P     = M->other;
  int*      emask = Patch_list (U, P);

  int i, k, m;
  int nrows[_MAX_NEL+1];

  /* Allocate space for a Matrix_SC to hold 'nel' matrices */

  Matrix_SC *new = (Matrix_SC*) calloc (nel, sizeof(Matrix_SC));

  /* If a previous Matrix_SC array exists, use it... */

  if (M->SC) {
    rmap   = M->SC->rowmap;
    for (k = 0; k < nel; k++) {
      new[k].nrows  = nrows[k] = M->SC[k].nrows;
      new[k].rowmap            = M->SC[k].rowmap;
    }
  }

  /* ...otherwise, just re-compute everything.  */
    
  else {
    for (k = 0; k < nel; k++) {

      /* Check to see if this is a child.  If so, get the list of *
       * mapped nodes from the element's qmap...                  */

      if (emask[k]) {
	Qmap  q  = Qmap_alloc (& U[k], M, M_SOLVE | M_MAP);
	nrows[k] = icount (q.nq, q.solve, 1);
	Qmap_free(q);
      }

      /* ...otherwise, compute it the "normal" way. */
      
      else {
	nrows[k] = 0;
	for (nrows[k] = i = 0; i < nb; i++)
	  if (bmap[k][i] < bdof) nrows[k]++;
      }
    }
  }
  
  nrows[nel] = isum (nel, nrows, 1);

  /* Allocate memory for the vectors... */

  rmap = (rmap) ? rmap : ivector(0, nrows[nel]);
  if (nrows[nel]) {
    F = dvector (0, ni * nrows[nel] - 1);
    for (k = m = 0; k < nel; k++)
      m += nrows[k]*nrows[k];
    A = dvector (0, m-1);
  }
  C = dvector (0, (ni * (ni+1) / 2) * nfam - 1);

  for (k = 0; k < nel; k++) {
    new[k].nrows  = nrows[k];
    new[k].rowmap = rmap  ; (rmap += nrows[k]);
    new[k].A_11   = A     ; (A    += nrows[k]*nrows[k]);
    new[k].A_12   = F     ; (F    += nrows[k]*ni      );
    new[k].A_22   = C     +  Family_get(U+k)->id * (ni * (ni+1) / 2);
  }

  /* Set up the initial values of "rowmap" */
  
  for (k = 0; k < nel; k++) {
    if (emask[k]) {
      Qmap q = Qmap_alloc (& U[k], M, M_SOLVE | M_MAP);
      for (i = m = 0; m < q.nq; m++)
	if (q.map[m] < bdof) 
	  new[k].rowmap[i++] = m;
      Qmap_free (q);
    } else {
      for (i = m = 0; m < nb; m++) {
	if (bmap[k][m] < bdof) 
	  new[k].rowmap[i++] = m;
      }
    }
  }

  free (emask);
  return new;
}

/* --------------------------------------------------------------------- * 
 * Compute the bandwidth of the BC-corrected boundary stiffness matrix.  *
 * Based on the connectivity stored in bmap and the solve mask.          *
 * --------------------------------------------------------------------- */

static int bandwidth (Element *U, BSystem *B)
{
  int  nb     =(U->nr + U->ns - 2) << 1,
       nel    = B->elements,
     **bmap   = B->bmap, 
       bdof   = B->bdof,
       bw     = 0;
  register int i, k, b_max, b_min;

  for (k = 0; k < nel; k++) {
    b_max  = -(b_min = UNSET);              /* reset counters */
    for (i = 0; i < nb; i++)
      if (bmap[k][i] < bdof) {              /* only check members of the   */
	b_min = MIN (b_min, bmap[k][i]);    /* boundary coefficient matrix */
	b_max = MAX (b_max, bmap[k][i]);
      }
    bw = MAX(bw, b_max - b_min);
  }

  return bw;
}

/* ------------------------------------------------------------------------- *
 * Post() - Global Assembly                                                  *
 *                                                                           *
 * Add elemental matrix contributions to the right spot in the packed global *
 * matrix.  Any entry outside the range [1..bdof][1..bw] is ignored.         *
 * ------------------------------------------------------------------------- */

static void Post (double Kij, BSystem *B, int row, int column)
{
  const int  bdof = B->bdof,
         bw   = B->bandwidth;
  void  *Hp   = B->Hp;
  int    h;

  /* Out of range? */

  if (row >= bdof || column >= bdof || (bw > 0 && column-row > bw))     
    return;   

  if (bw)		
    h = bw * (column + 1) + row;
  else
    h = ((column * (column + 1)) >> 1) + row;
  
  ((double*) Hp)[h] += Kij;

  return;
}

/* ------------------------------------------------------------------------- *
 * build_links                                                               *
 *                                                                           *
 * A logical list of nodes which are "unlinked" from the boundary system, as *
 * in a child edge's interior.                                               *
 * ------------------------------------------------------------------------- */

static int build_links 
  (Element *U, Bedge *Ubc, int bpts, int scatter[][_MAX_NB], int links[])
{
  Bedge *bc;
  int i;

  ifill (bpts, 1, links, 1);             /* All links are initially on */

  for (bc = Ubc; bc ; bc = bc->next)     /* Turn OFF all child edges */
    if (bc->type == 'S') {
      const int k      = bc->elmt->id,
            np     = bc->edge->np,
            bindex = bc->edge->bindex;
      
      for (i = 0; i < np; i++)
	links [scatter[k][bindex + i]] = 0;
    }
  
  for (bc = Ubc; bc ; bc = bc->next)     /* Turn ON all standard vertices */
    if (bc->type != 'S') {
      const int k      = bc->elmt->id,
            n      = bc->edge->np-1,
            bindex = bc->edge->bindex;
      
      links [scatter[k][bindex + 0]] = 1;
      links [scatter[k][bindex + n]] = 1;
    }
  
  return icount (bpts, links, 1);
}

