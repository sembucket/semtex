/*
 * Spectral Element Stiffness Matrix
 *
 * Copyright (c) 1994 Ronald Dean Henderson
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

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "speclib/matrix.h"

/* Private functions */

static void     Post         (double a,  BSystem *B, int row, int col);

static void     Memory       (Field *U, BSystem *B),
                Factor       (Field *U, BSystem *B),
                Assemble     (Field *U, BSystem *B, OperatorSC A),
                Condense     (Field *U, BSystem *B);

static Matrix_SC*
                alloc_SC     (Field *U, BSystem *B);
static Matrix_IP*
                build_geof   (Field *U, BSystem *B);
static double*  build_mass   (Field *U, BSystem *B);
static int**    build_bmap   (Field *U, BSystem *B, int [][_MAX_NB], int []);

static int      
    read_mortars (Field *, FILE *,        int [][_MAX_NB]),
    global_mask  (Field *, Bedge*,   int, int [][_MAX_NB], int []),
    bandwidth    (Field *, BSystem *);

/* Run mason to build the connectivity file */

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
 * Matrix_init() - Set up the boundary system                                *
 *                                                                           *
 * This function does the initial setup for the boundary system of a spec-   *
 * tral element mesh, and returns a pointer to a BSystem structure.  The     *
 * function computes:                                                        *
 *                                                                           *
 *        - boundary system parameters (size, bandwidth, ...)                *
 *                                                                           *
 *        - numbering system for the boundary nodes                          *
 *                                                                           *
 *        - pointer arrays for the elemental matrices                        *
 *                                                                           *
 *        - indexing arrays for packing/unpacking the element matrices       *
 *                                                                           *
 * Matrix_init should only be called once to create a new BSystem, which can *
 * then be given to Matrix_build() to build the stiffness matrix.  Boundary  *
 * conditions are used only to decide which points must be included in the   *
 * boundary stiffness matrix, not to assign actual values.                   *
 *                                                                           *
 * To create a family of matrices, i.e., the same essential boundary condi-  *
 * tions but a different Helmholtz constant, you would use a loop like the   *
 * following:                                                                *
 *                                                                           *
 *       N = 10;                                                             *
 *       B = Matrix_init (U, Ubc, "foo");                                    *
 *       for (n = 1; n < N; n++)  {                                          *
 *           B->constant = 1. + n*n;    // Set the constant = 1 + n**2       *
 *           Matrix_build (U, B);       // Compute the Helmholtz matrix      *
 *           Frame_save(U, B, n);       // Save as frame number n            *
 *       }                                                                   *
 *                                                                           *
 *       Frame_load (U, B, 5);          // Load frame 5 and solve the        * 
 *       Solve      (U, F, Ubc, B);     // Helmholtz problem on it           *
 *                                                                           *
 *                                                                           *
 * NB: Frame_save/load must be used to store and initialize a given BSystem  *
 *     before a call to the Helmholtz solver.                                *
 * ------------------------------------------------------------------------- */

BSystem *Matrix_alloc (Field *U, Bedge *Ubc, char *fname)
{
  int solve  [_MAX_NEL * _MAX_NB];
  int scatter[_MAX_NEL] [_MAX_NB], bpts, nsolve;

  const int nel = Field_count(U);
  BSystem*  sys = (BSystem*) calloc (1, sizeof(BSystem));

  /* Open the connectivity file */

  FILE* mor = mason(fname);
  if (mor == (FILE *) NULL)
    speclib_error
      ("Matrix_alloc: "
       "cann't open/create the connectivity (session.mor) file\n");
  
  bpts   = read_mortars (U, mor, scatter);
  nsolve = global_mask  (U, Ubc, bpts, scatter, solve);
  fclose(mor);

  /* Set the members of the new boundary system structure */
    
  sys -> elements  = nel;
  sys -> families  = Family_count(U);
  sys -> bpts      = bpts;
  sys -> ipts      = (U->nr - 2) * (U->ns - 2) * nel;
  sys -> bdof      = nsolve;
  sys -> ndof      = nsolve + sys->ipts;
  sys -> bmap      = build_bmap (U, sys, scatter, solve);
  sys -> CG        = build_geof (U, sys);
  sys -> mass      = build_mass (U, sys);
  sys -> massinv   = dvector    (0, bpts-1);
  sys -> bandwidth = bandwidth  (U, sys);
  
  dvrecp (bpts, sys->mass, 1, sys->massinv, 1);
  
  sys -> singular  = UNSET;   /* All taken care of somewhere else */
  sys -> constant  = 0.;
  sys -> condition = 0.;
  
  sys -> SC        = NULL;
  sys -> Hp        = NULL;
  
  /* If the paramater MWORDS is set, then keep the size  *
   * of the boundary matrix (bdof x (bw + 1)) under this *
   * limit.                                              */
  
  if (iparam("MWORDS")) {
    const int bw = CLAMP (sys->bandwidth, 1, iparam("MWORDS")/nsolve - 1);
    sys->bandwidth = bw;
  }

  return sys;
}

/* ------------------------------------------------------------------------- *
 *  Matrix Static Condensation                                               *
 *                                                                           * 
 *  This function computes the Helmholtz matrices for a given geometry       *
 *  and Helmholtz constant (lambda), packs the boundary matrix into a con-   *
 *  densed form, and then factors it into LU-form.                           *
 * ------------------------------------------------------------------------- */

void Matrix_build (Field *U, BSystem *B)
{
  Matrix_build_A (U, B, HelmholtzSC);
  return;
}

void Matrix_build_A (Field *U, BSystem *B, OperatorSC A)
{
  Memory        (U, B);       /* Allocate system memory */
  Assemble      (U, B, A);    /* Construct the stiffness matrix */
  Condense      (U, B);       /* Perform the static condensation */
  Factor        (U, B);       /* Factor into LU-form for the solvers */

  /* At this point, the elemental boundary matrices can be deleted if  *
   * they aren't going to be used by the semi-direct solver.           */

  return;
}

/* Stub... */

void Matrix_assemble (Field *U, BSystem *B, OperatorSC A)
{
  Assemble (U, B, A);
  return;
}

/* ------------------------------------------------------------------------- *
 * Assemble() - Global Matrix Assembly                                       *
 *                                                                           *
 * Construct the global stiffness matrix from the elemental matrices.        *
 * ------------------------------------------------------------------------  */

static void Assemble (Field *U, BSystem *B, OperatorSC A)
{
  const int nrns = U->nr * U->ns;
  const int nb   = (U->nr + U->ns - 2) << 1;
  const int ni   = nrns - nb;
  const int nel  = B->elements;
  
  /*  Allocate memory for the elemental matrices */

  double** A_bb  = dmatrix (0, nb-1, 0, nb-1);
  double** A_bi  = dmatrix (0, nb-1, 0, ni-1);
  double** A_ii  = dmatrix (0, ni-1, 0, ni-1);

  int k;

  Family_reset();
  for (k = 0; k < nel; k++) {
    int     nrows  = B-> SC[k].nrows;
    int*    rmap   = B-> SC[k].rowmap;
    double* Abnd   = B-> SC[k].A_11;
    double* Acup   = B-> SC[k].A_12;
    Family* fam    = Family_get(&U[k]);
    int i;

    /* Get the matrices */

    A (&U[k], B, A_bb, A_bi, fam->set ? NULL : A_ii);

    for (i = 0; i < nrows; i++, Acup += ni, Abnd += nrows) {
      dcopy (ni,    A_bi[rmap[i]], 1, Acup, 1);
      dgathr(nrows, A_bb[rmap[i]], rmap, Abnd);
    }
    
    /* Factor the interior matrix */

    if (!fam->set) {
      double* C = B-> SC[k].A_22;
      int     p = 0;
      int     j, info;
      
      for (i = 0; i < ni; i++)
	for (j = 0; j <= i; j++)
	  C[p++] = A_ii[i][j];

      dpptrf ('U', ni, C, info);
      Family_set (fam);
    }
  }

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

static void Condense (Field *U, BSystem *M)
{
  const int nel  = M->elements;
  const int ni   = M->ipts / nel;
  int**   bmap   = M->bmap;
  Matrix_SC *SC  = M->SC;
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

    for (i = 0; i < nrows; i++) rmap[i] = bmap[k][rmap[i]];

    
    /* Solve [ A_22 ] [ X ] = [ A_12 ] */

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

static void Memory (Field *U, BSystem *M)
{
  int       bw      = M->bandwidth;
  const int nsolve  = MAX (M->bdof,1);
  size_t    unit_H  = sizeof(double);

  int size_Hp = (bw > 0) ? nsolve * (bw + 1) : nsolve * (nsolve + 1) / 2;

  /* Allocate storage for the matrices.  The following   *
   * line allocates the majority of storage used by the  *
   * code when computing the stiffness matrices.         */

  M -> SC  = alloc_SC (U, M);
  M -> Hp  = (void*) calloc (size_Hp, unit_H);

  return;
}

static void Factor (Field *U, BSystem *M)
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

static int read_mortars (Field *U, FILE *fp, int scatter[][_MAX_NB])
{
  int   nel    = Field_count(U),
        global = 0;
  char  buf [BUFSIZ];
  Edge *edg;
  int   element, face, offset, *p;

  register int k, n, np;

  rewind (fp);

  for (k = 0; k < nel; k++) {
    for (edg = U[k].edges; edg; edg = edg->next) {
      p  = & scatter[k][edg->bindex];
      np = edg->np;
      for (n = 0; n < np; n++, p++) {

	while (fgets (buf, BUFSIZ, fp) && *buf == '#');  /* comments */
	
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

static int global_mask 
  (Field *U, Bedge *Ubc, int bpts, int scatter[][_MAX_NB], int solve[])
{
  Bedge *bc;
  int n;
  
  ifill (bpts, 1, solve, 1);

  /* Loop through the list of boundary points and turn off all of those *
   * corresponding to essential (Dirichlet) boundaries                  */

  for (bc = Ubc; bc; bc = bc->next) {
    if (BC_getType(bc->type)==DIRICHLET) {
      int np    = bc->edge->np;
      int *node = bc->edge->bindex + scatter[bc->elmt->id];
      for (n = 0; n < np; n++)
	solve [ *(node++) ] = 0;
    }
  }

  return icount (bpts, solve, 1);
}

/* Build the global numbering map */

static int **build_bmap
   (Field *U, BSystem *B, int scatter[][_MAX_NB], int solve[])
{
  const int nel  =  B->elements;
  const int bpts =  B->bpts;
  const int nb   = (U->nr + U->ns - 2) << 1;
  int   known    = B->bdof;
  int   unknown  = 0;
  int** bmap     =  imatrix (0, nel-1, 0, nb);
  int   bnodes [_MAX_NEL * _MAX_NB];

  register int k, n;

  for  (n = 0; n < bpts; n++) 
    bnodes[n] = (solve[n]) ? unknown++ : known++;
  
  /* Assign to the global bmap */

  for  (k = 0; k < nel; k++)
    for (n = 0; n <= nb; n++)
      bmap[k][n] = bnodes[scatter[k][n]];

  return bmap;
}

/*
 * Form the global mass matrix 
 */

static double *build_mass (Field *U, BSystem *M)
{
  const int bpts = M->bpts;
  const int nel  = M->elements;
  const int nb   = (U->nr + U->ns - 2) << 1;
  int**     bmap = M->bmap;

  double *mass, tmp [_MAX_NB];
  register int i, k;

  if (!(mass = (double*) calloc (bpts, sizeof(double))))
    speclib_error("build_mass: out of memory");

  for (k = 0; k < nel; k++) {
    dgathr (nb, *U[k].mass, U[k].emap, tmp);
    for (i = 0; i < nb; i++)
      mass [bmap[k][i]] += tmp[i];
  }
  
  return mass;
}

/*
 * Build the Inner Product Geometric Factors
 */

static Matrix_IP *build_geof (Field *U, BSystem *B)
{
  const int  nel = B->elements;
  Matrix_IP* G   = (Matrix_IP*) malloc (nel*sizeof(Matrix_IP));

  Family *fam;
  register int k;

  Family_reset();
  for (k = 0; k < nel; k++)
    if (!(fam = Family_get(& U[k]))->set) {
      geofac (& U[k], & G[k]);
      Family_set(fam);
    } else
      memcpy (& G[k], & G[fam->parent->id], sizeof(Matrix_IP));

  return G;
}

/*
 * Memory allocation for the elemental systems
 */

static Matrix_SC *alloc_SC (Field *U, BSystem *M)
{
  const int nel  = M->elements;
  const int nfam = M->families;
  const int bdof = M->bdof;
  const int nb   = (U->nr + U->ns - 2) << 1;
  const int ni   = M->ipts / nel;

  int**   bmap = M->bmap;
  int*    rmap = NULL;
  double* F    = NULL;
  double* A    = NULL;
  double* C    = NULL;

  int i, k, m;
  int nrows[_MAX_NEL+1];

  /* Allocate space for a Matrix_SC to hold "nel" sub matrices */

  Matrix_SC *sys = (Matrix_SC*) calloc (nel, sizeof(Matrix_SC));

  /* If a previous Matrix_SC array exists, use it... */

  if (M->SC) {
    rmap = M->SC->rowmap;
    for (k = 0; k < nel; k++) {
      sys[k].nrows  = nrows[k] = M->SC[k].nrows;
      sys[k].rowmap = M->SC[k].rowmap;
    }
  } 

  /* ...otherwise, just re-compute everything.  */

  else {
    for (k = 0; k < nel; k++) {
      nrows[k] = 0;
      for (i = nrows[k]; i < nb; i++)
	if (bmap[k][i] < bdof) nrows[k]++;
    }
  }

  nrows[nel] = isum (nel, nrows, 1);

  /* Allocate memory for the A-matrix, C-matrix and F-matrix */

  rmap = (rmap) ? rmap : ivector (0, nrows[nel]);
  if (nrows[nel]) {
    for (m = 0, k = 0; k < nel; k++)
      m += nrows[k]*nrows[k];
    A = dvector (0, m-1);
    F = dvector (0, ni * nrows[nel]- 1);
  }

  C = dvector (0, (ni * (ni+1) / 2) * nfam - 1);

  for (k = 0; k < nel; k++) {
    sys[k].nrows  =  nrows[k];
    sys[k].rowmap =  rmap ; (rmap += nrows[k]);
    sys[k].A_11   =  A    ; (A    += nrows[k]*nrows[k]);
    sys[k].A_12   =  F    ; (F    += nrows[k]*ni      );
    sys[k].A_22   =  C    +  Family_get(U+k)->id * (ni * (ni+1) / 2);
  }

  /* Set up the initial values of "rowmap" */
  
  for (k = 0; k < nel; k++) {
    for (i = 0, m = 0; m < nb; m++) 
      if (bmap[k][m] < bdof) 
	sys[k].rowmap[i++] = m;
  }

  return sys;
}


/* --------------------------------------------------------------------- * 
 * Compute the bandwidth of the BC-corrected boundary stiffness matrix.  *
 * Based on the connectivity stored in bmap and the solve mask.          *
 * --------------------------------------------------------------------- */

static int bandwidth (Field *U, BSystem *B)
{
  const int nb   = (U->nr + U->ns - 2) << 1;
  const int nel  = B->elements;
  const int bdof = B->bdof;

  int** bmap = B->bmap;
  int   bw   = 0;
  int   i, k;

  for (k = 0; k < nel; k++) {
    int b_min  = UNSET;
    int b_max  = -b_min;                    /* reset counters */
    for (i = 0; i < nb; i++)
      if (bmap[k][i] < bdof) {              /* only check members of the   */
	b_min = MIN (b_min, bmap[k][i]);    /* boundary coefficient matrix */
	b_max = MAX (b_max, bmap[k][i]);
      }
    bw = MAX (bw, b_max - b_min);
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
  const int bdof = B->bdof;
  const int bw   = B->bandwidth;

  void*  Hp = B->Hp;
  int    h;

  if (row >= bdof || column >= bdof || (bw > 0 && column-row > bw))     
    return;   /* Out of range */

  if (bw)		
    h = bw * (column + 1) + row;
  else
    h = ((column * (column + 1)) >> 1) + row;
  
  ((double*) Hp)[h] += Kij;

  return;
}

/* ------------------------------------------------------------------------ */

/* Compatability: */

void Field_davg (Field *u, BSystem *B) 
{
  const int nz = FIELD_NZ(u);
  int k;

  for (k = 0; k < nz; k++) {
    Frame_set_one (k, u);
    Matrix_davg   (B, u);
  }
}

void Field_dsum (Field *u, BSystem *B, double *local, double *global) {
  Matrix_dsum(B,u,local,global);
}

void Field_local (Field *u, BSystem *B, double *local, double *global) {
  Matrix_local(B,u,local,global);
}

void Matrix_davg (const BSystem *M, Field *u)
{
  const int bpts   = M->bpts;
  const int nel    = M->elements;
  const int nb     =(u->nr + u->ns - 2) << 1;
  int**     bmap   = M->bmap;

  int i, k;

  tempVector (ub, bpts);
  dzero(bpts, ub, 1);

  for(k = 0; k < nel; k++)  
    for (i = 0; i < nb; i++) 
      ub[bmap[k][i]] += 
	(*u[k].mass ) [u[k].emap[i]] * (*u[k].field) [u[k].emap[i]] ;
  
  dvmul(bpts, M->massinv, 1, ub, 1, ub, 1);

  for(k = 0; k < nel; k++) 
    for (i = 0; i < nb; i++)
      (*u[k].field)[u[k].emap[i]] = ub[bmap[k][i]];

  freeVector (ub);
}

/* Direct Stiffness Summation */

void Matrix_dsum (const BSystem *M, Field *u, double *local, double *global)
{
  const int nel  = M->elements;
  const int bpts = M->bpts;
  const int nb   =(u->nr + u->ns - 2) << 1;
  const int ni   =(u->nr - 2)*(u->ns - 2);
  const int nrns = u->nr * u->ns;
  double *ui = global + bpts;
  
  int i, k;

  memset(global, 0, bpts*sizeof(double));

  for (k = 0; k < nel; k++, local += nrns, ui += ni) {
    int  *bmap = M->bmap[k];
    int  *emap = u[k].emap;

    for (i = 0; i < nb; i++)                   /* Boundary nodes */
      global[bmap[i]] += local[emap[i]];

    dgathr(ni, local, emap+nb, ui);           /* Interior nodes */
  }
}

/* Copy global values to local format */

void Matrix_local (const BSystem *M, Field *u, double *local, double *global)
{
  const int nel  = M->elements;
  const int bpts = M->bpts;
  const int nb   =(u->nr + u->ns - 2) << 1;
  const int ni   =(u->nr - 2)*(u->ns - 2);
  const int nrns = u->nr * u->ns;
  double*   ui   = global + bpts;

  int i, k;

  for (k = 0; k < nel; k++, local += nrns, ui += ni) {
    int *bmap = M->bmap[k];
    int *emap = u[k].emap;

    for (i = 0; i < nb; i++)                   /* Boundary nodes */
      local[emap[i]] = global[bmap[i]];

    dscatr(ni, ui, emap+nb, local);            /* Interior nodes */
  }
}
