/*
 * FUNCTIONS FOR SWAPPING MATRIX SETS (frames)
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
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include "speclib.h"
#include "tree.h"

#include "veclib/veclib.h"

/* Structure for holding matrix file info */

enum StatusFlags { Unknown, OnFile, InMemory, PMAGIC = 0xefefef };

typedef struct ms_tag {      /* .......... Matrix Set ........... */
  int            id     ;    /* Matrix ID number                  */
  char           type   ;    /* Field type                        */
  int            status ;    /* current state of the matrix set   */
  int*           famids ;    /* list of element family ID's       */
  union {                    /* Storage area...                   */
    FILE         *fp    ;    /*   Set is on file                  */
    BSystem      *mp    ;    /*   Set is loaded in memory         */
  } storage;                 /*                                   */
} *MatrixSet;                /* ................................. */

/* Private functions */

static MatrixSet  mset_add   (int id, Field *U),
                  mset_get   (int id, Field *U);
extern void       mset_show  (MatrixSet M);
static int        mset_load  (MatrixSet M, BSystem *B),
                  mset_save  (MatrixSet M, BSystem *B),
                  mset_copy  (MatrixSet M, BSystem *B);

static void       cleanup    (BSystem *B);

/* External variables */

Tree* matTree              = 0;     /* matrix tree */
char  matdir[FILENAME_MAX] = "";    /* matrix directory name */

/* ------------------------------------------------------------------------ *
 * Frame_init() - Initialize a MatrixSet                                    *
 *                                                                          *
 * The following convention is used for naming matrix files:                *
 *                                                                          *
 *    File name:    name.TYPE.ID.mat     where TYPE  = [ u, p, ...]         *
 *                                             FRAME = [ 0, 1, ...]         *
 *                                                                          *
 * Returns: 0 if a matrix set was successfully loaded                       *
 *          1 otherwise                                                     *
 * ------------------------------------------------------------------------ */ 

int Frame_init (Field *U, BSystem *B, int id, char *name)
{
  FILE *fp;
  MatrixSet M;
  char fname[FILENAME_MAX];

  cleanup (B);     /* de-allocate memory, etc. */

  if (!matTree) matTree = create_tree (mset_show, free);
  if (!(M = mset_get (id, U))) M = mset_add (id, U);

  /* See if the matrix file exists */

  sprintf (fname, "%s%s.%c.%d.mat", matdir, name, U->type, id);
  if((M->storage.fp = fp = fopen (fname, "rb")) == (FILE*) NULL ) {
    fp  = option("keep") ? fopen (fname, "wb+") : tmpfile();
    if (fp) {                    
      M->storage.fp = fp;   /* File opened OK.               */
      return 1;             /* Return 1 and compute the set. */
    } else 
      speclib_err ("unable to create a new matrix file");
  }    



  /* The file exists and has been opened successfully.  Try to load it. */

  if (mset_load(M, B) == Unknown) {

    /* The load failed.  Remove the bad file and call Frame_init *
    // again to try and set up a new file.                       */

    if (unlink(fname) == -1)
      speclib_err ("unable to remove a bad matrix file");
    else
      return Frame_init (U, B, id, name);    /* Try to recompute */
    }
  

  /* If we make it here, then the file was opened and the matrix set *
  // succesfully loaded.  The preprocessor doesn't have to worry     * 
  // about this matrix set any more.  The call to Frame_save will    * 
  // save the MatrixSet in memory if "core" is active.               */

  Frame_save (U, B, id);

  return 0;
}

/*
 * Load a matrix set
 */

void Frame_load (Field *U, BSystem *B, int id)
{
  MatrixSet M = mset_get(id, U);

  if (!M) {
    sprintf (speclib_err_msg, "Frame_load: no matrix for %c frame %d",
	U->type, id);
    speclib_err(NULL);
  }

  switch (M->status) {
  case OnFile:
    mset_load (M, B);
    break;

  case InMemory:
    mset_copy (M, B);
    break;

  case Unknown:
  default:
    sprintf (speclib_err_msg, "Frame_load: "
	     "don't know how to load the matrix for frame %c.%d", U->type, id);
    speclib_err (NULL);
    break;
  }
  return;
}

/*
 * Save the current frame
 */

void Frame_save (Field *U, BSystem *B, int id)
{
  const      keep = option("keep"),
             core = option("core"),
             prep = option("prep");
  MatrixSet  M    = mset_get (id, U);

  switch (M->status) {
  case Unknown:             /* MatrixSet has been computed but not saved. */
    if ( keep ||            /* If either "keep" is on (default)...        */
	 prep || !core )    /* ...or we're pre-processing, save the file. */

      mset_save (M, B);                    /* ...changes status to OnFile */
    
    /* ..... Fall through to the next level .....*/

  case OnFile:                      /* Loaded from file but not saved yet */

    /* If currently running in pre-processor mode or swapping the matrix  */
    /* files from disk, break out of the switch.  Otherwise, the matrices */
    /* will be saved in memory for the remainder of this run.             */

    if (prep || !core) break;     

    /* The next section is the ONLY place status is set to InMemory */

  case InMemory:
    fclose (M->storage.fp);
    memcpy (M->storage.mp = (BSystem *) malloc (sizeof(BSystem)),	    
               B, sizeof(BSystem));             /* save the current setup */
    B->SC     = NULL;                           /* flag B for re-building */
    B->Hp     = NULL;                           /*                        */
    M->status = InMemory;                       /* change status of M     */
    break;
    
  default: 
    speclib_err ("??? logic problem in Frame_save ???");
    break;
  }
  return;
}
  
/*
 * load a matrix set from a file 
 */

#define READ(ptr,size,nobj,fp)\
   if (fread(ptr,size,nobj,fp) != nobj) goto loadError

static int mset_load (MatrixSet M, BSystem *B)
{
  register k;
  BSystem  tmp;
  FILE*    fp = M->storage.fp;
  int      nobj, nrows, nsolve, nknown, size_Hp, eipts;
  size_t   unit_H = sizeof(double);

  cleanup (B);    /* free memory and prep for the load  */
  rewind (fp);    /* rewind the file for multiple reads */
  

  /* Check # of elements, # of families, global and local size *
   * parameters, and finally the Helmholtz constant.  All of   *
   * these have to match for the matrix set on file to work.   */

  READ (&tmp, sizeof(BSystem), 1, fp);
  if   (memcmp(&tmp, B, 7 * sizeof(int)) || 
          tmp.constant * FLT_EPSILON < fabs(tmp.constant - B->constant))
             goto loadError;
  B->condition  = tmp.condition;

  nsolve  = MAX(B->bdof, 1);
  nknown  = MAX(B->bpts - nsolve, 1);
  eipts   = B->ipts / B->elements;
  size_Hp = B->bandwidth ? nsolve*(B->bandwidth+1) : nsolve*(nsolve+1)/2;


  /* Read in the information for the mapping */

  B->SC  = (Matrix_SC *) calloc (B->elements, sizeof(Matrix_SC));
  for (k = 0, nrows = nobj = 0; k < (*B).elements; k++) {
    READ   (&B->SC[k].nrows, sizeof(int), 1, fp);
    nrows += B->SC[k].nrows;
    nobj  += B->SC[k].nrows * B->SC[k].nrows;
  }

  /* Allocate memory */

  B->SC->rowmap = ivector (0, nrows-1);
  B->SC->A_11   = dvector (0, nobj - 1);
  B->SC->A_12   = dvector (0, eipts * nrows - 1);
  B->SC->A_22   = dvector (0, eipts * eipts * B->families - 1);
  B->Hp         = (void*) malloc (size_Hp * unit_H);


  /* Read the rest of the file */

  READ (B->SC->rowmap, sizeof(int), nrows, fp);    
  READ (B->SC->A_11, sizeof(double), nobj, fp);
  nobj = eipts * nrows;
  READ (B->SC->A_12, sizeof(double), nobj, fp);
  nobj = eipts * eipts * B->families;
  READ (B->SC->A_22, sizeof(double), nobj, fp);     

  READ (B->Hp, unit_H, size_Hp, fp);


  /* Align pointers for the static condensation matrices */

  for (k = 1; k < (*B).elements; k++) {
    nrows           = B->SC[k-1].nrows;
    B->SC[k].rowmap = B->SC[k-1].rowmap    + nrows;
    B->SC[k].A_11   = B->SC[k-1].A_11 + nrows * nrows;
    B->SC[k].A_12   = B->SC[k-1].A_12 + eipts * nrows;
    B->SC[k].A_22   = B->SC->A_22     + eipts * eipts * M->famids[k]; 
  }
  
  return M->status = OnFile;    /* Set the status flag for M */

 loadError:
  perror  ("speclib: Frame_load");
  cleanup (B);
  return M->status = Unknown;
}

#undef READ


/*
 * copy a matrix set
 */

static int mset_copy (MatrixSet M, BSystem *B)
{
  memcpy (B, M->storage.mp, sizeof(BSystem));
  return InMemory;
}

/*
 * Save a matrix set 
 */

#define WRITE(ptr,size,nobj,fp)\
  if (fwrite (ptr,size,nobj,fp) != nobj) goto saveError

static int mset_save (MatrixSet M, BSystem *B)
{
  const  nel     = B->elements,
         bw      = B->bandwidth,
         eipts   = B->ipts / nel,
         nsolve  = MAX(B->bdof, 1),
         nknown  = MAX(B->bpts - nsolve,1),
         size_Hp = (bw > 0) ? nsolve * (bw + 1) : nsolve * (nsolve + 1) / 2;
  size_t unit_H  = sizeof(double);
  FILE   *fp     = M->storage.fp;
  register k, nrows, nobj;

  WRITE (B, sizeof(BSystem), 1, fp);
  for   (k = 0, nrows = nobj = 0; k < nel; k++) {
    nrows += B->SC[k].nrows;
    nobj  += B->SC[k].nrows * B->SC[k].nrows;
    WRITE  (&B->SC[k].nrows, sizeof(int), 1, fp);
  }


  WRITE (B->SC->rowmap, sizeof(int), nrows, fp);
  WRITE (B->SC->A_11, sizeof(double), nobj, fp);
  nobj = eipts * nrows;
  WRITE (B->SC->A_12, sizeof(double), nobj, fp);
  nobj = eipts * eipts * B->families;
  WRITE (B->SC->A_22, sizeof(double), nobj, fp);

  WRITE (B->Hp, unit_H, size_Hp, fp);

  fflush(fp);
  return M->status = OnFile;
  
 saveError:
  if (!(option ("core") && option ("keep")))
    speclib_err ("unable to save a matrix file");
  return M->status = Unknown;
}

#undef WRITE

/* Show a MatrixSet */

void mset_show (MatrixSet M)
{
  printf ("Matrix %c.%d: status = ", M->type, M->id);
  switch (M->status) {
  case OnFile:
    puts ("on file");   break;
  case InMemory:
    puts ("in memory"); break;
  default:
    puts ("unknown");   break;
  }
  return;
}

/* Find a MatrixSet */

static MatrixSet mset_get (int id, Field *U)
{
  Node*     np;
  char      key[16];
  MatrixSet mp = NULL;

  sprintf (key, "%c%d", U->type, id);

  if (np = tree_search (matTree->root, key))
    mp = (MatrixSet) np->other;

  return mp;
}

/* Create a new MatrixSet */

static MatrixSet mset_add (int id, Field *U)
{
  int       nel = Field_count (U);
  char      key[16];
  Node*     np;
  MatrixSet new;
  register  k;

  sprintf (key, "%c%d", U->type, id);
  new = (MatrixSet) calloc(1, sizeof(struct ms_tag));

  new->famids   = ivector (0, nel-1);     /* initialize the MatrixSet */
  new->status   = Unknown;
  new->type     = U->type;
  new->id       = id;

  for (k = 0; k < nel; k++)                    /* get the family ID's */
    new->famids [k] = Family_get(U + k)->id;

  np = create_node (key);        /* Create a node to store the matrix */
  np->other = (void*) new;
  tree_insert (matTree, np);

  return new;
}   

/*
 * Clean up the memory associated with a matrix set
 */

static void cleanup (BSystem *B)
{
  if (B->SC) {              /* free the interior matrices */
    free (B->SC->A_11);
    free (B->SC->A_12);
    free (B->SC->A_22);
    free (B->SC->rowmap);
    free (B->SC);
  }

  if (B->Hp) free (B->Hp);  /* free the boundary matrix */

  B->SC  =  NULL;           /* set everybody back to NULL to re-allocate */
  B->Hp  =  NULL;           /* memory on the next pass through.          */

  return;
}

/*
 * Set a new frame pointer
 */

void Frame_set (int frame_number, int nfields, ...)
{
  va_list ap;
  int     frame_index, nel;
  Field  *U;
  register int i, k;

  if (nfields < 1)
    speclib_err ("Frame_set called with less than 1 frame");
  
  va_start(ap, nfields);
  for (i = 0; i < nfields; i++) {
    U           = (Field*) va_arg(ap, Field*);
    nel         = Field_count (U);
    frame_index = U->ns * nel * frame_number;

    for (k = 0; k < nel; k++) {
      U[k].field = U[k].base + frame_index;
      U[k].frame = frame_number;
    }
  }
  
  va_end(ap);
  return;
}

void Frame_set_one (int frame, Field *U)
{
  int nel   = Field_count(U);
  int index = U->ns * nel * frame;
  register k;

  for(k = 0; k < nel; k++) {
    U[k].field = U[k].base + index;
    U[k].frame = frame;
  }

  return;
}
