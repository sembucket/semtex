/*
 * Statistics
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "prism/prism.h"
#include "prism/stats.h"

/* ------------------------------------------------------------------------- */

static void init (stat_t *s) 
{
  const int ndim = DIM;
  const int nfld = s->nfld;
  const int nstr = (ndim + 1) * ndim >> 1;
  const int npri = nfld - nstr;

  int i;

  char fname[FILENAME_MAX];
  FILE *fp;

  sprintf(fname, "%s.avg", s->name);
  fp = fopen(fname,"r");

  /* If there is no average file then set everything to zero */

  if (fp==NULL) {
    Prism_log(0, "\nStatistics: initialized to zero\n");
    for (i = 0; i < nfld; i++)
      Field_scal (0., s->avg[i]);
  } 


  /* Otherwise, restore the accumulated fields */

  else {
    FieldFile *f = FieldFile_alloc();

    Prism_log(0, "\nStatistics: initialized from %s\n", fname);

    FieldFile_read(f,fp);
    
    /* Restore number of averages */

    s->navg = FIELDFILE_STEP(f);

    /* Restore accumulated fields */
    
    for (i = 0; i < nfld; i++)
      FieldFile_get(f, s->avg[i]);
    
    /* Transform primitive variables */

    for (i = 0; i < npri; i++)
      Transform (s->avg[i], *s->avg[i]->base, Fourier);

    FieldFile_free(f);
    fclose(fp);
  }
}

/* ------------------------------------------------------------------------- */

stat_t *stat_alloc (struct domain *d)
{
  stat_t *s = (stat_t*) calloc(1,sizeof(stat_t));

  const int ndim = DIM;
  const int npri = ndim+1;
  const int nstr = (ndim + 1) * ndim >> 1;
  const int nfld = npri + nstr;

  int i;

  /* Initialize */

  s->name = strdup(d->name);
  s->nfld = nfld;
  s->navg = 0;
  s->src  = (Field**) malloc(npri*sizeof(Field*));
  s->avg  = (Field**) malloc(nfld*sizeof(Field*));
  s->d    = d;

  /* Store pointers to the primitive variables */

  s->src[0]    = d->U;
  s->src[1]    = d->V;
  s->src[2]    = d->W;
  s->src[ndim] = d->P;   /* overwrites ndim=2 */

  /* Create fields for the averages */

  for (i = 0; i < npri; i++)             /* primitives */
    s->avg[i] = Field_dup(s->src[i]);

  for (i = 0; i < nstr; i++) {           /* stresses */
    s->avg[npri + i] = Field_dup(d->U);
    FIELD_TYPE(s->avg[npri+i]) = (char) 'A' + i;
  }

  init(s);

  return s;
}

/* ------------------------------------------------------------------------- */

void stat_free (stat_t *s)
{
  const int nfld = s->nfld;
  int i;

  for (i = 0; i < nfld; i++)
    Field_free(s->avg[i]);

  free(s->src);
  free(s->avg);
  free(s->name);
  free(s);
}

/* ------------------------------------------------------------------------- */

void stat_update (stat_t *s)
{
  const int ndim = DIM;
  const int nfld = s->nfld;
  const int nstr = (ndim + 1) * ndim >> 1;
  const int npri = nfld - nstr;

  int i;

  Field *work[3];

  /* Borrow workspace from level 0 of the multistep arrays to compute the    *
   * velocities in physical space for the Reynolds stress calculation.       */

  work[0] = (s->d->Uf[0]);   
  work[1] = (s->d->Vf[0]);
  work[2] = (s->d->Wf[0]);
  for (i = 0; i < ndim; i++)
    Transform (s->src[i], *work[i]->base, Physical);

  for (i = 0; i < nfld; i++)
    Field_scal ((double) s->navg, s->avg[i]);

  /* Accumulate primitive variables in Fourier space */

  for (i = 0; i < npri; i++)
    Field_axpy (1., s->src[i], s->avg[i]);

  /* Accumulate Reynolds stresses */

  Field_vvtp (work[0], work[0], s->avg[npri+0]);
  Field_vvtp (work[0], work[1], s->avg[npri+1]);
  Field_vvtp (work[1], work[1], s->avg[npri+2]);

  if (ndim == 3) {
    Field_vvtp (work[0], work[2], s->avg[npri+3]);
    Field_vvtp (work[1], work[2], s->avg[npri+4]);
    Field_vvtp (work[2], work[2], s->avg[npri+5]);
  }

  /* Rescale all fields to produce average quantities */
  
  s->navg++;

  for (i = 0; i < nfld; i++)
    Field_scal (1./s->navg, s->avg[i]);
}

/* ------------------------------------------------------------------------- */

void stat_write (stat_t *s)
{
  static char *routine = "stat_write";

  const int ndim = DIM;
  const int nfld = s->nfld;
  const int nstr = (ndim + 1) * ndim >> 1;
  const int npri = nfld - nstr;
  const int init = iparam("STEP");

  FieldFile *f   = FieldFile_alloc();
  FILE      *fp  = NULL;

  int i;

  Prism_log(0, "Writing stats file [step=%d]\n", s->navg);

  ROOTONLY {
    char fname[FILENAME_MAX];
    char buf  [BUFSIZ];

    sprintf (fname, "%s.avg", s->name);

    /* Save a copy of the current file (ignore errors) */

    sprintf (buf, "%s.bak", fname);
    unlink  (buf);
    if (!link(fname,buf))
      unlink(fname);
    
    if ((fp=fopen(fname,"w")) == NULL)
      Prism_error("%s: can't open statistics file", routine);
  }

  for (i = 0; i < npri; i++)
    Transform (s->avg[i], *s->avg[i]->base, Physical);

  for (i = 0; i < nfld; i++)
    FieldFile_put (f, s->avg[i]);

  for (i = 0; i < npri; i++)
    Transform (s->avg[i], *s->avg[i]->base, Fourier);

  FieldFile_setName(f, s->name);
  FIELDFILE_STEP   (f) = s->navg;
  FieldFile_write  (f, fp);
  FieldFile_free   (f);

  if (fp) fclose(fp);
}
