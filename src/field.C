/*****************************************************************************
 * FIELD.C:  routines to deal with whole fields.                             *
 *****************************************************************************/

static char
RCSid[] = "$Id$";

#include "Fem.h"
#include <time.h>

#define  readLine(fp)  fgets(buf, STR_MAX, (fp))





Element *copyField (      FieldKind  kind,
		    const Element   *E   ,
		          char       name)
/* ========================================================================= *
 * Use one field as basis for producing another with same geometric factors, *
 * but with new (uninitialized) value storage area.  New field solve mask    *
 * and bmap storage is also allocated and is initialized by copying the      *
 * input field, so may need to be reset on return.                           *
 *                                                                           *
 * NB: this version assumes the input Element list is contiguous.            *
 * ========================================================================= */
{
  char          routine[] = "copyField";
  register int  i, j;

  int      nel       = countElmts (E); /* Can be a subset of whole field. */
  int      np        = E->np;
  int      ntot      = nel * SQR (np)* SQR (np);
  Element *copy      = (Element *) calloc (nel, sizeof(Element));


  if (!copy) message (routine, "unable to allocate copied field", ERROR);

  memcpy (copy, E, nel*sizeof(Element));

  for (i=0; i<nel-1; i++)
    (copy + i) -> next = copy + i + 1;
  (copy + i) -> next = NULL;

  if (!(copy->value = (double **) malloc (np * nel * sizeof(double*))))
    message (routine, "can't allocate value pointers", ERROR);
  if (!(*copy->value = dvector (0, ntot-1)))
    message (routine, "can't allocate storage area",   ERROR);
  if (!(copy->solve  = ivector (0, nel * 4 * (np - 1))))
    message (routine, "can't allocate solve mask",     ERROR);
  if (!(copy->bmap   = ivector (0, nel * 4 * (np - 1))))
    message (routine, "can't allocate field bmap",     ERROR);

  memcpy (copy->solve, E->solve, nel * 4 * (np - 1) * sizeof (int)); 
  memcpy (copy->bmap , E->bmap,  nel * 4 * (np - 1) * sizeof (int));

  for (i=0; i<nel; i++) {
    (copy+i)->fldtype  = kind;
    (copy+i)->name     = name;
    (copy+i)->value    = copy->value + i*np;
    (copy+i)->value[0] = copy->value[0] + i*SQR(np);
    (copy+i)->solve    = copy->solve + i*4*(np - 1);
    (copy+i)->bmap     = copy->bmap  + i*4*(np - 1);
    for (j=0; j<np; j++)
      (copy+i)->value[j] = (copy+i)->value[0] + j*np;
  }

  return copy;
}





void setField (Element *E    ,
	       double   value)
/* ========================================================================= *
 * Set field storage area to value.                                          *
 * ========================================================================= */
{
  dfill (iparam("NEL")*SQR(E->np)*SQR(E->np), value, *E->value, 1);
}





void addField (Domain *D, Element *E)
/* ========================================================================= *
 * Add field E to the list held by Domain D.                                 *
 * ========================================================================= */
{
  int       i;
  Element **X = (Element **) calloc (D -> nField + 1, sizeof (Element *));

  
  for (i = 0; i < D -> nField; i++) X[i] = D -> u[i];
  X[i] = E;
  free (D -> u);
  D -> u = X;
  strncat (D -> fTag, &(E->name), 1);
  D -> nField++;
}





void writeFields (const Domain *D)
/* ========================================================================= *
 * Put all the problem variable fields on D->field_fp.                       *
 *                                                                           *
 * Write SEM-compatible form.                                                *
 * ========================================================================= */
{
  char  routine[] = "writeFields";
  char  buf[STR_MAX];
  char *hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields Written\n",
    "%-25s "    "Format\n"
  };
  int     nfields = D -> nField;
  int     nel     = D -> nEl;
  int     np      = (D -> u)[0] -> np;
  int     ntot    = SQR(np) * nel;
  FILE   *fp      = D -> field_fp;

  time_t  tp      = time ((time_t*) NULL);


  strftime (buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  
  fprintf (fp,  hdr_fmt[0], D -> name                              );
  fprintf (fp,  hdr_fmt[1], buf                                    );
  sprintf (buf, "%1d %1d %1d %1d", np, np, 1, nel                  );
  fprintf (fp,  hdr_fmt[2], buf                                    );
  fprintf (fp,  hdr_fmt[3], D -> step                              );
  fprintf (fp,  hdr_fmt[4], D -> time                              );
  fprintf (fp,  hdr_fmt[5], dparam("DELTAT")                       );
  fprintf (fp,  hdr_fmt[6], dparam("KINVIS")                       );
  fprintf (fp,  hdr_fmt[7], dparam("BETA"  )                       );
  fprintf (fp,  hdr_fmt[8], D->fTag                                );
  fprintf (fp,  hdr_fmt[9], (option("BINARY")) ? "binary" : "ASCII");

  if (option ("BINARY")) {
    register int  n;
    for (n = 0; n < nfields; n++)
      if (ntot != fwrite (*(D -> u)[n] -> value, sizeof(double), ntot, fp))
	message (routine, "unable to write to field file", ERROR);

  } else {
    register int  i, n;
    for (i = 0; i < ntot; i++) {
      for (n = 0; n < nfields; n++)
        if (fprintf (fp, "%#16g ", (D -> u)[n] -> value[0][i]) < 0)
	  message (routine, "unable to write to field file", ERROR);
      fputc ('\n', fp);
    }
  }

  fflush (fp);
}





void restart (Domain *D)
/* ========================================================================= *
 * If a restart file (name "D -> name".rst) can be found, use it to load     *
 * velocity fields.  Otherwise initialize velocity fields to zero.           *
 * ========================================================================= */
{
  char  routine[] = "restart";
  char  buf[STR_MAX];

  char  restartfile[FILENAME_MAX];
  FILE *fp;

  int   i, np, ns, nz, nel, ntot;


  strcat (strcpy (restartfile, D -> name), ".rst");

  if (fp = fopen (restartfile, "r")) {

    for (i = 0; i < 3; i++) readLine (fp);

    sscanf (buf, "%d %d %d %d", &np, &ns, &nz, &nel);
    if ((np != D -> u[0] -> np) || (ns != D -> u[0] -> np))
      message (routine, "element size mismatch", ERROR);
    if (nz != 1)
      message (routine, "number of z planes mismatch", ERROR);
    if (nel != D -> nEl)
      message (routine, "number of elements mismatch", ERROR);
    ntot = np * np * nz * nel;

    for (i = 3; i < 10; i++) readLine (fp);

    if (strstr (buf, "binary")) {
      for (i = 0; i < DIM; i++)
	if (ntot != fread (*(D -> u)[i] -> value, sizeof(double), ntot, fp))
	  message (routine, "unable to read field from file", ERROR);

    } else if (strstr (buf, "ASCII")) {
      register int  j, n;

      for (j = 0; j < ntot; j++) {
	if (!(fgets (buf, STR_MAX, fp)))
	  message (routine, "premature EOF", ERROR);
	for (n = 0; n < DIM; n++)
	  if (sscanf (buf, "%lf", (D -> u)[n] -> value[0] + j) < 1)
	    message (routine, "unable to read field from file", ERROR);      
      }
    } else
      message (routine, "can't find file format specifier", ERROR);

  } else {
    ntot = SQR (D -> u[0] -> np) * iparam ("NEL") * iparam ("NZ");
    for (i = 0; i < DIM; i++) Veclib::zero (ntot, *D -> u[i] -> value, 1);
  }
}





void  setUForce (Domain *D, Element **Us[DIM])
/* ========================================================================= *
 * On entry, intermediate velocity storage u^^ is in lowest levels of Us.    *
 * Multiply by -1.0 / (DELTAT * KINVIS) to create forcing for viscous step.  *
 * ========================================================================= */
{
  int     i;
  int     ntot  = D -> nEl * SQR (D -> u[0] -> np);
  double  alpha = -1.0 / (dparam ("DELTAT") * dparam ("KINVIS"));


  for (i = 0; i < DIM; i++) Blas::scal (ntot, alpha, *Us[i][0] -> value, 1);
}
