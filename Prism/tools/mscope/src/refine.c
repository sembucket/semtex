/*
 * refine -- Refine a mesh
 *
 * Once you've defined some error indicator on a mesh, the next step in
 * adaptive refinement is to generate a new mesh that "equalizes" the 
 * error.  
 *
 * This is a simple program to perform the adaptive refinement in a manual
 * way.  It reads an input file and mesh, asks you for a field to examine 
 * (the error indicator), and a tolerance value.  It then refines each 
 * element where the error indicator is higher than the tolerance to 
 * create a new mesh.  
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <assert.h>
#include <cubit.h>
#include "mscope.h"

/* Prototypes */

int main (int argc, char *argv[]);
void parse_args (int argc, char *argv[]);
void load_error (char* fname);

/* Private Data */

static char *prog  = "refine";
static char *usage = "session field";
static char *x11   = "x11 -n Refine -t Refine";
static char *dev;

/* Shared Data */

extern Session Geometry;

/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[])
{
  int  nref;
  char buf[BUFSIZ];

  parse_args (argc, argv);

  /* Load the mesh */

  cubit_init();
  sprintf(buf, "load %s", argv[1]);
  DoParse(buf);

  /* Open the graphics device and display the mesh */

  sprintf  (buf, "dev %s", dev);
  DoParse  (buf);
  DoErase  ();
  DoGrid   ();

  user_refine = mscope_refine;
  user_bc     = mscope_bc;
  user_prune  = NULL;
  user_perm   = NULL;

  /* Load the error indicator to operate on */
  
  load_error (argv[2]);

  /* Start the mesh refinement loop */

  do {
    double tol, refine_frac = 1., hmin = 0.01;

    Error_init (Geometry.U);
    
    printf ("Enter tolerance (<0 to stop): ");
    scanf  ("%lf", &tol);

    if (tol > 0.) {
      nref = Error_adapt (Geometry.U, tol, refine_frac, hmin);
      printf ("Refined %d elements\n", nref);
    } else
      nref = 0;

    DoErase();
    DoGrid ();

  } while (nref > 0);

  /* Save the mesh */

  sprintf (buf, "lock save refined.rea");
  DoParse (buf);

  return 0;
}

void parse_args (int argc, char *argv[])
{
  dev = x11;

  if (argc < 3) {
    fprintf (stderr, "%s: usage: %s\n", prog, usage);
    exit (-1);
  }

  return;
}

void load_error (char* fname)
{
  char  type;
  FILE* fld;
  FieldFile f;

  /* Open the field file */
  
  assert (fld = fopen (fname, "r"));
  memset (&f, '\0', sizeof(FieldFile));
  if (FieldFile_read (&f, fld) == EOF) {
    perror (prog);
    exit   (-1);
  }

  /* Select one field as an error indicator */

  printf ("Available fields: %s\n", f.type);
  printf ("Select one for your error indicator: ");
  scanf  ("%c", &type);

  Geometry.U->type  = type;
  Geometry.mesh->nz = 1;
  FieldFile_copy (&f, Geometry.U);

  fclose (fld);
  return;
}
