/*
 * SINTERP - Spectral Interpolation Program
 *
 * The following program takes a PRISM field file and interpolates it to a
 * new mesh.  The interpolation is spectral in (x,y), but the two files must
 * have the same number of frames.  
 *
 * The mesh being interpolated TO does not need to be a "spectral element"
 * mesh.
 *
 * Use -d to specify a specific dump number, otherwise it will interpolate
 * every dump in the file.
 *
 * $Id$ 
 * ---------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "sem_utils.h"

/* Externals */

char *prog   = "sinterp";
char *author = "Ron Henderson";
char *rcsid  = "$Revision$";
char *usage  = "usage: sinterp [options] -r session[.rea] input[.fld]\n";
char *help   = 
  "-d #    ... read the specified dump number\n"
  "-m file ... read the mesh from the named file\n"
  "-r file ... read the named session file\n"
  "-o file ... write output to the named file\n"
  "\n"
  "If the input file is specified as '-' it will be read from <stdin>.\n"
  "The mesh is read from <stdin> or from a file specified with '-m'.\n"
  "This is the mesh to interpolate TO, whereas the session file spec-\n"
  "ifies the mesh/field you are interpolating FROM.\n";
  
/* Prototypes */

void       parse_args (int argc, char **argv, FileList *f);
FieldFile* Interpolate(Mesh *mesh, FieldFile *orig, Element *master);
Element*   Setup      (FileList *f);

/* ---------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  FileList f;
  Mesh *mesh;
  Element *master;
  int dump;
  int status;

  FieldFile *field, *image;

  fprintf (stderr, "warning: %s may be unstable!\n", prog);

  speclib_init();

  field = FieldFile_alloc();

  parse_args (argc = generic_args (argc, argv, &f), argv, &f);

  dump   = iparam    ("DUMP-req");
  mesh   = Mesh_read (f.mesh.fp);
  master = Setup     (&f);

  /* Unless -d was used to specify a dump, loop through the file and *
   * interpolate each one.  Otherwise, only interpolate dump -d #    */

  if (dump == UNSET) {
    while (FieldFile_read (field, f.in.fp) != FIELDFILE_EOF) 
      {
	FieldFile_write(image = Interpolate (mesh, field, master), f.out.fp);
	FieldFile_free (image);

	if (option("verbose"))
	  fprintf (stderr,"%s: finished t = %g\n",prog,FIELDFILE_TIME(field));
      }
  } else {
    while (FieldFile_read (field, f.in.fp) != FIELDFILE_EOF && --dump)
           FieldFile_free (field);

    if (dump)
      fprintf (stderr, "%s: sorry, only %d dumps available\n", prog,
	       iparam("DUMP-req")-dump);
    else {
      FieldFile_write (image = Interpolate (mesh, field, master), f.out.fp);
      FieldFile_free  (image);
    }
  }

  FieldFile_free (field);

  speclib_exit();
  return 0;
}

/* --------------------------------------------------------------------- *
 * Interpolate() -- Interpolate a field to a new mesh                    *
 *                                                                       *
 * This function takes an input field and a mesh, and interpolates the   *
 * input field onto the mesh.  Any points outside the original mesh are  *
 * set to zero.                                                          *
 * --------------------------------------------------------------------- */

FieldFile *Interpolate (Mesh *mesh, FieldFile *field, Field *master)
{
  const int nr      = field->nr;
  const int ns      = field->ns;
  const int nrns    = nr * ns;
  const int ntot    = nr * ns * field->nel;
  const int nfields = strlen (field->type);
  const int nmesh   = mesh->nr * mesh->ns * mesh->nel;
  const int verbose = option("verbose");

  Field *U[_MAX_FIELDS];
  int i, j, k, kel, m, n, p;
  Probe *loc = NULL;
  int bad = 0;

  FieldFile *image = FieldFile_alloc();

  /* Check the sizes */

  if (field->nz != mesh->nz) 
    error_msg ("the field and mesh must have the same # of frames");

  /* Initialize the Fields */

  for (n = 0; n < nfields; n++) {
    U[n] = Field_aux (master, nr, ns, mesh->nz, FIELDFILE_TYPE(field,n));
    FieldFile_get (field, U[n]);
    FIELDFILE_TYPE(image,n) = FIELDFILE_TYPE(field,n);
  }

  /* Initialize the image field */

  memcpy (image, field, sizeof(Field));

  image->nr   = mesh->nr;
  image->ns   = mesh->ns;
  image->nz   = mesh->nz;
  image->nel  = mesh->nel;

  for (n = 0; n < nfields; n++)
    image->data[n] = dvector (0, nmesh * mesh->nz);
  for (i = 0; i < nmesh && loc == (Probe*) NULL; i++)
    loc = Probe_alloc (master, PROBE_XP, mesh->x[i], mesh->y[i]);
  if (loc == (Probe*) NULL)
    error_msg ("no points found in the mesh!");

  for (k = 0; k < (*mesh).nel; k++) {
    for (i = 0; i < (*mesh).ns; i++) {
      for (j = 0; j < (*mesh).nr; j++) {
	const int ip = (k * mesh->ns + i) * mesh->nr + j;

	if (! Probe_move(loc, mesh->x[ip], mesh->y[ip])) {
	  for (n = 0; n < nfields; n++) {
	    for (m = 0; m < (*mesh).nz; m++) {
	      Frame_set_one (m, U[n]);
	      image->data [n][ip + m*nmesh] = Probe_eval(loc, U[n]);
	    }
	  }
	  
	  /* Not found -- tally the error points */

	} else {
	  fprintf (stderr, "%s: the point (%g,%g) is not in the mesh\n",
		   prog, mesh->x[ip], mesh->y[ip]);
	  fprintf (stderr, "\t|r| = %g\n", hypot(mesh->x[ip], mesh->y[ip]));

	  for (n = 0; n < nfields; n++)
	    for (m = 0; m < (*mesh).nz; m++)
	      image->data[n][ip + m*nmesh] = 0.;
	  bad++;
	  continue;
	}
	if (verbose > 1)
	  fputc ('.', stderr);
      }
    }
    if (verbose > 0) 
      fputc     ('*', stderr);
  }
  

  if (verbose) {
    fputc   ('\n', stderr);
    fprintf (stderr, "Failed on %d of %d points [%3.0f%%]\n", bad, nmesh,
	     (float) bad / nmesh * 100.);
  }

  /* Clean up */

  for (n = 0; n < nfields; n++)
    Field_free (U[n]);

  return image;
}

/* --------------------------------------------------------------------- *
 * Setup() -- Initialize the master element array                        *
 *                                                                       *
 * Returns: pointer to an array of Elements                              *
 * --------------------------------------------------------------------- */

Element *Setup (FileList *f)
{
  Element *U;

  ReadParams(f->rea.fp);          /* First read the parameter list */

  if (iparam    ("NFRAMES") < iparam("NZ"))
      iparam_set("NFRAMES",   iparam("NZ"));

  U  = ReadMesh(f->rea.fp);        /* Generate the list of elements */
  
  return U;
}

/* --------------------------------------------------------------------- *
 * parse_args() -- Parse application arguments                           *
 *                                                                       *
 * This program only supports the generic utility arguments.             *
 * --------------------------------------------------------------------- */

void parse_args (int argc, char *argv[], FileList *f)
{
  char  c;
  int   i;
  char  fname[FILENAME_MAX];

  if (argc == 0) {
    fputs (usage, stderr);
    exit  (1);
  }

  while (--argc && (*argv)[0] == '-') {
    while (c = *++argv[0])                  /* more to parse... */
      switch (c) {
      default:
	fprintf(stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
    argv++;
  }
  
  /* check the rea file */

  if (!f->rea.fp)
    error_msg ("please use -r to specify the session (.rea) file");

  /* open the input file */

  if ((*argv)[0] == '-') {
    f->in.fp = stdin;
  } else {
    strcpy (fname, *argv);
    if ((f->in.fp = fopen(fname, "r")) == (FILE*) NULL) {
      sprintf(fname, "%s.fld", *argv);
      if ((f->in.fp = fopen(fname, "r")) == (FILE*) NULL) {
	fprintf(stderr, "%s: unable to open the input file -- %s or %s\n",
		prog, *argv, fname);
	exit(1);
      }
    }
    f->in.name = strdup(fname);
  }
  
  if (f->mesh.fp == stdin && f->in.fp == stdin) 
    error_msg ("either the mesh OR the input file can be read from <stdin>.");

  if (option("verbose")) 
    fprintf (stderr, "%s: in = %s, rea = %s, mesh = %s, out = %s\n", prog,
	     f->in.name   ? f->in.name   : "<stdin>",  f->rea.name,
	     f->mesh.name ? f->mesh.name : "<stdin>", 
	     f->out.name  ? f->out.name  : "<stdout>");

  return;
}
