/*
 * Convert a Prism field file to TECPLOT format
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "sem_utils.h"

void parse_args (int argc, char *argv[], FileList *f);

void TECPLOT_write (FieldFile *field, Mesh *mesh, FileList *f);
void TECPLOT_name  (FileList *f, FieldFile *field, char *name);

/* The following strings MUST be defined */

char *prog   = "sem2tec";
char *author = "Ron Henderson";
char *rcsid  = "$Revision$";
char *usage  = "usage: sem2tec [options] input[.fld]\n";
char *help   =
  "-n #    ... interpolate output to an N x N mesh\n"
  "-d #    ... extract the given dump from the field file\n"
  "-gll    ... output the mesh on the GLL-points\n"
  "-m file ... read the mesh from the named file\n"
  "-o file ... send output to the named file\n"
  "\n"
  "If the input file is specified as '-' it will be read from stdin. In this\n"
  "case you must use the -m option to include a mesh.\n"
  "\n"
  "By default, sem2tec interpolates the field file to an evenly-spaced   \n"
  "mesh with N x N points to improve the contour plotting in TECPLOT.    \n"
  "If you specify -gll then the output will be on the standard GLL-mesh. \n"
  "\n";

/* ----------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  FileList f;
  int dump;

  Mesh* mesh;
  FieldFile *ff;

  parse_args (argc = generic_args(argc, argv, &f), argv, &f);

  ff   = FieldFile_alloc();
  mesh = Mesh_read (f.mesh.fp);                /* Read the mesh */
  dump = iparam("DUMP-req");

  while (dump--) {
    if (FieldFile_read (ff, f.in.fp) == FIELDFILE_EOF) 
      error_msg ("ran out of field dumps!");
  }
  
  TECPLOT_write (ff, mesh, &f);            /* Write the TECPLOT file  */

  return 0;
}

/* ------------------------------------------------------------------------ *
 * TECPLOT_write() -- write out the tecplot file                            *
 *                                                                          *
 * This function does one of two things.  If an output file has been opened *
 * already, it simply writes the ASCII tecplot file there.  Otherwise, it   *
 * writes the ASCII version to a temporary file, runs preplot on it, and    *
 * then copies the preplot file to the "standard" output file (not stdout). *
 * ------------------------------------------------------------------------ */

void TECPLOT_write (FieldFile *f, Mesh *mesh, FileList *flist)
{
  FILE *fp;
  int   nr, ns, nz, nel, nrns, n, nfields;
  char  buf[BUFSIZ];

  char  tmp[FILENAME_MAX];

  const int N   = iparam("NORDER-req");
  const int uni = option("uniform");

  /* Check for compatibility between the mesh and field */

  if (mesh->nz != f->nz)
    error_msg ("mesh and fieldfile have different z-dimensions");

  /* Interpolate the mesh and the field */

  if (N != UNSET || uni != 0) {

    int nr = (N != UNSET) ? N : f->nr;
    int ns = (N != UNSET) ? N : f->ns;

    if (uni) {
      const double dx = 2./(nr-1.);
      const double dy = 2./(ns-1.);

      double xx[_MAX_NORDER];
      double yy[_MAX_NORDER];
      int i;

      for (i = 0; i < nr; i++)
	xx[i] = -1. + i*dx;
      for (i = 0; i < ns; i++)
	yy[i] = -1. + i*dy;

      FieldFile_interp  (f, nr, xx, ns, yy);
    } else 
      FieldFile_project (f, nr, ns);

    Mesh_interp (mesh, nr, ns, 0, uni);
  }

  
  if (flist->out.fp)            /* use the named output file */
    fp = flist->out.fp;
  else 
    if ((fp=fopen(tmpnam(tmp),"w+")) == (FILE *) NULL)
      error_msg ("unable to open a temporary file");
   
  nr      = f->nr;
  ns      = f->ns;
  nz      = f->nz;
  nel     = f->nel;
  nfields = FieldFile_getTypeList(f, buf);
  nrns    = nr * ns;

  if (nz == 1)
    fputs ("VARIABLES = X Y ", fp);
  else
    fputs ("VARIABLES = X Y Z ", fp);
  for (n = 0; n < nfields; n++)
    fprintf(fp, "%c ", toupper(buf[n]));
  fputc ('\n', fp);


  /* .......... 2-D Output File .......... */

  if (nz == 1) {
    int i, k, m;
    for (k = 0; k < nel; k++) {
      fprintf(fp, "ZONE T=\"Element %d\", I=%d, J=%d, F=POINT\n", 
	      k+1, nr, ns);
      for (i = 0; i < nrns; i++) {
	fprintf (fp, "%#14.7g %#14.7g ", 
		 (mesh->x)[k*nrns + i], (mesh->y)[k*nrns+i]);
	for (n = 0; n < nfields; n++)
	  fprintf(fp, "%#14.7g ", (f->data[n])[k*nrns + i]);
	fputc ('\n', fp);
      }
    }
  }


  /* .......... 3-D Output File .......... */

  else {
    int i, k, m;
    for (k = 0; k < nel; k++) {
      fprintf(fp, "ZONE T=\"Element %d\", I=%d, J=%d, K=%d, F=POINT\n", 
	      k+1, nr, ns, nz);
      for (m = 0; m < nz; m++) {
	for (i = 0; i < nrns; i++) {
	  fprintf (fp, "%#14.7g %#14.7g %#14.7g", 
		   (mesh->x)[k*nrns + i], (mesh->y)[k*nrns+i], (mesh->z)[m]);
	  for (n = 0; n < nfields; n++)
	    fprintf(fp, "%#14.7g ", (f->data[n])[(m * nel + k)*nrns + i]);
	  fputc ('\n', fp);
	}
      }
    }
  }

  fflush(fp);

  if (!(flist->out.fp)) {         /* Now decide if we need to run preplot. */
    char tec[FILENAME_MAX];
    TECPLOT_name (flist, f, tec);
    sprintf (buf, "preplot %s %s", tmp, tec);
    system  (buf);
    remove  (tmp);
  }

  return;
}

void TECPLOT_name (FileList *f, FieldFile *field, char *name)
{
  char *p;

  if (f->in.name) {                /* Named input file.  Use the file name */
    strcpy (name, f->in.name);     /* with ".plt" substituted for ".fld"   */
    if (p = strstr(name,".fld"))   /* if present, appended if not.         */
      strcpy (p, ".plt");
    else
      strcat (name, ".plt");
  } else                           /* No input file name.  Use the session */
    sprintf (name, "%s.plt", field->name);   /* name with ".plt" appended. */

  return;
}

/* ------------------------------------------------------------------------ *
 * parse_args() -- Parse appplication arguments                             *
 *                                                                          *
 * Sem2tec supports the standard arguments plus a -d # option to read a     *
 * dump other than the first one from a field file.                         *
 *                                                                          *
 * This program never writes to stdout.  If the -v option is turned on, it  *
 * passes the stdout of preplot on to the screen.                           *
 * ------------------------------------------------------------------------ */

void parse_args(int argc, char *argv[], FileList *f)
{
  char c, fname[FILENAME_MAX];
 
  if (argc == 0) {
    fputs (usage, stderr);
    exit  (1);
  }

  option_set ("uniform", 1);  /* To do the interpolation to an even mesh */

  while (--argc && (*argv)[0] == '-') {
    if (!strcmp("-gll", *argv))
      { option_set ("uniform", 0); argv++; continue; }
    while (c = *++argv[0])
      switch (c) {
      default:
	fprintf(stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
    argv++;
  }
  
  /* open the input file */

  if ((*argv)[0] == '-') {
    if (strlen (*argv) > 1) 
      { fputs (usage, stderr); exit (1); }
    if (f->mesh.fp == stdin) 
      error_msg("please use -m to specify the mesh file");
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

  /* Check for errors in the input arguments */

  if (iparam("NORDER-req") != UNSET)
    if (iparam("NORDER-req") != 0 && iparam("NORDER-req") < 3)
      error_msg ("the # of points should be 0 or > 2");
  if (iparam("DUMP-req") == UNSET) 
    iparam_set("DUMP-req", 1);


  /* If the output file has been specified via -o it will have its name *
   * defined already.  Otherwise, it's just set to the default.         */

  if (!f->out.name) f->out.fp = (FILE *) NULL;

  if (option("verbose")) {
    fprintf (stderr, "%s: in = %s, mesh = %s, out = %s", prog,
	     f->in.name   ? f->in.name   : "<stdin>",
	     f->mesh.name ? f->mesh.name : "<stdin>",
	     f->out.name  ? f->out.name  : "<preplot>");
    if (iparam("NORDER-req") != UNSET)
      fprintf (stderr, " @ %d x %d", 
	       iparam("NORDER-req"), iparam("NORDER-req"));
    fputc('\n', stderr);
  } else 

    /* Unless the verbose option is on, redirect stdout to /dev/null *
     * to prevent preplot from echoing the zone data.  This program  *
     * does not use stdout for anything else.                        */

    { freopen ("/dev/null", "w", stdout); }
      
  return;
}
