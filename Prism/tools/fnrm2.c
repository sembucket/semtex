/*
 * Compute the 2-norm of each component in a field file
 */

#include <stdio.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

static char *fld;

static char *prog  = "fnrm2";
static char *usage = "[options] input[.fld]";
static char *help  = 
"options:\n"
"none\n"
"\n"
"This program computes the 2-norm (Euclidean norm) of each component in a   \n"
"field file.\n";

main (int argc, char *argv[])
{
  FILE *fp;
  char fname[FILENAME_MAX];
  FieldFile *ff;

  speclib_init();

  /* parse the command line */

  parse_args(argc, argv);

  /* load the field file */

  strcpy (fname, fld);
  if (!(fp=fopen(fname,"r"))) {
    sprintf (fname, "%s.fld", fld);
    if (!(fp=fopen(fname,"r"))) {
      speclib_error ("can't open the field file");
    }
  }

  if (option("verbose")) fprintf (stderr, "%s: fld = %s\n", prog, fname);

  ff = FieldFile_alloc();

  FieldFile_read(ff,fp);
  report_nrm2   (ff);
  FieldFile_free(ff);

  fclose(fp);
  return 0;
}

int parse_args (int argc, char *argv[])
{
  int n;

  for (n = 1; n < argc; n++) {
    if (*argv[n]=='-') {
      const char c = argv[n][1];
      switch (c) {
      default:
	fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
    } else {
      fld = strdup(argv[n]);
    }
  }

  if (argc < 2) {
    fprintf (stderr, "usage: %s %s\n", prog, usage);
    exit(-1);
  }

  return 0;
}

int report_nrm2 (FieldFile *f)
{
  const int npts  = 
    FIELDFILE_NR(f)*FIELDFILE_NS(f)*FIELDFILE_NZ(f)*FIELDFILE_NEL(f);
  const int nflds = FIELDFILE_COUNT(f);

  int i;

  printf ("# file = %s\n", f->name);
  printf ("# npts = %d\n", npts);
  
  for (i = 0; i < nflds; i++) {
    const char   type = FIELDFILE_TYPE(f,i);
    const double nrm2 = dnrm2(npts, FIELDFILE_DATA(f,i), 1);
      
    printf ("%c = %g\n", type, nrm2);
  }

  return 0;
}
