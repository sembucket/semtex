/*
 * Compute the minimum and maximum of each component in a field file
 */

#include <stdio.h>
#include <limits.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

static char *fld;

static char *prog  = "fminmax";
static char *usage = "[options] input[.fld]";
static char *help  = 
"options:\n"
"none\n"
"\n"
"This program computes the min and max of each component in a field file.\n";

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
  report_minmax (ff);
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

int report_minmax (FieldFile *f)
{
  const int npts  = 
    FIELDFILE_NR(f)*FIELDFILE_NS(f)*FIELDFILE_NZ(f)*FIELDFILE_NEL(f);
  const int nflds = FIELDFILE_COUNT(f);

  int i, n;

  printf ("# file = %s\n", f->name);
  printf ("# npts = %d\n", npts);
  printf ("# type   min   max\n");

  for (n = 0; n < nflds; n++) {
    char  type = FIELDFILE_TYPE(f,n);
    double *u  = FIELDFILE_DATA(f,n);
    double min = FLT_MAX;
    double max = -min;

    for (i = 0; i < npts;i++) {
      if (u[i] < min)
	min = u[i];
      if (u[i] > max)
	max = u[i];
    }
    printf ("%c %#14.7g %#14.7g\n", type, min, max);
  }
  return 0;
}
