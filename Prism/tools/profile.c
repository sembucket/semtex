/*
 * Extract profiles from a FieldFile
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

static char *prog  = "profile";
static char *usage = "[options] -r session -p \"profile\" input";
static char *help  =
"options:\n"
"-o file ... output file\n"
"\n"
"If the input file is specified as '-' it will be read from <stdin>. Pro-  \n"
"files specified with -p (there should be at least one!) have the ollowing \n"
"format:\n"
"    \"profile\" = \"[n:] x0 y0 dx dy\"\n"
"where n = number of points, (x0,y0) is the beggining of the line, and     \n"
"(dx,dy) is the length of the line.";

static int verbose = 0;

/* static FILE *output = stdout; Changed hmb 2002 */
static FILE *output = NULL;
static FILE *rea    = NULL;

#define PROFILE_MAX  128
#define PROFILE_NP   128

typedef struct {
  int    np;
  double x0, dx;
  double y0, dy;
} profile_t;

static int profile_count = 0;
profile_t  profile[PROFILE_MAX];

/* ------------------------------------------------------------------------- */

void profile_add (char *spec)
{
  int  i  = profile_count;
  int  np = PROFILE_NP;
  char *s = strchr(spec,':');

  if (i == PROFILE_MAX) {
    fprintf (stderr, "%s: too many profiles!\n", prog);
    return;
  }

  if (s) {
    np   = atoi(spec);
    spec = ++s;
  }

  profile[i].np = np;
  profile[i].x0 = atof(strtok(spec," "));
  profile[i].y0 = atof(strtok(NULL," "));
  profile[i].dx = atof(strtok(NULL," "));
  profile[i].dy = atof(strtok(NULL," "));

  profile_count++;
}

/* ------------------------------------------------------------------------- */

void profile_compute (int nflds, Field *u[])
{
  int i, n;

  fputs("# profiles\n#\n# Fields: ", output);
  for (n = 0; n < nflds; n++)
    fprintf (output, "%c ", FIELD_TYPE(u[n]));
  fputs("\n#\n", output);

  for (i = 0; i < profile_count; i++) {
    const int    np = profile[i].np;
    const double x0 = profile[i].x0;
    const double dx = profile[i].dx/(profile[i].np-1.);
    const double y0 = profile[i].y0;
    const double dy = profile[i].dy/(profile[i].np-1.);

    Probe *p = Probe_alloc(u[0], PROBE_XP, x0, y0);
    int j;

    fprintf (output, "# profile: %d points, orig=(%g,%g), path=(%g,%g)\n",
	     np, x0, y0, profile[i].dx, profile[i].dy);

    for (j = 0; j < np; j++) {
      const double x = x0 + j*dx;
      const double y = y0 + j*dy;

      Probe_move (p,x,y);
      fprintf(output, "%g %g ", x, y);

      for (n = 0; n < nflds; n++) 
	fprintf(output,"%g ", Probe_eval(p,u[n]));
      fputc('\n',output);
    }

    Probe_free(p);
  }
}

/* ------------------------------------------------------------------------- */

void parse_args (int argc, char *argv[])
{
  int n;

  if (argc == 1) {
    fprintf (stderr, "usage: %s %s\n", prog, usage);
    exit (0);
  }

  for (n = 0; n < argc; n++) {
    if (*argv[n] == '-' && strlen(argv[n]) > 1) {
      const char c = argv[n][1];
      switch (c) {
      case 'h':
	fprintf (stderr, "usage: %s %s\n%s\n", prog, usage, help);
	exit(0);
	break;
      case 'v':
	verbose = 1;
	break;
      case 'p':
	profile_add(argv[++n]);
	break;
      case 'r':
	if (!(rea = fopen(argv[++n],"r"))) {
	  char fname[FILENAME_MAX];
	  sprintf (fname, "%s.rea", argv[n]);
	  if (!(rea = fopen(fname,"r"))) {
	    fprintf (stderr, "%s: unable to open %s or %s\n", prog, 
		     argv[n], fname);
	    exit(-1);
	  }
	}
	break;
      case 'o':
	if (!(output = fopen(argv[++n],"w"))) {
	  perror(prog);
	  exit (-1);
	}
	break;
      default:
	fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
    }
  }

  if (rea == NULL) {
    fprintf (stderr, "%s: please use -r to specify the .rea file\n", prog);
    exit (-1);
  }
}

/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[])
{
  FILE  *fp = stdin;
  Field *u[FIELDFILE_MAX];
  FieldFile *f;

  int n, nflds;

  output = stdout;

  speclib_init();
  parse_args(argc, argv);

  if (*argv[argc-1] != '-') {
    if (!(fp = fopen(argv[argc-1],"r"))) {
      char fname[FILENAME_MAX];
      sprintf (fname, "%s.fld", argv[argc-1]);
      if (!(fp = fopen(fname,"r"))) {
	fprintf (stderr, "%s: unable to open %s or %s\n", prog, 
		 argv[argc-1], fname);
	exit(-1);
      }
    }
  }

  /* Load the mesh and allocate fields */

  ReadParams(rea);
  u[0] = ReadMesh(rea);

  /* Load each field */

  f = FieldFile_alloc();
  while (FieldFile_read(f,fp) != FIELDFILE_EOF)
    continue;

  nflds = FIELDFILE_COUNT(f);
  for (n = 0; n < nflds; n++) {
    if (n) u[n] = Field_dup(u[0]);
    FIELD_TYPE(u[n]) = FIELDFILE_TYPE(f,n);
    FieldFile_get(f,u[n]);
  }

  profile_compute(nflds, u);

  for (n = 0; n < nflds; n++)
    Field_free(u[n]);

  FieldFile_free(f);
  speclib_exit();
  return 0;
}
