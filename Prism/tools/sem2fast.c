/* ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include "speclib/speclib.h"
#include "veclib/veclib.h"

/* Macros to access the FAST grid */

#define FAST_NPTS \
       (FAST_grid.nx * FAST_grid.ny * FAST_grid.nz)
#define FAST_XPNT(i,j,k) \
       (FAST_grid.x[FAST_grid.nx*((k)*FAST_grid.ny + (j)) + (i)])
#define FAST_YPNT(i,j,k) \
       (FAST_grid.y[FAST_grid.nx*((k)*FAST_grid.ny + (j)) + (i)])
#define FAST_ZPNT(i,j,k) \
       (FAST_grid.z[FAST_grid.nx*((k)*FAST_grid.ny + (j)) + (i)])

static struct {
  int nx;
  int ny;
  int nz;

  double *x;
  double *y;
  double *z;
} FAST_grid;

#define FAST_DATA(i,j,k,n) \
       (FAST_data.valu[n][FAST_grid.nx*((k)*FAST_grid.ny + (j)) + (i)])

static struct {
  int     count;
  char    type[_MAX_FIELDS];
  double *valu[_MAX_FIELDS];
} FAST_data;

static struct {
  char *name;
  FILE *fp;
} out;

static char *session;
static char *mesh;

static char *prog  = "sem2fast";
static char *usage = "[options] -m fastmesh -r session[.rea] input[.fld]";
static char *help  =
"-v                ... be verbose\n"
"-o file           ... send output to a named file\n"
"\n"
"This program interpolates a Prism field file onto a FAST grid.\n";

/* ------------------------------------------------------------------------- */

int FAST_load (FILE *fp)
{
  int n;
  int i;

  fscanf (fp, "%d", &FAST_grid.nx);
  fscanf (fp, "%d", &FAST_grid.ny);
  fscanf (fp, "%d", &FAST_grid.nz);

  n = FAST_grid.nx * FAST_grid.ny * FAST_grid.nz;

  FAST_grid.x = (double*) malloc(n*sizeof(double));
  FAST_grid.y = (double*) malloc(n*sizeof(double));
  FAST_grid.z = (double*) malloc(n*sizeof(double));
  
  for (i = 0; i < n; i++)
    fscanf (fp, "%lf", FAST_grid.x + i);
  for (i = 0; i < n; i++)
    fscanf (fp, "%lf", FAST_grid.y + i);
  for (i = 0; i < n; i++)
    fscanf (fp, "%lf", FAST_grid.z + i);

  return 0;
}

int FAST_interp (FieldFile *f, Field *template)
{
  int i, j, k, n;

  const int verbose = option("verbose");
  const int nflds   = FIELDFILE_COUNT(f);

  Field *u[_MAX_FIELDS];

  Probe *p = Probe_alloc
    (template, PROBE_XP, FAST_XPNT(0,0,0), FAST_YPNT(0,0,0));

  if (FAST_grid.nz != f->nz) {
    fprintf (stderr, 
	     "%s: z-dimension must match in FAST grid [%d] and field file "
	     "[%d]\n", prog, FAST_grid.nz, f->nz);
    exit(-1);
  }

  for (n = 0; n < nflds; n++) {
    u[n] = Field_dup(template);
    FIELD_TYPE(u[n]) = FIELDFILE_TYPE(f,n);
    FieldFile_get(f,u[n]);
  }

  FAST_data.count = nflds;
  for (n = 0; n < nflds; n++) {
    FAST_data.type[n] = FIELDFILE_TYPE(f,n);
    FAST_data.valu[n] = (double*) calloc(FAST_NPTS, sizeof(double));
  }

  for (j = 0; j < FAST_grid.ny; j++) {
    for (i = 0; i < FAST_grid.nx; i++) {
      const double x = FAST_XPNT(i,j,0);
      const double y = FAST_YPNT(i,j,0);

      if (Probe_move(p, x, y)) {
	for (n = 0; n < nflds; n++)
	  for (k = 0; k < FAST_grid.nz; k++) {
	    Field_setFrame(u[n], k);
	    FAST_DATA(i,j,k,n) = 0.;
	  }
	if (verbose) fputc('x',stderr);
      } else {
	for (n = 0; n < nflds; n++)
	  for (k = 0; k < FAST_grid.nz; k++) {
	    Field_setFrame(u[n], k);
	    FAST_DATA(i,j,k,n) = Probe_eval(p,u[n]);
	  }
	if (verbose) fputc('.',stderr);
      }
    }
  }

  if (verbose) fputc('\n',stderr);

  Probe_free(p);

  for (n = 0; n < nflds; n++)
    Field_free(u[n]);

  return 0;
}

int FAST_save (void)
{
  return 0;
}

/* ------------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  Field *u;
  FieldFile *ff;
  FILE *fp;
  char fname[FILENAME_MAX];

  speclib_init();

  parse_args(argc, argv);

  /* load the field file */

  strcpy(fname,argv[argc-1]);
  if (!(fp=fopen(fname,"r"))) {
    sprintf (fname, "%s.fld", argv[argc-1]);
    if (!(fp=fopen(fname,"r"))) {
      fprintf (stderr, "%s: can't open the input file -- %s or %s\n", 
	       prog, argv[argc-1], fname);
      exit(-1);
    }
  }

  ff = FieldFile_alloc();
  FieldFile_read(ff,fp);
  option_set("norder", ff->nr);

  /* Load the mesh */

  strcpy(fname,session);
  if (!(fp=fopen(fname,"r"))) {
    sprintf (fname, "%s.rea", session);
    if (!(fp=fopen(fname,"r"))) {
      fprintf (stderr, "%s: can't open the .rea file -- %s or %s\n", 
	       prog, session, fname);
      exit(-1);
    }
  }

  ReadParams(fp);
  u = ReadMesh(fp);

  /* Load the FAST grid */

  if (!(fp=fopen(mesh,"r"))) {
    fprintf (stderr, "%s: can't open the FAST grid -- %s\n", prog, mesh);
    exit(-1);
  }

  FAST_load  (fp);
  FAST_interp(ff,u);
  FAST_save  ();

  speclib_exit();

  return 0;
}

/* ------------------------------------------------------------------------- */

int parse_args (int argc, char *argv[])
{
  int n;

  for (n = 1; n < argc; n++) {
    if (*argv[n]=='-') {
      const char c = argv[n][1];

      switch (c) {
      case 'o':
	out.name = strdup(argv[++n]);
	break;
      case 'm':
	mesh = strdup(argv[++n]);
	break;
      case 'r':
	session = strdup(argv[++n]);
	break;
      case 'v':
	option_set("verbose", 1);
	break;
      case 'h':
	fprintf (stderr, "usage: %s %s\n%s", prog, usage, help);
	exit(0);
	break;
      default:
	fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
    }
  }

  if (argc < 6) {
    fprintf (stderr, "usage: %s %s\n", prog, usage);
    exit(-1);
  }

  if (!session) {
    fprintf (stderr, "%s: please use -r to specify a .rea file\n", prog);
    exit(-1);
  }

  if (!mesh) {
    fprintf (stderr, "%s: please use -m to specify a FAST grid\n", prog);
    exit(-1);
  }

  return 0;
}
