/*
 * xplane
 *
 * This program extracts cutting planes from a field file.  See further 
 * comments below.
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "speclib/speclib.h"
#include "veclib/veclib.h"

#define NPTS 64    /* default number of points */

typedef enum {     /* flags for coordinate axes */
  None = 0,
  X = 'x',
  Y = 'y',
  Z = 'z'
} AXIS;

static struct {
  int    dir    ;  /* orthogonal direction */
  int    npts[2];  /* resolution along each axis */
  double orig[2];  /* origin of the cutting plane */
  double size[2];  /* size of the cutting plane */
} xplane;

static struct {
  char *name;
  FILE *fp;
} out;

static char *fld;
static char *rea;
static char *session;

static char *prog  = "xplane";
static char *usage = "[options] -r session[.rea] input[.fld]";
static char *help  =
"options:\n"
"-xy \"x0,y0,dx,dy\" ... xy-cutting plane\n"
"-xz \"x0,z0,dx,dz\" ... xz-cutting plane\n"
"-yz \"y0,z0,dy,dz\" ... yz-cutting plane\n"
"-orig #           ... origin of the cutting plane along the orthogonal axis\n"
"-nx #             ... resolution along the x-axis\n"
"-ny #             ... resolution along the y-axis\n"
"-tec              ... write TECPLOT-formatted output\n"
"-sm               ... write SM-formatted binary output\n"
"-v                ... be verbose\n"
"-o file           ... send output to a named file\n"
"\n"
"This program extracts information along a cutting plane from a 2D or 3D    \n"
"field file.  The cutting plane is a rectangular region aligned with any two\n"
"of the three coordinate directions.  It is defined by four numbers: coord- \n"
"inates of its origin and lengths of its sides.  The two-letter combination \n"
"in the command line specification (-xy,-xz,-yz) specifies the alignment of \n"
"the cutting plane.  The option -orig moves it along the third orthogonal   \n"
"direction.  You can also change the default resolution of the cutting plane\n"
"in the (x,y)-directions.\n"
"\n"
"If you do not specify a cutting plane on the command line, the default is  \n"
"is an xy-plane that covers the entire computational domain.\n"
"\n"
"Output is a simple text file, a TECPLOT-formated file, or a set of binary  \n"
"SM image files (each file holds one field).\n"
"\n"
"Examples:\n"
"xplane -r mygrid myfield > mycut.xy\n"
"xplane -nx 256 -ny 64 -xy \"0,0,5,1\" -r mygrid myfield > mycut.xy\n"
"xplane -xz \"0,4,0,4\" -orig 1.5 -r mygrid myfield > mycut.xz\n";

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

  strcpy(fname,fld);
  if (!(fp=fopen(fname,"r"))) {
    sprintf (fname, "%s.fld", fld);
    if (!(fp=fopen(fname,"r"))) {
      fprintf (stderr, "%s: can't open the input file -- %s or %s\n", 
	       prog, fld, fname);
      exit(-1);
    }
  }

  if (option("verbose")) fprintf (stderr, "%s: fld = %s\n", prog, fname);

  ff = FieldFile_alloc();
  FieldFile_read(ff,fp);
  option_set("norder", ff->nr);

  /* Load the mesh */

  strcpy(fname,rea);
  if (!(fp=fopen(fname,"r"))) {
    sprintf (fname, "%s.rea", rea);
    if (!(fp=fopen(fname,"r"))) {
      fprintf (stderr, "%s: can't open the .rea file -- %s or %s\n", 
	       prog, rea, fname);
      exit(-1);
    }
  }

  if (option("verbose")) fprintf (stderr, "%s: rea = %s\n", prog, fname);

  ReadParams(fp);
  u = ReadMesh(fp);

  /* Set all defaults.  At the moment xz and yz cutting planes must span  *
   * the computational domain in the z-direction.                         */

  switch (xplane.dir) {
  case 'x':
    xplane.npts[0] = option("ny") ? option("ny") : NPTS;
    xplane.npts[1] = ff->nz;
    xplane.orig[1] = 0.;
    xplane.size[1] = dparam("LZ");
    break;
  case 'y':
    xplane.npts[0] = option("nx") ? option("nx") : NPTS;
    xplane.npts[1] = ff->nz;
    xplane.orig[1] = 0.;
    xplane.size[1] = dparam("LZ");
    break;
  default:
    xplane.npts[0] = option("nx") ? option("nx") : NPTS;
    xplane.npts[1] = option("ny") ? option("ny") : NPTS;
    if (xplane.dir != 'z') {
      xplane.dir = 'z';
      xplane_bbox(u);
    }
    break;
  }

  if (!out.name) {
    out.name = strdup(session);
    out.fp   = stdout;
  } else if (!option("sm")) {
    out.fp   = fopen(out.name,"w");
  }

  compute (ff,u);
  speclib_exit();

  return 0;
}

/* ------------------------------------------------------------------------- */

static char *root (char *s) {
  char *p = s + strlen(s)-1;
  while (p > s && *p != '.') p--;
  if (p != s) *p = '\0';
  return s;
}

int parse_args (int argc, char *argv[])
{
  int n;

  for (n = 1; n < argc; n++) {
    if (*argv[n]=='-') {
      const char c = argv[n][1];

      if (strcmp(argv[n],"-orig")==0) {
	dparam_set("origin", atof(argv[++n]));
	continue;
      }

      switch (c) {
      case 'x':
      case 'y': {
	char *dir  = argv[n];
	char *spec = argv[++n];
	xplane_init(dir,spec);
	break;
      }
      case 'n': {
	const char dir = argv[n][2];
	if (dir=='x')
	  option_set("nx", atoi(argv[++n]));
	else if (dir=='y')
	  option_set("ny", atoi(argv[++n]));
	else
	  fprintf (stderr, "%s: invalid direction: %c\n", prog, dir);
	break;
      }
      case 'o':
	out.name = strdup(argv[++n]);
	break;
      case 'r':
	rea = strdup(argv[++n]);
	break;
      case 's':
	option_set("sm", 1);
	break;
      case 't':
	option_set("tecplot", 1);
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
    } else {
      fld     = strdup(argv[n]);
      session = root(strdup(fld));
    }
  }

  if (argc < 4) {
    fprintf (stderr, "usage: %s %s\n", prog, usage);
    exit(-1);
  }

  if (!rea) {
    rea = (char*) malloc(strlen(session)+5);
    strcpy(rea,session);
    strcat(rea,".rea");
  }

  return 0;
}

/* ------------------------------------------------------------------------- */

int xplane_init (char *dir, const char *spec)
{
  char buf[BUFSIZ];
  char *p;
  double info[4];
  int n;

  n = 0;
  p = strtok(strcpy(buf,spec), ",");
  do {
    info[n++] = atof(p);
  } 
  while (p = strtok(NULL,","));

  if (n != 4) {
    fprintf (stderr, "%s: invalid specification -- %s\n", prog, spec);
    return -1;
  }

  if (strcmp(dir,"-yz")==0)
    xplane.dir = X;
  else if (strcmp(dir,"-xz")==0)
    xplane.dir = Y;
  else if (strcmp(dir,"-xy")==0)
    xplane.dir = Z;
  else {
    fprintf (stderr, "%s: invalid direction -- %s\n", prog, dir);
    return -1;
  }

  xplane.orig[0] = info[0];
  xplane.orig[1] = info[1];
  xplane.size[0] = info[2];
  xplane.size[1] = info[3];

  return 0;
}

/* ------------------------------------------------------------------------- */

int xplane_bbox (Field *u)
{
  Element *elmt;
  const int npts = FIELD_NR(u) * FIELD_NS(u);
  int i;

  xplane.orig[0] = xplane.orig[1] =  FLT_MAX;
  xplane.size[0] = xplane.size[1] = -FLT_MAX;

  for (elmt = Field_head(u); elmt; elmt = elmt->next) {
    for (i = 0; i < npts; i++) {
      xplane.orig[0] = MIN(xplane.orig[0], (*elmt->xmesh)[i]);
      xplane.orig[1] = MIN(xplane.orig[1], (*elmt->ymesh)[i]);
      xplane.size[0] = MAX(xplane.size[0], (*elmt->xmesh)[i]);
      xplane.size[1] = MAX(xplane.size[1], (*elmt->ymesh)[i]);
    }
  }

  xplane.size[0] -= xplane.orig[0];
  xplane.size[1] -= xplane.orig[1];

  return 0;
}

/* ------------------------------------------------------------------------- */

int compute (FieldFile *f, Field *u)
{
  const int nflds = FIELDFILE_COUNT(f);
  Field **uu = (Field**) calloc(nflds,sizeof(Field*));

  int n;

  /* extract data from the field file */

  for (n = 0; n < nflds; n++) {
    uu[n] = Field_dup(u);
    FIELD_TYPE(uu[n]) = FIELDFILE_TYPE(f,n);
    FieldFile_get(f,uu[n]);
  }
  
  switch (xplane.dir) {
  case Z:
    compute_xy(n,uu);
    break;
  case Y:
    compute_xz(n,uu);
    break;
  case X:
    compute_yz(n,uu);
    break;
  default:
    fprintf (stderr, "%s: cutting plane not implemented\n", prog);
    break;
  }

  for (n = 0; n < nflds; n++)
    Field_free(uu[n]);

  return 0;
}
  
#define XPLANE_DATA(data,n,i,j) \
(data[((n)*xplane.npts[1] + (j))*xplane.npts[0] + (i)])

int compute_xy (int nflds, Field *u[])
{
  int i, j, n;

  const int verbose = option("verbose");
  const int tecplot = option("tecplot");

  const double d0 = xplane.size[0] / (xplane.npts[0]-1.);
  const double d1 = xplane.size[1] / (xplane.npts[1]-1.);

  double x0, x1;

  double *data = (double*) 
    calloc(nflds*xplane.npts[0]*xplane.npts[1], sizeof(double));

  Probe *p = Probe_alloc(u[0], PROBE_XP, xplane.orig[0], xplane.orig[1]);

  if (option("verbose")) {
    fprintf(stderr, "%s: xy<%g,%g,%g,%g>, npts = %d x %d\n",
	    prog,
	    xplane.orig[0], xplane.orig[1],
	    xplane.size[0], xplane.size[1],
	    xplane.npts[0], xplane.npts[1]);
  }

  /* extract the data */

  for (j = 0; j < xplane.npts[1]; j++) {
    x1 = xplane.orig[1] + d1*(double)j;
    for (i = 0; i < xplane.npts[0]; i++) {
      x0 = xplane.orig[0] + d0*(double)i;
      if (Probe_move(p,x0,x1)) {
	if (verbose) fputc('x', stderr);
	continue;
      } else {
	for (n = 0; n < nflds; n++)
	  XPLANE_DATA(data,n,i,j) = Probe_eval(p,u[n]);
	if (verbose) fputc('.', stderr);
      }
    }
  }

  if (verbose) fputc('\n', stderr);

  output (nflds, u, data);

  Probe_free(p);
  free(data);
  return 0;
}


int compute_xz (int nflds, Field *u[])
{
  int i, j, n;

  const int verbose = option("verbose");
  const int tecplot = option("tecplot");

  const double d0 = xplane.size[0] / (xplane.npts[0]-1.);
  const double d1 = dparam("LZ")   / (xplane.npts[1]-1.);
  const double y  = dparam("origin");

  double x0, x1;

  double *data = (double*) 
    calloc(nflds*xplane.npts[0]*xplane.npts[1], sizeof(double));

  Probe *p = Probe_alloc(u[0], PROBE_XP, xplane.orig[0], y);

  if (option("verbose")) {
    fprintf(stderr, "%s: xz<%g,%g,%g,%g>, orig = %g, npts = %d x %d\n",
	    prog,
	    xplane.orig[0], xplane.orig[1],
	    xplane.size[0], xplane.size[1], y,
	    xplane.npts[0], xplane.npts[1]);
  }

  /* extract the data */

  for (i = 0; i < xplane.npts[0]; i++) {
    x0 = xplane.orig[0] + d0*(double)i;
    if (Probe_move(p,x0,y)) {
      if (verbose) fputc('x', stderr);
      continue;
    } else {
      for (j = 0; j < xplane.npts[1]; j++) {
	for (n = 0; n < nflds; n++) {
	  Field_setFrame(u[n],j);
	  XPLANE_DATA(data,n,i,j) = Probe_eval(p,u[n]);
	}
	if (verbose) fputc('.', stderr);
      }
    }
  }

  if (verbose) fputc('\n', stderr);

  output (nflds, u, data);

  Probe_free(p);
  free(data);
  return 0;
}

int compute_yz (int nflds, Field *u[])
{
  int i, j, n;

  const int verbose = option("verbose");
  const int tecplot = option("tecplot");

  const double d0 = xplane.size[0] / (xplane.npts[0]-1.);
  const double d1 = xplane.size[1] / (xplane.npts[1]-1.);
  const double x  = dparam("origin");

  double x0, x1;

  double *data = (double*) 
    calloc(nflds*xplane.npts[0]*xplane.npts[1], sizeof(double));

  Probe *p = Probe_alloc(u[0], PROBE_XP, x, xplane.orig[0]);

  if (option("verbose")) {
    fprintf(stderr, "%s: yz<%g,%g,%g,%g>, orig = %g, npts = %d x %d\n",
	    prog,
	    xplane.orig[0], xplane.orig[1],
	    xplane.size[0], xplane.size[1], x,
	    xplane.npts[0], xplane.npts[1]);
  }

  /* extract the data */

  for (i = 0; i < xplane.npts[0]; i++) {
    x0 = xplane.orig[0] + d0*(double)i;
    if (Probe_move(p,x,x0)) {
      if (verbose) fputc('x', stderr);
      continue;
    } else {
      for (j = 0; j < xplane.npts[1]; j++) {
	for (n = 0; n < nflds; n++) {
	  Field_setFrame(u[n],j);
	  XPLANE_DATA(data,n,i,j) = Probe_eval(p,u[n]);
	}
	if (verbose) fputc('.', stderr);
      }
    }
  }

  if (verbose) fputc('\n', stderr);

  output (nflds, u, data);

  Probe_free(p);
  free(data);
  return 0;
}

/* ------------------------------------------------------------------------- */

int output (int nflds, Field *u[], double *data) 
{
  if (option("sm"))
    return output_sm(nflds, u, data);
  else if (option("tecplot"))
    return output_tecplot(nflds, u, data);
  else
    return output_xplane(nflds, u, data);
}

int output_xplane (int nflds, Field *u[], double *data)
{
  int i, j, n;

  const double d0 = xplane.size[0] / (xplane.npts[0]-1.);
  const double d1 = xplane.size[1] / (xplane.npts[1]-1.);

  double x0, x1;

  fputs ("# cutting plane extracted from ", out.fp);
  fputs (session, out.fp);
  fputs ("\n# variables = ", out.fp);
  switch (xplane.dir) {
  case 'x':
    fputs ("y z", out.fp);
    break;
  case 'y':
    fputs ("x z", out.fp);
    break;
  case 'z':
    fputs ("x y", out.fp);
    break;
  }
  
  for (n = 0; n < nflds; n++)
    fprintf (out.fp, " %c", FIELD_TYPE(u[n]));
  fprintf (out.fp, "\n# npts = %d x %d\n",
	   xplane.npts[0], xplane.npts[1]);
    
  for (j = 0; j < xplane.npts[1]; j++) {
    x1 = xplane.orig[1] + d1*(double)j;
    for (i = 0; i < xplane.npts[0]; i++) {
      x0 = xplane.orig[0] + d0*(double)i;
      fprintf (out.fp, "%g %g", x0, x1);
      for (n = 0; n < nflds; n++)
	fprintf (out.fp, " %g", XPLANE_DATA(data,n,i,j));
      fputc ('\n', out.fp);
    }
  }
  return 0;
}

int output_tecplot (int nflds, Field *u[], double *data)
{
  int i, j, n;

  const double d0 = xplane.size[0] / (xplane.npts[0]-1.);
  const double d1 = xplane.size[1] / (xplane.npts[1]-1.);

  double x0, x1;

  fputs ("VARIABLES = ", out.fp);
  switch (xplane.dir) {
  case 'x':
    fputs ("Y Z", out.fp);
    break;
  case 'y':
    fputs ("X Z", out.fp);
    break;
  case 'z':
    fputs ("X Y", out.fp);
    break;
  }
  
  for (n = 0; n < nflds; n++)
    fprintf (out.fp, " %c", toupper(FIELD_TYPE(u[n])));
  fprintf (out.fp, "\nZONE T=\"x-plane\", I=%d, J=%d, F=POINT\n",
	   xplane.npts[0], xplane.npts[1]);
    
  for (j = 0; j < xplane.npts[1]; j++) {
    x1 = xplane.orig[1] + d1*(double)j;
    for (i = 0; i < xplane.npts[0]; i++) {
      x0 = xplane.orig[0] + d0*(double)i;
      fprintf (out.fp, "%g %g", x0, x1);
      for (n = 0; n < nflds; n++)
	fprintf (out.fp, " %g", XPLANE_DATA(data,n,i,j));
      fputc ('\n', out.fp);
    }
  }
  return 0;
}

int output_sm (int nflds, Field *u[], double *data)
{
  int i, j, n;
  char fname[FILENAME_MAX];
  char *base = root(strdup(fld));

  if (option("verbose")) {
    fprintf (stderr, "%s: out = %s.[", prog, base);
    for (n = 0; n < nflds; n++)
      fprintf (stderr, "%c", FIELD_TYPE(u[n]));
    fprintf (stderr, "]\n");
  }

  for (n = 0; n < nflds; n++) {
    sprintf (fname, "%s.%c", base, FIELD_TYPE(u[n]));
    out.fp = fopen(fname,"w");

    fwrite (&xplane.npts[0],sizeof(int),1,out.fp);
    fwrite (&xplane.npts[1],sizeof(int),1,out.fp);

    for (j = 0; j < xplane.npts[1]; j++) {
      for (i = 0; i < xplane.npts[0]; i++) {
	float tmp = XPLANE_DATA(data,n,i,j);
	fwrite(&tmp,sizeof(float),1,out.fp);
      }
    }
    fclose(out.fp);
  }

  free(base);
  return 0;
}
