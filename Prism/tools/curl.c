/*
 * curl -- compute the curl of a vector field
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>

#include "speclib/speclib.h"

static char *prog  = "curl";
static char *usage = "[options] -r session input > output";
static char *help  =
"options:\n"
"-h          print this help message\n"
"-v          print diagnostics\n"
"-o fname    send output to a named file\n"
"-2d         force calculation of vorticity for a 2d field\n"
"\n"
"This program computes the curl of a vector.  If the input is a 2D field   \n"
"(type = uv) the output is a scalar (type = t).  If the input is a 3D field\n"
"(type = uvw) the output is a vector (type = rst).  In either case, the    \n"
"vorticity is appended to the input data, and everything is written back   \n"
"out.";

static int verbose = 0;

static char *session = NULL;
/* static FILE *output  = stdout; Changed hmb Jan 2002 */
static FILE *output  = NULL;
static FILE *rea     = NULL;
static char *progid  = "$Id$";

static BSystem *matrix;

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
      case 'v': {
 	char vers[16], date[16];
	sscanf(progid,"%*s %*s %s %s", vers, date);
	fprintf (stderr, "%s: version %s built %s\n", prog, vers, date);
	verbose = 1;
	break;
      }
      case 'h':
	fprintf (stderr, "usage: %s %s\n%s\n", prog, usage, help);
	exit(0);
	break;
      case 'r': {
	char fname[FILENAME_MAX];
	if (!(rea = fopen(strcpy(fname,argv[++n]),"r"))) {
	  if (!(rea = fopen(strcat(fname,".rea"),"r"))) {
	    fprintf (stderr, "%s: unable to open %s or %s\n", prog, 
		     argv[n], fname);
	    exit(-1);
	  }
	}
	session = strdup(fname);
	break;
      }
      case 'o':
	if (!(output = fopen(argv[++n],"w"))) {
	  perror("curl");
	  exit (-1);
	}
	break;
      case '2':
	option_set("force2d",1);
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


/* Compute full 3D (x,y)-derivatives */

static void Field_DX (const Field *u, Field *du) 
{
  const int nz = FIELD_NZ(u);
  int k;

  for (k = 0; k < nz; k++) {
    Field_setFrameMulti (k, 2, u, du);
    Field_dx(u, du);
    Field_davg (du, matrix);
  }
}

static void Field_DY (const Field *u, Field *du) 
{
  const int nz = FIELD_NZ(u);
  int k;

  for (k = 0; k < nz; k++) {
    Field_setFrameMulti (k, 2, u, du);
    Field_dy(u, du);
    Field_davg (du, matrix);
  }
}

/* Compute the z-derivative of a Field using the FFT */

static void Field_DZ (const Field *u, Field *du) 
{
  const double beta = dparam("BETA");
  const int    n    = Field_frameSize (u);
  const int    m    = Field_frameCount(u);

  int i, j;

  Field_copy (u, du);
  Field_FFT (du, -1);

  /* Compute: \partial_z { u(x,y) exp[i q z] } := i q u(x,y) exp[i q z]     *
   *                                                                        *
   * The first 2*n data points are the real and complex components of mode  *
   * zero, so they are simply set to zero.                                  */

  for (i = 0; i < 2*n; i++)
    FIELD_FLAT(du,i) = 0.;

  /* The remaining components are i q (a + i b) = -(q b) + i (q a).  Recall *
   * that j = even are the real parts and j = odd are the imaginary parts.  */
  
  for (j = 2; j < m; j += 2) {
    const double q = (j/2.) * beta;
    for (i = 0; i < n; i++) {
      const double real = FIELD_FRAME(du,n,i,j);
      const double imag = FIELD_FRAME(du,n,i,j+1);

      FIELD_FRAME(du,n,i,j)   = -q * imag;
      FIELD_FRAME(du,n,i,j+1) =  q * real;
    }
  }

  /* transform back */

  Field_FFT (du, +1);
}

/* ------------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  FILE *fp = stdin;
  FieldFile *f;

  Field *u[3];  /* velocity vector */
  Field *q[3];  /* vorticity vector */
  Bedge *bc;
  
  int n;

  output = stdout;

  speclib_init();

  parse_args (argc, argv);

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

  /* Set parameters for z-direction, if applicable */

  if (dparam("LZ") > 0.)
    dparam_set ("BETA", scalar("2*PI/LZ"));
  else if (dparam("BETA") > 0.)
    dparam_set ("LZ", scalar("2*PI/BETA"));

  /* Allocate fields for the velocity and vorticity */

  u[0] = ReadMesh(rea);    FIELD_TYPE(u[0]) = 'u';
  u[1] = Field_dup(u[0]);  FIELD_TYPE(u[1]) = 'v';
  u[2] = Field_dup(u[0]);  FIELD_TYPE(u[2]) = 'w';
  q[0] = Field_dup(u[0]);  FIELD_TYPE(q[0]) = 'r';
  q[1] = Field_dup(u[0]);  FIELD_TYPE(q[1]) = 's';
  q[2] = Field_dup(u[0]);  FIELD_TYPE(q[2]) = 't';

  bc     = ReadBCs (rea, 0, u[0]);
  matrix = Matrix_alloc(u[0],bc,session);

  /* Load each field and compute the curl */

  f = FieldFile_alloc();

  while (FieldFile_read(f,fp) != FIELDFILE_EOF) {
    if (!FieldFile_checkType(f,'w') || option("force2d")) {
      for (n = 0; n < 2; n++)
	FieldFile_get(f, u[n]);

      Field_DX   (u[1], q[2]);
      Field_DY   (u[0], q[1]);
      Field_axpy (-1.,  q[1], q[2]);
      FieldFile_put (f, q[2]);

    } else {
      Field *tmp = Field_dup(u[0]);

      for (n = 0; n < 3; n++) 
	FieldFile_get(f, u[n]);

      Field_DY   (u[2], q[0]);
      Field_DZ   (u[1], tmp);
      Field_axpy (-1.0, tmp, q[0]);

      Field_DZ   (u[0], q[1]);
      Field_DX   (u[2], tmp);
      Field_axpy (-1.0, tmp, q[1]);

      Field_DX   (u[1], q[2]);
      Field_DY   (u[0], tmp);
      Field_axpy (-1.0, tmp, q[2]);

      for (n = 0; n < 3; n++)
	FieldFile_put(f, q[n]);

      Field_free (tmp);
    }

    FieldFile_write (f, output);
  }

  for (n = 0; n < 3; n++) {
    Field_free (u[n]);
    Field_free (q[n]);
  }

  FieldFile_free (f);
  speclib_exit();
  return 0;
}
