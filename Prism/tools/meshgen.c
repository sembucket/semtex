/*
 * meshgen --- spectral element mesh generator
 *
 * This program reads a Prism input file and generates a mesh from it.
 * It then writes ASCII output into either a user-supplied file or 
 * stdout.
 *
 * NOTES
 * -----
 * This program assumes you are working with Prism and thus the input file
 * must follow Prism's conventions.  3-D meshes are thus limited to (x,y) +
 * some number of Fourier modes defined on an evenly spaced mesh.
 *
 * EXAMPLE
 * -------
 * To generate a mesh for "sample.rea" you would run the following:
 *
 *     % meshgen -o sample.mesh sample.rea 
 * 
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ----------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "mesh.h"

/* Externals */

static FILE *fp_in  = 0;
/* static FILE *fp_out = stdout; Changed hmb Jan 2002 */

static FILE *fp_out = 0;


static int   nr, ns, nz;

static char *prog   = "meshgen";
static char *rcsid  = "$Revision$";
static char *usage  = "usage: meshgen [options] input[.rea]\n";
static char *help   = 
"options:\n"
"-n ... generate a mesh with the specified N-order\n"
"-o ... write output to the named file (instead of stdout)\n"
"-z ... generate a mesh with the given # of z-coordinates\n";

static void error_msg(char *message)
{ fprintf(stderr, "%s: %s\n", prog, message); exit(-1); }


/* parse command line arguments */

void parse_args(int argc, char *argv[])
{
  char session [FILENAME_MAX];
  char fname   [FILENAME_MAX];
  char c, *p;

  while (--argc && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {
      case 'h':
	fputs (usage, stderr);
	fputs (help , stderr);
	exit  (0);
	break;
      case 'n':
	if (*++argv[0])
	  nr = ns = atoi (*argv);
	else {
	  nr = ns = atoi (*++argv);
	  argc--;
	}
	(*argv)[1] = '\0';
	break;
      case 'o':
	if (*++argv[0])
	  strcpy(fname, *argv);
	else {
	  strcpy(fname, *++argv);
	  argc--;
	}

	if ((fp_out=fopen(fname, "w")) == (FILE*) NULL)
	  error_msg ("unable to open an output file");
	(*argv)[1] = '\0';
	break;
      case 'v': {
	double v;
	sscanf(rcsid, "%*s%lf", &v);
	fprintf(stderr, "%s v%#4.2lf -- by Ron Henderson\n", prog, v);
	option_set("verbose", 2);
	break;
      }
      case 'z':
	if (*++argv[0])
	  nz = atoi (*argv);
	else {
	  nz = atoi (*++argv);
	  argc -= 1;
	}
	(*argv)[1] = '\0';
	break;

      default:
	fprintf(stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }

  if (argc != 1) {
    fputs(usage, stderr);
    exit (1);
  }

  if (p = strstr(*argv, ".rea"))
    strncpy(session, *argv, p - *argv);
  else
    strcpy (session, *argv);

  /* open the input file */

  if ((fp_in=fopen(*argv,"r")) == (FILE*) NULL) {
    sprintf(fname, "%s.rea", session);
    if ((fp_in=fopen(fname,"r")) == (FILE*) NULL) {
      fprintf(stderr, "%s: unable to open %s or %s for reading\n",
	      prog, session, fname);
      exit(1);
    }
  }
}

/* ----------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  Field *u;
  Mesh  *mesh;

  fp_out = stdout;

  speclib_init();

  parse_args (argc, argv);       
  ReadParams (fp_in);

  /* Modify the parameter BETA in case LZ was set */

  if (dparam("LZ")>0.) 
    dparam_set("BETA", scalar("2*PI/LZ"));

  /* Set the resolution parameters for the mesh */

  if (nr != 0) iparam_set("NR", nr); else nr = iparam("NR");
  if (ns != 0) iparam_set("NS", ns); else ns = iparam("NS");
  if (nz != 0) iparam_set("NZ", nz); else nz = iparam("NZ");

  u = ReadMesh (fp_in);
  mesh = Mesh_alloc (nr, ns, nz, Field_count(u));

  Mesh_define (mesh, u);
  Mesh_write  (mesh, fp_out);

  fclose (fp_in);
  fclose (fp_out);

  return 0;
}
