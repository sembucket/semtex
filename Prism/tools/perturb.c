/*
 * perterb
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * This program creates a 3D initial condition by combining a 2D base flow
 * and a 2D or 3D perturbation field.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "speclib/speclib.h"
#include "veclib/veclib.h"

static char *base;
static char *perb;
static char *rea;

static struct {
  char *name;
  FILE *fp;
} out = { 
  NULL, 
  stdout 
};

static Field *U[2];
static Field *u[3];

static char *prog  = "perturb";
static char *usage = "[options] -r session[.rea] -b base[.fld] -p perb[.fld]";
static char *help  =
"options:\n"
"-eps #       specify the amplitude of the perturbation\n"
"-nz #        specify the number of spanwise points for 3D output\n"
"-o fname     specify the output file [default=stdout]\n"
"-v           be verbose\n"
"-h           print this help message\n"
"\n"
"This program creates a linear combination of a 2D base flow and a 2D or 3D\n"
"perturbation field.\n";

/* ------------------------------------------------------------------------- */

void parse_args (int argc, char *argv[])
{
  int n;

  for (n = 1; n < argc; n++) {
    if (*argv[n]=='-') {
      const char c = argv[n][1];
      switch(c) {
      case 'n':
	option_set("nz", atoi(argv[++n]));
	break;
      case 'o':
	out.name = strdup(argv[++n]);
	out.fp   = fopen(out.name,"w");
	break;
      case 'e':
	dparam_set("epsilon", atof(argv[++n]));
	break;
      case 'p':
	perb = strdup(argv[++n]);
	break;
      case 'b':
	base = strdup(argv[++n]);
	break;
      case 'r':
	rea = strdup(argv[++n]);
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

  if (rea == NULL) {
    fprintf(stderr, "%s: please use -r to specify the mesh\n", prog);
    exit(-1);
  }

  if (base == NULL) {
    fprintf(stderr, "%s: please use -b to specify the base flow\n", prog);
    exit(-1);
  }

  if (perb == NULL) {
    fprintf(stderr, "%s: please use -p to specify the perturbation\n", prog);
    exit(-1);
  }

}      

/* ------------------------------------------------------------------------- */

void get_mesh()
{
  char fname[FILENAME_MAX];
  FILE *fp;

  strcpy(fname,rea);

  if (!(fp = fopen(fname,"r"))) {
    strcat(fname,".rea");
    if (!(fp = fopen(fname,"r"))) {
      fprintf (stderr, "%s: unable to open %s or %s\n", prog, rea, fname);
      exit(-1);
    }
  }

  ReadParams(fp);
  U[0] = ReadMesh(fp);
  fclose(fp);

  if (option("verbose")) {
    fprintf (stderr, "%s: mesh size = [%d %d] x %d\n", prog,
	     FIELD_NR(U[0]), FIELD_NS(U[0]), FIELD_NELMT(U[0]));
  }
}

/* ------------------------------------------------------------------------- */

void get_base() 
{
  FieldFile *f = FieldFile_alloc();

  double norm[2];
  char fname[FILENAME_MAX];
  FILE *fp;
  int i;

  strcpy(fname,base);
  if (!(fp=fopen(fname,"r"))) {
    strcat(fname,".fld");
    if (!(fp=fopen(fname,"r"))) {
      fprintf (stderr,"%s: unable to open %s or %s\n", prog, base, fname);
      exit(-1);
    }
  }

  /* Allocate base flow field */

  U[1] = Field_dup(U[0]);
  for (i = 0; i < 2; i++)
    FIELD_TYPE(U[i]) = 'u'+i;

  FieldFile_read(f,fp);
  for (i = 0; i < 2; i++)
    FieldFile_get(f,U[i]);
  
  for (i = 0; i < 2; i++)
    norm[i] = Field_L2(U[i]);
  dparam_set("amp:base", dnrm2(2,norm,1));

  if (option("verbose")) 
    fprintf (stderr, "%s: base norm = %-14g [%g %g]\n", 
	     prog, dparam("amp:base"), norm[0], norm[1]);

  FieldFile_free(f);
  fclose(fp);
}

/* ------------------------------------------------------------------------- */

void get_perb() 
{
  FieldFile *f = FieldFile_alloc();

  double norm[3];
  char fname[FILENAME_MAX];
  FILE *fp;
  int i;

  strcpy(fname,perb);
  if (!(fp=fopen(fname,"r"))) {
    strcat(fname,".fld");
    if (!(fp=fopen(fname,"r"))) {
      fprintf (stderr,"%s: unable to open %s or %s\n", prog, base, fname);
      exit(-1);
    }
  }

  /* Allocate perturbation field */

  for (i = 0; i < 3; i++) {
    u[i] = Field_dup(U[0]);
    FIELD_TYPE(u[i]) = 'u'+i;
  }

  FieldFile_read   (f,fp);
  FieldFile_project(f,FIELD_NR(u[0]),FIELD_NS(u[0]));

  for (i = 0; i < 3; i++) {
    if (FieldFile_checkType(f,FIELD_TYPE(u[i])))
      FieldFile_get(f,u[i]);
    else 
      Field_scal (0., u[i]);
  }

  for (i = 0; i < 3; i++)
    norm[i] = Field_L2(u[i]);
  dparam_set("amp:perb", dnrm2(3,norm,1));

  if (option("verbose"))
    fprintf (stderr, "%s: perb norm = %-14g [%g %g %g]\n", 
	     prog, dparam("amp:perb"), norm[0], norm[1], norm[2]);
    
  FieldFile_free(f);
  fclose(fp);
}

/* ------------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  int nz, i, n, npts;
  double dz;
  double amp;
  FieldFile *f;

  speclib_init();

  dparam_set("epsilon", 0.01);

  parse_args (argc, argv);

  get_mesh();
  get_base();
  get_perb();

  if ((nz=option("nz"))==0)
    nz = 1;
  dz   = 2.*M_PI/nz;
  amp  = dparam("epsilon")*dparam("amp:base")/dparam("amp:perb");
  npts = FIELD_NR(U[0])*FIELD_NS(U[0])*FIELD_NELMT(U[0]);

  f = FieldFile_alloc();
  FieldFile_setSize(f, 3, FIELD_NR(U[0]),FIELD_NS(U[0]),nz,FIELD_NELMT(U[0]));
  for (n = 0; n < 3; n++)
    FIELDFILE_TYPE(f,n) = FIELD_TYPE(u[n]);

  if (option("verbose"))
    fprintf (stderr, "%s: creating 3D flow, nz=%d, amp=%g\n", prog, nz, amp);

  for (i = 0; i < npts; i++) {
    for (n = 0; n < nz; n++) {
      const double z = n*dz;

      FIELDFILE_DATA(f,0)[i+n*npts] = 
	amp*FIELD_FLAT(u[0],i)*cos(z) + FIELD_FLAT(U[0],i);
      
      FIELDFILE_DATA(f,1)[i+n*npts] = 
	amp*FIELD_FLAT(u[1],i)*cos(z) + FIELD_FLAT(U[1],i);
      
      FIELDFILE_DATA(f,2)[i+n*npts] = 
	amp*FIELD_FLAT(u[2],i)*sin(z);
    }
  }

  FieldFile_write(f,out.fp);
  FieldFile_free (f);

  return 0;
}
    
