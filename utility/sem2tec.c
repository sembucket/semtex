/*****************************************************************************
 * sem2tec: convert a SEM field file to Tecplot format.
 *
 * Usage: sem2tec [-h] [-o output] [-m mesh] [-n #] input[.fld]
 *
 * Based on the original code by Ron Henderson.
 *****************************************************************************/

static char
RCSid[] = "$Id$";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <femdef.h>
#include <alplib.h>
#include <femlib.h>

#define   MAXFIELDS 10

static char usage[] = 
  "usage: sem2tec [options] session[.fld]\n"
  "options:\n"
  "-h       ... print this message\n"
  "-o file  ... write output to the named file instead of running preplot\n" 
  "-m file  ... read the mesh from the named file (instead of stdin)\n"
  "-n #     ... evaluate the solution on an evenly-spaced mesh with N x N\n"
  "             points.  If N = 0, then no interpolation is done, i.e., the\n"
  "             output mesh will be on a standard GLL-spectral element mesh\n";

static FILE *fp_fld = 0,          /* default input files */
            *fp_msh = stdin;

static char *tecfile;             /* output file name */

static int   nr, ns, nz, nel, nfields, preplot_it = 1, np = 1;

static char    type[MAXFIELDS];
static double *data[MAXFIELDS], *x, *y;

static void    write_tec   (FILE*);
static void    parse_args  (int, char**);
static void    read_mesh   (FILE*);
static int     read_data   (FILE*);
static int     interpolate (void);
static double* do_interp   (double**, double**, double**, double**, double*);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  char  fname[L_tmpnam];
  char  buf  [BUFSIZ];
  FILE  *fp, *fp_tec;
  
  /* open a temporary file for the ascii tecplot format */

  if ((fp = fopen (tmpnam (fname),"w+")) == (FILE*) NULL) {
    fprintf (stderr, "sem2tec: unable to open a temporary file\n");
    exit    (EXIT_FAILURE);
  }

  parse_args  (argc, argv);

  read_mesh   (fp_msh);
  read_data   (fp_fld);
  interpolate ();
  write_tec   (fp);

  if (preplot_it) {
    sprintf (buf, "preplot %s %s", fname, tecfile);
    system  (buf);
    remove  (fname);
  } else {
    rewind (fp);
    fp_tec = fopen (tecfile, "w");
    while  (fgets(buf, BUFSIZ, fp)) fputs(buf, fp_tec);
    fclose (fp_tec);
    fclose (fp);
  }

  return EXIT_SUCCESS;
}


static void write_tec (FILE *fp)
{
  int i, j, k, nrns;

  nrns = nr * ns;

  fprintf(fp, "VARIABLES = X Y ");
  for (i = 0; i < nfields; i++)
    fprintf(fp, "%c ", toupper(type[i]));
  fprintf(fp, "\n");

  for (k = 0; k < nel; k++) {
    fprintf(fp, "ZONE T=\"Element %d\", I=%d, J=%d, F=POINT\n", k+1, nr, ns);
    for (i = 0; i < nrns; i++) {
      fprintf(fp, "%#14.7g %#14.7g ", x[k*nrns + i], y[k*nrns + i]);
      for (j = 0; j < nfields; j++)
	fprintf(fp, "%#14.7g ", data[j][k*nrns + i]);
      fprintf(fp, "\n");
    }
  }

  for (i = 0; i < nfields; i++) free (data[i]);
  fflush(fp);

  return;
}


static void parse_args (int    argc,
			char** argv)
{
  char c, fname[FILENAME_MAX];

  while (--argc && (*++argv)[0] == '-') 
    while (c = *++argv[0])
      switch (c) {
      case 'h':
	fputs (usage, stdout);
	exit  (EXIT_SUCCESS);
	break;
      case 'm':
	sprintf(fname, "%s", *++argv);
	if ((fp_msh = fopen(fname, "r")) == (FILE*) NULL) {
	  fprintf(stderr, "sem2tec: unable to open the mesh file -- %s\n", 
		  fname);
	  exit(EXIT_FAILURE);
	}
	argv[0] += strlen(*argv)-1; argc--;
	break;
      case 'o':
	tecfile = (char*) malloc(strlen(*++argv)+1);
	strcpy(tecfile, *argv);
	argv[0] += strlen(*argv)-1; argc--;
	preplot_it = 0;
	break;
      case 'n':
	if (*++argv[0])
	  np = atoi(*argv);
	else {
	  np = atoi(*++argv);
	  argc--;
	}
	(*argv)[1] = '\0';
	break;
      default:
	fprintf(stderr, "sem2tec: unknown option -- %c\n", c);
	break;
      }

  if (argc != 1) {
    fputs(usage, stderr);
    exit(EXIT_FAILURE);
  }

  /* open the input file */

  if ((fp_fld = fopen(*argv, "r")) == (FILE*) NULL) {
    sprintf(fname, "%s.fld", *argv);
    if ((fp_fld = fopen(fname, "r")) == (FILE*) NULL) {
      fprintf(stderr, "sem2tec: unable to open %s or %s\n", *argv, fname);
      exit(EXIT_FAILURE);
    }
  }

  /* get the name of the ouput file (if not supplied) */

  if (!tecfile) {
    int len = strlen(*argv);
    if (strcmp(*argv + len-4, ".fld") == 0) {
      tecfile = (char*) malloc(len+1);
      strncpy(tecfile, *argv, len-4);
    } else {
      tecfile = (char*) malloc(len+5);
      strcpy (tecfile, *argv);
    }
    strcat (tecfile, ".plt");
  }

  return;
}


static void read_mesh (FILE *fp)
{
  int  n, i;
  char buf [BUFSIZ];

  fgets(buf, BUFSIZ, fp);
  if (sscanf(buf, "%d%d%d%d", &nr, &ns, &nz, &nel) != 4) {
    fputs("error while reading mesh\n", stderr);
    exit (EXIT_FAILURE);
  }

  n  = nr * ns * nel;
  x  = dvector(0, n);
  y  = dvector(0, n);

  for (i = 0; i < n; i++) {
    fgets(buf, BUFSIZ, fp);
    if (sscanf(buf, "%lf%lf", x+i, y+i) != 2) {
      fputs("error while reading mesh data\n", stderr);
      exit (EXIT_FAILURE);
    }
  }

  return;
}


static int read_data (FILE *fp)
{
  int  i, n, ntot, machine, swab;
  char buf[BUFSIZ], *c;
  
  /* read the header down to the field list */

  for (n = 0; n < 9; n++) fgets(buf, BUFSIZ, fp);

  /* read the list of fields */

  n = 0;
  c = buf;
  nfields = 0;
  while (n++ < 25 && nfields < MAXFIELDS) 
    if (isalpha(*c++)) type[nfields++] = *(c-1);

  if (nfields == MAXFIELDS) {
    fprintf(stderr, "sem2tec: a maximum of %d fields may be converted.\n", 
	    MAXFIELDS);
    exit(EXIT_FAILURE);
  }

  /* allocate memory */

  ntot = nr * ns * nel;
  for (n = 0; n < nfields; n++)
    data[n] = (double*) malloc(ntot * sizeof(double));

  /* check the format */

  c = fgets(buf, BUFSIZ, fp); 
  while (isspace(*c)) c++;

  switch(tolower(*c)) {                     /* ascii or binary read */
  case 'a':
    for (i = 0; i < ntot; i++)
      for (n = 0; n < nfields; n++)
	if (fscanf(fp, "%lf", data[n] + i) < 0) {
	  fputs("sem2tec: an error has occured while reading\n", stderr);
	  exit (EXIT_FAILURE);
	}
    break;

  case 'b':
    machine  = iformat();

    if (     strstr (buf, "little") && machine == 0) swab = 1;
    else if (strstr (buf, "big"   ) && machine == 1) swab = 1;
    else                                             swab = 0;

    for (n = 0; n < nfields; n++) {
      if (fread(data[n], sizeof(double), ntot, fp) != ntot) {
	fputs("sem2tec: an error has occured while reading\n", stderr);
	exit (EXIT_FAILURE);
      }
      if (swab) dbrev(ntot, data[n], 1, data[n], 1);
    }
    break;

  default:
    fprintf(stderr, "sem2tec: unknown format flag -- %c\n", c);
    exit(EXIT_FAILURE);
    break;
  }

  return 1;
}


static int interpolate (void)
/* ------------------------------------------------------------------------- *
 * Interpolate from the GLL mesh to an evenly-spaced mesh.
 * ------------------------------------------------------------------------- */
{
  int      nrns = nr * ns;
  double   **imr, **itmr, **ims, **itms, *mesh_x, *mesh_y;
  double   *zr, *zs, *zm, *p;
  register k;

  switch (np) {
  case 0:              /* interpolation turned off */
    return;
    break;
    
  case 1:              /* no size specified ... use (NR|NS) */
    np = MAX (nr,ns);
    break;

  default:             /* size specified on the command line */
    break;
  }
  
  imr    = dmatrix(0, np-1, 0, nr-1);
  itmr   = dmatrix(0, nr-1, 0, np-1);
  ims    = dmatrix(0, np-1, 0, ns-1);
  itms   = dmatrix(0, ns-1, 0, np-1);
  zm     = dvector(0, np-1);

  /* compute the GLL-mesh and the new mesh */

#if 0
  coef  (nr);
  getops(nr, &zr, 0, 0, 0);
  coef  (ns);
  getops(ns, &zs, 0, 0, 0);

  zm[0] = -1.;
  zm[1] =  2./(np-1);
  dramp (np, zm, zm+1, zm, 1);

  /* compute the interpolation matrices */

  igllm (imr, itmr, zr, zm, nr, np);
  igllm (ims, itms, zs, zm, ns, np);
#endif

  dMeshOps (GLL, STD, nr, np, 0, &imr, &itmr, 0, 0);
  dMeshOps (GLL, STD, ns, np, 0, &ims, &itms, 0, 0);

  /* interpolate the mesh */

  mesh_x = do_interp (imr, itmr, ims, itms, x);
  mesh_y = do_interp (imr, itmr, ims, itms, y);

  for (k = 0; k < nfields; k++) {
    p = do_interp (imr, itmr, ims, itms, data[k]);
    free (data[k]);
    data[k] = p;
  }

  free (x); x = mesh_x;
  free (y); y = mesh_y;

  nr = ns = np;

  return;
}


static double *do_interp (double** imr ,
			  double** itmr,
			  double** ims ,
			  double** itms,
			  double*  data)
/* ------------------------------------------------------------------------- *
 * Wrapper for tensor-product interpolation.
 * ------------------------------------------------------------------------- */
{
  register int k;
  const int    nrns = nr * ns;
  const int    ntot = np * np;
  double       *new  = dvector (0, ntot * nel - 1),
               *tmp  = dvector (0, np*nr),
               *p    = new;

  for (k = 0; k < nel; k++, data += nrns, p += ntot) {
    dmxm (*ims, np,  data, ns, tmp, nr);
    dmxm ( tmp, np, *itmr, nr, p  , np);
  }

  free (tmp);

  return new;
}
