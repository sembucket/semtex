/*****************************************************************************
 * repeatz.c: read a field file, and output the number of repetitions
 * in the z direction indicated on the command line.  Field must be
 * binary format.  By default, output the original file, no repeats.
 * Adust Nz and Beta in header as appropriate.
 *
 * Copyright (c) 2002 Hugh Blackburn.
 *
 * USAGE
 * -----
 * repeatz [-h] [-n <rep>] [input[.fld]
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <femlib.h>
#include <femdef.h>
#include <alplib.h>

static void getargs (int, char**, FILE**, int*);

static char prog[] = "repeatz";
static const char *hdr_fmt[] = { 
  "%-25s "              "Session\n",
  "%-25s "              "Created\n",
  "%-5d%-5d%-5d%-10d "  "Nr, Ns, Nz, Elements\n",
  "%-25d "              "Step\n",
  "%-25.6g "            "Time\n",
  "%-25.6g "            "Time step\n",
  "%-25.6g "            "Kinvis\n",
  "%-25.6g "            "Beta\n",
  "%-25s "              "Fields written\n",
  "%-25s "              "Format\n"
};


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Wrapper.
 * ------------------------------------------------------------------------- */
{
  char   buf[STR_MAX], fields[STR_MAX], fmt[STR_MAX];
  int    i, j, n, np, nz, nel, nrep = 1;
  int    nfields, nplane, npts, swab;
  FILE   *fp_in = stdin, *fp_out = stdout;
  double *data, beta;

  getargs (argc, argv, &fp_in, &nrep);
  format  (fmt);

  while (fgets (buf, STR_MAX, fp_in)) { 

    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    
    if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    fprintf (fp_out, hdr_fmt[2], np, np, nz*nrep, nel);

    n = 4;
    while (n--) { fgets (buf, STR_MAX, fp_in); fputs (buf, fp_out); }

    fgets (buf, STR_MAX, fp_in);
    sscanf (buf, "%lf", &beta);
    beta /= nrep;
    fprintf (fp_out, hdr_fmt[7], beta);   

    fgets (buf, STR_MAX, fp_in); fputs (buf, fp_out);
    for (nfields = 0, i = 0; i < 25; i++) if (isalpha(buf[i])) nfields++;

    fgets (buf, STR_MAX, fp_in); fputs (buf, fp_out);
    if (!strstr(buf, "binary"))
      message (prog, "input file not binary format", ERROR);
    swab = (strstr (buf, "big") && strstr (fmt, "little")) || 
           (strstr (fmt, "big") && strstr (buf, "little"));

    /* -- Set sizes, allocate storage. */

    nplane = np * np * nel;
    npts   = nz * nplane;
    data   = dvector (0, npts - 1);
    
    /* -- Read and write all data fields. */

    for (i = 0; i < nfields; i++) {
      if (fread (data, sizeof (double), npts, fp_in) != npts)
	message (prog, "an error occured while reading", ERROR);
      if (swab) dbrev (npts, data, 1, data, 1);
      for (j = 0; j < nrep; j++) {
	if (fwrite (data, sizeof (double), npts, fp_out) != npts)
	  message (prog, "an error occured while writing", ERROR);
      }
    }
    
    freeDvector (data, 0);
  } 
  
  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     FILE** fp_in,
		     int*   nrep )
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] = "repeatz [-h] [-n <rep>] [input[.fld]\n";

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0]) *nrep = atoi (*argv);
      else           {*nrep = atoi (*++argv); argc--;}
      break;
    default:
      fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
      break;
    }
  
  if (argc == 1)
    if ((*fp_in = fopen(*argv, "r")) == (FILE*) NULL) {
      sprintf(fname, "%s.fld", *argv);
      if ((*fp_in = fopen(fname, "r")) == (FILE*) NULL) {
	fprintf(stderr, "%s: unable to open input file -- %s or %s\n",
		prog, *argv, fname);
	exit (EXIT_FAILURE);
      }
    }
}
