/*****************************************************************************
 * wavestress.c: from a 2D/complex modal data file, compute 2D
 * distributions of streamwise-averaged Reynolds stresses.  Data must
 * be binary format and contain only fields u v p.
 *
 * Copyright (c) 2011 <--> $Date$, Hugh Blackburn
 *
 * USAGE
 * -----
 * wavestress [-h] [input[.fld]
 *
 * INPUT FILE
 * ----------
 * Contains only fields uvwp and have N_Z = 2 (Real and Imaginary parts).
 *
 * OUTPUT FILE
 * -----------
 * Is a standard 2D/real (N_Z = 1) Reynolds stress file containing
 * uvpABC, with uvp set to zero and
 *
 * A = u.Re^2 + u.Im^2
 * B = u.Re*v.Re + u.im*v.Im
 * C = v.Re^2 + v.Im^2
 *****************************************************************************/

static char RCS[] = "$Id$";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

static void getargs (int, char**, FILE**);
static int  _index  (const char*, char);

static char prog[] = "wavestress";
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
  int    i, j, n, np, nz, nel, swab = 0;
  int    nfields, nplane;
  FILE   *fp_in = stdin, *fp_out = stdout;
  double **idata, **odata, *vcmpt1, *vcmpt2;

  getargs (argc, argv, &fp_in);
  format  (fmt);

  while (fgets (buf, STR_MAX, fp_in)) { 

    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    
    if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    if (nz != 2) {
      sprintf (fields, "input must have nz = 2 (here %1d)", nz);
      message (prog, fields, ERROR);
    }
    fprintf (fp_out, hdr_fmt[2], np, np, 1, nel);

    n = 6;
    while (--n) {
      fgets (buf, STR_MAX, fp_in);
      fputs (buf, fp_out);
    }

    fgets(fields, STR_MAX, fp_in);
    memset(fields+25, '\0', STR_MAX-25);
    for (nfields = 0, i = 0; i < 25; i++) if (isalpha(fields[i])) nfields++;
    if (!((nfields == 4) && (strstr (fields, "uvwp")))) 
	message (prog, "input must have only fields u v w p.", ERROR);
    fprintf (fp_out, hdr_fmt[8], "uvpABC");

    fgets (buf, STR_MAX, fp_in);
    for (i = 0; i < strlen (buf); i++) buf[i] = tolower (buf[i]);

    if (!strstr(buf, "binary"))
      message (prog, "input file not binary format", ERROR);
    if (!strstr (buf, "endian"))
      message (prog, "input field file in unknown binary format", WARNING);
    else
      swab = ((strstr (buf, "big") && strstr (fmt, "little")) ||
	      (strstr (fmt, "big") && strstr (buf, "little")) );
    sprintf (buf, "%s %s", "binary", fmt);
    fprintf (fp_out, hdr_fmt[9], buf);

    /* -- Set sizes, allocate storage. */

    nplane = np * np * nel;

    idata = dmatrix (0, 3, 0, nplane * 2);
    odata = dmatrix (0, 5, 0, nplane);
    
    /* -- Read in all data fields. */

    dzero (4*nplane, idata[0], 1);
    dzero (6*nplane, odata[0], 1);

    for (i = 0; i < nfields; i++) {
      if (fread (idata[i], sizeof (double), nplane * 2, fp_in) != nplane * 2)
	message (prog, "an error occured while reading", ERROR);
      if (swab) dbrev (nplane*2, idata[i], 1, idata[i], 1);
    }
    
    /* -- Compute A. */

    vcmpt1 = idata[0]; 		/* -- Real part of u. */
    vcmpt2 = idata[0] + nplane; /* -- Imag part of u. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[3], 1, odata[3], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[3], 1, odata[3], 1);

    /* -- Compute B . */

    vcmpt1 = idata[0];		/* -- Real part of u. */
    vcmpt2 = idata[1];		/* -- Real part of v. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[4], 1, odata[4], 1);

    vcmpt1 = idata[0] + nplane; /* -- Imag part of u. */
    vcmpt2 = idata[1] + nplane;	/* -- Imag part of v. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[4], 1, odata[4], 1);

    /* -- Compute C. */

    vcmpt1 = idata[1];		/* -- Real part of v. */
    vcmpt2 = idata[1] + nplane; /* -- Imag part of v. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[5], 1, odata[5], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[5], 1, odata[5], 1);

    /* -- Write out uvpABC in binary. */
    
    for (i = 0; i < 6; i++)
      if (fwrite (odata[i], sizeof (double), nplane, fp_out) != nplane)
	message (prog, "an error occured while writing", ERROR);

    freeDmatrix (idata, 0, 0);
    freeDmatrix (odata, 0, 0);
  } 
  
  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     FILE** fp_in)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] = "wavestress [-h] [input[.fld]\n";

  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
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


static int _index (const char* s,
		   char        c)
/* ------------------------------------------------------------------------- *
 * Return index of c in s, -1 if not found.
 * ------------------------------------------------------------------------- */
{
  int       i;
  const int len = strlen (s);

  for (i = 0; i < len; i++) if (s[i] == c) return i;

  return -1;
}
