/*****************************************************************************
 * CONVERT:  Format conversion program for PRISM-compatible field files.
 *
 * USAGE: convert [-h] [-o output] [input[.fld]
 *
 * Default action is to convert binary to ASCII files or vice versa.
 *
 * Output is either to stdout or an optional file argument.
 *
 * For binary input files, automatic conversion from written format to
 * machine's internal format is carried out before ASCII output, if possible.
 *
 * sample                      Session
 * Mon Apr 22 18:23:13 1991    Created
 * 9 9 1 8                     Nr, Ns, Nz, Nel
 * 50                          Step
 * 0.05                        Time
 * 0.001                       Time step
 * 0.025                       Kinvis
 * 1                           Beta-z
 * uvp                         Fields written
 * ASCII                       Format
 *
 * Other recognized formats are: ascii (for backwards compatibility)
 *                               binary, IEEE little-endian
 *                               binary, IEEE big-endian
 *****************************************************************************/

static char
rcsid[] = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


static char  prog[]  = "convert";
static FILE* fp_in   = stdin;
static FILE* fp_out  = stdout;
static int   swap    = 0;
static char  usage[] = "usage: convert [-h] [-o output] [input[.fld]]\n"
                       "options\n"
                      "-h ... print this message\n";

static void parse_args   (int, char**);
static void err_msg      (const char*);
static void a_to_b       (int, int, FILE*, FILE*);
static void b_to_a       (int, int, FILE*, FILE*);
static void dswap        (int, double*);
static int  count_fields (const char*);
static int  iformat      ();
static void format       (char*);
static int  swap_format  (const char*, const char*);


int main (int argc, char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine.
 * ------------------------------------------------------------------------- */
{
  char buf[BUFSIZ];
  char fmt[BUFSIZ];
  char *c;
  int  n, nr, ns, nz, nel;


  parse_args (argc, argv);

  while (fgets (buf, BUFSIZ, fp_in)) {

    swap = 0;

    /* -- Find size information. */

    n = 3;
    while (--n) {
      fputs (buf, fp_out);
      fgets (buf, BUFSIZ, fp_in);
    }

    if (sscanf(buf, "%d%d%d%d", &nr, &ns, &nz, &nel) != 4)
      err_msg ("unable to read the file size");               

    n = 7;
    while (--n) {
      fputs (buf, fp_out);
      fgets (buf, BUFSIZ, fp_in);
    }

    /* -- Find number of fields & format string. */
    
    n = count_fields(buf);
    fputs (buf, fp_out);
    fgets (buf, BUFSIZ, fp_in);
    c = buf;
    while (isspace (*c)) c++;
    format (fmt);

    switch (*c) {
    case 'a': case 'A':
      fprintf (fp_out, "binary, %-17s Format\n", "binary", fmt);
      a_to_b  (nr * ns * nz * nel, n, fp_in, fp_out);
      break;

    case 'b': case 'B':
      swap = swap_format (c, fmt); 
      fprintf (fp_out, "%-25s Format\n", "ASCII");
      b_to_a  (nr * ns * nz * nel, n, fp_in, fp_out);
      break;

    default:
      sprintf (buf, "unknown format flag -- %c", c);
      err_msg (buf);
      break;
    }
  }

  return EXIT_SUCCESS;
}


void a_to_b (int npts, int nfields, FILE *in, FILE *out)
/* ------------------------------------------------------------------------- *
 * Convert ASCII to binary.
 * ------------------------------------------------------------------------- */
{
  int    i, j;
  char   buf[128];
  double **data;

  /* -- Allocate memory. */

  data = (double**) malloc (nfields * sizeof(double*));
  for (i = 0; i < nfields; i++) 
    data[i] = (double*) malloc (npts*sizeof(double));

  /* -- Read the numbers from the file (ASCII). */

  for (j = 0; j < npts; j++) {
    for (i = 0; i < nfields; i++)
      if (fscanf (in, "%lf", data[i] + j) != 1) {
	sprintf (buf, "unable to read a number -- line %d, field %d\n",
		 j+1, i+1);
	err_msg (buf);
      }
    fgets (buf, 128, in);   /* -- Scan to newline. */
  }

  /* -- Byte-reverse numbers. */

  if (swap) for (i = 0; i < nfields; i++) dswap (npts, data[i]);

  /* -- Write binary output. */

  for (i = 0; i < nfields; i++)
    if (fwrite (data[i], sizeof(double), npts, out) != npts)
      err_msg ("an error has occured while writing");
  
  /* -- Free temporary storage. */

  for (i = 0; i < nfields; i++)
    free (data[i]);
  free (data);
}


void b_to_a (int npts, int nfields, FILE *in, FILE *out)
/* ------------------------------------------------------------------------- *
 * Convert binary to ASCII.
 * ------------------------------------------------------------------------- */
{
  int    i, j;
  double **data;

  /* -- Allocate memory. */

  data = (double**) malloc (nfields * sizeof(double*));
  for (i = 0; i < nfields; i++) 
    data[i] = (double*) malloc (npts*sizeof(double));

  /* -- Read numbers from file. */

  for (i = 0; i < nfields; i++)
    if (fread (data[i], sizeof (double), npts, in) != npts)
      err_msg ("an error has occured while reading");

  /* -- Byte-reverse numbers. */

  if (swap) for (i = 0; i < nfields; i++) dswap (npts, data[i]);

  /* -- Write ASCII output. */

  for (j = 0; j < npts; j++) {
    for (i = 0; i < nfields; i++)
      if (fprintf (out, "%#16.10g ", data[i][j]) < 0)
	err_msg ("an error has occured while writing");
    fputs ("\n", out);
  }
  
  /* -- Free temporary data storage. */

  for (i = 0; i < nfields; i++)
    free (data[i]);
  free (data);
}


static int count_fields (const char *s)
/* ------------------------------------------------------------------------- *
 * Count the number of field names in a string.
 * ------------------------------------------------------------------------- */
{
  int n = 0, i = 0;

  while (i++ < 25) if (isalpha (*s++)) n++;
  
  return n;
}


static void err_msg (const char *s)
/* ------------------------------------------------------------------------- *
 * Print an error message and die.
 * ------------------------------------------------------------------------- */
{
  fprintf (stderr, "%s: %s\n", prog, s);
  exit (EXIT_FAILURE);
}





static void parse_args (int argc, char *argv[])
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char  c;
  int   i;
  char  fname[FILENAME_MAX];

  while (--argc && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {
      case 'h':
	fputs (usage, stderr);
	exit  (EXIT_SUCCESS);
	break;
      case 'r':
	swap = 1;
	break;
      case 'o':
	if (*++argv[0])
	  strcpy (fname, *argv);
	else {
	  strcpy (fname, *++argv);
	  argc--;
	}
	if ((fp_out = fopen (fname,"w")) == (FILE*) NULL) {
	  fprintf (stderr, "convert: unable to open the output file -- %s\n", 
		   fname);
	  exit (EXIT_FAILURE);
	}
	*argv += strlen (*argv) - 1;
	break;

      default:
	fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }

  if (argc == 1)
    if ((fp_in = fopen (*argv, "r")) == (FILE*) NULL) {
      sprintf (fname, "%s.fld", *argv);
      if ((fp_in = fopen (fname, "r")) == (FILE*) NULL) {
	fprintf(stderr, "%s: unable to open input file -- %s or %s\n",
		prog, *argv, fname);
	exit (EXIT_FAILURE);
      }
    }

  return;
}


static void dswap (int n, double *x)
/* ------------------------------------------------------------------------ *
 * Byte-reversal routine.
 * ------------------------------------------------------------------------ */
{
  double *cx = x;
  char   *c  = (char*) x;
  register int i,j;

  for(i = 0; i < n; i++, cx++, c = (char*) cx)
    for(j = 0; j < 4; j++){
      char d = c[j];
      c[j]   = c[7-j];
      c[7-j] = d;
    }

  return;
}


static int iformat (void)
/* ------------------------------------------------------------------------- *
 * Return 1 if machine floating-point format is IEEE little-endian,
 * 0 if IEEE big-endian, -1 for unrecognized format.
 * ------------------------------------------------------------------------- */
{
  union { float f;  int i;    unsigned char c[4]; } v;
  union { double d; int i[2]; unsigned char c[8]; } u;
  int   reverse = (-1);
  u.d = 3;
  v.i = 3;
  if      (u.c[0] == 64 && u.c[1] == 8 && v.c[3] == 3) reverse = 0;
  else if (u.c[7] == 64 && u.c[6] == 8 && v.c[0] == 3) reverse = 1;

  return (reverse);
}


static void format (char* s)
/* ------------------------------------------------------------------------- *
 * Fill s with a string describing machine's floating-point storage format.
 * ------------------------------------------------------------------------- */
{
  switch (iformat ()) {
  case -1:
    sprintf (s, "unknown");
    break;
  case 1:
    sprintf (s, "IEEE little-endian");
    break;
  case 0: default:
    sprintf (s, "IEEE big-endian");
    break;
  }
}


static int swap_format (const char* input, const char* machine)
/* ------------------------------------------------------------------------- *
 * Compare strings describing binary input format and this machine's
 * internal storage format.
 *
 * Return:
 *   1 to flag byte reversal IEEE big-endian <--> IEEE little-endian,
 *   0 for no reversal.
 *
 * Unrecognized input formats cause abortion.
 * ------------------------------------------------------------------------- */
{
  if (!strstr (input, "IEEE")) {
    char s[BUFSIZ];
    sprintf (s, "unrecognized binary format \"%s\"", input);
    err_msg (s);
  }
  
  if (   (strstr (input, "little") && strstr (machine, "big"   ))
      || (strstr (input, "big"   ) && strstr (machine, "little"))) return 1;
  
  return 0;
}
