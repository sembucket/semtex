/*****************************************************************************
 * deltaweight.c: take roll forcing output from stressdiv and apply
 * delta-function-like restriction to the critical layer region (as
 * determined from wave base flow w velocity field).
 * 
 * Copyright (c) 2011 <--> $Date$, Hugh Blackburn
 *
 * USAGE
 * -----
 * deltaweight [-h] [-w <num>] wave[.bse] force[.fld]
 *   -h       ... print this message
 *   -w <num> ... give weight function width (chi). Default: 0.001.
 *
 * SYNOPSIS
 * --------

 * Read in all the (unweighted) forcing data, and as that happens also
 * reproduce its header info on stdout. Then deal with the wave base
 * flow: check it conforms (N_P, N_Z, NEL) with forcing data, read in
 * only its 'w' data, use that to create the delta function
 * approximation exp(-w*w/chi)/sqrt(PI*chi). Apply it as a multiplying
 * factor/weight to the forcing data, then put that out to stdout and
 * exit.
 *
 * INPUT FILES
 * -----------
 * Both input files must have N_Z = 1, same N_P and NEL.
 * wave.bse has fields uvwp
 * force.fld has fields uv
 *
 * OUTPUT FILE
 * -----------
 * Same format as force.fld.
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

static void getargs (int, char**, double*, FILE**, FILE**);
static int  _index  (const char*, char);

static char prog[] = "deltaweight";
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
  int    i, np, nz, nel, NP, NZ, NEL, swab = 0;
  int    nfields, nplane;
  FILE   *fp_wave = NULL, *fp_forc = NULL, *fp_out = stdout;
  double **data, chi = 0.001;

  getargs (argc, argv, &chi, &fp_wave, &fp_forc);
  format  (fmt);

  /* -- Force file. */

  fgets (buf, STR_MAX, fp_forc); fputs (buf, fp_out); 
  fgets (buf, STR_MAX, fp_forc); fputs (buf, fp_out);

  fgets (buf, STR_MAX, fp_forc);  
  if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
    message (prog, "unable to read force file size", ERROR);

  if (nz != 1) {
    sprintf (fields, "inputs must have nz = 1 (here %1d)", nz);
    message (prog, fields, ERROR);
  }
  fprintf (fp_out, hdr_fmt[2], np, np, 1, nel);

  i = 6; while (--i) { fgets (buf, STR_MAX, fp_forc); fputs (buf, fp_out); }

  fgets(fields, STR_MAX, fp_forc);
  memset(fields+25, '\0', STR_MAX-25);
  for (nfields = 0, i = 0; i < 25; i++) if (isalpha(fields[i])) nfields++;
  if (!((nfields == 2) && (strstr (fields, "uv")))) 
    message (prog, "force file must have only fields u v.", ERROR);
  fprintf (fp_out, hdr_fmt[8], "uv");

  fgets (buf, STR_MAX, fp_forc);
  for (i = 0; i < strlen (buf); i++) buf[i] = tolower (buf[i]);
  if (!strstr(buf, "binary"))
    message (prog, "input file not binary format", ERROR);
  if (!strstr (buf, "endian"))
    message (prog, "force field file in unknown binary format", WARNING);
  else
    swab = ((strstr (buf, "big") && strstr (fmt, "little")) ||
	    (strstr (fmt, "big") && strstr (buf, "little")) );
  sprintf (buf, "%s %s", "binary", fmt);
  fprintf (fp_out, hdr_fmt[9], buf);

  /* -- Set sizes, allocate storage for all fields. */

  nplane = np * np * nel;
  data = dmatrix (0, 3, 0, nplane); /* -- First 2 are uv, 3rd is (wave) w. */
  dzero (3*nplane, data[0], 1);
  
  /* -- Read in force data fields. */

  for (i = 0; i < 2; i++) {
    if (fread (data[i], sizeof (double), nplane, fp_forc) != nplane)
      message (prog, "an error occured while reading force data", ERROR);
    if (swab) dbrev (nplane, data[i], 1, data[i], 1);
  }

  /* -- Now deal with wave base flow. We only want 'w' velocity. */

  fgets (buf, STR_MAX, fp_wave);
  fgets (buf, STR_MAX, fp_wave);

  fgets (buf, STR_MAX, fp_wave);  
  if (sscanf (buf, "%d%*s%d%d", &NP, &NZ, &NEL) != 3)
    message (prog, "unable to read wave size", ERROR);

  if (nz != NZ || np != NP || nel != NEL) 
    message (prog, "wave base and force files have different size", ERROR);
  
  i = 6; while (--i) { fgets (buf, STR_MAX, fp_wave); }

  fgets(fields, STR_MAX, fp_wave);
  memset(fields+25, '\0', STR_MAX-25);
  for (nfields = 0, i = 0; i < 25; i++) if (isalpha(fields[i])) nfields++;
  if (!((nfields == 4) && (strstr (fields, "uvwp")))) 
    message (prog, "wave file must have only fields u v w p.", ERROR);

  fgets (buf, STR_MAX, fp_wave);
  for (i = 0; i < strlen (buf); i++) buf[i] = tolower (buf[i]);
  if (!strstr(buf, "binary"))
    message (prog, "wave base file not binary format", ERROR);
  if (!strstr (buf, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);
  else
    swab = ((strstr (buf, "big") && strstr (fmt, "little")) ||
	    (strstr (fmt, "big") && strstr (buf, "little")) );

  for (i = 0; i < 3; i++)
    if (fread (data[2], sizeof (double), nplane, fp_wave) != nplane)
      message (prog, "an error occured while reading wave data", ERROR);
  if (swab) dbrev (nplane, data[2], 1, data[2], 1);

  /* -- We are done with the inputs. */

  close (fp_forc);
  close (fp_wave);

  /* -- Compute delta function weighting approx exp(-w*w/chi)/sqrt(PI*chi). */

  dvmul (nplane, data[2], 1, data[2], 1, data[2], 1);
  dsmul (nplane, -1.0/chi, data[2], 1, data[2], 1);
  dvexp (nplane, data[2], 1, data[2], 1);
  dsmul (nplane, 1.0/sqrt (M_PI*chi), data[2], 1, data[2], 1);

  /* -- Multiply this into Reynolds stress forcing uv. */

  dvmul (nplane, data[2], 1, data[0], 1, data[0], 1);
  dvmul (nplane, data[2], 1, data[1], 1, data[1], 1);
  
  /* -- Write out uv in binary. */
    
  for (i = 0; i < 2; i++)
    if (fwrite (data[i], sizeof (double), nplane, fp_out) != nplane)
      message (prog, "an error occured while writing", ERROR);

  freeDmatrix (data, 0, 0);
  
  return EXIT_SUCCESS;
}


static void getargs (int     argc   ,
		     char**  argv   ,
		     double* chi    ,
		     FILE**  fp_wave,
		     FILE**  fp_forc)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] =
    "deltaweight [-h] [-w <num>] wave[.bse] force[.fld]\n"
    "   -h       ... print this message\n"
    "   -w <num> ... give weight function width (chi). Default: 0.001.\n";

  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'w':
      if (*++argv[0]) *chi = atof (*argv);
      else {*chi = atof (*++argv); argc--;}
      break;
    default:
      fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
      break;
    }

  if (argc != 2) message (prog, "need 2 input files", ERROR);

  *fp_wave = efopen (argv[0], "r");
  *fp_forc = efopen (argv[1], "r");
}
