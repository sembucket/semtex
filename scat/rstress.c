/*****************************************************************************
 * rstress.c:  Compute Reynolds stresses.
 *
 * NB: This version is specialised for use with files made by scat. See
 * chknames() below.
 *
 * Copyright (c) 1997 <--> $Date$, Hugh Blackburn
 *
 * --
 *
 * This file is part of Semtex.
 * 
 * Semtex is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * Semtex is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Semtex (see the file COPYING); if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 *
 * --
 *
 * SYNOPSIS
 * --------
 * Rstress deals with field-average files.  If called with just an
 * average file (e.g. session.avg) as input, it tries to compute
 * Reynolds stress components using the correlations and average
 * fields it is assumed to contain.  If called with an additional
 * field file, the average values contained in the average file are
 * subtracted from the field file (i.e. the field file is assumed to
 * contain instantaneous values from which the average is to be
 * subtracted).  Note that this also allows rstress to be used to
 * subtract (instantaneous) values in one field dump from another.
 *
 * Both average and field files are assumed to be in binary, double,
 * physical-space format.  The field files and average files (if both
 * present) must have the same element/Fourier structure.  Only the
 * first dump in each file is dealt with.
 *
 * Product terms are computed without dealiasing.
 * 
 * USAGE
 * -----
 * rstress [options] avg.file [field.file]
 * options:
 * -h ... print this message
 *
 * This version is spacialised to deal with averages that also contain
 * the scalr field 'c'. Here are the names of the correlation
 * variables:
 *
 * Average "velocity" fields are called: u, v, c (& w); Product
 * average fields are called: A, B, C, D, E, F (& G, H, I, J)
 *
 *                   / uu uv uc uw \     /  A  B  D  G \
 *                   | .  vv vc uv |  =  |  .  C  E  H |
 *                   | .  .  cc cw |     |  .  .  F  I |
 *                   \ .  .  .  ww /     \  .  .  .  J /
 *****************************************************************************/

static char RCS[] = "$Id$";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include <cfemdef.h>
#include <cveclib.h>

static char  prog[]    = "rstress";
static char* hdr_fmt[] = {	 /* -- Header output formatting. */
  "%-25s "             "Session\n",
  "%-25s "             "Created\n",
  "%-5d%-5d%-5d%-10d " "Nr, Ns, Nz, Elements\n",
  "%-25d "             "Step\n",
  "%-25.6g "           "Time\n",
  "%-25.6g "           "Time step\n",
  "%-25.6g "           "Kinvis\n",
  "%-25.6g "           "Beta\n",
  "%-25s "             "Fields written\n",
  "%-25s "             "Format\n"
};
typedef struct {		 /* -- Data information structure. */
  char     session [STR_MAX];
  char     creation[STR_MAX];
  int      np               ;
  int      nz               ;
  int      nel              ;
  int      step             ;
  double   time             ;
  double   timestep         ;
  double   kinvis           ;
  double   beta             ;
  char     field [STR_MAX]  ;
  char     format[STR_MAX]  ;
  double** data             ;
} Dump;


static int  _index    (const char*, char);
static void getargs   (int, char**, FILE**, FILE**);
static void getheader (FILE*, Dump*);
static void getdata   (FILE*, Dump*);
static void chknames  (const char*);
static void covary    (Dump*);
static void demean    (Dump*, Dump*);
static void printup   (FILE*, Dump*);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine.
 * ------------------------------------------------------------------------- */
{
  FILE *avgfile = 0, *fldfile = 0;
  Dump *heada   = 0, *headf   = 0;
  int  i, npts, ntot;

  getargs (argc, argv, &avgfile, &fldfile);

  /* -- Read in the average file. */

  if (fldfile) {		 /* -- Remove averages from fldfile. */
    heada = (Dump*) calloc (1, sizeof (Dump));
    headf = (Dump*) calloc (1, sizeof (Dump));

    getheader (avgfile, heada);
    getheader (fldfile, headf);

    if (heada -> np  != headf -> np  ||
	heada -> nz  != headf -> nz  ||
	heada -> nel != headf -> nel)
      message (prog, "structure of files don't match",           ERROR);

    for (i = 0; i < strlen(headf -> field); i++)
      if (!strchr (heada -> field, headf -> field[i]))
      message (prog, "average fields don't match dumped fields", ERROR);

    getdata   (avgfile, heada);
    getdata   (fldfile, headf);
    demean    (heada,   headf);
    printup   (stdout,  headf);

  } else {			 /* -- Reynolds stresses using avgfile. */
    heada = (Dump*) calloc (1, sizeof (Dump));

    getheader (avgfile, heada);
    chknames  (heada -> field);
    getdata   (avgfile, heada);
    covary    (heada);
    printup   (stdout, heada);
  }
    
  /* -- Printup. */

  return (EXIT_SUCCESS);
}


static void getargs (int    argc   ,
		     char** argv   ,
		     FILE** avgfile,
		     FILE** fldfile)
/* ------------------------------------------------------------------------- *
 * Parse command-line arguments.
 * ------------------------------------------------------------------------- */
{
  char usage[] = "usage: rstress [options] avg.file [field.file]\n"
                 "options:\n"
                 "  -h ... display this message\n";
  char err[STR_MAX], c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf (stderr, usage); exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (err, "illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  switch (argc) {
  case 1:
    *avgfile = efopen (argv[0], "r");
    break;
  case 2:
    *avgfile = efopen (argv[0], "r");
    *fldfile = efopen (argv[1], "r");
    break;
  default:
    fprintf (stderr, usage); exit (EXIT_FAILURE);
    break;
  }  
}


static void getheader (FILE* f,
		       Dump* h)
/* ------------------------------------------------------------------------- *
 * Find header information.
 * ------------------------------------------------------------------------- */
{
  char buf[STR_MAX];

  fgets  (h -> session,  STR_MAX, f);
  fgets  (h -> creation, STR_MAX, f);
  fscanf (f, "%d %*s %d %d", &h -> np, &h -> nz, &h -> nel);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%d", &h -> step);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> time);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> timestep);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> kinvis);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> beta);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%s", h -> field);
  fgets  (buf, STR_MAX, f);
  fgets  (h -> format, STR_MAX, f);

  if (!strstr (h -> format,      "binary"))
    message (prog, "input field file not in binary format",     ERROR);
  else if (!strstr (h -> format, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);
}


static void getdata (FILE* f,
		     Dump* h)
/* ------------------------------------------------------------------------- *
 * Find the number of fields, allocate storage & do binary read.
 * ------------------------------------------------------------------------- */
{
  char      localfmt[STR_MAX];
  int       i, swab;
  const int npts    = h -> np * h -> np * h -> nz * h -> nel;
  const int nfields = strlen (h -> field);
  const int ntot    = npts * nfields;

  h -> data = dmatrix (0, nfields - 1, 0, npts - 1);

  if (fread (h -> data[0], sizeof (double), ntot, f) != ntot) {
    sprintf (localfmt, "could not read fields: %s", h -> field);
    message (prog, localfmt, ERROR);
  }

  format (localfmt);
  swab = (strstr (h -> format, "big") && strstr (localfmt,    "little")) || 
         (strstr (localfmt,    "big") && strstr (h -> format, "little"));

  if (swab) dbrev (ntot, h -> data[0], 1, h -> data[0], 1);
}


static void chknames (const char* field)
/* ------------------------------------------------------------------------- *
 * Check that the names of the enclosed fields make sense for Reynolds
 * stress computations.  Computations could be 2D or 3D.
 *
 * Average "velocity" fields are called: u, v, c (& w);
 * Product average    fields are called: A, B, C, D, E, F (& G, H, I, J)
 *
 *                   / uu uv uc uw \     /  A  B  D  G \
 *                   | .  vv vc uv |  =  |  .  C  E  H |
 *                   | .  .  cc cw |     |  .  .  F  I |
 *                   \ .  .  .  ww /     \  .  .  .  J /
 * ------------------------------------------------------------------------- */
{
  char err[STR_MAX];

  if (!strstr (field, "uv") || !strstr (field, "c")) {
    sprintf (err, "field names (%s) should contain \"uv\" & \"c\"", field);
    message (prog, err, ERROR);
  }
  if (strstr (field, "w"))
    if (!strstr (field, "ABCDEFGHIJ")) {
      sprintf (err, "field names (%s) should contain \"ABCDEFGHIJ\"", field);
      message (prog, err, ERROR);
    }
  else
    if (!strstr (field, "ABCDEF")) {
      sprintf (err, "field names (%s) should contain \"ABCDEF\"", field);
      message (prog, err, ERROR);
    }
}


static int _index (const char* s, char c)
/* ------------------------------------------------------------------------- *
 * Return index of c in s, -1 if not found.
 * ------------------------------------------------------------------------- */
{
  int       i;
  const int len = strlen (s);

  for (i = 0; i < len; i++) if (s[i] == c) return i;

  return -1;
}


static void covary (Dump* h)
/* ------------------------------------------------------------------------- *
 * On input, h contains data for the average velocity fields and the
 * average product velocity fields.  Convert the average product
 * velocity fields to covariances by subtracting the products of the
 * average velocity fields.  This is simple since the data are already
 * in physical space and we don't dealias.
 * ------------------------------------------------------------------------- */
{
  int       i, j, k, m, n;
  const int npts = h -> np * h -> np * h -> nz * h -> nel;
  
  /* -- 2D. */

  i = _index (h -> field, 'u');
  j = _index (h -> field, 'v');
  k = _index (h -> field, 'c');

  m = _index (h -> field, 'A');
  dvvvtm (npts, h->data[m], 1, h->data[i], 1, h->data[i], 1, h->data[m], 1);
  m = _index (h -> field, 'B');
  dvvvtm (npts, h->data[m], 1, h->data[i], 1, h->data[j], 1, h->data[m], 1);
  m = _index (h -> field, 'C');
  dvvvtm (npts, h->data[m], 1, h->data[j], 1, h->data[j], 1, h->data[m], 1);
  m = _index (h -> field, 'D');
  dvvvtm (npts, h->data[m], 1, h->data[i], 1, h->data[k], 1, h->data[m], 1);
  m = _index (h -> field, 'E');
  dvvvtm (npts, h->data[m], 1, h->data[j], 1, h->data[k], 1, h->data[m], 1);
  m = _index (h -> field, 'F');
  dvvvtm (npts, h->data[m], 1, h->data[k], 1, h->data[k], 1, h->data[m], 1);

  if (!strstr (h -> field, "w")) return;

  /* -- 3D. */
  
  m = _index (h -> field, 'w');

  n = _index (h -> field, 'G');
  dvvvtm (npts, h->data[n], 1, h->data[i], 1, h->data[m], 1, h->data[n], 1);
  n = _index (h -> field, 'H');
  dvvvtm (npts, h->data[n], 1, h->data[j], 1, h->data[m], 1, h->data[n], 1);
  n = _index (h -> field, 'I');
  dvvvtm (npts, h->data[n], 1, h->data[k], 1, h->data[m], 1, h->data[n], 1);
  n = _index (h -> field, 'J');
  dvvvtm (npts, h->data[n], 1, h->data[m], 1, h->data[m], 1, h->data[n], 1);
}


static void demean (Dump* a,
		    Dump* f)
/* ------------------------------------------------------------------------- *
 * Subtract mean values in "a" from field values in "f".
 * ------------------------------------------------------------------------- */
{
  int       i, j;
  const int nfields = strlen (f -> field);
  const int npts    = a -> np * a -> np * a -> nz * a -> nel;

  for (i = 0; i < nfields; i++) {
    j = _index (a -> field, f -> field[i]);
    dvsub (npts, f -> data[i], 1, a -> data[j], 1, f -> data[i], 1);
  }
}


static void printup (FILE* f,
		     Dump* h)
/* ------------------------------------------------------------------------- *
 * Output the modified data.
 * ------------------------------------------------------------------------- */
{
  int       i;
  const int ntot = h -> np * h -> np * h -> nz * h -> nel;
  const int nfields = strlen (h -> field);
  char      s1[STR_MAX], s2[STR_MAX];
  time_t    tp = time (0);

  fprintf  (f, h -> session);
  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[2], h -> np, h -> np, h -> nz, h -> nel);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[3], h -> step);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[4], h -> time);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[5], h -> timestep);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[6], h -> kinvis);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[7], h -> beta);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[8], h -> field);
  fprintf  (f, s1);
  sprintf  (s2, "binary "); format (s2 + strlen (s2));
  sprintf  (s1, hdr_fmt[9], s2);
  fprintf  (f, s1);

  for (i = 0; i < nfields; i++) 
    if (fwrite (h -> data[i], sizeof (double), ntot, f) != ntot)
      message (prog, "error occurred while writing", ERROR);
}