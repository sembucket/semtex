/*****************************************************************************
 * RELV: add in reference frame velocities to field-file velocities.
 *       Optionally modify associated mesh file to reflect change of position.
 *
 * USAGE: relv [-h] [-g grid.msh grid.out] [-n select]
 *             [-o out.fld] [-p] [-s state.sta] [-V U V] [-v] [file.fld]
 *
 * Relv uses at least two files, one is an (ASCII) history file containing
 * reference frame positions and velocities in the X & Y directions (with
 * other state variables) and associated time step information, the other
 * is an field file.  It is assumed that the first two fields in this file
 * are X & Y direction velocity fields; the corresponding step number is
 * searched for in the history file and the appropriate reference frame
 * velocities are added to the velocity fields.  By default, the result is
 * written to stdout [optionally on a file named using the -o option].
 * NB: if the state file is NOT selected then zero values are substituted
 * for reference frame positions and velocities.
 *
 * An extension to these functions is the modification of an associated mesh
 * file (produced by meshpr) to reflect the changes in postion of the mesh
 * at the same times as for the data contained in the field file.
 * The resultant meshes are written (in concatenated form if appropriate)
 * on a file named by grid.out.  Relv may be used to alter the mesh position
 * ONLY, without changing velocities, by selecting option -p.
 *
 * By default, relv will read through an extended field file with many
 * dumps, converting as it goes.  Optionally [-n select] a particular dump
 * may be selected from the file (starting at 1).
 *
 * The following is a typical input file header:
 *
 * sample                      Session
 * Mon Apr 22 18:23:13 91      Created
 * 9 9 1 8                     Nr, Ns, Nz, Nelt
 * 50                          Step
 * 0.05                        Time
 * 0.001                       Time step
 * 0.025                       Kinvis
 * 1                           Beta-z
 * U V P                       Fields written
 * ascii                       Format
 *
 * The following is a typical line from a state-variable history file:
 *
 * step  tttt  xxxx  xvxv xaxa prxf vixf toxf yyyy yvyv yaya pryf viyf toyf
 * (step)
 *  where:
 *        step is step number,              tttt is time,
 *
 *        xxxx is X position,               yyyy is Y position,
 *        xvxv is X velocity,               yvyv is Y velocity,
 *        xaxa is X acceleration,           yaya is Y acceleration,
 *        prxf is X pressure force,         pryf is Y pressure force,
 *        vixf is X viscous  force,         viyf is Y viscous  force,
 *        toxf is X total    force,         toyf is Y total    force.
 *****************************************************************************/

static char
RCSid[] = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <femdef.h>
#include <alplib.h>

static void getargs      (int, char**, char**, char**, char**, char**,
			  char**, double*, double*, int*, int*, int*);
static void do_ascii     (FILE*, FILE*, int, int, double, double, int);
static void do_binary    (FILE*, FILE*, int, int, double, double, int, int);
static void find_step    (FILE*, int, double*, double*, double*, double*,
			  double*, double*, double*, double*);
static int  count_fields (char*);

static char prog[] = "relv";
static const char *hdr_fmt[] = { 
  "%-25s "    "Session\n",
  "%-25s "    "Created\n",
  "%-25s "    "Nr, Ns, Nz, Elements\n",
  "%-25d "    "Step\n",
  "%-25.6g "  "Time\n",
  "%-25.6g "  "Time step\n",
  "%-25.6g "  "Kinvis\n",
  "%-25.6g "  "Beta\n",
  "%-25s "    "Fields written\n",
  "%-25s "    "Format\n"
};


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  char    fmt[STR_MAX], *c;
  char    *field=0, *state=0, *outfile=0, *meshin=0, *meshout = 0;
  int     i, nr, ns, nz, nel, npts, nfields, step, his_step, swab = 0;
  int     nread = 0, selected = 1, nselect = 0, verbose = 0, position = 0;
  double  X, Y, Xp, Yp, Xv, Yv, Xa, Ya, Xf, Yf, U=0.0, V=0.0;
  FILE    *fps_in = 0, *fpf_in = 0, *fpf_out = 0, *fpm_in = 0, *fpm_out = 0;

  getargs (argc, argv, &field, &state, &outfile, &meshin, &meshout, 
	   &U, &V, &nselect, &verbose, &position);

  format (fmt);

  /* -- Open files. */

  if (field) {
    if ( !(fpf_in = fopen (field, "r")) ) {
      sprintf(buf, "%s.fld", field);
      if ( !(fpf_in = fopen (buf, "r")) ) {
	sprintf (buf, "unable to open input file -- %s or %s", field, buf);
	message (prog, buf, ERROR);
      }
    }
  } else
    fpf_in = stdin;

  fpf_out = (outfile) ? fopen (outfile, "w") : stdout;

  if (state)
    if ( !(fps_in = fopen (state, "r")) ) {
      sprintf (buf, "unable to open state file %s", state); 
      message (prog, buf, ERROR);
    }

  if (meshin) {
    if ( !(fpm_in = fopen (meshin, "r")) ) {
      sprintf (buf, "unable to open mesh input file %s", meshin);
      message (prog, buf, ERROR);
    }
    if ( !(fpm_out = fopen (meshout, "w")) ) {
      sprintf (buf, "unable to open mesh output file %s", meshout);
      message (prog, buf, ERROR);
    }
  }

  while (fgets (buf, STR_MAX, fpf_in)) {

    nread++;
    if (nselect) selected = nselect == nread;

    /* -- Find the number of field data. */

    i = 3;
    while (--i) {
      if (selected) fputs (buf, fpf_out);
      fgets (buf, STR_MAX, fpf_in);
    }

    if (sscanf (buf, "%d%d%d%d", &nr, &ns, &nz, &nel) != 4)
      message (prog, "unable to read the file size", ERROR);
    npts = nr * nr * nz * nel;
    if (selected) fputs (buf, fpf_out);

    /* --
     * Get step number, rewind state-history file & search for frame
     * velocities at this step.  Die if can't.  Action is over-ridden if no
     * state variable file was selected.                                    */
   
    fgets (buf, STR_MAX, fpf_in);
    if (selected) fputs (buf, fpf_out);
    
    if (state) {
      if (!sscanf (buf, "%d", &step)) 
	message (prog, "unable to read the step number", ERROR);
      find_step (fps_in, step, &Xp, &Yp, &Xv, &Yv, &Xa, &Ya, &Xf, &Yf);
      if (position) Xv = Yv = 0.0;
    } else
      Xp = Yp = Xv = Yv = Xa = Ya = Xf = Yf = 0.0;

    Xv -= U;
    Yv -= V;

    /* -- Reproduce unused strings, look for number of fields. */

    fgets (buf, STR_MAX, fpf_in);
    i = 5;
    while (--i) {
      if (selected) fputs (buf, fpf_out);
      fgets (buf, STR_MAX, fpf_in);
    }

    if ((nfields = count_fields(buf)) < 2)
      message (prog, "expected at least two fields", ERROR);
    if (selected) fputs (buf, fpf_out);

    /* -- Check file format "ASCII"/"ascii" or "BINARY"/"binary". */

    c = fgets (buf, STR_MAX, fpf_in);
    while (*c++ = tolower(*c));
    
    if (strstr (buf, "binary")) {
      if (!strstr (buf, "endian"))
	message (prog, "input field file in unknown binary format", WARNING);
      else {
	swab = (   (strstr (buf, "big") && strstr (fmt, "little"))
	        || (strstr (fmt, "big") && strstr (buf, "little")) );
      }
      sprintf (buf, "binary ");
      format  (buf + strlen(buf));
      if (selected) fprintf (fpf_out, hdr_fmt[9], buf);
      do_binary (fpf_in, fpf_out, npts, nfields, Xv, Yv, selected, swab);
      
    } else if (strstr(buf, "ascii")) {
      if (selected) fprintf (fpf_out, hdr_fmt[9], "ASCII");
      do_ascii  (fpf_in, fpf_out, npts, nfields, Xv, Yv, selected);

    } else
      message (prog, "got a bad file format message", ERROR);

    /* -- Now do the grid-file work if requested: mesh file is 2-col ASCII. */

    if (selected && fpm_in) {

      rewind (fpm_in);

      fgets (buf,  STR_MAX, fpm_in);
      fputs (buf,  fpm_out);

      while (fgets(buf, STR_MAX, fpm_in)) {
	sscanf (buf, "%lf %lf", &X, &Y);
	X += Xp;
	Y += Yp;
	fprintf(fpm_out, "%#14.7g %14.7g\n", X, Y);
      }
    }

    /* -- Put state variables on stderr. */

    if (verbose && selected) {
      fprintf(stderr, "X-STATE: %14g %14g %14g %14g\n", Xp, Xv, Xa, Xf);
      fprintf(stderr, "Y-STATE: %14g %14g %14g %14g\n", Yp, Yv, Ya, Yf);
    }

    if (nselect && selected) break;
  }

  if (nread < nselect) {
    sprintf (buf, "asked for dump number %1d, only read %1d", nselect, nread);
    message (prog, buf, ERROR);
  }
    
  return EXIT_SUCCESS;
}


static void do_ascii (FILE*  fpi     ,
		      FILE*  fpo     ,
		      int    npts    ,
		      int    nfields ,
		      double xv      ,
		      double yv      ,
		      int    selected)
/* ------------------------------------------------------------------------- *
 * Process ASCII-format input, produce ASCII-format output.
 * ------------------------------------------------------------------------- */
{
  register int i, j;
  double*      data = dvector (0, nfields - 1);
  
  /* -- Read the numbers from the file, add frame velocities, print up. */

  for (j=0; j<npts; j++) {
    for (i=0; i<nfields; i++)
      if (fscanf(fpi, "%lf", &data[i]) != 1) {
	sprintf (buf, "unable to read a number--line %d, field %d\n", j+1,i+1);
	message (prog, buf, ERROR);
      }
    fgets(buf, STR_MAX, fpi);

    data[0] += xv;
    data[1] += yv;

    if (selected) {
      for (i=0; i<nfields; i++)
	if (fprintf (fpo, "%#16.10g ", data[i]) < 0)
	  message (prog, "an error has occured while writing", ERROR);
      fputs("\n", fpo);
    }
  }

  freeDvector (data, 0);
}


static void do_binary (FILE*  fpi     ,
		       FILE*  fpo     ,
		       int    npts    ,
		       int    nfields ,
		       double xv      ,
		       double yv      ,
		       int    selected,
		       int    swab    )
/* ------------------------------------------------------------------------- *
 * Process binary-format input, produce binary-format output.
 * ------------------------------------------------------------------------- */
{
  register int i;
  double**     data = dmatrix (0, nfields - 1, 0, npts - 1);

  /* -- Read the numbers from the file, byte-swapping if needed. */

  for (i = 0; i < nfields; i++) {
    if (fread(data[i], sizeof(double), npts, fpi) != npts)
      message (prog, "an error has occured while reading", ERROR);
    if (swab) dbrev (npts, data[i], 1, data[i], 1);
  }

  /* -- Add frame velocities. */

  dsadd (npts, xv, data[0], 1, data[0], 1);
  dsadd (npts, yv, data[1], 1, data[1], 1);

  /* -- Write binary output. */

  if (selected)
    for (i = 0; i < nfields; i++)
      if (fwrite(data[i], sizeof(double), npts, fpo) != npts)
	message (prog, "an error has occured while writing", ERROR); 

  freeDmatrix (data, 0, 0);
}


static int count_fields (char* s)
/* ------------------------------------------------------------------------- *
 * Count the number of field names in a string.
 * ------------------------------------------------------------------------- */
{
  int n = 0, i = 0;

  while (i++ < 25) if (isalpha(*s++)) n++;

  return(n);
}


static void getargs (int     argc,
		     char**  argv,
		     char**  field,
		     char**  state,
		     char**  outfile,
		     char**  meshin,
		     char**  meshout,
		     double* U,
		     double* V,
		     int*    nselect,
		     int*    verbose,
		     int*    position)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c;
  static char *usage  =
    "relv [options] [file.fld]\n"
    "  [options] are:\n"
    "  -h                   ... print this message\n"
    "  -g grid.msh grid.out ... grid.msh from gridgen, output grid.out\n"
    "  -n select            ... deal with only this dump in field file\n"
    "  -o out.fld           ... named output field file\n"
    "  -p                   ... change only mesh positions, not velocities\n"
    "  -s state.sta         ... select a state-variable file (D: no motion)\n"
    "  -V U V               ... subtract off global velocities U & V too\n"
    "  -v                   ... verbose output: body position on stderr\n";

  if (argc < 3) fprintf (stderr, usage);

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf(stderr, usage);
      exit(0);
      break;
    case 'g':
      --argc;
      *meshin  = *++argv;
      --argc;
      *meshout = *++argv;
      break;
    case 'n':
      if (*++argv[0])
	*nselect = atoi(*argv);
      else {
	--argc;
	*nselect = atoi(*++argv);
      }
      break;      
    case 'o':
      if (*++argv[0])
	*outfile = *argv;
      else {
	--argc;
	*outfile = *++argv;
      }
      break;
    case 'p':
      *position = 1;
      break;
    case 's':
      if (*++argv[0])
	*state = *argv;
      else {
	--argc;
	*state = *++argv;
      }
      break;
    case 'V':
      --argc;
      *U = atof(*++argv);
      --argc;
      *V = atof(*++argv);
      break;
    case 'v':
      *verbose = 1;
      break;
    default:
      fprintf(stderr, "relv: illegal option: %c\n", c);
      break;
    }

  if (argc == 1) *field = *argv;
} 


static void find_step (FILE*  fp,
		       int    step,
		       double *Xp, double *Yp,
		       double *Xv, double *Yv,
		       double *Xa, double *Ya,
		       double *Xf, double *Yf)
/* ------------------------------------------------------------------------- *
 * Rewind file to find info for step number.
 * ------------------------------------------------------------------------- */
{
  int n, found = 0;

  rewind (fp);
  while (fgets (buf, STR_MAX, fp) && !found) {
    sscanf (buf, "%d %*s %lf %lf %lf %*s %*s %lf %lf %lf %lf %*s %*s %lf",
	    &n, Xp, Xv, Xa, Xf, Yp, Yv, Ya, Yf);
    found = (n == step);
  }

  if (!found) {
    sprintf (buf, "couldn't find step %d in state history file", step);
    message (prog, buf, ERROR);
  }
}
