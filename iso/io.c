/*****************************************************************************
 * io.c: I/O, miscellaneous routines.
 *
 * Copyright (C) 1992, 1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


void message (const char* routine,
	      const char* text   ,
	      int         level  )
/* ------------------------------------------------------------------------- *
 * A simple error handler.
 * ------------------------------------------------------------------------- */
{
  switch (level) {
  case WARNING: fprintf (stderr, "WARNING: %s: %s\n", routine, text); break;
  case ERROR:   fprintf (stderr, "ERROR: %s: %s\n",   routine, text); break;
  case REMARK:  fprintf (stdout, "%s: %s\n",          routine, text); break;
  default:
    fprintf (stderr, "bad error level in message: %d\n", level);
    exit (EXIT_FAILURE);
    break;
  }

  if (level == ERROR) exit (EXIT_FAILURE);
}


FILE* efopen (const char* file,
	      const char* mode)
/* ------------------------------------------------------------------------- *
 * fopen file, die if can't.   Legal modes are "r",  "w",  "a" and the
 * extended forms "r+", "w+", "a+": see the man page for more details.
 * ------------------------------------------------------------------------- */
{
  FILE* fp;
  char  s[STR_MAX];
  
  if ((fp = fopen (file, mode)) != NULL) return fp;

  sprintf (s, "can't open file %s mode %s\n", file, mode);
  message ("efopen", s, ERROR);
  return 0;
}


void readCVF (FILE*      fp,
	      CVF        Z )
/* ------------------------------------------------------------------------- *
 * Read the 3 components of Z.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "readCVF";
  const int  Npts = N * N * K;

  if (fread (&Z[1][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read first component data",  ERROR);
  if (fread (&Z[2][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read second component data", ERROR);
  if (fread (&Z[3][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read third component data",  ERROR);
}
  

void writeCVF (FILE*      fp,
	       const CVF  Z )
/* ------------------------------------------------------------------------- *
 * Write the 3 components of Z.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "writeCVF";
  const int  Npts = N * N * K;

  if (fwrite (&Z[1][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write first component data",  ERROR);
  if (fwrite (&Z[2][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write second component data", ERROR);
  if (fwrite (&Z[3][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write third component data",  ERROR);
}


void printParam (FILE*        fp,
		 const Param* H )
/* ------------------------------------------------------------------------- *
 * Output Param info (ASCII).
 * ------------------------------------------------------------------------- */
{
  fprintf (fp, "# Session name                : %s\n",  H -> session);
  fprintf (fp, "# Grid size                   : %1d\n", H -> ngrid  );
  fprintf (fp, "# Timestepping order          : %1d\n", H -> norder );
  fprintf (fp, "# Steps between field dumps   : %1d\n", H -> io_fld );
  fprintf (fp, "# Steps between history dumps : %1d\n", H -> io_his );
  fprintf (fp, "# Maximum number of steps     : %1d\n", H -> nstep  );
  fprintf (fp, "# Kinematic viscosity         : %g\n",  H -> kinvis );
  fprintf (fp, "# Time step                   : %g\n",  H -> dt     );
  fprintf (fp, "# Time                        : %g\n",  H -> time   );

  if (H -> chkpnt == TRUE)
    fprintf (fp, "# Checkpointing               : on\n\n");
  else
    fprintf (fp, "# Checkpointing               : off\n\n");

  fflush (fp);
}


void readParam (FILE*  fp,
		Param* I )
/* ------------------------------------------------------------------------- *
 * Read Param info from top of field file.  Must be matched to writeParam.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "readParam";
  char   s[STR_MAX], err[STR_MAX];
  int    ngrid;
  double time, kinvis;

  fgets  (s, STR_MAX, fp);
  fgets  (s, STR_MAX, fp);
  sscanf (s, "%d", &ngrid);
  if (!ispow2 (ngrid)) {
    sprintf (err, "value read for modes (%1d) not a power of 2", ngrid);
    message (routine, err, ERROR);
  }
  fgets  (s, STR_MAX, fp);
  fgets  (s, STR_MAX, fp);
  sscanf (s, "%lf",   &kinvis);
  fgets  (s, STR_MAX, fp);
  fgets  (s, STR_MAX, fp);
  sscanf (s, "%lf",   &time);

  I -> ngrid  = ngrid;
  I -> kinvis = kinvis;
  I -> time   = time;
}


void writeParam (FILE*        fp,
		 const Param* I )
/* ------------------------------------------------------------------------- *
 * Output Param field file header info .
 * ------------------------------------------------------------------------- */
{
  char routine[] = "writeParam";
  char s[STR_MAX], err[STR_MAX];
  char *hdr_fmt[]  = {
    "%-25s  Session\n",
    "%-25d  Grid size\n",
    "%-25d  Step\n",
    "%-25g  Kinematic viscosity\n",
    "%-25g  Time step\n",
    "%-25g  Time\n" };

  fprintf (fp, hdr_fmt[0], I -> session);
  fprintf (fp, hdr_fmt[1], I -> ngrid);
  fprintf (fp, hdr_fmt[2], I -> step);
  fprintf (fp, hdr_fmt[3], I -> kinvis);
  fprintf (fp, hdr_fmt[4], I -> dt);
  fprintf (fp, hdr_fmt[5], I -> time);
}


void startup (Param* I)
/* ------------------------------------------------------------------------- *
 * Open restart file and extract mesh size, open run-time files.
 * ------------------------------------------------------------------------- */
{
  char  s[STR_MAX];
  FILE* fp;

  strcat (strcpy (s, I -> session), ".rst");
  fp = efopen (s, "rb"); readParam (fp, I); fclose (fp);
  
  if   (I -> chkpnt) strcat (strcpy (s, I -> session), ".chk");
  else               strcat (strcpy (s, I -> session), ".fld");
  I -> fld_dmp = efopen (s, "w");

  strcat (strcpy (s, I -> session), ".his");
  I -> his_dmp = efopen (s, "w");
  fprintf (I -> his_dmp, "#     Time Mode         Energy\n");
  fprintf (I -> his_dmp, "# ----------------------------\n");
}


void restart (CVF    U,
	      Param* I)
/* ------------------------------------------------------------------------- *
 * Get initial conditions from restart file, check match with declared size.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "restart";
  char  s[STR_MAX], err[STR_MAX];
  int   ngrid  = I -> ngrid;	/* -- Save to over-ride restart file values. */
  real  kinvis = I -> kinvis;
  real  dt     = I -> dt;
  FILE* fp;
  
  strcat (strcpy (s, I -> session), ".rst");
  fp = efopen (s, "rb");

  readParam (fp, I);
  if (ngrid != I -> ngrid) {
    sprintf (err, "restart ngrid (%1d) clash (%1d) (NEVER HAPPEN)",
	     ngrid, I -> ngrid);
    message (routine, err, ERROR);
  }
  I -> kinvis = kinvis;
  I -> dt     = dt;

  readCVF (fp, U);
}


void analyze (CVF            U   ,
	      Param*         I   ,
	      const complex* Wtab)
/* ------------------------------------------------------------------------- *
 * Carry out periodic analysis of data.
 * ------------------------------------------------------------------------- */
{
  printf ("Step: %-8d Time: %-8g Energy: %-8g",
	  I -> step, I -> time, energyF (U));

#ifdef TG                /* -- Diagnostic for inviscid Taylor--Green vortex. */
  printf (" %-8g %-8g", rmsEns (U), Brachet (I -> time));
#endif

  printf ("\n");
}


void dump (const CVF  U,
	   Param*     I)
/* ------------------------------------------------------------------------- *
 * Write/append a field dump to output file, do history output.
 * ------------------------------------------------------------------------- */
{
  if (!(I -> step % I -> io_his)) {
    real* spec  = (real*) malloc (K * sizeof (real));
    int   i;

    energySpec (U, spec);

    for (i = 1; i < K; i++)	/* -- Ignore mode zero. */
      fprintf (I -> his_dmp, "%g %3d %g\n", I -> time, i, spec[i]);

    fflush (I -> his_dmp);
    fflush (stdout);

    free (spec);
  }

  if (!(I -> step % I -> io_fld) || (I -> step == I -> nstep)) {
    if (I -> chkpnt) {
      char s[STR_MAX], b[STR_MAX], c[STR_MAX];
      
      fclose (I -> fld_dmp);
      strcat (strcpy (s, I -> session), ".chk");
      strcat (strcpy (b, I -> session), ".chk.bak");

      sprintf (c, "mv ./%s ./%s", s, b);
      system  (c);
      I -> fld_dmp = efopen (s, "w");
    }
    writeParam (I -> fld_dmp, I);
    writeCVF   (I -> fld_dmp, U);
    fflush     (I -> fld_dmp);
  }
}


void cleanup (Param* I)
/* ------------------------------------------------------------------------- *
 * Final operations on field files.
 * ------------------------------------------------------------------------- */
{
  if (I -> chkpnt) {
    char s[STR_MAX], b[STR_MAX], c[STR_MAX];

    fclose (I -> fld_dmp);
    strcat (strcpy (s, I -> session), ".chk");
    strcat (strcpy (b, I -> session), ".fld");

    sprintf (c, "mv ./%s ./%s", s, b);
    system  (c);
    strcat  (s, ".bak");
    sprintf (c, "rm -f ./%s", s);
    system  (c);
  } else {
    fclose (I -> fld_dmp);
  }
  
  fclose (I -> his_dmp);
}


void format (char* fmt)
/* ------------------------------------------------------------------------- *
 * Describe binary format of this machine.
 * ------------------------------------------------------------------------- */
{
  union { float  f; int i;    unsigned char c[4]; } v;
  union { double d; int i[2]; unsigned char c[8]; } u;
  int reverse = (-1);
  u.d = 3;
  v.i = 3;
  if      (u.c[0]==64 && u.c[1]==8 && v.c[3]==3) reverse = 0;
  else if (u.c[7]==64 && u.c[6]==8 && v.c[0]==3) reverse = 1;

  switch (reverse) {
  case 1:
    strcpy (fmt, "IEEE little-endian");
    break;
  case -1:
    strcpy (fmt, "unknown");
    break;
  case 0: default:
    strcpy (fmt, "IEEE big-endian");
    break;
  }
}
