/*****************************************************************************
 * io.c: I/O, miscellaneous routines.
 *
 * : io.c,v 2.2 1995/11/23 08:14:45 hmb Exp hmb $
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


void  readCVF (FILE*      fp ,
	       CVF        Z  ,
	       const int* Dim)
/* ------------------------------------------------------------------------- *
 * Read the 3 components of Z.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "readCVF";
  const int  Npts = Dim[1] * Dim[2] * Dim[3];

  if (fread (&Z[1][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read first component data",  ERROR);
  if (fread (&Z[2][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read second component data", ERROR);
  if (fread (&Z[3][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read third component data",  ERROR);
}
  

void  writeCVF (FILE*      fp ,
		const CVF  Z  ,
		const int* Dim)
/* ------------------------------------------------------------------------- *
 * Write the 3 components of Z.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "writeCVF";
  const int  Npts = Dim[1] * Dim[2] * Dim[3];

  if (fwrite (&Z[1][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write first component data",  ERROR);
  if (fwrite (&Z[2][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write second component data", ERROR);
  if (fwrite (&Z[3][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write third component data",  ERROR);
}


void printParam (FILE*        fp  ,
		 const Param* H   ,
		 const char*  prog,
		 const char*  rev )
/* ------------------------------------------------------------------------- *
 * Output Param info (ASCII).
 * ------------------------------------------------------------------------- */
{
  if (prog) fprintf (fp, "%s; %s\n",  prog, rev);

  fprintf (fp, "Session name:                  %s\n",  H -> name    );
  fprintf (fp, "Grid size:                     %1d\n", H -> modes   );
  fprintf (fp, "Number of steps between dumps: %1d\n", H -> stepSave);
  fprintf (fp, "Maximum number of steps:       %1d\n", H -> stepMax );
  fprintf (fp, "Reynolds number:               %g\n",  H -> Re      );
  fprintf (fp, "Time step:                     %g\n",  H -> dt      );
  fprintf (fp, "Time:                          %g\n",  H -> time    );
}


void readParam (FILE*  fp,
		Param* I )
/* ------------------------------------------------------------------------- *
 * Read Param info (binary).  Must be matched to writeParam.
 * ------------------------------------------------------------------------- */
{
  char   routine[] = "readParam";
  char   s[STR_MAX], err[STR_MAX];
  int    modes;
  double time, Re;

  fread  (s, 1, sizeof (s), fp);
  fread  (s, 1, sizeof (s), fp);
  sscanf (s, "%d", &modes);
  if (!ispow2 (modes)) {
    sprintf (err, "value read for modes (%1d) not a power of 2", modes);
    message (routine, err, ERROR);
  }
  fread  (s, 1, sizeof (s), fp);
  fread  (s, 1, sizeof (s), fp);
  sscanf (s, "%lf", &Re);
  fread  (s, 1, sizeof (s), fp);
  fread  (s, 1, sizeof (s), fp);
  sscanf (s, "%lf", &time);

  I -> Re    = Re;
  I -> modes = modes;
  I -> time  = time;
}


void writeParam (FILE*        fp,
		 const Param* I )
/* ------------------------------------------------------------------------- *
 * Output Param info (binary).
 * ------------------------------------------------------------------------- */
{
  char routine[] = "writeParam";
  char s[STR_MAX], err[STR_MAX];
  char *hdr_fmt[]  = {
    "%-25s  Session\n",
    "%-25d  Fourier modes\n",
    "%-25d  Step\n",
    "%-25g  Reynolds number\n",
    "%-25g  Time step\n",
    "%-25g  Time\n" };

  memset  (s, '\0', STR_MAX);
  sprintf (s, hdr_fmt[0], I -> name);
  fwrite  (s, 1,  sizeof (s), fp);

  memset  (s, '\0', sizeof (s));
  sprintf (s, hdr_fmt[1], I -> modes);
  fwrite  (s, 1, sizeof (s), fp);
  
  memset  (s, '\0', sizeof (s));
  sprintf (s, hdr_fmt[2], I -> step); 
  fwrite  (s, 1, sizeof (s), fp);

  memset  (s, '\0', sizeof (s));
  sprintf (s, hdr_fmt[3], I -> Re); 
  fwrite  (s, 1, sizeof (s), fp);

  memset  (s, '\0', sizeof (s));
  sprintf (s, hdr_fmt[4], I -> dt); 
  fwrite  (s, 1, sizeof (s), fp);

  memset  (s, '\0', sizeof (s));
  sprintf (s, hdr_fmt[5], I -> time); 
  fwrite  (s, 1, sizeof (s), fp);
}


void startup (FILE*       fp      ,
	      Param*      I       ,
	      const char* session ,
	      const int   chkpoint)
/* ------------------------------------------------------------------------- *
 * Read runtime directives from fp.
 * ------------------------------------------------------------------------- */
{
  char routine[] = "startup";
  char s[STR_MAX], err[STR_MAX];

  /* -- Save session name, create output file. */

  strncpy (I -> name, session, STR_MAX - 5);

  if   (chkpoint) strcat (strcpy (s, I -> name), ".chk");
  else            strcat (strcpy (s, I -> name), ".fld");
  I -> output = efopen (s, "w");

  /* -- Strip title. */

  fgets (s, STR_MAX, fp);

  /* -- Spatial resolution. */

  fgets (s, STR_MAX, fp);
  if (strstr (s, "NGRID")) {
    sscanf (s, "%d", &I -> modes);
    if (!ispow2 (I -> modes) && I -> modes < 4) {
      sprintf (err, "need NGRID to be a power of 2 and >= 4, read: %s", s);
      message (routine, err, ERROR);
    }
  } else {
    sprintf (err, "expected NGRID, read: %s", s);
    message (routine, err, ERROR);
  }

  /* -- Time step. */

  fgets (s, STR_MAX, fp);
  if (strstr (s, "DT")) {
    double d_t;
    sscanf (s, "%lf", &d_t);
    I -> dt = d_t;
  } else {
    sprintf (err, "expected DT, read: %s", s);
    message (routine, err, ERROR);
  }

  /* -- Field dump interval. */

  fgets (s, STR_MAX, fp);
  if (strstr (s, "IO_FLD")) {
    sscanf (s, "%d", &I -> stepSave);
  } else {
    sprintf (err, "expected IO_FLD, read: %s", s);
    message (routine, err, ERROR);
  }
  
  /* -- Number of steps. */

  fgets (s, STR_MAX, fp);
  if (strstr (s, "STEPS")) {
    sscanf (s, "%d", &I -> stepMax);
  } else {
    sprintf (err, "expected STEPS, read: %s", s);
    message (routine, err, ERROR);
  }

  /* -- Reynolds number. */

  fgets (s, STR_MAX, fp);
  if (strstr (s, "RE")) {
    double  re;
    sscanf (s, "%lf", &re);
    I -> Re = re;
  } else {
    sprintf (err, "expected RE, read: %s", s);
    message (routine, err, ERROR);
  }

  I -> time = 0.0;
}


void initialize (CVF        U  ,
		 Param*     I  ,
		 const int* Dim)
/* ------------------------------------------------------------------------- *
 * Get initial conditions from restart file, check match with declared size.
 * ------------------------------------------------------------------------- */
{
  char  routine[STR_MAX] = "initialize";
  char  s[STR_MAX], err[STR_MAX];
  FILE* fp;
  int   modes = I -> modes;
  real  Re    = I -> Re;
  real   dt    = I -> dt;
  
  strcat (strcpy (s, I -> name), ".rst");
  fp = efopen (s, "rb");

  readParam (fp, I);
  if (modes != I -> modes) {
    sprintf (err, "restart modes (%1d) clashes with session declaration (%1d)",
	     modes, I -> modes);
    message (routine, err, ERROR);
  }
  I -> Re = Re;
  I -> dt = dt;

  readCVF   (fp, U, Dim);
}


void analyze (CVF            U   ,
	      Param*         I   ,
	      const complex* Wtab,
	      const int*     Dim )
/* ------------------------------------------------------------------------- *
 * Carry out periodic analysis of data.
 * ------------------------------------------------------------------------- */
{
  printf ("Step: %-8d Time: %-8g Energy: %-8g",
	  I -> step, I -> time, energyF (U, Dim));

#ifdef TG                /* -- Diagnostic for inviscid Taylor--Green vortex. */
  printf (" %-8g %-8g", rmsEns (U, Dim), Brachet (I -> time));
#endif

  printf ("\n");
}


void dump (const CVF  U     ,
	   Param*     I     ,
	   const int  chkpnt,
	   const int* Dim   )
/* ------------------------------------------------------------------------- *
 * Write/append a field dump to output file.
 * ------------------------------------------------------------------------- */
{
  if (I -> step % I -> stepSave) return;

  if (chkpnt) {
    char s[STR_MAX], b[STR_MAX], c[STR_MAX];

    fclose (I -> output);
    strcat (strcpy (s, I -> name), ".chk");
    strcat (strcpy (b, I -> name), ".chk.bak");

    sprintf (c, "mv ./%s ./%s", s, b);
    system  (c);
    I -> output = efopen (s, "w");
  }
  writeParam (I -> output, I);
  writeCVF   (I -> output, U, Dim);
  fflush     (I -> output);
}


void cleanup (Param*    I     ,
	      const int chkpnt)
/* ------------------------------------------------------------------------- *
 * Final operations on field files.
 * ------------------------------------------------------------------------- */
{
  if (chkpnt) {
    char s[STR_MAX], b[STR_MAX], c[STR_MAX];

    fclose (I -> output);
    strcat (strcpy (s, I -> name), ".chk");
    strcat (strcpy (b, I -> name), ".fld");

    sprintf (c, "mv ./%s ./%s", s, b);
    system  (c);
    strcat  (s, ".bak");
    sprintf (c, "rm -f ./%s", s);
    system  (c);
  } else {
    fclose (I -> output);
  }
}
