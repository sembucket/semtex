/*****************************************************************************
 * io.c: I/O, miscellaneous routines.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


void  message (const char* routine, const char* text, int level)
/* ------------------------------------------------------------------------- *
 * A simple error handler.
 * ------------------------------------------------------------------------- */
{
  switch (level) {
  case WARNING:
    fprintf (stderr, "WARNING: %s: %s\n", routine, text); 
    break;
  case ERROR:
    fprintf (stderr, "ERROR: %s: %s\n", routine, text); 
    break;
  case REMARK:
    fprintf (stdout, "%s: %s\n", routine, text);
    break;
  default:
    fprintf (stderr, "bad error level in message: %d\n", level);
    exit (EXIT_FAILURE);
    break;
  }

  if (level == ERROR) exit (EXIT_FAILURE);
}


FILE*  efopen (char* file, char* mode)
/* ------------------------------------------------------------------------- *
 * fopen file, die if can't.   Legal modes are "r",  "w",  "a" and the
 * extended forms "r+", "w+", "a+": see the man page for more details.
 * ------------------------------------------------------------------------- */
{
  FILE*  fp;
  char   s[STRMAX];
  
  if ((fp = fopen (file, mode)) != NULL) return fp;

  sprintf (s, "can't open file %s mode %s\n", file, mode);
  message ("efopen", s, ERROR);
}


void  read_start_file (FILE*  fp, header*  Info)
/* ------------------------------------------------------------------------- *
 * Start-up file is a normal ASCII file, with information about simulation
 * parameters.  We store it all in Run-info.  Do some elementary checks.
 * If it all looks OK, we put in a magic number to identify (binary) restart
 * files.  Note some precision-dependent sscanfs.
 * ------------------------------------------------------------------------- */
{
  char    routine[] = "read_start_file";
  string  line;
  int     i;

  fgets (Info -> Title, STRMAX, fp);
  i = STRMAX;
  while (!(isprint (Info -> Title[i-1])) && i>0) i--;
  Info->Title[i] = '\0';
  fprintf (stdout, "Title: %s\n", Info -> Title);

  i = 0;
  fgets (line, STRMAX, fp);
  while (isspace (line[i])) i++;
  while (!(isspace (line[i])) && isprint(line[i])) i++;
  line[i] = '\0';
  sscanf  (line, "%s", Info -> IC_File);
  fprintf (stdout, "IC file:             %s\n", Info -> IC_File);

  fgets   (line, STRMAX, fp);
  sscanf  (line, "%d", &Info -> N_Grid);
  fprintf (stdout, "Grid size:           %d\n", Info -> N_Grid);

  fgets (line, STRMAX, fp);
  if (sizeof (real) == sizeof (float)) sscanf (line, "%f",  &Info -> Delta_T);
  else                                 sscanf (line, "%lf", &Info -> Delta_T);
  fprintf (stdout, "Time step:           %g\n", Info -> Delta_T);

  fgets   (line, STRMAX, fp);
  sscanf  (line, "%d", &Info -> N_Save);
  fprintf (stdout, "Save increment:      %d\n", Info -> N_Save);

  fgets   (line, STRMAX, fp);
  sscanf  (line, "%d", &Info -> Max_Step);
  fprintf (stdout, "Maximum steps:       %d\n", Info -> Max_Step);

  fgets (line, STRMAX, fp);
  if (sizeof (real) == sizeof (float)) sscanf (line, "%f",  &Info -> K_Visc);
  else                                 sscanf (line, "%lf", &Info -> K_Visc);
  fprintf (stdout, "Kinematic viscosity: %g\n", Info -> K_Visc);

  if (!(ispow2 (Info -> N_Grid)))
    message (routine, "Number of modes must be a power of 2",           ERROR);
  if (Info -> Delta_T <= 0.0)
    message (routine, "Need a positive timestep",                       ERROR);
  if (Info -> N_Save > Info -> Max_Step)
    message (routine, "Save increment exceeds maximum number of steps", ERROR);
  if (Info -> K_Visc < 0.0)
    message (routine, "Need a non-negative viscosity",                  ERROR);

  Info -> Magic = MAGIC;
}
  

void  read_field (FILE* fp, CVF Z, const int  Npts)
/* ------------------------------------------------------------------------- *
 * Read the 3 components of Z.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "read_field";

  if (fread (&Z[1][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read first component data",  ERROR);
  if (fread (&Z[2][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read second component data", ERROR);
  if (fread (&Z[3][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to read third component data",  ERROR);
}
  

void write_field (FILE* fp, CVF Z, const int Npts)
/* ------------------------------------------------------------------------- *
 * Write the 3 components of Z.
 * ------------------------------------------------------------------------- */
{
  const char routine[] = "write_field";

  if (fwrite (&Z[1][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write first component data",  ERROR);
  if (fwrite (&Z[2][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write second component data", ERROR);
  if (fwrite (&Z[3][0][0][0], sizeof (complex), Npts, fp) != Npts)
    message (routine, "Unable to write third component data",  ERROR);
}


void make_file_name (const header* Info, string  Restart_file)
/* ------------------------------------------------------------------------- *
 * Restart file names are generated by appending a numeric tag to the end of
 * the initial condition filename given in Info.  The tag number is the
 * number of timesteps divided by the save interval.
 * ------------------------------------------------------------------------- */
{
  string tag;

  sprintf (tag, "%d", Info -> N_Step / Info -> N_Save); 
  Restart_file = strcpy (Restart_file, Info -> IC_File);
  Restart_file = strcat (Restart_file, ".");
  Restart_file = strcat (Restart_file, tag);
}


void write_restart (FILE*          fp      , 
		    const header*  Info, 
		    const CVF      Vel     ,
		    const CVF      G_old   ,
		    const int      Npts    )
/* ------------------------------------------------------------------------- *
 * Save a binary restart-file.
 * ------------------------------------------------------------------------- */
{
  write_header (fp, Info);
  write_field  (fp, Vel,   Npts);
  write_field  (fp, G_old, Npts);
}


void  print_header (FILE* fp, const header* H)
/* ------------------------------------------------------------------------- *
 * Output header info (ASCII).
 * ------------------------------------------------------------------------- */
{
  fprintf (fp, "Simulation name:                  %s\n",  H -> Title   );
  fprintf (fp, "Root for file names:              %s\n",  H -> IC_File );
  fprintf (fp, "Grid size:                        %1d\n", H -> N_Grid  );
  fprintf (fp, "Time step:                        %g\n",  H -> Delta_T );
  fprintf (fp, "Number of steps between restarts: %1d\n", H -> N_Save  );
  fprintf (fp, "Maximum number of steps:          %1d\n", H -> Max_Step);
  fprintf (fp, "Number of steps so far:           %1d\n", H -> N_Step  );
  fprintf (fp, "Kinematic viscosity:              %g\n",  H -> K_Visc  );
}


void  read_header (FILE* fp, header* H)
/* ------------------------------------------------------------------------- *
 * Read header info (binary).
 * ------------------------------------------------------------------------- */
{
  char routine[] = "read_header";

  if (fread (H, 1, sizeof (header), fp) != sizeof (header))
    message (routine, "Can't get header from file", ERROR);

  if (H -> Magic != MAGIC)
    message (routine, "input file not an ISO field file", ERROR);   
}


void  write_header (FILE* fp, const header* H)
/* ------------------------------------------------------------------------- *
 * Output header info (binary).
 * ------------------------------------------------------------------------- */
{
  if (fwrite (H, 1, sizeof (header), fp) != sizeof (header))
    message ("write_header", "Can't put header on file", ERROR);
}
