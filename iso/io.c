/*===========================================================================
 * RCS Information:
 * ----------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 *===========================================================================*/

#include "globals.h"


int ispow2(int num)
/*===========================================================================*/
/* Is num a strictly positive integer power of two? Yes ==> 1, No ==> 0      */
/*===========================================================================*/
{
  while (num % 2 == 0 && num > 2) num /= 2;

  if (num != 2) return 0; else return 1;
}


FILE *efopen(char *file, char *mode)
/*===========================================================================*/
/* fopen file, die if can't.   Legal modes are "r",  "w",  "a" and the       */
/* extended forms "r+", "w+", "a+": see the man page for more details.       */
/*===========================================================================*/
{
  FILE *fp;
  
  if ((fp = fopen(file,mode)) != NULL)
    return fp;
  fprintf(stderr, "can't open file %s mode %s\n", file, mode);
  exit(1);
  /* NOTREACHED */
}

void error(/* message */ string error_text)
/*===========================================================================*/
/* Error handler...abort.                                                    */
/*===========================================================================*/
{
  (void) fprintf(stderr, "iso: Run-time error...aborting\n");
  (void) fprintf(stderr, "%s\n", error_text);
  exit(1);
}


void getargs(int argc, char *argv[], int *Start, string Input_file)
/*===========================================================================*/
/* Process command-line arguments.                                           */
/* Usage: iso -s <start_file> || -r <restart_file>                           */
/*===========================================================================*/
{
  if (argc == 3) {
    if (argv[1][0] == '-' && argv[1][1] == 'r') {
      *Start = FALSE;
      (void) strcpy(Input_file, argv[2]);
    } else if (argv[1][0] == '-' && argv[1][1] == 's') {
      (void) strcpy(Input_file, argv[2]);
    } else {
      (void) fprintf(stderr, "iso: bad flag %s\n", argv[1]);
      error("Usage: iso -s <start_file> || -r <restart_file>");
    }
  } else {
    (void) fprintf(stderr, "iso: arg count\n");
    error("Usage: iso -s <start_file> || -r <restart_file>");
  }
}


void read_start_file(/* from   */ FILE   *fp,
                     /* return */ header *Run_info)
/*===========================================================================*/
/* Start-up file is a normal ASCII file, with information about simulation   */
/* parameters.  We store it all in Run-info.  Do some elementary checks.     */
/* If it all looks OK, we put in a magic number to identify (binary) restart */
/* files.                                                                    */
/*===========================================================================*/
{
  string line;
  int    i;


  (void) fgets(Run_info->Title, MAXSTR, fp);
  i = MAXSTR;
  while (!(isprint(Run_info->Title[i-1])) && i>0) i--;
  Run_info->Title[i] = '\0';
  (void) fprintf(stdout, "Title: %s\n", Run_info->Title);

  i = 0;
  (void) fgets(line, MAXSTR, fp);
  while (isspace(line[i])) i++;
  while (!(isspace(line[i])) && isprint(line[i])) i++;
  line[i] = '\0';
  (void) sscanf(line, "%s", Run_info->IC_File);
  (void) fprintf(stdout, "IC file:             %s\n", Run_info->IC_File);

  (void) fgets(line, MAXSTR, fp);
  (void) sscanf(line, "%d", &Run_info->N_Grid);
  (void) fprintf(stdout, "Grid size:           %d\n", Run_info->N_Grid);

  (void) fgets(line, MAXSTR, fp);
  (void) sscanf(line, "%f", &Run_info->Delta_T);
  (void) fprintf(stdout, "Time step:           %.2e\n", Run_info->Delta_T);

  (void) fgets(line, MAXSTR, fp);
  (void) sscanf(line, "%d", &Run_info->N_Save);
  (void) fprintf(stdout, "Save increment:      %d\n", Run_info->N_Save);

  (void) fgets(line, MAXSTR, fp);
  (void) sscanf(line, "%d", &Run_info->Max_Step);
  (void) fprintf(stdout, "Maximum steps:       %d\n", Run_info->Max_Step);

  (void) fgets(line, MAXSTR, fp);
  (void) sscanf(line, "%f", &Run_info->K_Visc);
  (void) fprintf(stdout, "Kinematic viscosity: %.2e\n", Run_info->K_Visc);


  if (!(ispow2(Run_info->N_Grid)))
    error("Number of modes must be a power of 2");
  if (Run_info->Delta_T <= 0.0)
    error("Need a positive timestep");
  if (Run_info->N_Save > Run_info->Max_Step)
    error("Save increment must be less than the maximum number of steps");
  if (Run_info->K_Visc < 0.0)
    error("Need a non-negative viscosity");

  Run_info->Magic = MAGIC;
}
  

void read_components(/* from  */ FILE                  *fp,
		     /* get   */ complex_vector_field  Z,
		     /* using */ int                   Npts)
/*===========================================================================*/
/* Read the 3 vector components of Z for each location in Fourier space.     */
/*===========================================================================*/
{
  if (fread((char *)&Z[1][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to read first component data");
  if (fread((char *)&Z[2][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to read second component data");
  if (fread((char *)&Z[3][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to read third component data");
}


void make_file_name(/* using  */ header Run_info,
                    /* return */ string Restart_file)
/*===========================================================================*/
/* Restart file names are generated by appending a numeric tag to the end of */
/* the initial condition filename given in Run_info.  The tag number is the  */
/* number of timesteps divided by the save interval.                         */
/*===========================================================================*/
{
  string tag;

  (void) sprintf(tag, "%d", Run_info.N_Step / Run_info.N_Save); 
  Restart_file = strcpy(Restart_file, Run_info.IC_File);
  Restart_file = strcat(Restart_file, ".");
  Restart_file = strcat(Restart_file, tag);
}


void write_restart(/* on    */ FILE                  *fp, 
		   /* put   */ header                Run_info, 
		               complex_vector_field  Vel,
		               complex_vector_field  G_old,
                   /* using */ int                   Npts)
/*===========================================================================*/
/* Save a binary restart-file.                                               */
/*===========================================================================*/
{
  if (fwrite((char *)&Run_info, 1, sizeof(header), fp) != sizeof(header))
    error("Couldn't put header on restart file");

  if (fwrite((char *)&Vel[1][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to write first velocity component data on restart file");
  if (fwrite((char *)&Vel[2][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to write second velocity component data on restart file");
  if (fwrite((char *)&Vel[3][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to write third velocity component data on restart file");

  if (fwrite((char *)&G_old[1][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to write first G_old component data on restart file");
  if (fwrite((char *)&G_old[2][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to write second G_old component data on restart file");
  if (fwrite((char *)&G_old[3][0][0][0],sizeof(complex), Npts, fp) != Npts)
    error("Unable to write third G_old component data on restart file");
}


  
