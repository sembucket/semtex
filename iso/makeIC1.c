/*****************************************************************************
 * MakeIC1: generate a set of initial condition Fourier coefficients for DNS.
 *
 * Usage: makeIC1 -n <cubesize> -o <outfilename> [-s <seed>]
 *
 * Method: in this version, the Fourier coefficients are generated with a
 * "top-hat" spectrum: i.e., uniform over a range of wavenumber magnitude.
 * If K = <cubesize> / 2 (Nyquist wavenumber), then we give a uniform magni-
 * tude to E(k) over K/3 <= k <= 2K/3, where k = sqrt(k1^2 + k2^2 + k3^2).
 * The method used to generate the various components of the velocity
 * Fourier coefficients is that discussed in detail by Rogallo (1981).
 * We generate the positive k3 part of the spectrum only, since we make use
 * of the fact that the Fourier coefficients are conjugate-symmetric in the
 * simulation.
 *
 * The optional command-line argument -s <seed> can be invoked with a
 * positive integer (e.g., -s 1) to generate different random components.
 * The IC components are normalized such that q^2 = <UiUi>/2 = 1.
 * The maximum velocity component is found (in physical space), and the
 * maximum timestep allowed by the CFL condition Delta_X/Max_Vel is printed
 * on standard output.  The grid-size is calculated with the assumption that
 * the length of the sides of the computational box is 2PI.
 *
 * Files: write the data to a file (named <outfile>) which contains a header
 * structure, of which only the size of the cube is of any interest (it is
 * checked later by the simulation programme) and then the complex data for
 * the three velocity components in turn.
 *
 * $Id$
 *
 * ------------------------------------------------------------------------- */

#include "iso.h"


int main (int argc, char *argv[])
{
  FILE                  *fp;
  string                filename;
  int                   c, i, argnr, cubesize, Npts;
  int                   paramerr = FALSE, seed = -1;
  CVF  IC;
  real**      head;
  int*               Dimension;
  complex*               Wtab;
  header                IC_info;
  real                 Max_Vel;


  /* -- Process command-line arguments. */

  if (argc > 4 && argc < 8) {
    argnr = 1;
    do {
     if (argv[argnr][0] == '-') {
       if (argv[argnr][1] == 'n') {
	 argnr++;
	 cubesize = atoi(argv[argnr]);
	 if (!ispow2(cubesize)) {
	   paramerr = TRUE;
	   fprintf(stderr, "size must be power of 2\n");
	 }
       } else if (argv[argnr][1] == 'o') {
	 argnr++;
	 (void) strcpy(filename, argv[argnr]);
       } else if (argv[argnr][1] == 's') {
	 argnr++;
	 seed = atoi(argv[argnr]);
       } else {
	 paramerr = TRUE;
	 fprintf(stderr, "unrecognized flag: %s\n", argv[argnr]);
       }
     } else {
       paramerr = TRUE;
       fprintf(stderr, "bad flag: %s\n", argv[argnr]);
     }
     argnr++;
   } while (!(paramerr) && argnr < argc);
  } else {
    paramerr = TRUE;
    fprintf(stderr, "arg count\n");
  }
  
  if (paramerr) {
    fprintf(stderr,
	    "Usage: makeIC -n <cubesize> -o <outfilename> [-s <seed>]\n");
    exit (1);
  }

  /* -- Allocate storage of IC components, zero all locations. */

  Dimension    = ivect (1, 3);
  Dimension[1] = (Dimension[2] = cubesize);
  Dimension[3] = cubesize / 2;

  Npts = Dimension[1] * Dimension[2] * Dimension[3];
  
  head = cfield (Dimension, &IC);
  Wtab = cvect  (0, Dimension[3]-1);

  for (c = 1; c <= 3; c++)
    for (i = 0; i < 2*Npts; i++)
      head[c][i] = 0.0;
 
  /* -- Generate initial condition Fourier coefficients. */

  tophat (Dimension, IC, seed);

  /* -- Normalize the velocities to get q^2 = 1. 
   *    Find out the largest velocity component in physical space. */

  preFFT (Dimension[3], Wtab);
  Max_Vel = normalize (Dimension, Wtab, head, IC);

  /* -- Output to file. */

  IC_info.N_Grid = cubesize;
  fp = efopen (filename, "w");
  if (fwrite (&IC_info, 1, sizeof (header), fp) != sizeof (header)) {
    fprintf (stderr, "makeIC1: couldn't output header onto file\n");
    exit (EXIT_FAILURE);
  }
  
  for (c = 1; c <= 3; c++)
    if (fwrite (&IC[c][0][0][0], sizeof (complex), Npts, fp) != Npts) {
      fprintf (stderr, "Unable to write component %d data on IC file", c);
      exit (EXIT_FAILURE);
    }

  fclose (fp);

  /* -- Calculate the CFL condition timestep. */

  printf ("%.3e\n", 2.0*M_PI / (cubesize * Max_Vel));


/* Check IC energy ... <UiUi>/2:
   printf("q^2: %f\n", energy (Dimension, IC));
*/
  return EXIT_SUCCESS;
}


  
