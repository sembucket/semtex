/*****************************************************************************
 * MakeIC1: generate a set of initial condition Fourier coefficients for DNS.
 *
 * Usage: makeIC1 -k <peak> -n <cubesize> -o <outfilename> [-s <seed>]
 *
 * Files: write the data to a file (named <outfile>) which contains a header
 * structure, of which only the size of the cube is of any interest (it is
 * checked later by the simulation programme) and then the complex data for
 * the three velocity components in turn.
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include "iso.h"

int N, K, FourKon3;		/* -- Global grid size variables. */


int main (int    argc,
	  char** argv)
{
  FILE*    fp;
  char     filename[STR_MAX];
  int      c, i, argnr;
  int      paramerr = FALSE;
  int      seed = -1;
  CVF      IC;
  CF       work;
  real**   head;
  complex* Wtab;
  Param*   Info = (Param*) calloc (1, sizeof (Param));
  real     peaK = 4.0, energy, lambda;

  /* -- Process command-line arguments. */

  if (argc > 4 && argc < 8) {
    argnr = 1;
    do {
      if (argv[argnr][0] == '-') {
	if (argv[argnr][1] == 'n') {
	  argnr++;
	  Info -> ngrid = atoi (argv[argnr]);
	  if (!ispow2 (Info -> ngrid)) {
	    paramerr = TRUE;
	    fprintf (stderr, "size must be power of 2\n");
	  }
	} else if (argv[argnr][1] == 'k') {
	  argnr++;
	  peaK = atof (argv[argnr]);
	} else if (argv[argnr][1] == 'o') {
	  argnr++;
	  strcpy (filename, argv[argnr]);
	} else if (argv[argnr][1] == 's') {
	  argnr++;
	  seed = atoi (argv[argnr]);
	} else {
	  paramerr = TRUE;
	  fprintf (stderr, "unrecognized flag: %s\n", argv[argnr]);
	}
      } else {
	paramerr = TRUE;
	fprintf (stderr, "bad flag: %s\n", argv[argnr]);
      }
      argnr++;
    } while (!(paramerr) && argnr < argc);
  } else {
    paramerr = TRUE;
    fprintf (stderr, "arg count\n");
  }
  
  if (paramerr) {
    fprintf (stderr,
      "Usage: makeIC1 -k <peak> -n <cubesize> -o <outfilename> [-s <seed>]\n");
    exit (EXIT_FAILURE);
  }

  Info -> dt     = 0.01;
  Info -> step   = 0;
  Info -> kinvis = 0.01;

  strcpy ((Info -> session =
	   malloc (sizeof (int) * strlen ("Decaying spectrum") + 1)),
	  "Decaying spectrum");

  /* -- Allocate storage of IC components, zero all locations. */

  N        = Info -> ngrid;
  K        = N / 2;
  FourKon3 = (4 * K) / 3;

  head = cfield (&IC);
  cbox (0, N-1, 0, N-1, 0, K-1, &work);

  Wtab = cvector (0, K - 1);
  preFFT (Wtab, K);  

  /* -- Generate divergence-free initial condition Fourier coefficients. */

  randomise (IC, &seed, work);
#ifdef DEBUG
  printf ("Check energy: k=0.5<UiUi>: %g\n", energyF (IC));
#endif

  ispectrum (IC, 1.0, peaK, unknown01);
#ifdef DEBUG
  printf ("Check energy: k=0.5<UiUi>: %g\n", energyF (IC));
#endif

  normaliseVF (IC);

  energy = energyF (IC);
  lambda = microF  (IC, work);

  printf ("Check energy: k=0.5<UiUi>: %g\n", energy);
  printf ("Taylor microscale        : %g\n", lambda);

  /* -- Output to file. */

  fp = efopen (filename, "w");
  writeParam  (fp, Info);
  writeCVF    (fp, IC);
  fclose      (fp);

  return EXIT_SUCCESS;
}


  
