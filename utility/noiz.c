/* ************************************************************************* *
 * NOIZ:  add a random Gaussian perturbation to a velocity field.            *
 *                                                                           *
 * USAGE: noiz [-h] [-o output] [-p perturb] [input[.fld]                    *
 *                                                                           *
 * SYNOPSIS:                                                                 *
 * Noiz reads a field file and adds a gaussian-distributed random variable   *
 * of specified standard deviation to each velocity datum.  Fields may be in * 
 * ASCII or binary format, output is in same format.                         *
 *                                                                           *
 * NOTES:                                                                    *
 * This is a simple adaption of Ron's CONVERT program.                       *
 * Default value of perturbation is 1.2x10^-7.                               *
 *                                                                           *
 * The following is a typical header:                                        *
 *                                                                           *
 * sample                      Session                                       *
 * Mon Apr 22 18:23:13 91      Created                                       *
 * 9 9 1 8                     Nr, Ns, Nz, Nelt                              *
 * 50                          Step                                          *
 * 0.05                        Time                                          *
 * 0.001                       Time step                                     *
 * 0.025                       Kinvis                                        *
 * 1                           Beta-z                                        *
 * U V P                       Fields written                                *
 * ascii                       Format                                        *
 *                                                                           *
 * ************************************************************************* */

/*------------------*
 * RCS Information: *
 *------------------*/
static char
  RCSid[] = "$Id$";


#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <nrutil.h>


#define IA    16807
#define IM    2147483647
#define AM    (1.0/IM)
#define IQ    127773
#define IR    2836
#define NTAB  32
#define NDIV  (1+(IM-1)/NTAB)
#define EPS1  1.2e-7		/* single precision */
#define EPS2  1.2e-13		/* double precision */
#define RNMX  (1.0-EPS2)


void   parse_args   (int, char **, FILE **, FILE **, double *);
void   a_to_a       (int, int, FILE *, FILE *,  char *, double);
void   b_to_b       (int, int, FILE *, FILE *,  char *, double);
void   dswap        (int, double *);
double gasdev       (long *);   
int    count_fields (char *);




main(int argc, char *argv[])
/* ========================================================================= *
 * Wrapper.                                                                  *
 * ========================================================================= */
{
  char   buf[BUFSIZ], fields[BUFSIZ], *c;
  int    n, nr, ns, nz, nel;
  FILE  *fp_in   = stdin,
        *fp_out  = stdout;
  double pert = EPS1;
  

  parse_args(argc, argv, &fp_in, &fp_out, &pert);

  pert *= pert;	       /* we will need a perturbation variance later */

  while (fgets(buf, BUFSIZ, fp_in)) { 

    n = 3;
    while (--n) {
      fputs(buf, fp_out);
      fgets(buf, BUFSIZ, fp_in);
    }

    if (sscanf(buf, "%d%d%d%d", &nr, &ns, &nz, &nel) != 4)
      message("noiz", "unable to read the file size", ERROR);               

    n = 6;
    while (--n) {
      fputs(buf, fp_out);
      fgets(buf, BUFSIZ, fp_in);
    }
    fputs(buf, fp_out);   

    fgets(fields, BUFSIZ, fp_in);
    
    n = count_fields(fields);
    fputs (fields, fp_out);
    fgets (buf, BUFSIZ, fp_in);
    fputs (buf, fp_out);

    c = buf;
    while (isspace(*c)) c++;

    switch (*c) {
    case 'a': case 'A':
      a_to_a (nr * ns * nz * nel, n, fp_in, fp_out, fields, pert);
      break;

    case 'b': case 'B':
      b_to_b (nr * ns * nz * nel, n, fp_in, fp_out, fields, pert);
      break;

    default:
      sprintf (buf, "unknown format flag -- %c", c);
      message ("noiz", buf, ERROR);
      break;
    }

  } 

  return;
}





void a_to_a(int    npts,
	    int    nfields, 
	    FILE  *in, 
	    FILE  *out, 
	    char  *fields, 
	    double pert)
/* ========================================================================= *
 * ASCII input.                                                              *
 * ========================================================================= */
{
  int    i, j;
  char   buf[128];
  double datum;
  long   seed = 0;


  for (j = 0; j < npts; j++) {
    for (i = 0; i < nfields; i++)
      if (fscanf(in, "%lf", &datum) == 1) {
	switch (fields[i]) {
	case 'u': case 'v': case 'w' :
	  datum += pert*gasdev(&seed);
	default:
	  if (fprintf(out, "%#16.10g ", datum) < 1)
	    message("noiz", "an error has occured while writing", ERROR);
	  break;
	} 
      } else {
	sprintf(buf, "unable to read a number -- line %d, field %d\n",
		j+1, i+1);
	message("noiz", buf, ERROR);
      }
    fgets(buf, BUFSIZ, in);
    fprintf(out, "\n");
  }
}





void b_to_b(int    npts,
	    int    nfields, 
	    FILE  *in, 
	    FILE  *out, 
	    char  *fields, 
	    double pert)
/*===========================================================================*
 * Convert binary to ASCII.                                                  *
 *===========================================================================*/
{
  int      i, j;
  double **data;
  long     seed = 0;


  data = dmatrix(0, nfields-1, 0, npts-1);

  for (i = 0; i < nfields; i++) {
    if (fread(data[i], sizeof(double), npts, in) != npts)
      message("noiz", "an error has occured while reading", ERROR);
    switch (fields[i]) {
    case 'u': case 'v': case 'w':
      for (j = 0; j < npts; j++)
	data[i][j] += pert * gasdev(&seed);
    default: break;
    }
  }

  if (fwrite(data[0], sizeof(double), nfields*npts, out) != nfields*npts)
    message("noiz",  "an error has occured while writing", ERROR);
  
  free_dmatrix(data, 0, 0);
  }





int count_fields(char *s)
/* ========================================================================= *
 * Count the number of field names in a string.                              *
 * ========================================================================= */
{
  int n = 0, i = 0;

  while (i++ < 25) if (isalpha(*s++)) n++;
  
  return n;
}





void parse_args(int argc, char *argv[], FILE **fp_in, FILE **fp_out, double *p)
/* ========================================================================= *
 * Parse command line arguments.                                             *
 * ========================================================================= */
{
  char  c;
  int   i;
  char  fname[FILENAME_MAX];
  char *usage   = "usage: convert [-h] [-o output] [-p perturb] "
                       "[input[.fld]]\n"
                       "options:\n"
                       "-h         ... print this help message\n"
                       "-o output  ... write to named file\n"
                       "-p perturb ... standard deviation of perutrbation\n";


  while (--argc && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {
      case 'h':
	fputs(usage, stderr);
	exit (0);
	break;
      case 'o':
	if (*++argv[0])
	  strcpy(fname, *argv);
	else {
	  strcpy(fname, *++argv);
	  argc--;
	}
	if ((*fp_out = fopen(fname,"w")) == (FILE*) NULL) {
	  fprintf(stderr, "convert: unable to open the output file -- %s\n", 
		  fname);
	  exit(1);
	}
	*argv += strlen(*argv)-1;
	break;
      case 'p':
	if (*++argv[0])
	  *p = atof(*argv);
	else {
	  *p = atof(*++argv);
	  argc--;
	}
	*argv += strlen(*argv)-1;
	break;
      default:
	fprintf(stderr, "convert: unknown option -- %c\n", c);
	break;
      }

  if (argc == 1)
    if ((*fp_in = fopen(*argv, "r")) == (FILE*) NULL) {
      sprintf(fname, "%s.fld", *argv);
      if ((*fp_in = fopen(fname, "r")) == (FILE*) NULL) {
	fprintf(stderr, "convert: unable to open input file -- %s or %s\n",
		*argv, fname);
	exit(1);
      }
    }

  return;
}




double ran1(long *idum)
/* ========================================================================= *
 * Generate IUD random variates on (0, 1).  Numerical Recipes.               *
 * ========================================================================= */
{
  int         j;
  long        k;
  static long iy = 0;
  static long iv[NTAB];
  double      temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum)/IQ;
      *idum = IA * (*idum - k*IQ) - IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j]= *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k*IQ) - IR*k;
  if (*idum < 0) *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp = AM * iy) > RNMX) return RNMX;
  else return temp;
}
 


double gasdev(long *idum)
/* ========================================================================= *
 * Generate normally distributed deviate with zero mean & unit variance.     *
 * Numerical Recipes.                                                        *
 * ========================================================================= */
{
  static int    iset = 0;
  static double gset;
  double        fac, r, v1, v2;
    
  if  (iset == 0) {
    do {
      v1= 2.0 * ran1(idum) - 1.0;
      v2= 2.0 * ran1(idum) - 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1.0);
    fac = sqrt(-2.0 * log(r) / r);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}
