/*
 * Crossing Frequency 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXVALS    64     /* Maximum number of values */
#define MAXCYCLES 256     /* Maximum number of cycles */

/* Prototypes */

int  next_line  (FILE *, double []);
void parse_args (int argc, char *argv[]);

/* External variables */

int    col = 1;      /* Column to use as the dependent variable */
double val = 0.;     /* Threshold value for crossing frequency  */         

main (int argc, char *argv[])
{
  int nx  = 0;      /* Number of crossing over the full series */
  int nc  = 0;      /* Number of cycles (even crossings)       */
  int ncol, n;

  double a[MAXVALS], x[MAXVALS], freq[MAXCYCLES];

  parse_args (argc, argv);

  puts ("#\n# Crossing points\n#");

  /* Scan through the file to find each point in time where the   *
   * dependent variable crosses the threshold value.  Print the   *
   * values of all variables at these points.                     */

  ncol = next_line (stdin, a);  
  while (next_cross(stdin, a, x, val, col, ncol)) {

    for (n = 0; n < ncol; n++)
      printf ("%#10.6g ", x[n]);
    putchar ('\n');

    /* Every even crossing is a full cycle.  Compute the average cycle  *
     * frequency based on the time origin and number of crossings.  For *
     * now just store the time of each crossing.                        */
   
    if (!(nx++ & 1))               /* Increment crossing count */
      freq [nc++] = x[0];          /* Increment cycle count    */
  }


  /* Summarize the cycle frequencies */

  puts ("#\n"
        "# Cycle   T   Average    f   Average\n#");

  for (n = 1;  n < nc; n++) {
    const double T  = (freq[n]-freq[n-1]),
                 Ta = (freq[n]-freq[0])/n;
    printf ("%3d %#10.6g %#10.6g %#10.6g %#10.6g\n", n, T, Ta, 1./T, 1./Ta);
  }

  return 0;
}

/* Find the next crossing */

int next_cross 
  (FILE *fp, double a[], double x[], double val, int col, int N)
{
  double ta, tb, va, vb, alpha;
  double b[MAXVALS];
  int crossed = 0;
  int n;

  while (!crossed && next_line(fp, b) == N) {
    
    ta = a[0];   va = a[col];
    tb = b[0];   vb = b[col];
    
    if ( ( va <= val && vb > val ) ||
         ( va >= val && vb < val ) ) 
      {
	alpha = (val - va) / (vb - va);
	for (n = 0; n < N; n++)
	  x [n] = a[n] + alpha * (b[n] - a[n]);
	crossed = 1;
      }

    memcpy (a, b, N*sizeof(double));
  }

  return crossed;
}

/* Read the next line of input */

int next_line (FILE *fp, double vbuf [MAXVALS])
{
  char buf [BUFSIZ], *p;
  int  n;

  /* Load the next line */

  do 
    p = fgets (buf, BUFSIZ, fp); 
  while
    (p != NULL && *buf == '#');

  if (!p) return -1;   /* Reached the end of the file */

  /* Read values into vbuf */

  for (n = 0, p = strtok(buf, " "); n < MAXVALS && p != NULL; 
       n++,   p = strtok(NULL, " \n"))
    {
      vbuf [n] = atof (p);
    }

  return n;
}

/* Parse command line arguements */

void parse_args (int argc, char *argv[])
{
  char c;

  while (--argc && (*++argv)[0] ==  '-')
    switch (c = *++argv[0]) {
    case 'c':
      col = atoi (*++argv) - 1;
      argc--;
      printf ("# Dependent variable in column %d\n", col+1);
      break;
    case 'x':
      val = atof (*++argv);
      argc--;
      printf ("# Crossing value is %g\n", val);
      break;
    case 'h':
      fputs ("usage: xfreq [-c column] [-x value] [file]\n", stderr);
      exit  (0);
      break;
    default:
      fprintf (stderr, "xfreq: uknown option -- %c\n", c);
      break;
    }

  return;
}
