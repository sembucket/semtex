/*****************************************************************************
 * demean.c: a minimal utility to remove mean values from columnar inputs.
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define COLMAX 32
#define STRMAX 1024

int main ()
{
  char   buf  [STRMAX], sep[] = " \t", *tok;
  FILE*  store[COLMAX];
  double mean [COLMAX];
  int    numcols, numrows, i, j;
  double datum;

  gets (buf);
  tok     = strtok (buf, sep);
  numcols = 1;
  numrows = 1;
  mean[0] = atof (tok);
  while (tok = strtok (NULL, sep)) mean[numcols++] = atof (tok);
  
  for (i = 0; i < numcols; i++) {
    store[i] = tmpfile();
    fwrite (mean + i, sizeof (double), 1, store[i]);
  }

  while (gets (buf)) {
    mean[0] += datum = atof (strtok (buf, sep));
    fwrite (&datum, sizeof (double), 1, store[0]);
    for (i = 1; i < numcols; i++) {
      mean[i] += datum = atof (strtok (NULL, sep)); 
      fwrite (&datum, sizeof (double), 1, store[i]);
    }
    numrows++;
  }

  for (i = 0; i < numcols; i++) {
    mean[i] /= numrows;
    rewind (store[i]);
  }

  for (j = 0; j < numrows; j++) {
    for (i = 0; i < numcols; i++) {
      fread (&datum, sizeof (double), 1, store[i]);
      datum -= mean[i];
      fprintf (stdout, " %.12e", datum);
    }
    fprintf (stdout, "\n");
  }

  return EXIT_SUCCESS;
}
