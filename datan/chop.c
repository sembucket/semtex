/*****************************************************************************
 * chop: read an input file and reproduce a specified number of lines on
 * standard output.
 *
 * Usage: chop [-h] [-s nstart] [-n nlines] [-S nskip] [-B nblock] [file]
 *
 * The first two command line arguments specify the first line of the
 * input file to reproduce, and the number of subsequent lines.  The
 * third argument gives a skip between lines of output that are
 * reproduced.
 *
 * Can be used as a filter.  If nlines not specified, read through
 * until EOF.  Lines are assumed to be BUFSIZ characters long at most.
 * If nblock is set, then data are output in blocks of nblock lines,
 * with nskip lines omitted between each block.
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>


static void chop (FILE* fp, int s, int n, int S, int B);


int main (int    argc, 
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Wrapper for chop(), which does the work.  Here we do administration.
 * ------------------------------------------------------------------------- */
{
  static char usage[] = 
    "usage: chop [options] [input]\n"
    "  options:\n"
    "  -h         ... display this message\n"
    "  -n <lines> ... reproduce this many lines of file [Default: to EOF]\n"
    "  -s <line>  ... start at this line number              [Default: 1]\n"
    "  -S <num>   ... skip <num> lines between each output   [Default: 1]\n"
    "  -B <num>   ... output blocks size <num>, skip between [Default: 0]\n";
  int   i, nstart = 1, nlines = 0, nskip = 1, nblock = 0;
  FILE* fp;
  
  while (--argc && **++argv == '-') {
    switch (*++argv[0]) {
    case 'h':
      fprintf (stderr, usage);
      exit    (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0])
	nlines = atoi (*argv);
      else {
	--argc;
	nlines = atoi (*++argv);
      }
      break;
    case 's':
      if (*++argv[0])
	nstart = atoi (*argv);
      else {
	--argc;
	nstart = atoi (*++argv);
      }
      break;
    case 'S':
      if (*++argv[0])
	nskip = atoi (*argv);
      else {
	--argc;
	nskip = atoi (*++argv);
      }
      break;
    case 'B':
      if (*++argv[0])
	nblock = atoi (*argv);
      else {
	--argc;
	nblock = atoi (*++argv);
      }
      break;
    default:
      fprintf (stderr, "chop: unknown arg %s\n", *argv);
      exit    (EXIT_FAILURE);
    }
  }

  if (argc) {			/* -- Input from list of named files. */
    for (i = 0; i < argc; i++)
      if (!(fp = fopen (argv[i], "r"))) {
	fprintf (stderr, "chop: can't open %s\n",argv[i]);
	fprintf (stderr, usage);
	exit    (EXIT_FAILURE);
      } else {
	chop    (fp, nstart, nlines, nskip, nblock);
	fclose  (fp);
      }

  } else {			/* -- Input from stdin. */
    chop  (stdin, nstart, nlines, nskip, nblock);
    while ((i = getchar()) != EOF);
  }
  
  return EXIT_SUCCESS;
}


static void chop (FILE* file  , 
		  int   nstart,
		  int   nline ,
		  int   nskip ,
		  int   nblock)
/* ------------------------------------------------------------------------- *
 * This does the real work.
 * ------------------------------------------------------------------------- */
{
  char buf[BUFSIZ], *OK;
  int  n = 0;

  if (nstart > 1)
    while (nstart > 1 && (OK = fgets (buf, BUFSIZ, file))) {
      if (!OK) {
	fprintf (stderr, "chop: Reached EOF before start line.\n");
	exit    (EXIT_FAILURE); 
      }
      nstart--;
    }
  
  if (nblock) {
    const int quant = nblock + nskip;

    while (fgets (buf, BUFSIZ, file)) {
      if ((n % quant) < nblock) fputs (buf, stdout);
      n++;
    }

  } else

    if (nline) {
      while (nline && (OK = fgets (buf, BUFSIZ, file))) {
	if (!OK) {
	  fprintf (stderr, "chop: Reached EOF prematurely.\n") ;
	  exit    (EXIT_FAILURE);
	}
	if (n % nskip == 0) {
	  fputs (buf, stdout);
	  nline--;
	}
	n++;
      }

    } else
      while (fgets (buf, BUFSIZ, file)) {
	if (n % nskip == 0) fputs (buf, stdout);
	n++;
      }
}




