/*****************************************************************************
 * chop: read an input file and reproduce a specified number of lines on
 * standard output.
 *
 * Usage: chop [-h] [-s startline] [-n number of lines] [file]
 *
 * The two command line arguments specify the first line of the input file
 * to reproduce, and the number of subsequent lines.  Can be used as a
 * filter.  If number of lines not specified, read through until EOF.
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

static void chop (FILE* fp, int s, int n);


int main (int argc, char *argv[])
/* ------------------------------------------------------------------------- *
 * Wrapper for chop(), which does the work.  Here we do administration.
 * ------------------------------------------------------------------------- */
{
  static char usage[] = 
    "usage: chop [options] [input]\n"
    "[options]:\n"
    "-h          ... display this message\n"
    "-n <lines>  ... reproduce this many lines of file  [Default: to EOF]\n"
    "-s <line>   ... start at this line number          [Default: 1]\n";
  int   i, start=1, nlines=0;
  FILE *fp;
  
  while (argc > 1 && argv[1][0] == '-') {
    switch (argv[1][1]) {
    case 's':		             /* -s: start at line of next argv.      */
      argc--; argv++;
      start = atoi(argv[1]);
      break;
    case 'n':		             /* -n: reproduce n lines of input file. */
      argc--; argv++;
      nlines = atoi(argv[1]);
      break;
    case 'h':
      fprintf(stderr, usage);
      exit(0);
      break;
    default:
      fprintf(stderr, "%s: unknown arg %s\n", argv[0], argv[1]);
      exit(1);
    }
    argc--;
    argv++;
  }

  if (argc == 1) {	      /* no file was specified as input, read stdin */
    chop(stdin, start, nlines);
    while((i = getchar()) != EOF);			/* mop up remainder */
  } else
    for (i = 1; i < argc; i++)
      if ((fp=fopen(argv[i], "r")) == NULL) {
	fprintf(stderr, "chop: can't open %s\n",argv[i]);
	fprintf(stderr, usage);
	exit(1);
      } else {
	chop(fp, start, nlines);
	fclose(fp);
      }
  
  return EXIT_SUCCESS;
}


static void chop (FILE *fp, int s, int n)
/* ------------------------------------------------------------------------- *
 * This does the real work.
 * ------------------------------------------------------------------------- */
{
  int c;
  
  if (s > 1) {
    while (s > 1 && (c = getc(fp)) != EOF) if (c=='\n') s--;
    if (c == EOF) {
      fprintf(stderr,"chop: Reached EOF before start line.\n");
      exit(1); 
    }
  }

  if (n > 0) {
    while ((c = getc(fp)) != EOF && n > 0) {
      putchar(c);
      if (c == '\n') n--;
    }
    if (n != 0) {
      fprintf(stderr,"chop: Reached EOF prematurely.\n") ;
      exit(1);
    }
  }
  else while ((c = getc(fp)) != EOF) putchar(c);     /* Read thru until EOF. */
}
