/*
 * ff -- field filter
 *
 * $Revision$
 *
 * Author: R. D. Henderson
 *
 * This is a tool for converting Prism field files.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

#include "speclib/speclib.h"
#include "speclib/fieldfile.h"

static char *prog  = "ff";
static char *usage = "[options] input > output";
static char  help[] = 
"options:\n"
"-h          ... print this help message\n"
"-v          ... print diagnostics\n"
"-o fname    ... send output to a named file\n"
"-d #        ... only write the specified field dump number\n"
"-t #        ... reset the value of \"Time\" in the file\n"
"-n #        ... interpolate data to a different polynomial order\n"
"-m #        ... interpolate data to a different number of Fourier modes\n"
"-64         ... write the output in 64-bit binary\n"
"-32         ... write the output in 32-bit binary\n"
"-ascii      ... write the output in text (ascii) format\n"
"-swap       ... byte-swap the data\n"
"\n"
"This program applies varies filters to a FieldFile, including data convers-\n"
"ions and simple types of intepolation.  It converts Prism field files so   \n"
"that binary files can be ported to various platforms.  You can convert     \n"
"files to 32 bit format to reduce the amount of data that needs to be       \n"
"You can convert files to text (ascii) format to export them to platforms   \n"
"that aren't binary compatible.\n"
"\n"
"To convert a file as part of a data processing pipeline, specify the input \n"
"file name as '-' and it will be read from stdin.\n"
"\n"
"The program converts all dumps in a file by default.  Use -d to convert a  \n"
"specific dump, i.e. -d 1 for the first, -d 2 for the second.  Use -d -1 to \n"
"convert only the last dump in the file.\n";

static int   verbose = 0;
static int   dump    = 0;
static int   format  = BINARY;
static float time    = -1.;
/* static FILE* output  = stdout; Changed hmb Jan 2002 */
static FILE* output  = NULL;

/* ------------------------------------------------------------------------- */

void parse_args (int argc, char *argv[])
{
  int n;

  if (argc == 1) {
    fprintf (stderr, "usage: %s %s\n", prog, usage);
    exit (0);
  }

  for (n = 1; n < argc; n++) {
    if (*argv[n] == '-' && strlen(argv[n]) > 1) {
      const char c = argv[n][1];
      switch (c) {
      case 'h':
	fprintf (stderr, "usage: %s %s\n%s\n", prog, usage, help);
	exit(0);
      case 'd':
	dump = atoi(argv[++n]);
	break;
      case 't':
	time = atof(argv[++n]);
	option_set("time",1);
	break;
      case 'm':
	option_set("nmodes", atoi(argv[++n]));
	break;
      case 'n':
	option_set("norder", atoi(argv[++n]));
	break;
      case 'v':
	verbose = 1;
	break;
      case '6':
	format = BINARY;
	break;
      case '3':
	format = PACKED;
	break;
      case 'a':
	format = ASCII;
	break;
      case 's':
	option_set("swap",1);
	break;
      case 'o':
	if (!(output = fopen(argv[++n],"w"))) {
	  perror("ff");
	  exit (-1);
	}
	break;
      default:
	fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
    }
  }
}

void DoubleReverse(int len, double *buf)
{
  union f {
    double x;
    char c[8];
  } *a,b;
  int i;
  a=(union f *)buf;
  for(i=0;i<len;i++) {
    b.c[0]=a->c[7];
    b.c[1]=a->c[6];
    b.c[2]=a->c[5];
    b.c[3]=a->c[4];
    b.c[4]=a->c[3];
    b.c[5]=a->c[2];
    b.c[6]=a->c[1];
    b.c[7]=a->c[0];
    a->x=b.x;
    a++;
  }
}

void swap (FieldFile *f)
{
  const int n   = FIELDFILE_COUNT(f);
  const int len = 
    FIELDFILE_NR(f)*FIELDFILE_NS(f)*FIELDFILE_NZ(f)*FIELDFILE_NEL(f);
  int i;
  for (i = 0; i < n; i++)
    DoubleReverse (len, FIELDFILE_DATA(f,i));
}

/* ------------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  FILE *fp  = stdin;
  int ndump = 0;
  FieldFile *f;

  output  = stdout;

  speclib_init(); 

  parse_args (argc, argv);

  if (*argv[argc-1] != '-') {
    if (!(fp = fopen(argv[argc-1],"r"))) {
      char fname[FILENAME_MAX];
      sprintf ("%s.fld", argv[argc-1]);
      if (!(fp = fopen(fname,"r")))
	fprintf (stderr, "%s: unable to open %s or %s\n", prog,
		 argv[argc-1], fname);
    }
  }

  f = FieldFile_alloc();

  while (FieldFile_read(f,fp) != FIELDFILE_EOF) {
    ndump++;

    if (option("swap")) 
      swap(f);

    if (option("nmodes"))
      FieldFile_projectz (f, option("nmodes"));
    if (option("norder")) 
      FieldFile_project (f, option("norder"), option("norder"));
    
    switch (format) {
    case BINARY:
      FieldFile_setFormat(f,BINARY);
      break;
    case PACKED:
      FieldFile_setFormat(f,PACKED);
      break;
    case ASCII:
      FieldFile_setFormat(f,ASCII);
      break;
    default:
      break;
    }

    if (option("time"))
      FIELDFILE_TIME(f) = time;
    
    if (dump == 0 || dump == ndump) {
      if (verbose)
	fprintf (stderr, "%s: writing dump %d\n", prog, ndump);
      FieldFile_write(f, output);
    }
  }

  if (dump == -1) {
    if (verbose)
      fprintf (stderr, "%s: writing dump %d\n", prog, ndump);
    FieldFile_write(f, output);
  }

  FieldFile_free(f);
  return 0;
}

