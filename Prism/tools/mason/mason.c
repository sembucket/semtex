/*
 * MASON -- Mortar Space Mapping Program
 *
 * This program builds the connectivity information for a spectral element
 * mesh.  The elements can be either three- (triangles) or four-sided
 * (quadrilaterals), based on the flag -tri.  The default is for four-sided
 * elements.
 *
 * The input is a Prism mesh file, normally called a ".rea" file.  Output is 
 * an ASCII file containing element/face numbering information.
 *
 * The program performs three types of bandwidth optimization based on the -O 
 * flag.  The default is a single pass of the Reverse Cuthill-McKee algorithm.
 * With -O2 it will search all boundary nodes to find the best "root" for the
 * numbering algorithm.  The -O3 option runs a greedy algorithm that can some-
 * times outperform RCM for meshes that contain internal "holes".
 *
 *
 * RCS Information
 * ------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "mason.h"

char *prog  = "mason";
char *author= "Ron Henderson";
char *usage = "usage: mason [options] file[.rea]\n";
char *help  = "options:\n"
  "-h      ... print this help string\n"
  "-v      ... be verbose about things\n"
  "-n #    ... set the N-order to #\n"
  "-o file ... send output to the named file instead of stdout\n"
  "-tri    ... build connectivity info for a triangular mesh\n"
  "-rand   ... randomize the numbering scheme (worst case)\n"
  "-number ... generate patch numbering info\n"
  "-O0     ... turn off all optimizations\n"
  "-O1,-O  ... Reverse Cuthill-McKee algorithm in single-pass mode [default]\n"
  "-O2     ... RCM w/searching for the best starting node\n"
  "-O3     ... Greedy algorithm w/searching\n";

char *rcsid = "$Revision$";

int oplevel = 1;    /* Requested optimization level */
int opstart = 1;    /* Requested starting element   */
int rflag   = 0;    /* Randomize number system ?    */
int number  = 0;    /* Numbering flag for patching  */

typedef struct {
  char  *name;
  FILE  *fp;
} File;

static struct { 
  File   in ;
  File   out;
} FileList = 
/* Changed hmb Jan 2002 
  { { "<stdin>" , stdin  }, 
  { "<stdout>", stdout } }; */
  { { "<stdin>" , NULL  }, 
    { "<stdout>", NULL } };

int Patch_number (File);

/* ------------------------------------------------------------------------ */

int main (int argc, char *argv[])
{
  Domain *omega;

  FileList.in.fp  = stdin;
  FileList.out.fp = stdout;

  omega = parse_args (argc, argv);

  if (number) Patch_number (FileList.in);

  loadParams (omega, FileList.in.fp);  /* ........ Setup ......... */
  loadMesh   (omega, FileList.in.fp);
  loadBCs    (omega, FileList.in.fp);

  linkEdges  (omega);                  /* ...... Processing ...... */
  linkVertex (omega);
  linkPatches(omega);
  linkNodes  (omega);

  output     (omega, FileList.out.fp); /* ........ Output ........ */

  if (omega->verbose) {
    showPatch (omega);
    showEdges (omega);
  }

  return 0;
}

/* ------------------------------------------------------------------------ *
 * parse_args() -- Parse application arguments                              *
 *                                                                          *
 * See the usage/help strings for supported options.                        *
 * ------------------------------------------------------------------------ */

Domain *parse_args (int argc, char *argv[])
{
  Domain *new;
  char  c, fname[FILENAME_MAX], buf[BUFSIZ];
  
  new = (Domain *) calloc (1, sizeof(Domain));

  while (--argc && (*++argv)[0] == '-') {  /* start parsing options ... */

    if (strcmp (*argv+1,"tri") == 0) {     /* check for the -tri option */
      new->tri = 1;
      continue;
    }

    if (strcmp (*argv+1,"rand") == 0) {    /* check for the -rand option */
      rflag   =  1;
      oplevel = -1;
      continue;
    }

    if (strcmp (*argv+1,"number") == 0) {  /* generate patch numbering */
      number = 1;
      continue;
    }

    switch (c = *++argv[0]) {
    case 'h':                               /* print the help string */
      fputs(usage, stderr);
      fputs(help,  stderr);
      exit (0);
      break;

    case 'n':
      new->norder = (*++argv[0]) ? atoi (*argv) : (argc--, atoi(*++argv));
      break;

    case 'O':
      oplevel = (*++argv[0]) ? atoi (*argv) : 1;
      break;

    case 'o':                         /* direct the output to a file */
      if (*++argv[0])
	strcpy (fname, *argv);
      else {
	strcpy (fname, *++argv);
	argc--;
      }

      if (fname[0] != '-') {
	if ((FileList.out.fp = fopen(fname,"w")) == (FILE*) NULL) {
	  fprintf(stderr, "%s: unable to open the output file -- %s\n", 
		  prog, fname);
	  exit(1);
	}
	FileList.out.name = strdup(fname);
      } else 
	FileList.out.fp   = stdout;
      break;

    case 'v': {               /* be verbose and echo the version */
      float v;
      sscanf (rcsid, "%*s%f", &v);
      fprintf(stderr, "%s v%g -- by %s\n", prog, v, author);
      new->verbose++;
      break;
    }
      
    default: {
      sprintf (buf, "unknown option -- %c", c);
      error_msg(buf);
      break;
    }
    }
  }

  /* Open the input file */

  if (argc) {
    if (strstr (*argv, ".rea"))
      strcpy (fname, *argv);
    else
      sprintf(fname, "%s.rea", *argv);

    if (!(FileList.in.fp = fopen (fname, "r"))) {
      if (!(FileList.in.fp = fopen (*argv,"r"))) {
	sprintf (buf, "unable to open the input file -- %s or %s", 
		 *argv, fname);
	error_msg(buf);
      } else
	strcpy (fname, *argv);
    }
    FileList.in.name = strdup (fname);
  } else {
    fputs (usage, stderr);
    exit  (1);
  }


  /* Open the output (.mor) file */

  if (!FileList.out.fp) {
    char *p;
    strcpy (fname, FileList.in.name);
    if ((p = strstr(fname,".rea")))
      strcpy (p, ".mor");
    else 
      strcat (fname, ".mor");
    if (!(FileList.out.fp = fopen(fname,"w"))) {
      sprintf (buf, "unable to open the output file -- %s\n", fname);
      error_msg(buf);
    }
    FileList.out.name = strdup (fname);
  }


  if (new->verbose) {
    fprintf (stderr, "%s: in = %s, out = %s", prog,
	     FileList.in.name, FileList.out.name);
    if (new->norder)
      fprintf (stderr, ", N-order = %d", new->norder);
    if (oplevel)
      fprintf (stderr, "\n\tOptimization Level = %d, starting element = %d", 
	       oplevel, opstart);
    putc ('\n', stderr);
  }

  return new;
}

static 
char   _error_buf [128];
char   *error_buf = _error_buf;

void error_msg (char *msg)
{
  fprintf (stderr, "%s: %s\n", prog, msg);
  exit (1);
}

#include <unistd.h>

int Patch_number (File orig)
{
  char *tmp = tmpnam(NULL);
  char  buf [BUFSIZ];
  FILE *fp;

  static char* ident = "mason: number";

  if ((fp = fopen(tmp, "w+")) == (FILE*) NULL)
    perror(ident);

  /* Copy the input file to "tmp" */

  rewind (orig.fp);
  while (fgets (buf, BUFSIZ, orig.fp))
    fputs (buf, fp);
  rewind (fp);
  fclose (orig.fp);

  sprintf (buf, "number %s > /dev/null", tmp);
  if (system(buf) != 0)
    perror(ident);

  sprintf (buf, "sed -f rea.sed < %s > %s", tmp, orig.name);
  if (system(buf) != 0)
    perror (ident);

  fclose  (fp);
  unlink  (tmp);
  unlink  ("rea.sed");

  if ((orig.fp = fopen (orig.name, "r")) == (FILE*) NULL)
    perror (ident);

  return 0;
}
