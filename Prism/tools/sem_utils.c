/*
 * Generic utility functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "sem_utils.h"

/* Externals */

static char *genargs = 
  "options:\n"
  "-h      ... print this help message\n"
  "-v      ... be verbose about things\n";

/* ---------------------------------------------------------------------- *
 * generic_args() -- Generic command line arguments for all utilities     *
 *                                                                        *
 * This function processes the command line arguments common to all of    *
 * the utility programs.  Any arguments not processed by this function    *
 * are returned in an updated argv[].  It NEVER processes the last arg-   *
 * ument since that's reserved for the input file name.                   *
 *                                                                        *
 * The generic options are the following:                                 *
 *                                                                        *
 *     -h      ...  print a help message                                  *
 *     -v      ...  be verbose about things                               *
 *     -n #    ...  specify an N-order to interpolate to                  *
 *     -z #    ...  specify a number of NZ-planes to interpolate to       *
 *     -o file ...  send output to the named file                         *
 *     -m file ...  read the mesh from the named file                     *
 *     -r file ...  read the session (.rea) from the named file           *
 *                                                                        *
 * Your application doesn't have to support all of these options, but it  *
 * should list the ones it does in the external character array "help",   *
 * except for the -h and -v options (they're automatic).  If you're ap-   *
 * plication DOES support any of these features, it MUST use these argu-  *
 * ments.  Likewise, you shouldn't use these symbols for anything else    *
 * since that makes things more confusing for everybody.                  *
 *                                                                        *
 * Return value: number of arguments passed through to the application    *
 * ---------------------------------------------------------------------- */

int generic_args (int argc, char *argv[], FileList *f)
{
  int   appargc = 0;
  char *appargv [MAXARGS], **orig;
  char  c, fname[FILENAME_MAX];
  int i;
  
  /* Set up the default values for each file pointer */

  memset (f, '\0', sizeof(FileList));

  f -> in.fp    =  stdin;
  f -> out.fp   =  stdout;
  f -> mesh.fp  =  stdin;
  
  orig = argv;          /* save the orignal argument vector    */
  manager_init();       /* initialize the symbol table manager */

  /* Mark the following parameters with the UNSET flag. Otherwise, *
   * the parser will generate an undefined variable error.         */

  iparam_set("NORDER-req", UNSET);
  iparam_set("NZ-req",     UNSET);
  iparam_set("DUMP-req",   UNSET);

  while (--argc && appargc < MAXARGS)       /* start parsing options ... */

    if ((*++argv)[0] == '-' && strlen(*argv) > 1) {
      switch (c = *++argv[0]) {
      case 'h':                               /* print the help string */
	fputs(usage,   stderr);
	fputs(genargs, stderr);
	fputs(help,    stderr);
	exit (0);
	break;
      case 'n':                        /* Queue an NR|NS-interpolation */
	if (*++argv[0]) 
	  iparam_set("NORDER-req", atoi(*argv));
	else {
	  iparam_set("NORDER-req", atoi(*++argv));
	  argc--;
	}
	break;
      case 'd':                                /* Queue a dump to read */
	if (*++argv[0]) 
	  iparam_set("DUMP-req", atoi(*argv));
	else {
	  iparam_set("DUMP-req", atoi(*++argv));
	  argc--;
	}
	break;
      case 'z':                           /* Queue an NZ interpolation */
	if (*++argv[0]) 
	  iparam_set("NZ-req", atoi(*argv));
	else {
	  iparam_set("NZ-req", atoi(*++argv));
	  argc--;
	}
	break;

      case 'm':                           /* read the mesh from a file */
	if (*++argv[0])
	  strcpy(fname, *argv);
	else {
	  strcpy(fname, *++argv);
	  argc--;
	}
	if (!(f->mesh.fp = fopen(fname,"r")))
	  if (!(f->mesh.fp = fopen(strcat (fname, ".mesh"), "r"))) {
	    fprintf(stderr, "%s: unable to open the mesh file -- %s or %s\n", 
		    prog, *argv, fname);
	    exit(1);
	  }
	f->mesh.name = strdup(fname);
	break;

      case 'r':                           /* read the session file */
	if (*++argv[0])
	  strcpy(fname, *argv);
	else {
	  strcpy(fname, *++argv);
	  argc--;
	}
	if (!(f->rea.fp = fopen(fname,"r")))
	  if (!(f->rea.fp = fopen(strcat(fname,".rea"),"r"))) {
	    fprintf(stderr,"%s: unable to open the session file -- %s or %s\n",
		    prog, *argv, fname);
	    exit(1);
	  }
	f->rea.name = strdup(fname);

	/* Try to open the connectivity file.  This file is usually *
	 * optional, so the application is responsible for error    *
	 * checking.                                                */
	
	strcpy(fname + strlen(fname)-3, "mor");
	f->mor.fp   = fopen (fname,"r");
	f->mor.name = strdup(fname);
	break;

      case 'o':                         /* direct the output to a file */
	if (*++argv[0])
	  strcpy(fname, *argv);
	else {
	  strcpy(fname, *++argv);
	  argc--;
	}
	if ((f->out.fp = fopen(fname,"w")) == (FILE*) NULL) {
	  fprintf(stderr, "%s: unable to open the output file -- %s\n", 
		  prog, fname);
	  exit(1);
	}
	f->out.name = strdup(fname);
	break;

      case 'v': {               /* be verbose and echo the version */
	float v;
	sscanf (rcsid, "%*s%f", &v);
	fprintf(stderr, "%s v%g -- by %s\n", prog, v, author);
	if (isdigit (*++argv[0]))
	  option_set ("verbose", atoi(*argv));
	else
	  option_set("verbose", 1);
	break;
      }

      default:                  /* transfer it to the un-processed list */
	sprintf 
	  ( appargv[appargc++] = (char *) malloc(strlen(*argv)+2),
	    "-%s", *argv );
	break;
      }

      /* Current argv[] is not an option, so just copy it to the *
       * application list to be processed later.                 */
      
    } else {
      strcpy (appargv[appargc++] = (char *) malloc (strlen(*argv)+1), *argv);
    }
  
  /* Now assign the application's arguments */

  while ((i = argc++) < appargc) orig[i] = appargv[i];
  
  return appargc;
}

void error_msg (char *msg)
{
  fprintf (stderr, "%s: %s\n", prog, msg);
  exit (1);
}


