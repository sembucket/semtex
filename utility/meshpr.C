/*****************************************************************************
 * MESHPR.C:  utility to generate mesh nodes from mesh description file.
 *
 * Usage: meshpr [options] [file]
 *   options:
 *   -h   ... display this message
 *   -v   ... set verbose output
 *   -n N ... override element order to be N
 *
 *****************************************************************************/

// $Id$


#include "Fem.h"


static char prog[] = "meshpr";


static void getargs (int, char **, char **, int *, int*);


int main (int argc, char **argv)
// ---------------------------------------------------------------------------
// From ASCII file named on command line or on stdin, generate mesh node
// information and print up on stdout.
// ---------------------------------------------------------------------------
{
  char      *session = 0;
  int        verb    = 0,
             np      = 0;
  ifstream   file;
  Mesh*      M = new Mesh;
  Field*     F;

  initialize();

  getargs   (argc, argv, &session, &verb, &np);
  setIparam ("VERBOSE", verb);

  file.open (session);
  if (file.bad()) message (prog, "couldn't open session file", ERROR);

  seekBlock   (file, "parameter");
  readOptions (file);     
  readIparams (file);     
  readFparams (file);     
  endBlock    (file);

  seekBlock (file, "mesh");
  file >> *M;
  endBlock (file);

  if (np) setIparam ("N_POLY", np);
  M -> connectSC (iparam ("N_POLY"));

  F = new Field (*M, iparam ("N_POLY"));
  Field::printMesh (F);

  return EXIT_SUCCESS;
}





static void getargs(int     argc    , 
		    char  **argv    ,
		    char  **session ,
		    int    *verb    ,
		    int    *np      )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: meshpr [options] session\n"
                 "options:\n"
                 "  -h   ... display this message\n"
                 "  -v   ... set verbose output\n"
		 "  -n N ... override number of element knots to be N\n";
  char err[StrMax];
  char c;


  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      for (*verb = 1; **++argv == 'v'; ++*verb);
      break;
    case 'n':
      if (*++argv[0])
	*np = atoi (*argv);
      else {
	--argc;
	*np = atoi (*++argv);
      }
      break;
    default:
      sprintf (err, "%s: illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if   (argc == 1) *session = *argv;
  else             message (prog, "must provide session file", ERROR);
}
