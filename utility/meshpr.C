///////////////////////////////////////////////////////////////////////////////
// meshpr.C:  utility to generate mesh nodes from mesh description file.
//
// Copyright (c) 1995--1999 Hugh Blackburn
//
// Prism-compatible output.
//
// Usage: meshpr [options] file
//   options:
//   -h   ... display this message
//   -c   ... disable checking of mesh connectivity
//   -v   ... set verbose output
//   -u   ... set uniform spacing [Default: GLL]
//   -n N ... override element order to be N
//   -z N ... override number of planes to be N
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

#include <femdef.h>
#include <Femlib.h>
#include <Utility.h>
#include <Mesh.h>

static char prog[] = "meshpr";
static void getargs (int, char**, char*&, integer&,
		     integer&, integer&, integer&, integer&);


int main (int     argc,
	  char**  argv)
// ---------------------------------------------------------------------------
// From FEML file named on command line, generate mesh knot
// information and print up on standard output.
// ---------------------------------------------------------------------------
{
  // -- Set defaults & parse command line.

  char*   session = 0;
  integer verb    = 0,
          check   = 1,
          np      = 0,
          nz      = 0,
          basis   = GLL;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, verb, check, np, nz, basis);

  // -- Set up to read from file, initialize Femlib parsing.

  FEML feml (session);

  if (verb)            Femlib::value ("VERBOSE", verb);
  if   (np)            Femlib::value ("N_POLY",  np  );
  else  np = (integer) Femlib::value ("N_POLY"       );
  if   (nz)            Femlib::value ("N_Z",     nz  );
  else  nz = (integer) Femlib::value ("N_Z"          );

  // -- Build mesh from session file information.

  Mesh M (&feml, check);

  // -- Generate mesh knots and print up.

  const integer    NEL  = M.nEl();
  const integer    NTOT = np * np;
  register integer ID, j;
  vector<real>     x (np*np), y (np*np);
  const real*      z;

  cout << np << " " << np << " " << nz << " " << NEL << " NR NS NZ NEL"<< endl;

  Femlib::mesh (basis, basis, np, np, &z, 0, 0, 0, 0);

  // -- Print out x-y mesh.

  for (ID = 0; ID < NEL; ID++) {
    M.meshElmt (ID, np, z, x(), y());

    for (j = 0; j < NTOT; j++)
      cout << setw (15) << x (j) << setw (15) << y (j) << endl;
  }
  
  // -- Print out z-mesh.

  if (nz > 1) {
    const real dz = Femlib::value ("TWOPI/BETA") / nz;
    for (j = 0; j <= nz; j++)
      cout << setw (15) << j * dz << endl;
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (integer  argc   , 
		     char**   argv   ,
		     char*&   session,
		     integer& verb   ,
		     integer& check  ,
		     integer& np     ,
		     integer& nz     ,
		     integer& basis  )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: meshpr [options] session\n"
                 "options:\n"
                 "  -h   ... display this message\n"
                 "  -c   ... disable checking of mesh connectivity\n"
                 "  -v   ... set verbose output\n"
                 "  -u   ... set uniform spacing [Default: GLL]\n"
		 "  -n N ... override number of element knots to be N\n"
                 "  -z N ... override number of planes to be N\n";
  char err[StrMax], c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      for (verb = 1; *++argv[0] == 'v'; verb++);
      break;
    case 'c':
      check = 0;
      break;
    case 'u':
      basis = STD;
      break;
    case 'n':
      if (*++argv[0])
	np = atoi (*argv);
      else {
	--argc;
	np = atoi (*++argv);
      }
      break;
    case 'z':
      if (*++argv[0])
	nz = atoi (*argv);
      else {
	--argc;
	nz = atoi (*++argv);
      }
      break;
    default:
      sprintf (err, "%s: illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if   (argc == 1) session = *argv;
  else             message (prog, "must provide session file", ERROR);
}
