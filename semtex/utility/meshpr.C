///////////////////////////////////////////////////////////////////////////////
// meshpr.C:  utility to generate mesh nodes from mesh description file.
//
// Prism-compatible output.
//
// Usage: meshpr [options] file
//   options:
//   -h   ... display this message
//   -v   ... set verbose output
//   -u   ... set uniform spacing [Default: GLL]
//   -n N ... override element order to be N
//
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

#include <femdef.h>
#include <Femlib.h>
#include <Utility.h>
#include <Mesh.h>

static char prog[] = "meshpr";
static void getargs (integer, char**, char*&, integer&, integer&, integer&);


integer main (integer argc,
	      char**  argv)
// ---------------------------------------------------------------------------
// From FEML file named on command line, generate mesh knot
// information and print up on standard output.
// ---------------------------------------------------------------------------
{
  // -- Set defaults & parse command line.

  char*   session = 0;
  integer verb    = 0,
          np      = 0,
          basis   = GLL;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, verb, np, basis);

  // -- Set up to read from file, initialize Femlib parsing.

  FEML feml (session);

  if (verb)            Femlib::value ("VERBOSE", verb);
  if   (np)            Femlib::value ("N_POLY", np);
  else  np = (integer) Femlib::value ("N_POLY");

  // -- Build mesh from session file information.

  Mesh M (feml);

  // -- Generate mesh knots and print up.

  const integer    NEL  = M.nEl();
  const integer    NTOT = np * np;
  const integer    NZ   = (integer) Femlib::value ("N_Z");
  register integer ID, j;
  vector<real>     x (np*np), y (np*np);
  const real*      z;

  cout << np << " " << np << " " << NZ << " " << NEL << " NR NS NZ NEL"<< endl;

  Femlib::mesh (basis, basis, np, np, &z, 0, 0, 0, 0);

  // -- Print out x-y mesh.

  for (ID = 0; ID < NEL; ID++) {
    M.meshElmt (ID, np, z, x(), y());

    for (j = 0; j < NTOT; j++)
      cout << setw (15) << x (j) << setw (15) << y (j) << endl;
  }
  
  // -- Print out z-mesh.

  if (NZ > 1) {
    const real dz = Femlib::value ("TWOPI/BETA") / NZ;
    for (j = 0; j <= NZ; j++)
      cout << setw (15) << j * dz << endl;
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (integer  argc   , 
		     char**   argv   ,
		     char*&   session,
		     integer& verb   ,
		     integer& np     ,
		     integer& basis  )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: meshpr [options] session\n"
                 "options:\n"
                 "  -h   ... display this message\n"
                 "  -v   ... set verbose output\n"
                 "  -u   ... set uniform spacing [Default: GLL]\n"
		 "  -n N ... override number of element knots to be N\n";
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
    default:
      sprintf (err, "%s: illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if   (argc == 1) session = *argv;
  else             message (prog, "must provide session file", ERROR);
}
