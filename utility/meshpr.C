///////////////////////////////////////////////////////////////////////////////
// meshpr.C:  utility to generate mesh nodes from mesh description file.
//
// Copyright (c) 1995--1999 Hugh Blackburn
//
// Prism-compatible output.
//
// Usage: meshpr [options] file
//   options:
//   -h       ... display this message
//   -c       ... disable checking of mesh connectivity
//   -v       ... set verbose output
//   -u       ... set uniform spacing [Default: GLL]
//   -3       ... produce 3D mesh output: Np*Np*Nz*Nel*(x y z)
//   -n <num> ... override element order to be num
//   -z <num> ... override number of planes to be num
//   -b <num> ... override wavenumber beta to be <num> (3D)
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
static void getargs (int, char**, char*&, integer&, integer&,
		     integer&, integer&, integer&, integer&, real&);


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
          threed  = 0,
          np      = 0,
          nz      = 0,
          basis   = GLL;
  real    beta    = -1.;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, verb, check, np, nz, threed, basis, beta);

  // -- Set up to read from file, initialize Femlib parsing.

  FEML feml (session);

  if (verb)            Femlib::value ("VERBOSE", verb);
  if   (np)            Femlib::value ("N_POLY",  np  );
  else  np = (integer) Femlib::value ("N_POLY"       );
  if   (nz)            Femlib::value ("N_Z",     nz  );
  else  nz = (integer) Femlib::value ("N_Z"          );

  if (nz > 1 && beta > 0.0) Femlib::value ("BETA", beta);

  // -- Build mesh from session file information.

  Mesh M (&feml, check);

  // -- Generate mesh knots and print up.

  const integer    NEL  = M.nEl();
  const integer    NTOT = np * np;
  const real       dz   = Femlib::value ("TWOPI/BETA") / nz;
  register integer ID, j, k;
  vector<real>     x (np*np), y (np*np);
  const real*      zero;
  real             z;

  if (!threed)
    cout
      << np  << " "
      << np  << " "
      << nz  << " "
      << NEL << " NR NS NZ NEL"<< endl;

  Femlib::mesh (basis, basis, np, np, &zero, 0, 0, 0, 0);

  if (threed) {

    // -- Print out x, y, z for every mesh location, in planes.

    nz = (nz > 1) ? nz : 0;
    for (k = 0; k <= nz; k++) {
      z = k * dz;
      for (ID = 0; ID < NEL; ID++) {
	M.meshElmt (ID, np, zero, x(), y());
	for (j = 0; j < NTOT; j++)
	  cout << x(j) << '\t' << y(j) << '\t' << z << endl;
      }
    }

  } else {
   
    // -- Print out x-y mesh.

    for (ID = 0; ID < NEL; ID++) {
      M.meshElmt (ID, np, zero, x(), y());
      for (j = 0; j < NTOT; j++)
	cout << setw(15) << x(j) << setw(15) << y(j) << endl;
    }
  
    // -- Print out z-mesh.
    
    if (nz > 1) for (j = 0; j <= nz; j++) cout << setw(15) << j * dz << endl;
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
		     integer& threed ,
		     integer& basis  ,
		     real&    beta   )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: meshpr [options] session\n"
    "options:\n"
    "  -h       ... display this message\n"
    "  -c       ... disable checking of mesh connectivity\n"
    "  -v       ... set verbose output\n"
    "  -u       ... set uniform spacing [Default: GLL]\n"
    "  -3       ... produce 3D mesh output: Np*Np*Nz*Nel*(x y z)\n"
    "  -n <num> ... override number of element knots to be num\n"
    "  -z <num> ... override number of planes to be num\n"
    "  -b <num> ... override wavenumber beta to be <num> (3D)\n";
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
    case 'b':
      if (*++argv[0]) beta = atof (*argv);
      else { --argc;  beta = atof (*++argv); }
      break;
    case 'c':
      check = 0;
      break;
    case 'u':
      basis = STD;
      break;
    case '3':
      threed = 1;
      break;
    case 'n':
      if (*++argv[0]) np = atoi (*argv);
      else { --argc;  np = atoi (*++argv); }
      break;
    case 'z':
      if (*++argv[0]) nz = atoi (*argv);
      else { --argc;  nz = atoi (*++argv); }
      break;
    default:
      sprintf (err, "%s: illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if   (argc == 1) session = *argv;
  else             message (prog, "must provide session file", ERROR);
}
