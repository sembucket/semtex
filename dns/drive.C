//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element DNS for incompressible flows.
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// USAGE:
// -----
// dns [options] session
//   options:
//   -h       ... print usage prompt
//   -i[i]    ... use iterative solver for viscous [and pressure] steps
//   -v[v...] ... increase verbosity level
//   -chk     ... checkpoint field dumps
//   -O <num> ... set numbering scheme optimisation level, Default=3
//
// AUTHOR:
// ------
// Hugh Blackburn
// CSIRO
// PO Box 56
// Highett, Vic 3190
// Australia
// hugh.blackburn@csiro.au
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "dns_h"

static char prog[] = "dns";
static void getargs    (int, char**, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&);

void integrateNS (Domain*, DNSAnalyser*, Flowrate*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char*            session;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;
  DNSAnalyser*     analyst;
  Flowrate*        discharge = 0;
  
  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session);

  preprocess (session, file, mesh, elmt, bman, bsys, domain);

  analyst = new DNSAnalyser (domain, file);
  
  if (file -> seek ("CUT")) discharge = new Flowrate (domain, file);

  domain -> restart();

  ROOTONLY domain -> report();
  
  integrateNS (domain, analyst, discharge);

  Femlib::finalize();

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -i[i]    ... use iterative solver for viscous [& pressure] steps\n"
    "  -v[v...] ... increase verbosity level\n"
    "  -chk     ... checkpoint field dumps\n"
    "  -O <num> ... set numbering scheme optimisation level, Default=3\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      do
	Femlib::value ("ITERATIVE",
		       static_cast<int>(Femlib::value ("ITERATIVE") + 1));
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE",
		       static_cast<int>(Femlib::value ("VERBOSE") + 1));
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	Femlib::value ("CHKPOINT", 
		       static_cast<int>(1));
      } else {
	fprintf (stdout, usage, prog);
	exit (EXIT_FAILURE);	  
      }
      break;
    case 'O': {
      int level;
      if (*++argv[0]) level = atoi (*argv);
      else {level = atoi (*++argv); argc--;}
      level = clamp (level, -1, 3);
      Femlib::value ("NUMOPTLEVEL", level);
      break;
    }
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
  
  if   (argc != 1) message (routine, "no session definition file", ERROR);
  else             session = *argv;
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			BoundarySys*&     bsys   ,
			Domain*&          domain )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed in order of creation.
// ---------------------------------------------------------------------------
{
  const integer      verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  integer            i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_POLY");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.

  VERBOSE cout << "Building boundary condition manager ..." << endl;

  bman = new BCmgr (file, elmt);

  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.

  VERBOSE cout << "Building domain ..." << endl;

  domain = new Domain (file, elmt, bman);

  VERBOSE cout << "done" << endl;
}
