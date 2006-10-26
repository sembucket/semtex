//////////////////////////////////////////////////////////////////////////////
// traction.C: Compute tractions on wall boundaries from field file.
//
// Copyright (c) 2006 <--> $Date$, Hugh Blackburn
//
// USAGE:
// -----
// traction session [file]
//
// Essentially this carries out the same computation as is done during
// execution of dns, but as a standalone utility.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>

static char prog[] = "traction";
static void getargs    (int, char**, char*&, istream*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&);
static bool getDump    (Domain*, istream&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  Geometry::CoordSys         system;
  char*                      session;
  istream*                   file;
  FEML*                      F;
  Mesh*                      M;
  BCmgr*                     B;
  Domain*                    D;
  vector<Element*>           E;

  Femlib::initialize (&argc, &argv);

  // -- Read command line.

  getargs (argc, argv, session, file);

  // -- Set up domain.

  preprocess (session, F, M, E, B, D);

  // -- Loop over all dumps in field file, compute and print traction.

  while (getDump (D, *file)) {

  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int        argc   ,
		     char**     argv   ,
		     char*&     session,
		     istream*&  file   )
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session [file]"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -v[vv..] ... increase verbosity level\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE", Femlib::ivalue ("VERBOSE") + 1);
      while (*++argv[0] == 'v');
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  switch (argc) {
  case 1:
    session = argv[0];
    file = &cin;
    break;
  case 2:
    session = argv[0];
    file = new ifstream (argv[1]);
    if (file -> bad()) {
      cerr << usage;
      sprintf (buf, "unable to open field file: %s", argv[1]);
      message (prog, buf, ERROR);
    }
    break;
  default:
    cerr << usage;
    message (prog, "session file not supplied", ERROR);
    break;
  }  
}


static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			Domain*&          domain )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed in order of creation.
// ---------------------------------------------------------------------------
{
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  int_t              i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
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


static bool getDump (Domain*    D   ,
		     istream&   dump)
// ---------------------------------------------------------------------------
// Read next set of field dumps from file.
// ---------------------------------------------------------------------------
{
  dump >> *D;
  return dump.good ();
}
