//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element DNS for incompressible flows.
//
// Copyright (C) 1994, 1999  Hugh Blackburn.
//
// USAGE:
// -----
// dns [options] session
//   options:
//   -h       ... print usage prompt
//   -i[i]    ... use iterative solver for viscous [and pressure] steps
//   -v[v...] ... increase verbosity level
//   -chk     ... checkpoint field dumps
//
// AUTHOR:
// ------
// Hugh Blackburn
// CSIRO Division of Building, Construction and Engineering
// P.O. Box 56
// Highett, Vic 3190
// Australia
// hugh.blackburn@dbce.csiro.au
//
// RCSid[] = "$Id$";
//
//////////////////////////////////////////////////////////////////////////////

#include <stab.h>

static char prog[] = "stab";
//static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs    (int, char**, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&);

void NavierStokes (Domain*, STABAnalyser*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  //  set_new_handler (&memExhaust);
#if !defined(__DECCXX)
  ios::sync_with_stdio();
#endif

  char             *session;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;
  STABAnalyser*    adjunct;

  char             host[256];
  integer          pid;
  
  // Get process id and hostname
  pid = getpid();
  gethostname(host, 256);

  cout << "host = " << host << endl;
  cout << "pid  = " << pid << endl;

  Femlib::initialize (&argc, &argv);

  // -- Set default values for eigensystem routines.

  Femlib::value ("DUMP"      , 1 );  // -- toggle field dump for each NS call
  Femlib::value ("IO_STEP"   , 1 );  // -- 0 = no step output
  Femlib::value ("IO_EIG"    , 2 );  // -- # of eigenvalues to output
  Femlib::value ("MAX_BASE_F", 50);  // -- max # of base files
  Femlib::value ("KDIM"      , 20);  // -- dimension of krylov matrix
  Femlib::value ("NVEC"      , 5 );  // -- number of eigenvectors to test
  Femlib::value ("IO_KRY"    , 10);  // -- how often to write a krylov file

  getargs (argc, argv, session);

  preprocess (session, file, mesh, elmt, bman, bsys, domain);
 
  adjunct = new STABAnalyser (domain, file);

  domain -> restart();

  domain -> loadbase();

  domain -> report();
  
  NavierStokes (domain, adjunct);

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
  const char routine[] = "getargs";
  char       buf[StrMax];
  char       usage[]   = "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h        ... print this message\n"
    "  -i[i]     ... use iterative solver for viscous [& pressure] steps\n"
    "  -v[v...]  ... increase verbosity level\n"
    "  -chk      ... checkpoint field dumps\n"
    "  -d        ... diagnostic test\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      do
	Femlib::value ("ITERATIVE", (integer) Femlib::value ("ITERATIVE") + 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE",   (integer) Femlib::value ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	Femlib::value ("CHKPOINT",  (integer) 1);
      } else {
	fprintf (stdout, usage, prog);
	exit (EXIT_FAILURE);	  
      }
      break;
    case 'd':
      cout << "Enabling Diagnostic Test" << endl;
      break;
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
  const integer      verbose = (integer) Femlib::value ("VERBOSE");
  Geometry::CoordSys space;
  const real*        z;
  integer            i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   = mesh -> nEl();
  np    =  (integer) Femlib::value ("N_POLY");
  nz    =  (integer) Femlib::value ("N_Z");
  space = ((integer) Femlib::value ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, mesh, z, np);

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
