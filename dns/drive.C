//////////////////////////////////////////////////////////////////////////////
// drive.C
//
// SYNOPSIS:
// --------
// Control spectral element DNS for incompressible flows.
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
// hmb@dbce.csiro.au
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>
#include <new.h>

static char prog[] = "dns";
static void memExhaust   () { message ("new", "free store exhausted", ERROR); }
static void getargs      (integer, char**, char*&);
       void NavierStokes (Domain*, Analyser*);


integer main (integer argc,
	      char**  argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#if !defined(__DECCXX)
  ios::sync_with_stdio();
#endif

  Geometry::CoordSys system;
  char      *session, fields[StrMax];
  integer   np, nz, nel;
  FEML*     F;
  Mesh*     M;
  BCmgr*    B;
  Domain*   D;
  Analyser* A;
  
  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session);
  
  F = new FEML  (session);
  M = new Mesh  (*F);
  B = new BCmgr (*F);

  nel    = M -> nEl();  
  np     =  (integer) Femlib::value ("N_POLY");
  nz     =  (integer) Femlib::value ("N_Z"   );
  system = ((integer) Femlib::value ("CYLINDRICAL") ) ?
                                Geometry::Cylindrical : Geometry::Cartesian;  

  Geometry::set (np, nz, nel, system);
  if   (nz > 1) strcpy (fields, "uvwp");
  else          strcpy (fields, "uvp");

  D = new Domain   (*F, *M, *B, fields, session);
  A = new Analyser (*D, *F);

  D -> initialize();
  ROOTONLY D -> report();
  
  NavierStokes (D, A);

  Femlib::finalize();

  return EXIT_SUCCESS;
}


static void getargs (integer argc   ,
		     char**  argv   ,
		     char*&  session)
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
    "  -chk      ... checkpoint field dumps\n";
 
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
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
  
  if   (argc != 1) message (routine, "no session definition file", ERROR);
  else             session = *argv;
}
