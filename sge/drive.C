//////////////////////////////////////////////////////////////////////////////
// drive.C
//
// SYNOPSIS:
// --------
// Control spectral element unsteady incompressible flow solver for
// simulation of advective/diffusive transport in coated capillary tube.
// Cartesian coordinates only.
//
// USAGE:
// -----
// chroma [options] session
//   options:
//   -h       ... print usage prompt
//   -i       ... use iterative solver for diffusion step.
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

#include <chroma.h>
#include <new.h>

static char  prog[] = "chroma";
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }
static void  getargs    (int, char**, char*&);
static real* gasdata    (const char*, const int, const int);
       void  transport  (Domain*, RunInfo*, const real*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#if !defined(__DECCXX)
  ios::sync_with_stdio();
#endif

  char*    session;
  int      np, nz, nel;
  FEML*    F;
  Mesh*    M;
  BCmgr*   B;
  Domain*  D;
  RunInfo* I;
  real*    G;
  
  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session);

  F = new FEML  (session);
  M = new Mesh  (*F);
  B = new BCmgr (*F);

  nel = M -> nEl();  
  np  = (int) Femlib::value ("N_POLY");
  nz  = (int) Femlib::value ("N_Z");
  if (nz < 2) message (prog, "require at least 2 data planes", ERROR);
  
  Geometry::set (np, nz, nel, Geometry::Cartesian);

  I = new RunInfo (*D, *F);
  D = new Domain  (*F, *M, *B, "c", session);
  G = gasdata     (session, np, nel);

  D -> initialize(); ROOTONLY D -> report();

  transport (D, I, G);

  Femlib::finalize ();

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
  char buf[StrMax];
  char usage[]   =
    "Usage: %s [options] session-file\n"
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
	Femlib::value ("ITERATIVE", (int) Femlib::value ("ITERATIVE") + 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE",   (int) Femlib::value ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	Femlib::value ("CHKPOINT", 1);
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
  
  if   (argc != 1) message (prog, "no session definition file", ERROR);
  else             session = *argv;
}


static real* gasdata (const char* session,
		      const int   np     ,
		      const int   nel    )
// ---------------------------------------------------------------------------
// Load the (w-only) 2D gas phase velocity data from file
// "session.gas" This is a standard 2D SEM field file.  Check that it
// has the correct element order and number of elements matches
// declaration NEL_GAS in session file.  This must not be greater than
// Geometry::nEl.  Remaining values are set to zero (i.e. the velocity
// in the stationary phase).  NB the implication is that the first
// NEL_GAS elements in all session files correspond to the gas phase.
// ---------------------------------------------------------------------------
{
  char routine[] = "gasdata", buf[StrMax], err[StrMax];
  int       i, NP, NZ, NEL, swab;
  const int nProc = Geometry::nProc();
  const int nelG  = (int) Femlib::value ("NEL_GAS");
  real*     w     = new real [np * np * nel];

  ROOTONLY {
    strcat (strcpy (buf, session), ".gas");
    ifstream file  (buf);

    if (!file) {
      sprintf (err, "file with gas phase velocity data, %s, not found", buf);
      message (routine, err, ERROR);
    }

    file.getline(buf, StrMax);
    if (!strstr (buf, "Session")) {
      sprintf (err, "%s.gas doesn't look like a field file", session);
      message (routine, err, ERROR);
    }

    file.getline(buf, StrMax).getline(buf, StrMax);
    istrstream (buf, strlen(buf)) >> NP >> NP >> NZ >> NEL;
    
    if (NP != np) {
      sprintf (err, "size mismatch, np: %1d/%1d", np, NP);
      message (routine, err, ERROR);
    }

    if (nelG != NEL) {
      sprintf (err, "size mismatch: nel (%1d) != NEL_GAS (%1d)", NEL, nelG);
      message (routine, err, ERROR);
    }

    file.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);
    file.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);
    file.getline(buf, StrMax);

    Veclib::describeFormat (err);

    if (!strstr (buf, "binary"))
      message (routine, "input field file not in binary format", ERROR);
  
    if (!strstr (buf, "endian"))
      message (routine, "input field file in unknown binary format", WARNING);
    else {
      swab = ((strstr (buf, "big") && strstr (err, "little")) ||
	      (strstr (err, "big") && strstr (buf, "little")) );
    }

    file.read ((char*) w, np * np * nelG * sizeof (real));
    if (swab) Veclib::brev (np * np * nelG, w, 1, w, 1);
    Veclib::zero (np * np * (nel - NEL), w + np * np * NEL, 1);
  
    if (nProc > 1)
      for (i = 1; i < nProc; i++)
	Femlib::send (w, np * np * nel, i);
  } else			// -- This must be a parallel run.
    Femlib::recv (w, np * np * nel, 0);
}


