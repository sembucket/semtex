#include <stab.h>
#include <krylov.h>

static const char prog[] = "floquet";
static void getargs    (int, char**, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&, Krylov*&);
void NavierStokes (Domain*, STABAnalyser*);
void Initialize(Krylov*, Domain* D, STABAnalyser*);


// --------------------------------------------------------------------------
// ------------------- MAIN -------------------------------------------------
// --------------------------------------------------------------------------

int main (int    argc,
	  char** argv)
{
  char             *session;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;
  STABAnalyser*    adjunct;
  Krylov*          krylov;
  integer          converged = 0;
  real             norm;
  int              kdim;


  Femlib::initialize(&argc, &argv);
  getargs(argc, argv, session);

  preprocess (session, file, mesh, elmt, bman, bsys, domain, krylov);
 
  adjunct = new STABAnalyser (domain, file);

  // domain restart file op is part of initialise
  // only required if no krylov file.

  domain -> loadbase();

  domain -> report();

  kdim = krylov -> dim();

  if (krylov -> restart()) Initialize(krylov, domain, adjunct);

  // krylov matrix is now fully filled
  // either by a restart file or Initialize procedure.

  while (!converged) {

    norm = krylov -> norm(kdim);
    krylov -> scale(norm);
    krylov -> roll();
    krylov -> setDomain (domain, kdim-1 );

    // Setup and call navier stokes operation.
    domain -> step = 0;
    NavierStokes(domain, adjunct);

    krylov -> getDomain (domain, kdim ); 



    converged = 1;  //evtest.
    // -- Get subspace eigenvalues, test for convergence.

    // Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    // EV_small (Tseq, ntot, kdim, zvec, wr, wi, resnorm, verbose); 

    //converged = EV_test (itrn, kdim, zvec, wr, wi, resnorm, evtol, nvec,
    //			 Total_step, restart, domain->time);
    }

  krylov->dump();
}


void Initialize(Krylov* K, Domain* D, STABAnalyser* adjunct)
  // will initialise matrix if no restart file.
  // Krylov matrix is already randomized upon creation.
  // if domain restart field exists, use it!

{
  int  i;

  cout << "-- Constructing New Krylov matrix." << endl;

  D -> restart(); 
  K -> getDomain(D, 0);

  K -> normalise(0);

  // -- Fill initial Krylov sequence.
  for (i = 1; i <= K ->dim(); i++) {

    cout << "-- constructing vector " << i << " of " << K->dim() << endl;

    // copy column[i-1] across to domain structure 
    K -> setDomain(D, i-1);

    // call to Linear NS routine.
    D -> step = 0;
    NavierStokes (D, adjunct);
    //Total_step += D->step;

    K -> getDomain(D,i);
 
    //Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    //EV_small (Tseq, ntot, i, zvec, wr, wi, resnorm, verbose);
    //EV_test  (i, i, zvec, wr, wi, resnorm, evtol, i, Total_step,
    //      restart, domain->time);
  }

  cout << "-- matrix construction ...done" << endl;
}


// --------------------------------------------------------------------------
// --------------- GETARGS --------------------------------------------------
// --------------------------------------------------------------------------

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
    "  -v[v...]  ... increase verbosity level\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      do
	Femlib::value ("VERBOSE",   (integer) Femlib::value ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
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


// --------------------------------------------------------------------------
// --------------- PREPROCESS -----------------------------------------------
// --------------------------------------------------------------------------

static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			BoundarySys*&     bsys   ,
			Domain*&          domain ,
			Krylov*&          krylov )

  /* Create objects needed for execution, given the session file name.
     They are listed in order of creation. */

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

  // -- Build the krylov matrix

  VERBOSE cout << "Building Krylov Matrix ..." << endl;

  krylov = new Krylov(file                 , 
		      (int)Femlib::value("KDIM"),
		      (int)Femlib::value("NVEC"),
		      domain->nField()*Geometry::planeSize());

  VERBOSE cout << "done" << endl;
}

