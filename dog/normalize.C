///////////////////////////////////////////////////////////////////////////////
// normalize.C: normalize/rescale an eigenmode, either on the basis of
// kinetic energy or pressure.
//
// Copyright (c) 2011 <--> $Date$, Hugh Blackburn
//
// Synopsis:
// --------
// normalize [-h] [-p] [-s scale] session [session.fld]
//
// Description: 
// ----------- 
// 
// Given a velocity + pressure field, normalize it so that it has either
//
// \sqrt[(\int u . u dA)/A] = 1 (kinetic energy)
//
// or
// 
// \sqrt[(\int p^2 dA)/A] = 1 (pressure)
// 
// where A = \int dA is the 2D area of the domain. The mode has to be
// either real (nz = 1) or complex (nz = 2) and this is taken from the
// field file, not the session file. The session file has to
// correspond with the field file in respect of N_P and number of
// elements. If the session file indicates CYLINDRICAL then the
// integrands are multiplied by radius.  The integration is
// approximated by GLL quadrature at the same order as the field file
// (and no dealiasing of product terms is carried out).
//
// Following the normalization, optionally scale by a given
// multiplicative factor (supplied on command line).
//
// Take only the first dump in the field file. 
// 
// Note that this code is an update of "normalise.C", which uses the
// discrete 2-norm of the velocity vector, rather then the L_2 norm.
// 
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>
#include <data2df.h>

static char  prog[]  = "normalize";
static int_t verbose = 0;
static void  getargs  (int, char**, char*&, char*&, bool&, real_t&);
static int_t getDump  (istream&, vector<AuxField*>&, Header&, 
		       vector<Element*>&, const int_t, const int_t);
static bool  doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char               *session = 0, *dump = 0;
  istream            *fldfile;
  Header             hdr;
  int_t              i, NP, NEL, NCMP;
  real_t             norm, scale = 1.0;
  FEML*              F;
  Mesh*              M;
  bool               pressure = false;
  vector<Element*>   Esys;
  vector<AuxField*>  u;

  // -- Initialize.

  Femlib::initialize (&argc, &argv);
  getargs            (argc, argv, session, dump, pressure, scale);

  if (dump) {
    fldfile = new ifstream (dump);
    if (fldfile -> bad()) message (prog, "no field file", ERROR);
  } else fldfile = &cin;

  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh (F);
  NEL = M -> nEl();  
  NP  = Femlib::ivalue ("N_P");

  Esys.resize (NEL);
  for (i = 0; i < NEL; i++) Esys[i] = new Element (i, NP, M);
  
  // -- Load field file, reset nz, get number of velocity components;

  NCMP = getDump (*fldfile, u, hdr, Esys, NP, NEL);

 // -- Calculate norm.

  if (pressure)
    norm = sqrt (2.0 * u[NCMP] -> mode_L2 (0));
  else {
    for (norm = 0.0, i = 0; i < NCMP; i++) norm += u[i] -> mode_L2 (0);
    norm = sqrt (2 * norm);
  }

 // -- Normalisation + scaling.

  for (i = 0; i <= NCMP; i++) *u[i] *= scale/norm;

  // -- Output.

  cout << hdr;
  for (i = 0; i <= NCMP; i++) cout << *u[i];

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int     argc    ,
		     char**  argv    ,
		     char*&  session ,
		     char*&  dump    ,
		     bool&   pressure,
		     real_t& scale   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: normalize [options] session [session.fld]\n"
    "options:\n"
    "-h       ... print this message\n"
    "-p       ... mode normalization on pressure rather than velocity\n"
    "-s scale ... multiplicative scale applied after normalization\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'p':
      pressure = true;
      break;
    case 's':
      if   (*++argv[0]) scale = atof (*argv);
      else { --argc;    scale = atof (*++argv); }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if      (argc == 1)   session = argv[0];
  else if (argc == 2) { session = argv[0]; dump = argv[1]; }
  else                  message (prog, usage, ERROR);
}


static int_t getDump (istream&           file,
		      vector<AuxField*>& u   ,
		      Header&            hdr ,
		      vector<Element*>&  Esys,
		      const int_t        np  ,
		      const int_t        nel )
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// Check we have either uvp or uvwp, return number of velocity components.
//
// Warning: SIDE EFFECT: nz (which must be 1 or 2 in field file) is reset.
// ---------------------------------------------------------------------------
{
  char    buf[StrMax], fields[StrMax];
  int_t   i, swab, nf;
  real_t* alloc;

  file >> hdr;

  if (hdr.nz < 1 || hdr.nz > 2)
    message (prog, "nz in field file must be 1 or 2", ERROR);

  if (np != hdr.nr || nel != hdr.nel)
    message (prog, "size of dump mismatch with session file", ERROR);

  nf = strlen (hdr.flds);
  Femlib::ivalue ("N_Z", hdr.nz);
  Geometry::set  (hdr.nel, nf - 1);

  // -- Check we have (only) uvp or uvwp.

  if (hdr.flds[0] != 'u')
    message (prog, "fields must start with u", ERROR);

  if (!(strstr ("uvp",  hdr.flds) || strstr ("uvwp", hdr.flds)))
    message (prog, "fields must start with either uvp or uvwp", ERROR);

  if (strstr ("uvp", hdr.flds) && nf != 3)
    message (prog, "2-component field file must have only u v p", ERROR);

  if (strstr ("uvwp", hdr.flds) && nf != 4)
    message (prog, "3-component field file must have only u v w p", ERROR);

  // -- Arrange for byte-swapping if required.

  swab = doSwap (hdr.frmt);

  // -- Create AuxFields.

  u.resize (nf);
  for (i = 0; i < nf; i++) {
    alloc = new real_t [Geometry::nTotProc()];
    u[i]  = new AuxField (alloc, hdr.nz, Esys, fields[i]);
  }

  // -- Read binary field data.

  for (i = 0; i < nf; i++) {
    file >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return nf - 1;		// -- Number of velocity components.
}


static bool doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping is required to make sense of binary input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}
