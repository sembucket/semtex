///////////////////////////////////////////////////////////////////////////////
// integral.C: return the domain integral of all fields in dump file.
//
// Copyright (c) 1999 Hugh Blackburn
//
// Synopsis:
// --------
// integral [-h] [-v] session [file]
//
// Description:
// -----------
// Read in file, print up area of domain.  If 3D perform Fourier
// transform to get mean value into plane zero for each field.  Then
// return integral for each scalar field.  For 3D, values are
// multiplied by domain length, to produce volume integrals.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>

static char    prog[]  = "integral";
static integer verbose = 0;
static void    getargs  (int, char**, char*&, char*&);
static integer getDump  (ifstream&, vector<AuxField*>&, vector<Element*>&,
			 const integer, const integer, const integer);
static integer doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char              *session, *dump;
  ifstream          fldfile;
  integer           NP, NZ,  NEL;
  integer           np, nel, ntot, i;
  real              Lz, Area = 0.;
  const real        *z;
  FEML*             F;
  Mesh*             M;
  vector<Element*>  Esys;
  vector<AuxField*> u;

  // -- Initialize.

  Femlib::initialize (&argc, &argv);
  getargs            (argc, argv, session, dump);
  cout.precision     (8);

  fldfile.open (dump, ios::in);
  if (!fldfile) message (prog, "no field file", ERROR);

  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh (F);

  NEL = M -> nEl();  
  NP  = (integer)  Femlib::value ("N_POLY");
  NZ  = (integer)  Femlib::value ("N_Z"   );
  Lz  = (NZ > 1) ? Femlib::value ("TWOPI / BETA") : 1.;
  
  Geometry::set (NP, NZ, NEL, Geometry::Cartesian);
  Femlib::mesh  (GLL, GLL, NP, NP, &z, 0, 0, 0, 0);
  Esys.resize   (NEL);

  for (i = 0; i < NEL; i++) {
    Esys[i] = new Element (i, M, z, NP);
    Area   += Esys[i] -> area();
  }
  cout << Area << endl;
  
  // -- Load field file, interpolate within it.

  while (getDump (fldfile, u, Esys, NP, NZ, NEL)) {
    for (i = 0; i < u.size(); i++) {
      u[i] -> transform (+1);
      cout << u[i] -> name() << ": " << Lz * u[i] -> integral (0) << endl;
    }
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session,
		     char*& dump   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: integral [options] session [dump]\n"
    "options:\n"
    "-h ... print this message\n"
    "-v ... verbose output\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc != 2) message (prog, usage, ERROR);
  else { session = argv[0]; dump = argv[1]; }

  if (argc == 1) {
  } else {
    message (prog, usage, ERROR);
  }
}


static integer getDump (ifstream&          file,
			vector<AuxField*>& u   ,
			vector<Element*>&  Esys,
			const integer      np  ,
			const integer      nz  ,
			const integer      nel )
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  char    buf[StrMax], fields[StrMax];
  integer i, swab, nf, npnew, nznew, nelnew;
  real*   alloc;

  if (file.getline(buf, StrMax).eof()) return 0;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> npnew >> nznew >> nznew >> nelnew;
  file.getline (buf, StrMax);
  
  if (np != npnew || nz != nznew || nel != nelnew)
    message (prog, "size of dump mismatch with session file", ERROR);

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  nf = strlen  (fields);
  file.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);

  // -- Create AuxFields on first pass.

  if (u.size() == 0) {
    u.resize (nf);
    for (i = 0; i < nf; i++) {
      alloc = new real [Geometry::nTotProc()];
      u[i]  = new AuxField (alloc, nz, Esys, fields[i]);
    }
  } else if (u.size() != nf) 
    message (prog, "number of fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (i = 0; i < nf; i++) {
    file >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return file.good();
}


static integer doSwap (const char* ffmt)
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
