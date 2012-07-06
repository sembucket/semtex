///////////////////////////////////////////////////////////////////////////////
// forc.C: make 2D roll forcing from streak and wave solutions using
// asymptotics from Hall & Sherwin 2010.
//
// Copyright (c) 2010 <--> $Date$, Hugh Blackburn
//
// Synopsis:
// --------
// forc [-h] [-v] [-p] session streak.fld wave.fld
//
// Description: 
// -----------
//
// From the streak u velocity component field (2D/1C/real) and the
// wave pressure field (2D/1C/complex), compute the (real 2D/2C/real)
// forcing for the roll (v, w) velocity field. Output semtex field.
// 
// Streak file is the same as the base flow for the wave system: u is
// the third velocity component, and N_Z=1.  Wave file is an eigenmode
// from the wave stability problem, where there are three velocity
// components, then pressure, N_Z=2.
//
// NB: The session file needs to declare the constants 
//    BETA (which is the equivalent of alpha, streamwise wavenumber),
//    RHO which is the amplitude of the wave system,
//    CHI which is the scale of the spatial delta function parameter.
//
// Geometric singularities and vector divide operations can do nasty
// things here.  Search the output and set any nan values to zero.
//
// Flag -p forces output to be pressure mode magnitude rather than forcing.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>
#include <data2df.h>

static char  prog[]  = "forc";
static int_t verbose = 0;
static void  getargs  (int, char**, char*&, bool&, istream*&, istream*&);
static bool  getDump  (istream&, Header&, AuxField&, 
		       const char, const int_t, const int_t, const int_t);
static bool  needSwap (const char*);
static int_t _index   (const char*, char);
static void  d_z      (AuxField*, AuxField*, AuxField*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char               *session;
  istream            *streakfile, *wavefile;
  Header             hdr;
  int_t              np, nel, nz, ntot1D, ntot2D, i;
  real_t             Area = 0.0, integral;
  FEML*              F;
  Mesh*              M;
  BCmgr*             B;
  BoundarySys*       bsys;
  vector<Element*>   Esys;
  Field              *u;        // -- Make 1 field so we can do smoothing.
  AuxField           *P;	// -- This is the only complex AuxField.
  AuxField           *dP2, *ddP2, *dPr, *dPi;
  AuxField           *u2, *u_z, *u_y, *u_yy, *u_zz, *u_yz;
  AuxField           *f_z, *f_zz, *tmp, *Delta, *Delta_z, *lambda, *a, *a_z;
  AuxField           *J, *K, *f1, *f2, *F1, *F2, *delta;
  real_t*            work;
  bool               presmag = false;

  // -- Initialize.

  Femlib::initialize (&argc, &argv);
  getargs            (argc, argv, session, presmag, streakfile, wavefile);

  // -- Set up 2D mesh information.
  
  F = new FEML  (session);
  M = new Mesh  (F);


  const real_t n_0   = Femlib::value ("TWOPI*(2./3.)^(2./3.)*gamma(1.-2./3.)");
  const real_t alpha = Femlib::value ("BETA");
  const real_t rho   = Femlib::value ("RHO");
  const real_t chi   = Femlib::value ("CHI");
  
  nel = M -> nEl();  
  np  = Femlib::ivalue ("N_P");

  Esys.resize (nel);

  Geometry::set (np, nz = 2, nel, Geometry::Cartesian);
  ntot2D = Geometry::nTotProc();

  Geometry::set (np, nz = 1, nel, Geometry::Cartesian);
  ntot1D = Geometry::nTotProc();

  for (i = 0; i < nel; i++) {
    Esys[i] = new Element (i, np, M);
    Area   += Esys[i] -> area();
  }
  VERBOSE cerr << "Domain area: " << Area << endl;

  work = new real_t [ntot1D];

  // -- Create the initial storage to be pulled out of supplied files.
  // -- u is set up as a Field (rather than AuxField) so we can smooth.

  B    = new BCmgr       (F, Esys);
  bsys = new BoundarySys (B, Esys, 'u');

  P = new AuxField (      new real_t [ntot2D], nz = 2, Esys, 'P');
  u = new Field    (bsys, new real_t [ntot1D], nz = 1, Esys, 'u');

  // -- Read in these data.

  getDump (*streakfile, hdr, *u, 'w', np, nel, 1);
  getDump (*wavefile,   hdr, *P, 'p', np, nel, 2);

  // -- Normalize P.

  integral = 2.0 * P -> mode_L2 (0);
  VERBOSE cerr << "Pressure norm: " << integral << endl;
  *P /= sqrt(integral);
  integral = 2.0 * P -> mode_L2 (0);
  VERBOSE cerr << "Normalized pressure norm: " << integral << endl;

  // -- Derived values.

  dPr     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  dPi     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  dP2     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  ddP2    = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  u_z     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  u_y     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  u_yy    = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  u_zz    = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  u_yz    = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  f_z     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  f_zz    = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  Delta   = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  Delta_z = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  a       = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  a_z     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  tmp     = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  J       = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  K       = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  f1      = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  f2      = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  delta   = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  u2      = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  F1      = new AuxField (new real_t [ntot1D], nz = 1, Esys);
  F2      = new AuxField (new real_t [ntot1D], nz = 1, Esys);

  // -- Streak velocity derivatives.

  (*u_z  = *u) . gradient (0);
  (*u_y  = *u) . gradient (1);

  (*u_zz = *u_z) . gradient (0);
  (*u_yy = *u_y) . gradient (1);
  (*u_yz = *u_y) . gradient (0); // -- OR (*u_yz = *u_z) . gradient (1);

  lambda = u_y;			 // -- An alias.

  // -- f_z and f_zz.

  *f_z  = *u_z;
  *f_z *= -1.0;
  *f_z /= *u_y;
  
  *f_zz  = *u_yz;
  *f_zz += *u_yz;
  *f_zz *= *f_z;
  *f_zz += *u_zz;
  *tmp   = *u_yy;
  *tmp  *= *f_z;
  *tmp  *= *f_z;
  *f_zz += *tmp;
  *f_zz /= *u_y;
  *f_zz *= -1.0;

  // -- Find |dP/dz|^2.

  P   -> getPlane (0, work);
  dPr -> setPlane (0, work);

  d_z (dPr, f_z, tmp);

  P   -> getPlane (1, work);
  dPi -> setPlane (0, work);

  d_z (dPi, f_z, tmp);

  *dPr *= *dPr;
  *dPi *= *dPi;

  *dP2   = *dPr;
  *dP2 +=  *dPi;

  *ddP2 = *dP2;

  d_z (ddP2, f_z, tmp);
  
  // -- More derived values.

  *Delta  = *f_z;
  *Delta *= *f_z;
  *Delta +=  1.0;

  *Delta_z = *Delta;

  d_z (Delta_z, f_z, tmp);

  *a  = *lambda;
  *a /= *Delta;
  *a *= alpha;

  *a_z = *a;
  
  d_z (a_z, f_z, tmp);

  *K  = *dP2;
  *K *= *f_zz;
  *K *= -1.;

  ((*J   = *a_z    ) /= *a    ) *= -5./3.;
  ((*tmp = *Delta_z) /= *Delta) *= -7./2.;
  *J += *tmp;
  *J *= *dP2;
  *J += *ddP2;

  (*tmp = *Delta) . pow (-5.);
  (*a)            . pow (-5./3.);
  *tmp *= *a;
  *tmp *= n_0 * rho * rho;

  *J *= *tmp;
  *K *= *tmp;

  // -- These are smooth forcing fields (without spatial delta function).
  
  *f1   = *K;
  *tmp  = *J;
  *tmp *= *f_z;
  *tmp *= *Delta;
  *f1  -= *tmp;
  *f1  *= *lambda;

  *f2   = *K;
  *f2  *= *f_z;
  *tmp  = *J;
  *tmp *= *Delta;
  *f2  += *tmp;
  *f2  *= *lambda;
  *f2  *= -1.;

  // -- Final clean up (a hack): set any NaN values to zero.

  f1 -> zeroNaN();
  f2 -> zeroNaN();

  // -- Make spatial delta function.

  (*u2 = *u) *= *u;
  *delta = *u2;
  *delta /= -chi;
  delta -> exp ();
  *delta /= Femlib::value ("sqrt(PI*CHI)");

  // -- Ready to output, but first smooth outcomes. Note "reversed" ordering.
  
  F1 -> times (*f1, *delta);
  F2 -> times (*f2, *delta);

  u -> smooth (F1);
  u -> smooth (F2);

  sprintf (hdr.flds, "uv");
  hdr.nz = 1;
  cout << hdr;

  if (presmag) {
    P   -> getPlane (0, work); dPr -> setPlane (0, work);
    P   -> getPlane (1, work); dPi -> setPlane (0, work);
    *dPr *= *dPr;
    *dPi *= *dPi;
    *dPr += *dPi;

    cout << *dPr;
    cout << *dP2;
  } else {
    cout << *F2;
    cout << *F1;
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc   ,
		     char**    argv   ,
		     char*&    session,
		     bool&     presmag,
		     istream*& streak ,
		     istream*& wave   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: forc [options] session streak.fld wave.fld\n"
    "options:\n"
    "  -h ... print this message\n"
    "  -v ... verbose\n"
    "  -p ... output (normalised) pressure magnitudes |P|^2 and |dP/dz|^2\n";
  char *name;

 
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
      presmag = true;
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 3) { 
    session = argv[0]; 

    streak = new ifstream (argv[1]);
    if (streak -> bad()) message (prog, "bad streak file?", ERROR);

    wave   = new ifstream (argv[2]);
    if (wave   -> bad()) message (prog, "bad wave file?",   ERROR);
    
  } else
    message (prog, usage, ERROR);
}


static bool getDump (istream&    file ,
		     Header&     hdr  ,
		     AuxField&   field, 
		     const char  name ,
		     const int_t NP   ,
		     const int_t NEL  ,
		     const int_t NZ   )
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.  Check
// that the data has conforming storage. Input variable name is the
// appropriate field in storage. File is already open.
// ---------------------------------------------------------------------------
{
  char   buf[StrMax], fields[StrMax];
  int_t  i, idx, nf, npnew, nznew, nelnew;

  file >> hdr;
  
  if (NP != hdr.nr || NZ != hdr.nz || NEL != hdr.nel)
    message (prog, "size of dump mismatch with session file", ERROR);

  idx = _index (hdr.flds, name);
  if (idx < 0) {
    sprintf (buf, "could not find field %c in list %s", name, hdr.flds);
    message (prog, buf, ERROR);
  }

  VERBOSE cerr << "Extracting field: " << hdr.flds[idx] << endl;

  // -- Read binary field data.

  for (i = 0; i <= idx; i++) file >> field;

  if (needSwap (hdr.frmt)) field . reverse();

  return file.good();
}


static bool needSwap (const char* ffmt)
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


static int_t _index (const char* s, char c)
/* ------------------------------------------------------------------------- *
 * Return index of c in s, -1 if not found. Starts at zero.
 * ------------------------------------------------------------------------- */
{
  int_t       i;
  const int_t len = strlen (s);

  for (i = 0; i < len; i++) if (s[i] == c) return i;

  return -1;
}


static void d_z (AuxField* inout, AuxField* f_z, AuxField* work)
// ---------------------------------------------------------------------------
// Take z derivative along interface: d()/dz = d()/dy * f_z + d()/dz.
// ---------------------------------------------------------------------------
{
  ((*work = *inout) . gradient (1)) *= *f_z;

  ((*inout) . gradient (0)) += *work;
}

