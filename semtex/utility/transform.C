///////////////////////////////////////////////////////////////////////////////
// transform.C: carry out Fourier and/or 2D polynomial transform of data.
//
// Copyright (c) 1999-2004 Hugh Blackburn
//
// USAGE
// -----
// transform [options] [file]
// options:
// -h       ... print this message.
// -i       ... invert transform.
// -l       ... polynomial transform is Legendre       [Default: modal]
// -P||F||B ... Carry out DPT (P), DFT (F) or both (B) [Default: both]
// 
// If file is not present, read from standard input.  Write to
// standard output.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem_h>

class Field2DF
// ============================================================================
// Canonical field class, each np X np element is defined on [-1,1] X [-1, 1].
// Data are arranged element-ordered in 2D planes to create a 3D scalar field.
// ============================================================================
{
friend istream& operator >> (istream&, Field2DF&);
friend ostream& operator << (ostream&, Field2DF&);

public:
  Field2DF  (const int nP, const int nZ, const int nEl,
	     const char Name='\0');
  ~Field2DF () { delete data; delete plane; }

  char getName () { return name; }

  Field2DF& operator = (const Field2DF&);
  Field2DF& operator = (const real);

  Field2DF& DFT1D (const int);
  Field2DF& DPT2D (const int, const char);

  Field2DF& reverse   ();
  
private:
  const char name;
  const int  np, nz, nel, np2;
  int        nplane, ntot;
  real*      data;
  real**     plane;
};


Field2DF::Field2DF (const int  nP  ,
		    const int  nZ  ,
		    const int  nEl ,
		    const char Name) :

		    name       (Name   ),
                    np         (nP     ),
		    nz         (nZ     ),
		    nel        (nEl    ),
		    np2        (np * np)
// ---------------------------------------------------------------------------
// Field2DF constructor. 
// ---------------------------------------------------------------------------
{
  register int i;
  
  nplane = np * np * nel;
  if (nplane & 1) nplane++;
  ntot   = nplane * nz;

  data  = new real  [ntot];
  plane = new real* [nz];

  for (i = 0; i < nz; i++) plane [i] = data + i * nplane;
  Veclib::zero (ntot, data, 1);
}


Field2DF& Field2DF::DFT1D (const int sign)
// ---------------------------------------------------------------------------
// Carry out discrete Fourier transformation in z direction.
// ---------------------------------------------------------------------------
{
  if (nz > 2) Femlib::DFTr (data, nz, nplane, sign);

  return *this;
}


Field2DF& Field2DF::DPT2D (const int sign, const char basis)
// ---------------------------------------------------------------------------
// Carry out 2D discrete polynomial transform (element-by-element) on planes.
// ---------------------------------------------------------------------------
{
  int          i;
  vector<real> work (nplane);
  const real   *Fu, *Ft, *Bu, *Bt;

  if (basis == 'l')
    Femlib::legTran (np, &Fu, &Ft, &Bu, &Bt, 0, 0);
  else
    Femlib::modTran (np, &Fu, &Ft, &Bu, &Bt, 0, 0);

  if (sign == FORWARD)
    for (i = 0; i < nz; i++)
      Femlib::tpr2d (plane[i], plane[i], &work[0], Fu, Ft, np, np, nel);
  else
    for (i = 0; i < nz; i++)
      Femlib::tpr2d (plane[i], plane[i], &work[0], Bu, Bt, np, np, nel);

  return *this;
}


Field2DF& Field2DF::operator = (const Field2DF& rhs)
// ---------------------------------------------------------------------------
// If the two fields conform, copy rhs's data storage to lhs.
//
// Otherwise perform projection/interpolation of rhs's data area to lhs.
// Interpolation ASSUMES THAT FOURIER TRANSFORMATION HAS ALREADY OCCURRED
// in z direction if rhs is 3D.  Truncation of Fourier modes occurs if this
// Field2DF has less modes than rhs (to avoid aliasing).
// ---------------------------------------------------------------------------
{
  if (rhs.nel != nel)
    message ("Field2DF::operator =", "fields can't conform", ERROR);

  if (rhs.np == np && rhs.nz == nz)
    Veclib::copy (ntot, rhs.data, 1, data, 1);

  else {			// -- Perform projection.

    register int  i, k;
    register real *LHS, *RHS;
    const real    *IN, *IT;
    const int     nzm = min (rhs.nz, nz);
    vector<real>  work (rhs.np * np);
    real*         tmp = &work[0];

    //    Femlib::mesh (GLL, GLL, rhs.np, np, 0, &IN, &IT, 0, 0);
    Femlib::projection (&IN, &IT, rhs.np, GLL, 0.0, 0.0, np, GLL, 0.0, 0.0);

    for (k = 0; k < nzm; k++) {	// -- 2D planar projections.
      LHS = plane[k];
      RHS = rhs.plane[k];

      if (rhs.np == np)
	Veclib::copy (nplane, RHS, 1, LHS, 1);
      else
	for (i = 0; i < nel; i++, LHS += np2, RHS += rhs.np2) {
	  Blas::mxm (IN, np, RHS, rhs.np, tmp, rhs.np);
	  Blas::mxm (tmp, np, IT, rhs.np, LHS,     np);
	}
    }

    if ((i = nz - rhs.nz) > 0) // -- Zero pad for Fourier projections.
      Veclib::zero (i * nplane, data + rhs.ntot, 1);
  }

  return *this;
}


Field2DF& Field2DF::operator = (const real val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (ntot,      data, 1);
  else              Veclib::fill (ntot, val, data, 1);

  return *this;
}


Field2DF& Field2DF::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (ntot, data, 1, data, 1);

  return *this;
}


ostream& operator << (ostream&  strm,
		      Field2DF& F   )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  int i;
  
  for (i = 0; i < F.nz; i++)
    strm.write ((char*) F.plane[i], F.np * F.np * F.nel * sizeof (real));

  return strm;
}


istream& operator >> (istream&  strm,
		      Field2DF& F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.
// ---------------------------------------------------------------------------
{
  int i;
  
  for (i = 0; i < F.nz; i++)
    strm.read ((char*) F.plane[i], F.np * F.np * F.nel * sizeof (real));

  return strm;
}


static char prog[] = "transform";
static void getargs  (int, char**, int&, char&, char&, istream*&);
static int  getDump  (istream&, ostream&, vector<Field2DF*>&);
static void loadName (const vector<Field2DF*>&, char*);
static int  doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int               i, dir = FORWARD;
  char              type   = 'B', basis = 'm';
  istream*          input;
  vector<Field2DF*> u;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, dir, type, basis, input);
  
  while (getDump (*input, cout, u))
    for (i = 0; i < u.size(); i++) {
      if (type == 'P' || type == 'B') u[i] -> DPT2D (dir, basis);
      if (type == 'F' || type == 'B') u[i] -> DFT1D (dir);
      cout << *u[i];
    }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     int&      dir  ,
		     char&     type ,
		     char&     basis,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: transform [options] [file]\n"
    "options:\n"
    "-h ... print this message\n"
    "-P ... Discrete Polynomial Transform (2D)\n"
    "-F ... Discrete Fourier    Transform (1D)\n"
    "-B ... do both DPT & DFT [Default]\n"
    "-i ... carry out inverse transform instead\n"
    "-l ... use Legendre basis functions instead of modal expansions\n";
    
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      dir = INVERSE;
      break;
    case 'l':
      basis = 'l';
      break;
    case 'P':
      type = 'P';
      break;
    case 'F':
      type = 'F';
      break;
    case 'B':
      type = 'B';
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open input file", ERROR);
  } else input = &cin;
}


static void loadName (const vector<Field2DF*>& u,
		      char*                    s)
// --------------------------------------------------------------------------
// Load a string containing the names of fields.
// ---------------------------------------------------------------------------
{
  int i, N = u.size();

  for (i = 0; i < N; i++) s[i] = u[i] -> getName();
  s[N] = '\0';
}


static int doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping of input is required to make sense of input.
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


static int getDump (istream&           ifile,
		    ostream&           ofile,
		    vector<Field2DF*>& u    )
// ---------------------------------------------------------------------------
// Read next set of field dumps from ifile, put headers on ofile.
// ---------------------------------------------------------------------------
{
  static char* hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields written\n",
    "%-25s "    "Format\n"
  };
  char buf[StrMax], fmt[StrMax], fields[StrMax];
  int  i, j, swab, nf, np, nz, nel;

  if (ifile.getline(buf, StrMax).eof()) return 0;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  ofile << buf << endl;
  ifile.getline (buf, StrMax);
  ofile << buf << endl;

  // -- I/O Numerical description of field sizes.

  ifile >> np >> nz >> nz >> nel;
  ifile.getline (buf, StrMax);
  
  sprintf (fmt, "%1d %1d %1d %1d", np, np, nz, nel);
  sprintf (buf, hdr_fmt[2], fmt);
  ofile << buf;

  for (i = 0; i < 5; i++) {
   ifile.getline (buf, StrMax);
   ofile << buf << endl;
  }

  // -- I/O field names.

  ifile >> fields;
  nf = strlen (fields);
  for (j = 0, i = 0; i < nf; i++) fmt[j++] = fields[i];
  fmt[j] = '\0';
  sprintf (buf, hdr_fmt[8], fmt);
  ofile << buf;
  ifile.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  ifile.getline (buf, StrMax);

  swab = doSwap (buf);

  sprintf (buf, "binary ");
  Veclib::describeFormat (buf + strlen (buf));
  sprintf (fmt, hdr_fmt[9], buf);
  ofile << fmt;

  if (u.size() != nf) {
    u.resize (nf);
    for (i = 0; i < nf; i++) u[i] = new Field2DF (np, nz, nel, fields[i]);
  }

  for (i = 0; i < nf; i++) {
    ifile >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return ifile.good();
}
