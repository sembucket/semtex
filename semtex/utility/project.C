///////////////////////////////////////////////////////////////////////////////
// project.C:  Project solution files to different interpolation orders.
//
// Copyright (c) 1996 Hugh Blackburn
//
// SYNOPSIS
// --------
// Process sem field file, project to new interpolation order on the
// same mesh.  Each dump in file is expected to be the same size.
// Also it is assumed that the field file represents a vector field
// dump, so that if the input file has N space dimensions, the first N
// fields represent vector components.  Input file must be binary
// format.
//
// USAGE
// -----
// project [options] [file]
// options:
// -h       ... print this message.
// -n <num> ... project elements to num x num.
// -z <num> ... project to <num> planes in the homogeneous direction.
// -u       ... project elements to uniform internal grid [Default: GLL].
// -U       ... project from uniform grid to GLL.
// 
// If file is not present, read from standard input.  Write to
// standard output.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>

static int uniform = 0;


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

  Field2DF& transform (const int);
  Field2DF& reverse   ();
  
private:
  const char    name;
  const int np, nz, nel, np2;
  int       nplane, ntot;
  real*         data;
  real**        plane;
};


Field2DF::Field2DF (const int  nP  ,
		    const int  nZ  ,
		    const int  nEl ,
		    const char Name) :

		    name       (Name),
                    np         (nP  ),
		    nz         (nZ  ),
		    nel        (nEl ),
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


Field2DF& Field2DF::transform (const int sign)
// ---------------------------------------------------------------------------
// Carry out Fourier transformation in z direction.
// ---------------------------------------------------------------------------
{
  if (nz  > 2) Femlib::DFTr (data, nz, nplane, sign);
  if (sign == INVERSE && nz == 2)
    Veclib::copy (nplane, plane[0], 1, plane[1], 1);

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

    register     int i, k;
    register     real    *LHS, *RHS;
    const        real       **IN, **IT;
    const int    nzm = min (rhs.nz, nz);
    vector<real> work (rhs.np * np);
    real*        tmp = work();

    if      (uniform == +1)
      Femlib::mesh (GLL, STD, rhs.np, np, 0, &IN, &IT, 0, 0);
    else if (uniform == -1) 
      Femlib::mesh (STD, GLL, rhs.np, np, 0, &IN, &IT, 0, 0);
    else
      Femlib::mesh (GLL, GLL, rhs.np, np, 0, &IN, &IT, 0, 0);

    for (k = 0; k < nzm; k++) {	// -- 2D planar projections.
      LHS = plane[k];
      RHS = rhs.plane[k];

      if (rhs.np == np)
	Veclib::copy (nplane, RHS, 1, LHS, 1);
      else
	for (i = 0; i < nel; i++, LHS += np2, RHS += rhs.np2) {
	  Blas::mxm (*IN, np, RHS, rhs.np, tmp, rhs.np);
	  Blas::mxm (tmp, np, *IT, rhs.np, LHS,     np);
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


static char prog[] = "project";
static void getargs  (int, char**, int&, int&, int&);
static int  getDump  (istream&, ostream&, vector<Field2DF*>&,
		      int&, int&, int&, int&, int&);
static void loadName (const vector<Field2DF*>&, char*);
static int  doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char              fields[StrMax];
  int               i, j, nEl, fInc;
  int               nPnew = 0, nZnew = 0, keepW = 0;
  vector<Field2DF*> Uold, Unew;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, nPnew, nZnew, keepW);
  
  while (getDump (cin, cout, Uold, nPnew, nZnew, nEl, keepW, fInc)) {

    Unew.setSize (Uold.getSize() + fInc);

    switch (fInc) {

    case 0:			// -- Same number of fields out as in.
      for (i = 0; i < Uold.getSize(); i++) {
	Unew[i] = new Field2DF (nPnew, nZnew, nEl, Uold[i] -> getName());
       *Unew[i] = *Uold[i];
      }
      break;

    case 1:			// -- Add a new blank field called 'w'.
      for (i = 0; i < Uold.getSize(); i++) {
	 Unew[i] = new Field2DF (nPnew, nZnew, nEl, Uold[i] -> getName());
	*Unew[i] = *Uold[i];
      }
       Unew[i] = new Field2DF (nPnew, nZnew, nEl, 'w');
      *Unew[i] = 0.0;
      break;

    case -1:			// -- Delete field called 'w'.
      loadName (Uold, fields);
      if (!strchr (fields, 'w'))
	message (prog, "conflict: 3D-->2D but no field called 'w'", ERROR);
      for (j = 0, i = 0; i < Uold.getSize(); i++) {
	if (Uold[i] -> getName() == 'w') continue;
	 Unew[j] = new Field2DF (nPnew, nZnew, nEl, Uold[i] -> getName());
	*Unew[j] = *Uold[i];
	j++;
      }
      break;

    default:			// -- Whoops.
      message (prog, "unrecognized conversion", ERROR);
      break;
    }

    for (i = 0; i < Unew.getSize(); i++) {
      Unew[i] -> transform (-1);
      cout << *Unew[i];
    }
  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     int&   np   ,
		     int&   nz   ,
		     int&   keepW)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: project [options] [file]\n"
    "  options:\n"
    "  -h       ... print this message\n"
    "  -n <num> ... 2D projection onto num X num\n"
    "  -z <num> ... 3D projection onto num planes\n"
    "  -w       ... Retain w components in 3D-->2D proj'n [Default: delete]\n"
    "  -u       ... project to uniform grid from GLL\n"
    "  -U       ... project from uniform grid to GLL\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0])
	np = atoi (*argv);
      else {
	--argc;
	np = atoi (*++argv);
      }
      break;
    case 'z':
      if (*++argv[0])
	nz = atoi (*argv);
      else {
	--argc;
	nz = atoi (*++argv);
      }
      break;
    case 'u':
      uniform =  1;
      break;
    case 'U':
      uniform = -1;
      break;
    case 'w':
      keepW   = 1;
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    ifstream* inputfile = new ifstream (*argv);
    if (inputfile -> good())
      cin = *inputfile;
    else
      message (prog, "unable to open input file", ERROR);
  }
}


static void loadName (const vector<Field2DF*>& u,
		      char*                    s)
// --------------------------------------------------------------------------
// Load a string containing the names of fields.
// ---------------------------------------------------------------------------
{
  int i, N = u.getSize();

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
		    vector<Field2DF*>& u    ,
		    int&               npnew,
		    int&               nznew,
		    int&               nel  ,
		    int&               keepW,
		    int&               finc )
// ---------------------------------------------------------------------------
// Read next set of field dumps from ifile, put headers on ofile.
//
// Convert to Fourier space.
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
  int  i, j, swab, nf, np, nz;

  if (ifile.getline(buf, StrMax).eof()) return 0;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  ofile << buf << endl;
  ifile.getline (buf, StrMax);
  ofile << buf << endl;

  // -- I/O Numerical description of field sizes.

  ifile >> np >> nz >> nz >> nel;
  ifile.getline (buf, StrMax);
  
  if (!npnew) npnew = np;
  if (!nznew) nznew = nz;

  sprintf (fmt, "%1d %1d %1d %1d", npnew, npnew, nznew, nel);
  sprintf (buf, hdr_fmt[2], fmt);
  cout << buf;

  if      (nz >  1 && nznew == 1 && !keepW) finc = -1;
  else if (nz == 1 && nznew >  1)           finc = +1;
  else                                      finc =  0;

  for (i = 0; i < 5; i++) {
   ifile.getline (buf, StrMax);
   ofile << buf << endl;
  }

  // -- I/O field names.

  ifile >> fields;
  nf = strlen (fields);
  if (finc == 1 && strchr (fields, 'w')) finc = 0;
  for (j = 0, i = 0; i < nf; i++) {
    if (finc == -1 && fields[i] == 'w') continue;
    fmt[j++] = fields[i];
  }
  if (finc == 1) fmt[j++] = 'w';
  fmt[j] = '\0';
  sprintf (buf, hdr_fmt[8], fmt);
  cout << buf;
  ifile.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  ifile.getline (buf, StrMax);

  swab = doSwap (buf);

  sprintf (buf, "binary ");
  Veclib::describeFormat (buf + strlen (buf));
  sprintf (fmt, hdr_fmt[9], buf);
  cout << fmt;

  if (u.getSize() != nf) {
    u.setSize (nf);
    for (i = 0; i < nf; i++) u[i] = new Field2DF (np, nz, nel, fields[i]);
  }

  for (i = 0; i < nf; i++) {
    ifile >> *u[i];
    if (swab) u[i] -> reverse();
    u[i] -> transform (+1);
  }

  return ifile.good();
}


