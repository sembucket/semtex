///////////////////////////////////////////////////////////////////////////////
// probe.C: extract results from a field file at a set of 3D points.
//
// Synopsis:
// --------
// probe [-h] [-m file] -s session dump [file]
//
// Description:
// -----------
// Points must have three coordinates, and any number of them can be given,
// e.g.
//         1.32461       0.514135      1.00
//         1.31102       0.509459     -0.25
//            ..             ..         ..
//
// If the x--y values of a point cannot be located in the (2D) mesh,
// it is ignored, but if the z values falls outside [0, TWOPI/BETA],
// Fourier interpolation will be applied on assumption of periodicity.
//
// Points can either be supplied on standard input or in a named file.
// The field file must be in binary format, and the value of BETA in the
// session file will override the value given in the field file header.
//
// Output is always ASCII format.  Each line of output contains the values
// for the fields in the file in columns, in the order they were written
// to the field file.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <new.h>
#include <time.h>
#include <Sem.h>
#include <Stack.h>

static char prog[] = "probe";
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }

static void getargs    (int, char**, char*&, char*&, char*&);
static void loadPoints (ifstream&, int&, vector<Point*>&);
static void findPoints (vector<Point*>&, vector<Element*>&, vector<Element*>&,
			vector<real>&,   vector<real>&);
static int  getDump    (ifstream&, vector<AuxField*>&, vector<Element*>&,
			const int, const int, const int);
static int  doSwap     (const char*);
static void Finterp    (vector<AuxField*>&, const Point*, const Element*,
			const real, const real, const int, real*, real*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char              *session, *dump, *points = 0;
  int               NP, NZ,  NEL;
  int               i, j, k, nf, ntot, doff = 0, boff = 0;
  ifstream          fldfile, pntfile;
  FEML*             F;
  Mesh*             M;
  const real*       knot;
  vector<real>      r, s, work, data;
  vector<Point*>    point;
  vector<Element*>  elmt;
  vector<Element*>  Esys;
  vector<AuxField*> u;

  // -- Initialize.

  set_new_handler (&memExhaust);
  Femlib::prep    ();
  getargs         (argc, argv, session, dump, points);

  fldfile.open (dump, ios::in);
  if (!fldfile) message (prog, "no field file", ERROR);

  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh   (*F);

  NEL = M -> nEl();  
  NP  = (int) Femlib::value ("N_POLY");
  NZ  = (int) Femlib::value ("N_Z"   );
  
  Geometry::set (NP, NZ, NEL, Geometry::Cartesian);
  Femlib::mesh  (GLL, GLL, NP, NP, &knot, 0, 0, 0, 0);
  Esys.setSize  (NEL);

  for (k = 0; k < NEL; k++) {
    Esys[k] = new Element (k, *M, knot, NP, doff, boff);
    doff   += Esys[k] -> nTot();
    boff   += Esys[k] -> nExt();
  }

  // -- Set up FFT work areas.

  work.setSize  (3*NZ + 15);
  Femlib::rffti (NZ, work() + NZ);
  
  // -- Construct the list of points, then find them in the Mesh.

  if   (points) pntfile.open   (points, ios::in);
  else          pntfile.attach (0);

  ntot = 0;

  loadPoints (pntfile, ntot, point);
  findPoints (point, Esys, elmt, r, s);

  // -- Load field file, interpolate within it.

  cout.precision (8);

  while (getDump (fldfile, u, Esys, NP, NZ, NEL)) {

    data.setSize (nf = u.getSize());

    for (i = 0; i < ntot; i++)
      if (elmt[i]) {
	Finterp (u, point[i], elmt[i], r[i], s[i], NZ, work(), data()); 
	for (j = 0; j < nf; j++) cout << setw(15) << data[j];
	cout << endl;
      }
  }
  
  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session,
		     char*& dump   ,
		     char*& points )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: probe [options] -s session dump\n"
    "  options:\n"
    "  -h      ... print this message\n"
    "  -m file ... name file of point data [Default: stdin]\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'm':
      if (*++argv[0])
	points = *argv;
      else {
	--argc;
	points = *++argv;
      }
      break;
    case 's':
      if (*++argv[0])
	session = *argv;
      else {
	--argc;
	session = *++argv;
      }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static void loadPoints (ifstream&       pfile,
			int&            ntot ,
			vector<Point*>& point)
// ---------------------------------------------------------------------------
// Load data which describe location of points.
// ---------------------------------------------------------------------------
{
  int           num = 0;
  real          x, y, z;
  Point*        datum;
  Stack<Point*> data;

  while (pfile >> x >> y >> z) {
    datum = new Point;
    datum -> x = x;
    datum -> y = y;
    datum -> z = z;
    data.push (datum);
    num++;
  }

  ntot = num;
  point.setSize (ntot);

  while (num--) point[num] = data.pop();
}


static void findPoints (vector<Point*>&   point,
			vector<Element*>& Esys ,
			vector<Element*>& elmt ,
			vector<real>&     rloc ,
			vector<real>&     sloc )
// ---------------------------------------------------------------------------
// Locate points within elements, set Element pointer & r--s locations.
// ---------------------------------------------------------------------------
{
  int       i, k;
  real      x, y, z, r, s;
  const int NEL   = Esys .getSize();
  const int NPT   = point.getSize();
  const int guess = 1;

  elmt.setSize (NPT);
  rloc.setSize (NPT);
  sloc.setSize (NPT);

  elmt = 0;

  cerr.precision (8);

  for (i = 0; i < NPT; i++) {
    x = point[i] -> x;
    y = point[i] -> y;
    for (k = 0; k < NEL; k++) {
      r = s = 0.0;
      if (Esys[k] -> locate (x, y, r, s, guess)) {
	elmt[i] = Esys[k];
	rloc[i] = r;
	sloc[i] = s;
	break;
      }
    }

    if (!elmt[i])
      cerr << "point ("
	   << setw(15) << x << ","
	   << setw(15) << y << ","
	   << setw(15) << z << ") is not in the mesh" << endl;
  }
}


static int getDump (ifstream&          file,
		    vector<AuxField*>& u   ,
		    vector<Element*>&  Esys,
		    const int          np  ,
		    const int          nz  ,
		    const int          nel )
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  char buf[StrMax], fields[StrMax];
  int  i, swab, nf, npnew, nznew, nelnew;

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

  if (u.getSize() == 0) {
    u.setSize (nf);
    for (i = 0; i < nf; i++) u[i] = new AuxField (Esys, fields[i]);
  } else if (u.getSize() != nf) 
    message (prog, "number of fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (i = 0; i < nf; i++) {
    file >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return file.good();
}


static int doSwap (const char* ffmt)
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


static void Finterp (vector<AuxField*>& u   ,
		     const Point*       P   ,
		     const Element*     E   ,
		     const real         r   , 
		     const real         s   ,
		     const int          NZ  ,
		     real*              work,
		     real*              data)
// ---------------------------------------------------------------------------
// Carry out 2DxFourier interpolation.
// ---------------------------------------------------------------------------
{
  register int  i, k, Re, Im;
  register real phase;
  const int     NF    = u.getSize();
  const int     NZH   = NZ >> 1;
  const int     NHM   = NZH - 1;
  const real    betaZ = P -> z * Femlib::value("BETA");
  const real*   Wtab  = work + NZ;


  for (i = 0; i < NF; i++)	// -- For each field.

    if (NZ == 1)

      // -- Just 2D.

      data[i] = u[i] -> probe (E, r, s, 0);

    else {

      // -- 2D interpolation.
      
      for (k = 0; k < NZ; k++) work[k] = u[i] -> probe (E, r, s, k);

      // -- Fourier interpolation.

      Femlib::rfftf (NZ, work, Wtab);
      Blas::scal    (NZ, 2.0/NZ, work, 1);

      if (NZ & 1) {

	data[i] = 0.5 * work[0];
	for (k = 1; k < NZH; k++) {
	  Im       = k  + k;
	  Re       = Im - 1;
	  phase    = k * betaZ;
	  data[i] += work[Re] * cos (phase) - work[Im] * sin (phase);
	}

      } else {

	data[i] = 0.5 * work[0];
	for (k = 1; k < NHM; k++) {
	  Im       = k  + k;
	  Re       = Im - 1;
	  phase    = k * betaZ;
	  data[i] += work[Re] * cos (phase) - work[Im] * sin (phase);
	}
	data[i] += 0.5 * work[NZ - 1] * cos (NZH * betaZ);

      }
    }
}
