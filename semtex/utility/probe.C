///////////////////////////////////////////////////////////////////////////////
// probe.C: extract results from a field file at a set of 3D points.
//
// Synopsis
// --------
// Probe has three different user interfaces --- the internal
// mechanics of data extraction are the same in all cases, but the way
// the points are defined and data are output differ for each interface.
//
// Common information
// ------------------
// If the x--y values of a point cannot be located in the (2D) mesh,
// it is ignored, but if the z values falls outside [0, TWOPI/BETA],
// Fourier interpolation will be applied on assumption of periodicity.
//
// The field file must be in binary format, and the value of BETA in the
// session file will override the value given in the field file header.
//
// Interface 1: Extract data at set of points
// -----------
//
// Usage: probe [-h] [-p file] -s session dump [file]
//
// This is the most general form.  Extract data at a specified set
// of 3D points, given either on standard input or in a named file.
// Points must have three coordinates, and any number of them can be given,
// e.g.
//         1.32461       0.514135      1.00
//         1.31102       0.509459     -0.25
//            ..             ..         ..
//
// Points can either be supplied on standard input or in a named file.
//
// Output is always ASCII format.  Each line of output contains the values
// for the fields in the file in columns, in the order they were written
// to the field file.
//
// Interface 2: Extract data along a straight line
// -----------
//
// Usage: probeline [-h] -p "[n:]x0,y0,z0,dx,dy,dz" -s session dump
//
// Extract data along the straight line defined by the parameters.
// The number of points extracted is specified as <num>.  Output
// format is the same as for probe.
//
// Interface 3: Extract data at 2D array of points on x-y, x-z or y-z planes
// -----------
//
// Usage: probeplane [-h] [options] -s session dump
// options:
// -xy "x0,y0,dx,dy" ... xy-cutting plane
// -xz "x0,z0,dx,dz" ... xz-cutting plane
// -yz "y0,z0,dy,dz" ... yz-cutting plane
// -orig #           ... origin of the cutting plane along orthogonal axis
// -nx #             ... resolution along the x-axis
// -ny #             ... resolution along the y-axis
// -tec              ... write TECPLOT-formatted ASCII output
// -sm               ... write SM-formatted binary output
//
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <new.h>
#include <time.h>
#include <Sem.h>

static const int NPTS = 64;	// -- Default number of points for line/plane.

typedef enum {			// -- Flags for coordinate axes.
  None = 0,
  X    = 'x',
  Y    = 'y',
  Z    = 'z'
} AXIS;

static char *prog;
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }

static void getargs     (int, char**, char*&, char*&, char*&, char*&, char*&);
static int  loadPoints  (ifstream&, vector<Point*>&);
static int  linePoints  (vector<Point*>&);
static int  planePoints (vector<Point*>&, Mesh*);
static void findPoints  (vector<Point*>&, vector<Element*>&,
			 vector<Element*>&, vector<real>&, vector<real>&);
static int  getDump     (ifstream&, vector<AuxField*>&, vector<Element*>&,
			 const integer, const integer, const integer);
static void putData     (const char*, const char*, const char*,
			 int, vector<AuxField*>&,
			 vector<Element*>&, vector<Point*>&, matrix<real>&);
static void Finterp     (vector<AuxField*>&, const Point*, const Element*,
			 const real, const real, const integer, real*, real*);
static int  doSwap      (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char              *session, *dump, *format;
  char              *interface = 0, *points = 0;
  integer           NP, NZ,  NEL;
  integer           i, j, k, nf, ntot = 0, doff = 0, boff = 0;
  ifstream          fldfile, pntfile;
  FEML*             F;
  Mesh*             M;
  Point             lo, hi;
  const real*       knot;
  vector<real>      r, s, work, datum;
  vector<Point*>    point;
  vector<Element*>  elmt;
  vector<Element*>  Esys;
  vector<AuxField*> u;
  matrix<real>      data;

  // -- Initialize.

  prog = *argv;
  set_new_handler    (&memExhaust);
  Femlib::initialize (&argc, &argv);

  // -- Set defaults for probeplane interface.

  Femlib::value ("SIZED"  , None);
  Femlib::value ("NX"     , NPTS);
  Femlib::value ("NY"     , NPTS);
  Femlib::value ("ORTHO"  , Z   );
  Femlib::value ("OFFSET" , 0.0 );
  Femlib::value ("X_MIN"  , 0.0 );
  Femlib::value ("Y_MIN"  , 0.0 );
  Femlib::value ("Z_MIN"  , 0.0 );
  Femlib::value ("X_DELTA", 0.0 );
  Femlib::value ("Y_DELTA", 0.0 );
  Femlib::value ("Z_DELTA", 0.0 );

  // -- Parse command line.

  getargs (argc, argv, interface, format, session, dump, points);

  // -- Check presence of field file before proceeding.

  fldfile.open (dump, ios::in);
  if (!fldfile) message (prog, "no field file", ERROR);
  
  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh   (*F);

  NEL = M -> nEl();  
  NP  = (integer) Femlib::value ("N_POLY");
  NZ  = (integer) Femlib::value ("N_Z"   );
  
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

  // -- Construct list of points.

  if (strcmp (interface, "probe") == 0) {
    if   (points) pntfile.open   (points, ios::in);
    else          pntfile.attach (0);
    ntot = loadPoints (pntfile, point);
  } 
  else if (strcmp (interface, "probeline")  == 0) ntot = linePoints  (point);
  else if (strcmp (interface, "probeplane") == 0) ntot = planePoints (point,M);

  // -- Locate points in the Mesh.

  findPoints (point, Esys, elmt, r, s);

  // -- Load field file.

  if (!(getDump (fldfile, u, Esys, NP, NZ, NEL)))
    message (prog, "no data extracted", ERROR);

  datum.setSize (nf = u.getSize());
  data.setSize  (ntot, nf);

  // -- Interpolate within it.

  for (i = 0; i < ntot; i++)
    if (elmt[i]) {
      Finterp (u, point[i], elmt[i], r[i], s[i], NZ, work(), datum()); 
      for (j = 0; j < nf; j++) data (i, j) = datum (j);
    } else
      for (j = 0; j < nf; j++) data (i, j) = 0.0;

  // -- Output collected data.

  putData (session, interface, format, ntot, u, elmt, point, data);
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (integer argc     ,
		     char**  argv     ,
		     char*&  interface,
		     char*&  format   ,
		     char*&  session  ,
		     char*&  dump     ,
		     char*&  points   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  format = new char [16];
  strcpy (format, "free");	// -- Default output format.
  interface = *argv;

  if (strcmp (interface, "probe") == 0) {

    char usage[] =
      "Usage: probe [options] -s session dump\n"
      "  options:\n"
      "  -h      ... print this message\n"
      "  -p file ... name file of point data [Default: stdin]\n";

    while (--argc  && **++argv == '-')
      switch (*++argv[0]) {
      case 'h':
	cout << usage;
	exit (EXIT_SUCCESS);
	break;
      case 'p':
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

  } else if (strcmp (interface, "probeline") == 0) {

    char usage[] =
      "Usage: probeline [-h] -p \"[n:]x0,y0,z0,dx,dy,dz\" -s session dump\n";
    char *tok, *pspec;
    int set = 0;

    while (--argc && **++argv == '-')
      switch (*++argv[0]) {
      case 'h':
	cout << usage;
	exit (EXIT_SUCCESS);
	break;
      case 'p':
	if (*++argv[0])
	  pspec = *argv;
	else {
	  --argc;
	  pspec = *++argv;
	}
	if (strchr (pspec, ':')) {
	  if   (tok = strtok (pspec, ":")) Femlib::value ("NPTS", atof (tok));
	  else                             Femlib::value ("NPTS", NPTS);
	  if (tok = strtok (0, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    set = 1;
	  } else {
	    message (prog, "couldn't x0 parse number from string", ERROR);
	  }
	} else {
	  Femlib::value ("NPTS", NPTS);
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    set = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	}
	while (tok = strtok (0, ","))
	  switch (++set) {
	  case 2:
	    Femlib::value ("Y_MIN",   atof (tok));
	    break;
	  case 3:
	    Femlib::value ("Z_MIN",   atof (tok));
	    break;
	  case 4:
	    Femlib::value ("X_DELTA", atof (tok));
	    break;
	  case 5:
	    Femlib::value ("Y_DELTA", atof (tok));
	    break;
	  case 6:
	    Femlib::value ("Z_DELTA", atof (tok));
	    break;
	  default:
	    message (prog, "too many numbers in point string", ERROR);
	    break;
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
    if (set != 6) {
      message (prog, "wrong number of parameters to line string", ERROR);
    }

  } else if (strcmp (interface, "probeplane") == 0) {

    char *tok, *pspec, *usage =
      "Usage: probeplane [options] -s session dump\n"
      "options:\n"
      "-h                ... print this message\n"
      "-xy \"x0,y0,dx,dy\" ... xy-cutting plane\n"
      "-xz \"x0,z0,dx,dz\" ... xz-cutting plane\n"
      "-yz \"y0,z0,dy,dz\" ... yz-cutting plane\n"
      "-orig #           ... origin of the cutting plane along ortho axis\n"
      "-nx #             ... resolution along the x-axis\n"
      "-ny #             ... resolution along the y-axis\n"
      "-tec              ... write TECPLOT-formatted ASCII output\n"
      "-sm               ... write SM-formatted binary output\n";
    int nset = 0;

    while (--argc && **++argv == '-')
      switch (*++argv[0]) {
      case 'h':
	cout << usage;
	exit (EXIT_SUCCESS);
	break;
      case 'n':
	switch (argv[0][1]) {
	case 'x':
	  --argc;
	  Femlib::value ("NX", atof (*++argv));
	  break;
	case 'y':
	  --argc;
	  Femlib::value ("NY", atof (*++argv));
	  break;
	default:
	  message (prog, "can only specify nx or ny", ERROR);
	  break;
	}
      case 'o':
	--argc;
	Femlib::value ("OFFSET", atof (*++argv));
	break;
      case 'p':
	if (*++argv[0])
	  points = *argv;
	else {
	  --argc;
	  points = *++argv;
	}
	break;
      case 's':
	if (argv[0][1] == 'm') strcpy (format, "sm");
	else { --argc; session = *++argv; }
	break;
      case 't':
	strcpy (format, "tecplot");
	break;
      case 'x':
	switch (argv[0][1]) {
	case 'y':
	  Femlib::value ("ORTHO", Z);
	  --argc; pspec = *++argv;
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    nset = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	  while (tok = strtok (0, ","))
	    switch (++nset) {
	    case 2:
	      Femlib::value ("Y_MIN",   atof (tok));
	      break;
	    case 3:
	      Femlib::value ("X_DELTA", atof (tok));
	      break;
	    case 4:
	      Femlib::value ("Y_DELTA", atof (tok));
	      break;
	    default:
	      message (prog, "too many numbers in point string", ERROR);
	      break;
	    }
	  if (nset != 4)
	    message (prog, "need 4 parameters for cutting plane", ERROR);
	  Femlib::value ("SIZED", 1.0);
	  break;
	case 'z':
	  Femlib::value ("ORTHO", Y);
	  --argc; pspec = *++argv;
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("X_MIN", atof (tok));
	    nset = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	  while (tok = strtok (0, ","))
	    switch (++nset) {
	    case 2:
	      Femlib::value ("Z_MIN",   atof (tok));
	      break;
	    case 3:
	      Femlib::value ("X_DELTA", atof (tok));
	      break;
	    case 4:
	      Femlib::value ("Z_DELTA", atof (tok));
	      break;
	    default:
	      message (prog, "too many numbers in point string", ERROR);
	      break;
	    }
	  if (nset != 4)
	    message (prog, "need 4 parameters for cutting plane", ERROR);
	  Femlib::value ("SIZED", 1.0);
	  break;
	default:
	  message (prog, "extents can be xy, xz, or yz", ERROR);
	  break;
	}
	break;
      case 'y':
	switch (argv[0][1]) {
	case 'z':
	  Femlib::value ("ORTHO", X);
	  --argc; pspec = *++argv;
	  if (tok = strtok (pspec, ",")) {
	    Femlib::value ("y_MIN", atof (tok));
	    nset = 1;
	  } else {
	    message (prog, "couldn't parse number x0 from string", ERROR);
	  }
	  while (tok = strtok (0, ","))
	    switch (++nset) {
	    case 2:
	      Femlib::value ("Z_MIN",   atof (tok));
	      break;
	    case 3:
	      Femlib::value ("Y_DELTA", atof (tok));
	      break;
	    case 4:
	      Femlib::value ("Z_DELTA", atof (tok));
	      break;
	    default:
	      message (prog, "too many numbers in point string", ERROR);
	      break;
	    }
	  if (nset != 4)
	    message (prog, "need 4 parameters for cutting plane", ERROR);
	  Femlib::value ("SIZED", 1.0);
	  break;
	default:
	  message (prog, "extents can be xy, xz, or yz", ERROR);
	  break;
	}
	break;
      default:
	cerr << usage;
	exit (EXIT_FAILURE);
	break;
      }

  }

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static int loadPoints (ifstream&       pfile,
                       vector<Point*>& point)
// ---------------------------------------------------------------------------
// Probe point input for the "probe" interface, from file.
// ---------------------------------------------------------------------------
{
  int           ntot, num = 0;
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

  return ntot;
}


static int linePoints (vector<Point*>& point)
// ---------------------------------------------------------------------------
// Probe point generation for the "probeline" interface.
// ---------------------------------------------------------------------------
{
  int  i, ntot = (int) Femlib::value ("NPTS");
  real xmin = Femlib::value ("X_MIN");
  real ymin = Femlib::value ("Y_MIN");
  real zmin = Femlib::value ("Z_MIN");
  real dx   = (ntot == 1) ? 0.0 : Femlib::value ("X_DELTA") / (ntot - 1.0);
  real dy   = (ntot == 1) ? 0.0 : Femlib::value ("Y_DELTA") / (ntot - 1.0);
  real dz   = (ntot == 1) ? 0.0 : Femlib::value ("Z_DELTA") / (ntot - 1.0);

  point.setSize (ntot);

  for (i = 0; i < ntot; i++) {
    point[i] = new Point;
    point[i] -> x = xmin + i * dx;
    point[i] -> y = ymin + i * dy;
    point[i] -> z = zmin + i * dz;
  }

  return ntot;
}


static int planePoints (vector<Point*>& point,
			Mesh*           mesh )
// ---------------------------------------------------------------------------
// Probe point generation for the "probeplane" interface.
// ---------------------------------------------------------------------------
{
  int        i, j, k;
  const int  nx     = (int) Femlib::value ("NX");
  const int  ny     = (int) Femlib::value ("NY");
  const int  ortho  = (int) Femlib::value ("ORTHO");
  const real offset =       Femlib::value ("OFFSET");
  const int  ntot   = nx * ny;
  real       x0, y0, z0, dx, dy, dz;
  Point*     p;

  point.setSize (ntot);

  switch (ortho) {
  case X:
    x0 = offset;
    y0 = Femlib::value ("Y_MIN");
    z0 = Femlib::value ("Z_MIN");
    dy = Femlib::value ("Y_DELTA") / (nx - 1.0);
    dz = Femlib::value ("Z_DELTA") / (ny - 1.0);
    for (k = 0, j = 0; j < ny; j++)
      for (i = 0; i < nx; i++, k++) {
	point[k] = p = new Point;
	p -> x = x0;
	p -> y = y0 + i * dy;
	p -> z = z0 + j * dz;
      }
    break;
  case Y:
    x0 = Femlib::value ("X_MIN");
    y0 = offset;
    z0 = Femlib::value ("Z_MIN");
    dx = Femlib::value ("X_DELTA") / (nx - 1.0);
    dz = Femlib::value ("Z_DELTA") / (ny - 1.0);
    for (k = 0, j = 0; j < ny; j++)
      for (i = 0; i < nx; i++, k++) {
	point[k] = p = new Point;
	p -> x = x0 + i * dx;
	p -> y = y0;
	p -> z = z0 + j * dz;
      }
    break;
  case Z:
    if (!((int) Femlib::value ("SIZED"))) {
      Point lo, hi;
      mesh -> extent (lo, hi);
      Femlib::value ("X_MIN", lo.x);
      Femlib::value ("Y_MIN", lo.y);
      Femlib::value ("X_DELTA", hi.x - lo.x);
      Femlib::value ("Y_DELTA", hi.y - lo.y);
    }
    x0 = Femlib::value ("X_MIN");
    y0 = Femlib::value ("Y_MIN");
    dx = Femlib::value ("X_DELTA") / (nx - 1.0);
    dy = Femlib::value ("Y_DELTA") / (ny - 1.0);
    z0 = offset;
    for (k = 0, j = 0; j < ny; j++)
      for (i = 0; i < nx; i++, k++) {
	point[k] = p = new Point;
	p -> x = x0 + i * dx;
	p -> y = y0 + j * dy;
	p -> z = z0;
      }
    break;
  default:
    break;
  }
  
  return ntot;
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
  integer       i, k;
  real          x, y, z, r, s;
  const integer NEL   = Esys .getSize();
  const integer NPT   = point.getSize();
  const integer guess = 1;

  elmt.setSize (NPT);
  rloc.setSize (NPT);
  sloc.setSize (NPT);

  elmt = 0;

  cerr.precision (8);

  for (i = 0; i < NPT; i++) {
    x = point[i] -> x;
    y = point[i] -> y;
    z = point[i] -> z;
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
		     const integer      NZ  ,
		     real*              work,
		     real*              data)
// ---------------------------------------------------------------------------
// Carry out 2DxFourier interpolation.
// ---------------------------------------------------------------------------
{
  register integer i, k, Re, Im;
  register real    phase;
  const integer    NF    = u.getSize();
  const integer    NZH   = NZ >> 1;
  const integer    NHM   = NZH - 1;
  const real       betaZ = P -> z * Femlib::value("BETA");
  const real*      Wtab  = work + NZ;

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


static void putData (const char*        session  ,
		     const char*        interface,
		     const char*        format   ,
		     int                ntot     ,
		     vector<AuxField*>& u        ,
		     vector<Element*>&  elmt     ,
		     vector<Point*>&    point    ,
		     matrix<real>&      data     )
// ---------------------------------------------------------------------------
// Handle all the different output formats, according to the chosen interface.
// ---------------------------------------------------------------------------
{
  int i, j, k, n, nf = u.getSize();

  if (strstr (format, "free")) {
    cout.precision (6);
    for (i = 0; i < ntot; i++) {
      if (elmt[i]) {
	cout << setw (5) << i + 1 << " " 
	     << setw(12) << point[i] -> x << " " 
	     << setw(12) << point[i] -> y << " " 
	     << setw(12) << point[i] -> z;
	for (j = 0; j < nf; j++)
	  cout << setw(15) << data (i, j);
	cout << endl;
      }
    }
    return;

  } else if (strstr (format, "sm") && strstr (interface, "probeplane")) {

    char      fname[StrMax];
    const int nx = (int) Femlib::value ("NX");
    const int ny = (int) Femlib::value ("NY");
    ofstream  out;

    for (n = 0; n < nf; n++) {
      sprintf (fname, "%s.%c", session, u[n] -> name());
      out.open  (fname);
      out.write ((char*) &nx, sizeof (int));
      out.write ((char*) &ny, sizeof (int));
      
      for (k = 0, j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++, k++) {
	  float tmp = (float) data (k, n);
	  out.write ((char*) &tmp, sizeof (float));
	}
      }
      out.close();
    }

  }
}
