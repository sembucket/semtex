///////////////////////////////////////////////////////////////////////////////
// interp.C: interpolate results from a field file onto a set of 2D points.
//
// Copyright (c) 1997--1999 Hugh Blackburn
//
// Synopsis:
// --------
// interp [-h] [-v] [-m file] -s session dump [file]
//
// Description:
// -----------
// Interpolation is 2D and for 3D fields the data will be output in
// plane-by-plane order.  Output is always ASCII format.  Each line
// of output contains the values for the fields in the file in columns,
// in the order they were written to the field file.  The field file
// must be in binary format.
//
// The set of points can either be in the form output from meshpr, e.g.
// 12 12 4 422 NR NS NZ NEL
//         1.32461       0.514135
//         1.31102       0.509459
//            ..             ..
// in which case the output data will be in field dump format with header,
// or the input can be an (unstructured) set without a header, e.g.
//         1.32461       0.514135
//         1.31102       0.509459
//            ..             ..
// in which case the output has a matching lack of structure.
//
// If a point cannot be located in the mesh, zero values are output for
// that point location.  Points can either be supplied on standard input
// or in a named file.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <new.h>
#include <time.h>
#include <Sem.h>

static char prog[]  = "interp";
static integer  verbose = 0;
static integer  nreport = 100;

static void  memExhaust () { message ("new", "free store exhausted", ERROR); }

static void    getargs    (integer, char**, char*&, char*&, char*&);
static void    loadPoints (ifstream&, integer&, integer&,
			   integer&, vector<Point*>&);
static void    findPoints (vector<Point*>&, vector<Element*>&,
			   vector<Element*>&, vector<real>&, vector<real>&);
static integer getDump    (ifstream&, vector<AuxField*>&, vector<Element*>&,
			   const integer, const integer, const integer,
			   integer&, real&, real&, real&, real&);
static void    putHeader  (const char*, const vector<AuxField*>&,
			   const integer, const integer, const integer,
			   const integer, const real, const real,
			   const real, const real);
static integer doSwap     (const char*);
static void    loadName   (const vector<AuxField*>&, char*);


integer main (integer argc,
	      char**  argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char              *session, *dump, *points = 0;
  integer           NP, NZ,  NEL;
  integer           np, nel, ntot;
  integer           i, j, k, nf, step;
  ifstream          fldfile, pntfile;
  FEML*             F;
  Mesh*             M;
  real              c, time, timestep, kinvis, beta;
  const real*       z;
  vector<real>      r, s;
  vector<Point*>    point;
  vector<Element*>  elmt;
  vector<Element*>  Esys;
  vector<AuxField*> u;

  // -- Initialize.

  set_new_handler    (&memExhaust);
  Femlib::initialize (&argc, &argv);
  getargs            (argc, argv, session, dump, points);

  fldfile.open (dump, ios::in);
  if (!fldfile) message (prog, "no field file", ERROR);

  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh (F);

  NEL = M -> nEl();  
  NP  = (integer) Femlib::value ("N_POLY");
  NZ  = (integer) Femlib::value ("N_Z"   );
  
  Geometry::set (NP, NZ, NEL, Geometry::Cartesian);
  Femlib::mesh  (GLL, GLL, NP, NP, &z, 0, 0, 0, 0);
  Esys.setSize  (NEL);

  for (k = 0; k < NEL; k++) Esys[k] = new Element (k, M, z, NP);
  
  // -- Construct the list of points, then find them in the Mesh.

  if   (points) pntfile.open   (points, ios::in);
  else          pntfile.attach (0);

  np = nel = ntot = 0;

  loadPoints (pntfile, np, nel, ntot, point);
  findPoints (point, Esys, elmt, r, s);

  // -- Load field file, interpolate within it.

  cout.precision (8);
  while (getDump (fldfile, u, Esys, NP, NZ, NEL,
		  step, time, timestep, kinvis, beta)) {

    if (np) putHeader (session, u,  np, NZ, nel,
		       step, time, timestep, kinvis, beta);

    nf = u.getSize();
    for (k = 0; k < NZ; k++)
      for (i = 0; i < ntot; i++) {
	for (j = 0; j < nf; j++) {
	  if   (elmt[i]) c = u[j] -> probe (elmt[i], r[i], s[i], k);
	  else           c = 0.0;
	  cout << setw(15) <<  c;
	}
	if (verbose && !((i + 1)% nreport))
	  cerr 
	    << "interp: plane " << k + 1 << ", "
	    << i + 1 << " points interpolated" << endl;
	cout << endl;
      }
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (integer argc   ,
		     char**  argv   ,
		     char*&  session,
		     char*&  dump   ,
		     char*&  points )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: interp [options] -s session dump\n"
    "  options:\n"
    "  -h      ... print this message\n"
    "  -m file ... name file of point data [Default: stdin]\n"
    "  -v      ... verbose output\n";
 
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
    case 'v':
      verbose = 1;
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
			integer&        np   ,
			integer&        nel  ,
			integer&        ntot ,
			vector<Point*>& point)
// ---------------------------------------------------------------------------
// Load data which describe location of points.
// ---------------------------------------------------------------------------
{
  char          buf[StrMax];
  integer       nz = 0, num = 0;
  real          x, y;
  Point*        datum;
  Stack<Point*> data;

  pfile.getline (buf, StrMax);
  if (strstr (buf, "NEL")) {	// -- This is a structured set of points.
    sscanf (buf, "%d %d %d %d", &np, &np, &nz, &nel);
    ntot    = np * np * nel;
    nreport = np * np;

  } else {
    if   (sizeof (real) == sizeof (double)) sscanf (buf, "%lf %lf", &x, &y);
    else                                    sscanf (buf, "%f  %f",  &x, &y);

    datum = new Point;
    datum -> x = x;
    datum -> y = y;
    data.push (datum);
    num++;
  }

  while (pfile >> x >> y) {
    datum = new Point;
    datum -> x = x;
    datum -> y = y;
    data.push (datum);
    num++;
    if (nz && num == ntot) break;
  }

  if (np && num != ntot) {
    sprintf (buf, "No. of points (%1d) mismatches declaration (%1d)",num,ntot);
    message (prog, buf, ERROR);
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
  integer       i, k, kold;
  real          x, y, r, s;
  const integer guess = 1;
  const integer NEL   = Esys .getSize();
  const integer NPT   = point.getSize();

  elmt.setSize (NPT);
  rloc.setSize (NPT);
  sloc.setSize (NPT);

  elmt = 0;

  cerr.precision (8);

  kold = 0;
  for (i = 0; i < NPT; i++) {
    x = point[i] -> x;
    y = point[i] -> y;
    if (Esys[kold] -> locate (x, y, r, s, guess)) {
      elmt[i] = Esys[kold];
      rloc[i] = r;
      sloc[i] = s;
    } else {
      for (k = 0; k < NEL; k++) {
	if (k == kold) continue;
	if (Esys[k] -> locate (x, y, r, s, guess)) {
	  elmt[i] = Esys[k];
	  rloc[i] = r;
	  sloc[i] = s;
	  kold    = k;
	  break;
	}
      }
    }

    if (!elmt[i])
      cerr << "interp: point (" << setw(15) << x << ","
	   << setw(15) << y << ") is not in the mesh" << endl;
    if (verbose && !((i + 1)% nreport))
	cerr << "interp: " << i + 1 << " points found" << endl;
  }
  if (verbose) cerr << endl;
}


static integer getDump (ifstream&          file,
			vector<AuxField*>& u   ,
			vector<Element*>&  Esys,
			const integer      np  ,
			const integer      nz  ,
			const integer      nel ,
			integer&           step,
			real&              time,
			real&              tstp,
			real&              kinv,
			real&              beta)
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

  file >> step;
  file.getline (buf, StrMax);

  file >> time;
  file.getline (buf, StrMax);

  file >> tstp;
  file.getline (buf, StrMax);

  file >> kinv;
  file.getline (buf, StrMax);

  file >> beta;
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
    for (i = 0; i < nf; i++)
      u[i] = new AuxField (new real[Geometry::nTotal()], nz, Esys, fields[i]);
  } else if (u.getSize() != nf) 
    message (prog, "number of fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (i = 0; i < nf; i++) {
    file >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return file.good();
}


static void putHeader (const char*              session,
		       const vector<AuxField*>& u      ,
		       const integer            np     ,
		       const integer            nz     ,
		       const integer            nel    ,
		       const integer            step   ,
		       const real               time   ,
		       const real               tstp   ,
		       const real               kinv   ,
		       const real               beta   )
// ---------------------------------------------------------------------------
// Write header information for semtex/prism/nekton compatible file on cout.
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
  char   buf[StrMax], fmt[StrMax], fields[StrMax];
  time_t tp (::time (0));

  sprintf (buf, hdr_fmt[0], session);
  cout << buf;

  strftime (fmt, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (buf, hdr_fmt[1], fmt);
  cout << buf;

  sprintf (fmt, "%1d %1d %1d %1d", np, np, nz, nel);
  sprintf (buf, hdr_fmt[2], fmt);
  cout << buf;

  sprintf (buf, hdr_fmt[3], step);
  cout << buf;

  sprintf (buf, hdr_fmt[4], time);
  cout << buf;

  sprintf (buf, hdr_fmt[5], tstp);
  cout << buf;

  sprintf (buf, hdr_fmt[6], kinv);
  cout << buf;

  sprintf (buf, hdr_fmt[7], beta);
  cout << buf;

  loadName (u, fields);
  sprintf  (buf, hdr_fmt[8], fields);
  cout << buf;
  
  sprintf (buf, hdr_fmt[9], "ASCII");
  cout << buf;
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


static void loadName (const vector<AuxField*>& u,
		      char*                    s)
// --------------------------------------------------------------------------
// Load a string containing the names of fields.
// ---------------------------------------------------------------------------
{
  integer i, N = u.getSize();

  for (i = 0; i < N; i++) s[i] = u[i] -> name();
  s[N] = '\0';
}
