///////////////////////////////////////////////////////////////////////////////
// Input SEM mesh and field data.
//
// Copyright (C) 1999 Hugh Blackburn
// 
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sview.h>

static int  _swab;		// -- File-scope flag for byte-swapped input.
static int  iformat ();
static void format  (char*);
static void dbrev   (const int, const double*, double*);


Sem* loadMesh (const char* fname)
// ---------------------------------------------------------------------------
// Mesh header is of form:
//
// 9 9 32 206 NR NS NZ NEL
//
// then follows (x,y) locations of each element's 2D mesh points in
// row-major order, ASCII format.  In all local implementations,
// NR=NS=np.  Finally the NZ + 1 z locations of the 3D mesh are
// supplied.
//
// At present, each element is assumed the same order in each
// direction, and of equal orders also in the r & s directions.  This
// would be fairly easy to change here, and should carry through OK to
// the rest of the code.
//
// NB: the mean values of extent in the x, y, & z directions are subtracted
// from the mesh locations.
// ---------------------------------------------------------------------------
{
  const char      routine[] = "loadMesh";
  char            buf[StrMax], err[StrMax];
  int             nr, np, np2, nz, nel, ntot;
  float           x, y, z;
  float           xavg, yavg, zavg;
  float           xmin, xmax, ymin, ymax, zmin, zmax;
  ifstream        mfile (fname);
  Sem*            M = new Sem;
  register int    i, j, k;
  register float  *xg, *yg, *zg, r, cz, sz;

  if (!mfile) message (routine, "couldn't open mesh file", ERROR);
  
  mfile >> nr >> np >> nz >> nel;
  mfile.getline (buf, StrMax);

  if (!strstr (buf, "NR NS NZ NEL")) {
    sprintf (err, "mesh header line should include NR NS NZ NEL: %s", buf);
    message (routine, err, ERROR);
  }
  if (nr != np) {
    sprintf (err, "Element NR, NS orders must be equal: %1d: %1d", nr, np);
    message (routine, err, ERROR);
  }
  if (nz < 2) {
    sprintf (err, "Mesh mush be 3D (NZ >= 2): %1d", nz);
    message (routine, err, ERROR);
  }

  ++nz;				// -- Mesh will show full length of domain.
  M -> nrep = 1;		// -- Special default.
 
  M -> nel  = nel;
  M -> idim = new int [nel];
  M -> jdim = new int [nel];
  M -> kdim = new int [nel];

  for (i = 0; i < nel; i++) {
    M -> idim[i] = np;
    M -> jdim[i] = np;
    M -> kdim[i] = nz;
  }

  M -> xgrid = new float* [nel];
  M -> ygrid = new float* [nel];
  M -> zgrid = new float* [nel];

  ntot = np * np * nz * nel;

  xg = M -> xgrid[0] = new float [ntot];
  yg = M -> ygrid[0] = new float [ntot];
  zg = M -> zgrid[0] = new float [ntot];

  ntot = np * np * nz;

  for (i = 1; i < nel; i++) {
    M -> xgrid[i] = M -> xgrid[0] + i * ntot;
    M -> ygrid[i] = M -> ygrid[0] + i * ntot;
    M -> zgrid[i] = M -> zgrid[0] + i * ntot;
  }

  np2  = np * np;
  ntot = np2;

  for (i = 0; i < nel; i++)
    for (j = 0; j < np2; j++) {
      mfile >> x >> y;
      M -> xgrid[i][j] = x;
      M -> ygrid[i][j] = y;
      for (k = 1; k < nz; k++) {
	M -> xgrid[i][j + k * ntot] = x;
	M -> ygrid[i][j + k * ntot] = y;
      }
    }

  ntot = np * np;

  for (k = 0; k < nz; k++) {
    mfile >> z;
    for (i = 0; i < nel; i++)
      for (j = 0; j < np2; j++)
	M -> zgrid[i][j + k * ntot] = z;
  }

  ntot = np * np * nz * nel;

  // -- Convert to cylindrical coords if required.

  if (State.cylind)
    for (i = 0; i < ntot; i++) {
      r     = zg[i];
      cz    = cos (r);
      sz    = sin (r);
      r     = yg[i];
      yg[i] = r * cz;
      zg[i] = r * sz;
    }

  // -- Set maxima & minima.

  xmin =  1.0e35;
  xmax = -xmin;
  ymin =  xmin;
  ymax =  xmax;
  zmin =  xmin;
  zmax =  xmax;

  for (i = 0; i < ntot; i++) {
    x = xg[i];
    y = yg[i];
    z = zg[i];
    xmin = MIN (xmin, x);
    xmax = MAX (xmax, x);
    ymin = MIN (ymin, y);
    ymax = MAX (ymax, y);
    zmin = MIN (zmin, z);
    zmax = MAX (zmax, z);
  }

  xavg = 0.5 * (xmin + xmax);
  yavg = 0.5 * (ymin + ymax);
  zavg = 0.5 * (zmin + zmax);

  for (i = 0; i < ntot; i++) {
    xg[i] -= xavg;
    yg[i] -= yavg;
    zg[i] -= zavg;
  }

  State.xmin = xmin-xavg;
  State.xmax = xmax-xavg;
  State.ymax = ymin-yavg;
  State.ymax = ymax-yavg;
  State.zmin = zmin-zavg;
  State.zmax = zmax-zavg;
  State.xavg = xavg;
  State.yavg = yavg;
  State.zavg = zavg;

  State.length = MAX (hypot (xmax, ymax), hypot (xmax, zmax));
  State.length = MAX (State.length,       hypot (ymax, zmax));
  State.length *= 2.0;

  mfile.close();

  return M;
}


Data* setFields (const char* fname)
// ---------------------------------------------------------------------------
// Data is just a placeholder to keep a record of the number of fields,
// their names and offsets in a SEM (binary) data file.
// ---------------------------------------------------------------------------
{
  char  routine[] = "setFields";
  const int NP  = Mesh -> idim[0];
  const int NZ  = Mesh -> kdim[0] - 1;
  const int NEL = Mesh -> nel;
  int       i, nr, ns, nz, nel, ntot;
  char      buf[StrMax], err[StrMax];
  double*   tmp;
  Data*     D = new Data;

  D -> file.open (fname, ios::in);
  if (!D -> file) message (routine, "couldn't open field file", ERROR);

  D -> file.getline(buf, StrMax);
  
  if (!strstr (buf, "Session")) {
    sprintf (err, "not a valid field file; first line is: %s", buf);
    message (routine, err, ERROR);
  }

  D -> file.getline(buf, StrMax).getline(buf, StrMax);
  istringstream (buf) >> nr >> ns >> nz >> nel;

  if (nr  != ns ) message (routine, "element NR != NS",      ERROR);
  if (nr  != NP ) message (routine, "element size mismatch", ERROR);
  if (nz  != NZ ) message (routine, "no. z planes mismatch", ERROR);
  if (nel != NEL) message (routine, "no. elements mismatch", ERROR);
  
  tmp          = new double [ntot = nr * ns * nz * nel];
  D -> elmt    = new float* [nel];
  D -> elmt[0] = new float  [nr * ns * (nz + 1) * nel];
  for (i = 1; i < nel; i++)
    D -> elmt[i] = D -> elmt[0] + i * nr * ns * (nz + 1);

  D -> file.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);
  D -> file.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);

  i = 0;
  while (isalpha (buf[i])) { D -> name[i] = buf[i]; i++; }
  D -> name[i] = '\0';
  D -> nfields = strlen (D -> name);

  if (D -> nfields > FldMax) {
    sprintf (err, "no. of fields (%1d) exceeds maximum (%1d)",
	     D -> nfields, FldMax);
    message (routine, err, ERROR);
  }

  D -> file.getline(buf, StrMax);

  format (err);

  if (!strstr (buf, "binary"))
    message (routine, "input field file not in binary format", ERROR);
  
  if (!strstr (buf, "endian"))
    message (routine, "input field file in unknown binary format", WARNING);
  else {
    _swab = ((strstr (buf, "big") && strstr (err, "little")) ||
	     (strstr (err, "big") && strstr (buf, "little")) );
  }

  cout << "-- Mesh  : "
       << nel << " elements, " 
       << nr << " X " << ns << " X " << nz << endl;
  cout << "-- Fields: ";

  for (i = 0; i < D -> nfields; i++) {
    cout << D -> name[i];
    D -> fldPosn[i] = D -> file.tellg();
    D -> file.read ((char*) tmp, ntot * sizeof (double));
  }
  if (D -> file)
    cout << ", OK" << endl;
  else {
    cout << endl;
    message (routine, "problem reading binary data", ERROR);
  }

  // -- Set default field to the first available.

  delete [] tmp;

  loadData (D, D -> name[0]);

  return D;
}


int loadData (Data* F   ,
	      char  name)
// ---------------------------------------------------------------------------
// Load the temporary data array with named field, reading from file.
// The data structure used for interpolation (element-ordered) differs
// from that in the file (plane-ordered).  Thus data are first read into
// a temporary array, then reordered.  Byte swapping may also be required.
//
// There is a slight complication in that the domain is padded in the z 
// direction (allowing for periodicity) so that the last plane of data
// is a copy of the first plane.
// ---------------------------------------------------------------------------
{
  char routine[] = "loadData", err[StrMax];
  int             found = 0;
  const int       nzp  = Mesh -> kdim[0];
  const int       nz   = nzp - 1;
  const int       nel  = Mesh -> nel;
  const int       np2  = sqr (Mesh -> idim[0]);
  const int       ntot = np2 * nz  * nel;
  const int       ntop = np2 * nzp * nel;
  const int       skip = np2 * nel;
  register int    i, j, k, Linv1, Linv2;
  register float  *E, dmin =  1.0e35, dmax = -dmin, datum;
  register double *tmp;

  for (i = 0; !found && i < F -> nfields; i++)
    if (name == F -> name[i]) found = 1;

  if (!found) {
    sprintf (err, "field '%c' not found, no field set", name);
    message (routine, err, WARNING);
    F -> file.clear ();
    ((ifstream&)(F -> file)).seekg (0);
    return 0;
  }

  tmp = new double [ntop];

  // -- Reposition file and extract chosen field into temporary store.

  F -> file.clear ();
  F -> file.seekg (F -> fldPosn[i - 1]);
  F -> file.read  ((char*) tmp, ntot * sizeof (double));

  if (_swab) dbrev (ntot, tmp, tmp);
  
  // -- Copy last ((nz + 1)th) plane from first one.

  for (i = 0; i < skip; i++) tmp[ntot + i] = tmp[i];

  // -- Re-order data: [nel][nz][ns][nr] <-- [nz][nel][ns][nr].

  for (k = 0; k < nel; k++) {
    E = F -> elmt[k];
    for (j = 0; j < nzp; j++) {
      Linv1 = j * np2;
      Linv2 = j * skip + k * np2;
      for (i = 0; i < np2; i++) {
	datum = E[Linv1 + i] = tmp[Linv2 + i];
	dmin  = MIN (dmin, datum);
	dmax  = MAX (dmax, datum);
      }
    }
  }

  F -> current = name;

  cout << "Field '"   << name
       << "': min = " << dmin
       << ", max = "  << dmax << endl;

  delete [] tmp;

  return 1;
}


static int iformat ()
/* ------------------------------------------------------------------------- *
 * Return 1 if machine floating-point format is IEEE little-endian,
 * 0 if IEEE big-endian, -1 for unrecognized format.
 * ------------------------------------------------------------------------- */
{
  union { float  f; int i;    unsigned char c[4]; } v;
  union { double d; int i[2]; unsigned char c[8]; } u;
  int reverse = (-1);
  u.d = 3;
  v.i = 3;
  if      (u.c[0] == 64 && u.c[1] == 8 && v.c[3] == 3) reverse = 0;
  else if (u.c[7] == 64 && u.c[6] == 8 && v.c[0] == 3) reverse = 1;

  return (reverse);
}


static void format (char* s)
/* ------------------------------------------------------------------------- *
 * Fill s with a string describing machine's floating-point storage format.
 * ------------------------------------------------------------------------- */
{
  switch (iformat ()) {
  case -1:          sprintf (s, "unknown");            break;
  case  1:          sprintf (s, "IEEE little-endian"); break;
  case  0: default: sprintf (s, "IEEE big-endian");    break;
  }
}


static void dbrev (const int     n,
		   const double* x,
		         double* y)
// ---------------------------------------------------------------------------
// Reverse order of bytes in double vector x.
// ---------------------------------------------------------------------------
{
  register char *cx, *cy, d;
  register int  i;

  for (i = 0; i < n; i++) {
    cx = (char*) (x + i);
    cy = (char*) (y + i);
    d  = cx[0]; cy[0] = cx[7]; cy[7] = d;
    d  = cx[1]; cy[1] = cx[6]; cy[6] = d;
    d  = cx[2]; cy[2] = cx[5]; cy[5] = d;
    d  = cx[3]; cy[3] = cx[4]; cy[4] = d;
  }
}


void loadPnts (const char* pfile)
// ---------------------------------------------------------------------------
// If a file of point data can be found, open it and get point data from it.
// Intially put the data in a stack, then a global vector (Point).
//
// Points format:
// id time x y z
// 
// Initialise value entry to zero for each point.
// Convert locations from cylindrical to Cartesian if required.
// ---------------------------------------------------------------------------
{
  if (!pfile) return;

  using namespace std;

  stack< Pnt*, vector<Pnt*> > pstack;

  Pnt*     datum;
  int      i, N, id;
  double   r, cz, sz;
  float    time, x, y, z, val=0;
  ifstream file (pfile);

  if (!file) message ("loadPnts", "cannot open file", ERROR);

  while (file >> id >> time >> x >> y >> z) {

    if (State.cylind) 
      { r = y; cz = cos (z); sz = sin (z); y = r * cz; z  = r * sz; }

    x -= State.xavg;
    y -= State.yavg;
    z -= State.zavg;

    datum          = new Pnt;
    datum -> id    = id;
    datum -> time  = time;
    datum -> value = val;
    datum -> x     = x;
    datum -> y     = y;
    datum -> z     = z;

    pstack.push (datum);
  }

  Point.resize (N = pstack.size());
  
  for (i = N; i; i--) {
    Point[i-1] = pstack.top();
    pstack.pop();
  }
}

