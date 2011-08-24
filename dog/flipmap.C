///////////////////////////////////////////////////////////////////////////////
// flipmap.C: generate a list of index pairs for a symmetry-related transform.
//
// Copyright (c) 2002 <--> $Date$, Hugh Blackburn
//
// USAGE
// flipmap [options] [file]
// options:
//   -h             ... print this message
//   -x || -y || -d ... reflection symmetry to be used to generate map
//   -t <num>       ... set positional tolerance to num, default 6e-7 = EPSSP
//
// FILES
// Input file is a semtex/prism ASCII mesh file, generated e.g. by meshpr.
// Output is ASCII, with three header lines followed by a list of index pairs:
//   8 8 1 108 NR NS NZ NEL  # -- matches input file header
//   x                       # -- reflection symmetry generator
//   6768                    # -- NFLIP, no. of pairs to follow < NR*NS*NEL
//   11   4001
//   12   4000
//   13   3999
//   ... etc. (NFLIP lines)
//
// The mapping between the mesh points may not be unique; we take the
// first mapping found for each case (reflect positive->negative,
// *and* negative->positive).  This way a single gather will do the
// reflection, and leave no holes.  The nmap does not correspond to a
// global numbering scheme, it is simply the number of points in the
// mesh.  The numbers in the two lists are the correspondences between
// the indices of reflected points in a flat element-by-element
// ordering of the mesh points.
//
// To carry out the exchange of data, say all the NTOT=NR*NR*NEL data
// are first copied from array org to array tmp.  Exchange by
// gather-scatter:
//
//   for (i = 0; i < NFLIP; i++) org[neg[i]] = tmp[pos[i]];
//
// where pos and neg are the indices in the two lists.  To complete a
// symmetry operation may also require a negation (say if the data are
// one component of a vector).  This has to be determined
// independently.
//
// NB: -x means the reflection is in the x direction, i.e. the reflection
// occurs about the y axis!
//
// Symmetry generator 'd' means a double reflection, i.e. in both x and y.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>

static char prog[] = "flipmap";

const real_t EPS = EPSSP;	// -- Default positional tolerance.

static void  getargs  (int, char**, char&, real_t&, ifstream&);
static int_t header   (ifstream&);
static void  loadmesh (ifstream&, const int_t, const char,
		       vector<real_t>&, vector<real_t>&);
static void  findmap  (const char, const int_t, const real_t tol,
		       const vector<real_t>&, const vector<real_t>&,
		       vector<int_t>&, vector<int_t>&);
static void  printup  (const char, const int_t,
		       const vector<int_t>&, const vector<int_t>&);


int main (int    argc,
          char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  ifstream       file;
  char           generator;
  vector<real_t> x, y;
  vector<int_t>  orig, flip;
  int_t          npts;
  real_t         tol = EPS;

  getargs (argc, argv, generator, tol, file);

  npts = header (file);

  loadmesh (file, npts, generator, x, y);
  
  findmap  (generator, npts, tol, x, y, orig, flip);
  printup  (generator, npts, orig, flip);

  return EXIT_SUCCESS;
}


static void getargs (int       argc,
		     char**    argv,
		     char&     gen ,
		     real_t&   tol ,
		     ifstream& file)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: flipmap [-h] -x || -y || -d [-t <num>] [meshfile]\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 't':
      if (*++argv[0]) tol = atof (  *argv);
      else { --argc;  tol = atof (*++argv); }
      break;
    case 'x': gen = 'x'; break;
    case 'y': gen = 'y'; break;
    case 'd': gen = 'd'; break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (!(gen == 'x' || gen == 'y' || gen == 'd')) {
    cerr << prog << ": must specify -x, -y or -d" << endl;
    exit (EXIT_FAILURE);
  }

  if (argc == 1)
    file.open (argv[0], ios::in);
  else {
    cerr << usage;
    exit (EXIT_FAILURE);
  }

  if (!file) {
    cerr << prog << ": unable to open file" << endl;
    exit (EXIT_FAILURE);
  }
}


static int_t header (ifstream& file)
// ---------------------------------------------------------------------------
// If file header indicates this is a mesh file, output header and
// return number of points. Else die.
// ---------------------------------------------------------------------------
{
  int_t nr, ns, nz, nel;
  char  buf[StrMax];
  
  file >> nr >> ns >> nz >> nel;
  file.getline (buf, StrMax);

  if (!strstr (buf, "NR NS NZ NEL")) {
    cerr << prog
	 << ": mesh header line should include NR NS NZ NEL: " 
	 << buf << endl;
    exit (EXIT_FAILURE);
  }
  
  cout << nr << " " << ns << " " << 1 << " " << nel << " NR NS NZ NEL" << endl;

  return nr * ns * nel;
}


static void loadmesh (ifstream&       file,
		      const int_t     npts,
		      const char      gen ,
		      vector<real_t>& x   ,
		      vector<real_t>& y   )
// ---------------------------------------------------------------------------
// Load the vectors x and y.
// ---------------------------------------------------------------------------
{
  int_t i;

  x.resize (npts);
  y.resize (npts);

  for (i = 0; file && i < npts; i++) file >> x[i] >> y[i];
}


static void findmap (const char            gen ,
		     const int_t           npts,
		     const real_t          tol ,
		     const vector<real_t>& x   ,
		     const vector<real_t>& y   ,
		     vector<int_t>&        orig,
		     vector<int_t>&        flip)
// ---------------------------------------------------------------------------
// This is where the mapping gets constructed. Order npts*npts operation.
// Take the first available mapping index for each point.
// ---------------------------------------------------------------------------
{
  int_t i, j, k = 0;
  bool  found;

  orig.resize (npts);
  flip.resize (npts);

  if (gen == 'x')
    for (i = 0; i < npts; i++) {
      for (found = false, j = 0; !found && j < npts; j++)
	if (fabs(x[i] + x[j]) < tol && fabs(y[i] - y[j]) < tol) {
	  orig[k]   = i;
	  flip[k++] = j;
	  found     = true;
	}
      if (!found) cerr << "Warning: mesh point "
		       << x[i] << ",\t" << y[i] << "\tnot mirrored" << endl;
    }
  else if (gen == 'y')
    for (i = 0; i < npts; i++) {
      for (found = false, j = 0; !found && j < npts; j++) {
	if (fabs(y[i] + y[j]) < tol && fabs(x[i] - x[j]) < tol) {
	  orig[k]   = i;
	  flip[k++] = j;
	  found     = true;
	}
      }
      if (!found) cerr << "Warning: mesh point "
		       << x[i] << ",\t" << y[i] << "\tnot mirrored" << endl;
    }
  else				// -- 'd'.
    for (i = 0; i < npts; i++) {
      for (found = false, j = 0; !found && j < npts; j++)
	if (fabs(x[i] + x[j]) < tol && fabs(y[i] + y[j]) < tol) {
	  orig[k]   = i;
	  flip[k++] = j;
	  found     = true;
	}
      if (!found) cerr << "Warning: mesh point "
		       << x[i] << ",\t" << y[i] << "\tnot mirrored" << endl;
    }

  if (k != npts) {
    cerr << prog
	 << ": number of maps found, " << k
	 << ", not equal to number needed, " << npts << endl;
    exit (EXIT_FAILURE);
  }
}


static void printup (const char           gen ,
		     const int_t          npts,
		     const vector<int_t>& orig,
		     const vector<int_t>& flip)
// ---------------------------------------------------------------------------
// Print the map on standard output.
// ---------------------------------------------------------------------------
{
  int_t i;

  cout << gen  << endl;
  cout << npts << endl;
  
  for (i = 0; i < npts; i++)
    cout << orig[i] << '\t' << flip[i] << endl;
}
