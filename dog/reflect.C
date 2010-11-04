///////////////////////////////////////////////////////////////////////////////
// reflect.C: reflect a field dump defined on a half-mesh onto a full
// mesh.  Apply sign change to appropriate velocity component.  Built
// from semtex/utility/data2df_template.C.  See also dog/symmetrize.C
//
// Copyright (c) 2010 <--> $Date$, Hugh Blackburn
//
// USAGE
// -----
// reflect [options] -m mapfile fieldfile
// options:
// -h       ... print this message.
//
// If fieldfile is not present, read from standard input.  Write to
// standard output.  Note that the size of the output file will be
// about twice as large as the input file and the associated fields
// will have twice as many elements.  It is assumed that the mapfile
// has a structure appropriate to the task of the reflection mapping
// on the full domain.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>
#include <data2df.h>

static char prog[] = "reflect";
static void getargs  (int, char**, istream*&, istream*&);
static void loadmap  (Header&, istream&, char&,vector<int_t>&,vector<int_t>&);
static void halfread (istream*, Data2DF*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  int_t            i, j;
  istream          *input = 0, *mapping = 0;
  Header           h;
  vector<Data2DF*> u;
  Data2DF*         tmp;
  char             generator = 'x';
  vector<int_t>    positive, negative;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, mapping, input);

  *input >> h;
  h.nel += h.nel;		// -- Double the number of elements.
  loadmap (h, *mapping, generator, positive, negative);

  cout << h;

  u.resize (h.nFields());
  tmp = new Data2DF (h.nr, h.nz, h.nel, 0);

  for (i = 0; i < h.nFields(); i++) {

    const int_t ntot = u[i] -> _ntot;

    // - Create u[i].

    u[i] = new Data2DF (h.nr, h.nz, h.nel, h.flds[i]);

    // -- Read in u[i] and byte-swap if necessary.

    halfread (input, u[i]);	// -- Input file only has 1/2 the elements.
    if (h.swab()) u[i] -> reverse();

    // -- Copy the refkection of the other half of the data field,
    // -- negating appropriate velcity vector component.

    if (generator == 'y') {
      if (u[i] -> getName() == 'u' ||
	  u[i] -> getName() == 'w' ||
	  u[i] -> getName() == 'p') {
	*tmp = *u[i];
	u[i] -> reflect2D (positive, negative);
	for (j = 0; j < ntot; j++)
	  u[i]->_data[j] += (fabs(u[i]->_data[j] > EPSDP)) ? tmp->_data[j]:0.0;
      } else if (u[i] -> getName() == 'v') {
	*tmp = *u[i];
	u[i] -> reflect2D (positive, negative);
	for (j = 0; j < ntot; j++)
	  u[i]->_data[j] -= (fabs(u[i]->_data[j]) > EPSDP) ? tmp->_data[j]:0.0;
      }
    } else if (generator == 'x') {
      if (u[i] -> getName() == 'v' ||
	  u[i] -> getName() == 'w' ||
	  u[i] -> getName() == 'p') {
	*tmp = *u[i];
	u[i] -> reflect2D (positive, negative);
	for (j = 0; j < ntot; j++)
	  u[i]->_data[j] += (fabs(u[i]->_data[j] > EPSDP)) ? tmp->_data[j]:0.0;
      } else if (u[i] -> getName() == 'u') {
	*tmp = *u[i];
	u[i] -> reflect2D (positive, negative);
	for (j = 0; j < ntot; j++)
	  u[i]->_data[j] -= (fabs(u[i]->_data[j]) > EPSDP) ? tmp->_data[j]:0.0;
      }
    }

    // -- Output u[i];

    cout << *u[i];
  }

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     istream*& mapfl,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.  Based on dog/reflect.C.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: reflect [options] -m mapfile [file]\n"
    "options:\n"
    "-h       ... print this message\n";
  
  if (argc < 3) {
    cerr << usage ;
    exit (EXIT_FAILURE);
  }
    
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'm':
      --argc; ++argv;
      mapfl = new ifstream (*argv);
      if (mapfl -> bad()) message (prog, "unable to open map file", ERROR);
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


static void loadmap (Header&        headr    ,
		     istream&       file     ,
		     char&          generator,
		     vector<int_t>& positive ,
		     vector<int_t>& negative )
// ---------------------------------------------------------------------------
// Load symmetry mapping information from file.
// ---------------------------------------------------------------------------
{
  const int_t np  = headr.nr;
  const int_t nel = headr.nel;
  char        buf[StrMax], err[StrMax];
  int_t       i, NR, NS, NEL, NMAP;
  
  if (!file) {
    sprintf (err, "cannot find map file %s", buf);
    message (prog, err, ERROR);
  }

  file >> NR >> NS >> NEL >> NEL;
  file.ignore (StrMax, '\n');

  if (NR != np || NS != np || NEL != nel)
    message (prog, "map file doesn't conform with session file", ERROR);
  file >> generator;
  if (!(generator == 'x' || generator == 'y'))
    message (prog, "symmetry generator must be either 'x' or 'y'", ERROR);
  
  file >> NMAP;

  positive.resize (NMAP);
  negative.resize (NMAP);

  for (i = 0; i < NMAP; i++) file >> positive[i] >> negative[i];

  if (!file)
    message (prog, "bad (premature end of?) map file", ERROR);

}


static void halfread (istream* file,
		      Data2DF* dat )
// ---------------------------------------------------------------------------
// This is like the input operator for Data2DF but we only read in the
// first half of the data -- since the input data are defined for a
// domain which is only half as big as the output. Leave zero
// everything else (set by constructor).
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t ntot = dat -> _np * dat -> _np * (dat -> _nel >> 1);
  
  for (i = 0; i < dat -> _nz; i++)
    file -> read (reinterpret_cast<char*> (dat -> _plane[i]),
		  ntot * sizeof (real_t));
}
