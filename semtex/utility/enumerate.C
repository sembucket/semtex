///////////////////////////////////////////////////////////////////////////////
// enumerate.C:  utility to generate mesh numbering from mesh description file.
//
// Usage: enumerate [options] file
//   options:
//   -h       ... display this message
//   -v       ... set verbose output
//   -n N     ... override element order to be N
//   -O [0-3] ... set level of bandwidth optimization
//
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <stdlib.h>
#include <limits.h>
#include <iomanip.h>
#include <femdef.h>
#include <Femlib.h>
#include <Utility.h>
#include <Mesh.h>
#include <Veclib.h>
#include <List.h>

class Nsys {
friend void printup (vector<char>&, vector<Nsys*>&, const int);
public:
  Nsys (char, vector<int>&, vector<int>&, const int);

  int  match    (vector<int>&);
  void addField (char);

private:
  int          nel;
  int          nglobal;
  int          nbndry;
  int          nsolve;
  int          nbandw;
  int          optlev;
  vector<char> fields;
  vector<int>  bndmap;
  vector<int>  bndmsk;

  int  sortGid         (int*, int*);
  void renumber        (const int);
  int  buildAdjncy     (List<int>*) const;
  void fillAdjncy      (List<int>*, int*, int*, const int) const;
  void connectivSC     (List<int>*, const int*, const int*, const int) const;
  int  globalBandwidth () const;
  int  bandwidthSC     (const int*, const int*, const int) const;
};

static char prog[] =  "enumerate";
static void getargs   (int, char**, char*&, int&, int&, int&);
static void getfields (FEML&, vector<char>&);
static void printup   (vector<char>&, vector<Nsys*>&, const int);

static const int FldMax = 16;


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Determine, from BCs section of FEML file, list of fields for which
// numbering schemes are to be constructed.
//
// Generate a BC mask and initial numbering scheme for first field, using
// Mesh class routines.  Optimize numbering scheme according to selected level.
//
// For each succeeding field, first generate a BC mask and, if it matches
// a mask previously generated, add the field's name to the previous field's
// name vector but take no further action.  Otherwise, generate and optimize
// a new numbering system.
//
// Print up the masks and numbering schemes on cout.
// ---------------------------------------------------------------------------
{
  char* session = 0;
  int   verb    = 0,
        np      = 0,
        opt     = 0;

  getargs (argc, argv, session, verb, np, opt);

  FEML feml (session);

                   Femlib::value ("OPTIMIZE", opt);
  if (verb)        Femlib::value ("VERBOSE", verb);
  if   (np)        Femlib::value ("N_POLY", np);
  else  np = (int) Femlib::value ("N_POLY");

  vector<char> field;

  getfields (feml, field);

  Mesh          M (feml);
  vector<Nsys*> S (field.getSize());

  int         i, j, k = 0, found;
  const  int  NEL  = M.nEl();
  const  int  NTOT = 4 * NEL * (np - 1);
  vector<int> btog (NTOT);
  vector<int> mask (NTOT);

  M.buildMask (np, field[0], mask());
  M.buildMap  (np, btog());

  S[k++] = new Nsys (field[0], btog, mask, NEL);

  for (i = 1; i < field.getSize(); i++) {
    M.buildMask (np, field[i], mask());
    found = 0;
    for (j = 0; !found && j < k; j++)
      if (found = S[j] -> match (mask)) 
	S[j] -> addField (field[i]);
    if (!found) {
      M.buildMap (np, btog());
      S[k++] = new Nsys (field[i], btog, mask, NEL);
    }
  }

  printup (field, S, k);

  return EXIT_SUCCESS;
}


static void getargs (int    argc   , 
		     char** argv   ,
		     char*& session,
		     int&   verb   ,
		     int&   np     ,
		     int&   opt    )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: enumerate [options] session\n"
                 "options:\n"
                 "  -h       ... display this message\n"
                 "  -v       ... set verbose output\n"
		 "  -n N     ... override number of element knots to be N\n"
		 "  -O [0-3] ... bandwidth optimization level [Default: 1]\n";
  char err[StrMax];
  char c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      for (verb = 1; *++argv[0] == 'v'; verb++);
      break;
    case 'n':
      if (*++argv[0])
	np = atoi (*argv);
      else {
	--argc;
	np = atoi (*++argv);
      }
      break;
    case 'O':
      if (*++argv[0])
	opt = atoi (*argv);
      else {
	--argc;
	opt = atoi (*++argv);
      }
      break;
    default:
      sprintf (err, "%s: illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if   (argc == 1) session = *argv;
  else             message (prog, "must provide session file", ERROR);
}


static void getfields (FEML&         feml ,
		       vector<char>& field)
// ---------------------------------------------------------------------------
// The default action is to set up the list of fields according to the
// names found in the 'BCS' section of the FEML file.
//
// If no 'BCS' section is located (valid in the case of all periodic
// boundaries), the default list of variables is set to either "uvp"
// or "uvwp" depending on the value of N_Z: for a 2D/Fourier scheme,
// N_Z effectively decides between a 3 (N_Z == 1) or 4-variable 
// (N_Z > 1) solution.
// ---------------------------------------------------------------------------
{
  if (feml.seek ("BCS")) {
    int  j, id, nbcs;
    char fieldc, groupc, nextc, tag[StrMax];

    feml.attribute ("BCS", "NUMBER");
    
    while ((nextc = feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');

    feml.stream() >> id >> groupc >> nbcs;
    
    field.setSize (nbcs);

    for (j = 0; j < nbcs; j++) {
      feml.stream() >> tag >> fieldc;
      field[j] = fieldc;
      feml.stream().ignore (StrMax, '\n');
    }

  } else {
    if ((int) Femlib::value ("N_Z") > 1) {
      field.setSize (4);
      field[0] = 'u';
      field[1] = 'v';
      field[2] = 'w';
      field[3] = 'p';
    } else {
      field.setSize (3);
      field[0] = 'u';
      field[1] = 'v';
      field[2] = 'p';
    }
  }
};


static void printup (vector<char>&  F   ,
		     vector<Nsys*>& S   ,
		     const int      nSys)
// ---------------------------------------------------------------------------
// print up summary info followed by map & mask for eack system.
// ---------------------------------------------------------------------------
{
  register int i, j, k, side, soff;
  const    int nedge = S[0] -> nbndry / (4 * S[0] -> nel);
  
  cout << "# FIELDS         : " << F() << endl;

  cout << "# MATCHING       :";
  for (j = 0; j < nSys; j++) {
    cout << setw (12);
    cout << S[j] -> fields();
  }
  cout << endl;

  cout << "# NEL            :";
  for (j = 0; j < nSys; j++) {
    cout << setw (12);
    cout << S[j] -> nel;
  }
  cout << endl;

  cout << "# NBOUNDARY      :";
  for (j = 0; j < nSys; j++) {
    cout << setw (12);
    cout << S[j] -> nbndry;
  }
  cout << endl;

  cout << "# NGLOBAL        :";
  for (j = 0; j < nSys; j++) {
    cout << setw (12);
    cout << S[j] -> nglobal;
  }
  cout << endl;

  cout << "# NSOLVE         :";
  for (j = 0; j < nSys; j++) {
    cout << setw (12);
    cout << S[j] -> nsolve;
  }
  cout << endl;

  cout << "# OPTIMIZATION   :";
  for (j = 0; j < nSys; j++) {
    cout << setw (12);
    cout << S[j] -> optlev;
  }
  cout << endl;

  cout << "# BANDWIDTH      :";
  for (j = 0; j < nSys; j++) {
    cout << setw (12);
    cout << S[j] -> nbandw;
  }
  cout << endl;

  cout << "#" << endl;
  cout << "# elmt  side offst";
  for (j = 0; j < nSys; j++) cout << "  bmap  mask";
  cout << endl;

  for (i = 0, k = 1; k <= S[0] -> nel; k++)
    for (side = 1; side <= 4; side++)
      for (soff = 0; soff < nedge; soff++, i++) {
	cout << setw (6) << k << setw (6) << side << setw (6) << soff;
	for (j = 0; j < nSys; j++)
	  cout 
	    << setw (6) << S[j] -> bndmap (i)
	      << setw (6) << S[j] -> bndmsk (i);
	cout << endl;
      }
}


Nsys::Nsys (char         name,
	    vector<int>& map ,
	    vector<int>& mask,
	    const int    nEl )
// ---------------------------------------------------------------------------
// Constructor also carries out bandwidth optimization task.
// ---------------------------------------------------------------------------
{
  fields.setSize (FldMax);
  memset (fields(), '\0', FldMax);
  fields[0] = name;
  nel       = nEl;
  bndmap    = map;
  bndmsk    = mask;
  nbndry    = bndmap.getSize();
  optlev    = (int) Femlib::value ("OPTIMIZE");
  nglobal   = bndmap (Veclib::imax (nbndry, bndmap(), 1)) + 1;
  nsolve    = sortGid (bndmap(), bndmsk());

  renumber (optlev);
  nbandw = globalBandwidth ();
}


int Nsys::match (vector<int>& test)
// ---------------------------------------------------------------------------
// Return true if test and bndmsk match.
// ---------------------------------------------------------------------------
{
  if (test.getSize() != bndmsk.getSize()) return 0;

  return Veclib::same (bndmsk.getSize(), bndmsk(), 1, test(), 1);
}


void Nsys::addField (char name)
// ---------------------------------------------------------------------------
// Add a new field which matches the present one.
// ---------------------------------------------------------------------------
{
  int k = 0;

  while (fields (k)) k++;

  if   (k == FldMax) message (prog, "too many fields", ERROR);
  else               fields  (k) = name;
}


static int cmp1 (const void *a,
		 const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare first element (global node number) of two arrays.
// ---------------------------------------------------------------------------
{ return ((int *)a)[0] - ((int *)b)[0]; }


static int cmp2 (const void *a,
		 const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare second element (solve mask) of two arrays.
// ---------------------------------------------------------------------------
{ return ((int *)a)[1] - ((int *)b)[1]; }



int Nsys::sortGid (int* bmap,
		   int* bmsk)
// ---------------------------------------------------------------------------
// Global node numbers get sorted to place essential-BC nodes last:
// this simplifies the later partition of global matrices.
//
// The non-essential type node numbers can be further sorted to
// optimize global matrix bandwidths, but this is not done here.
//
// A globally-numbered table (reOrder) is constructed, each entry of which is
// two ints: the global node number and the essential BC mask value (0/1).
// This is then partitioned into "unknown" and "known" node numbers, with the
// "known" numbers last in the table.  Each partition is then sorted into
// ascending node number order, and the information is used to create
// a new element boundary-to-global numbering scheme in btog.
//
// Return number of global nodes at which solution is not set by essential BCs.
// ---------------------------------------------------------------------------
{
  vector<int> work (nbndry + 3 * nglobal);
  int         *bsave, *tmp, *reOrder;
  int         unknowns;

  bsave   = work();
  tmp     = bsave + nbndry;
  reOrder = tmp  + nglobal;

  Veclib::copy  (nbndry,  bmap, 1, bsave, 1);
  Veclib::scatr (nbndry , bmsk, bsave, tmp);
  Veclib::ramp  (nglobal, 0,    1, reOrder,     2);
  Veclib::copy  (nglobal, tmp,  1, reOrder + 1, 2);
  
  unknowns = nglobal - Veclib::count (nglobal, tmp, 1);

  if (unknowns < nglobal) {

    // -- Partition into "unknown" nodes & "essential BC" nodes.

    qsort (reOrder,                nglobal,            2 * sizeof (int), cmp2);

    // -- Sort each partition into ascending node number order.

    qsort (reOrder,                unknowns,           2 * sizeof (int), cmp1);
    qsort (reOrder + 2 * unknowns, nglobal - unknowns, 2 * sizeof (int), cmp1);

    // -- Reload new gids.

    Veclib::copy  (nglobal, reOrder, 2,  tmp, 1);
    Veclib::ramp  (nglobal, 0,    1, reOrder, 1);
    Veclib::scatr (nglobal, reOrder, tmp, reOrder + nglobal);
    Veclib::gathr (nbndry , reOrder + nglobal, bsave, bmap);
  }

  return unknowns;
}


void Nsys::renumber (const int optl)
// ---------------------------------------------------------------------------
// From the initial ordering specified in bndmap, use RCM to generate a
// reduced-bandwidth numbering scheme.  Reload into bndmap.
//
// Different optimization levels are allowed:
//
// 0: Do nothing (no renumbering).
// 1: Use FNROOT (trial root = 1) to find a pseudo-peripheral root node,
//    pass result to RCM for Reverse Cuthill McKee reordering.  Default level.
// 2: Use FNROOT to generate pseudo-peripheral nodes, but with trial roots
//    in steps of 10, up to n_solve.  Choose root to minimize global bandwidth.
// 3: Do not use FNROOT.  Try all unknown node numbers as trial roots for RCM.
//    Choose root to minimize global bandwidth.
//
// Reference:
//    A. George and J. W-H. Liu
//    Computer Solution of Large Sparse Positive Definite Systems
//    Prentice-Hall (1981)
// ---------------------------------------------------------------------------
{
  if (!optl) return;

  const int verb = (int) Femlib::value ("VERBOSE");

  if (verb)
    cout << "-- Bandwidth optimization (" << optl
      << "), Field '" << fields() << "'";

  register int i;
  int          root, nlvl;

  // -- Build node adjacency tables.
  
  List<int>* adjncyList = new List<int> [nsolve];
  const int  tabSize    = buildAdjncy (adjncyList);

  // -- Allocate memory.

  vector<int> work(tabSize + 1 + 4 * nsolve + 1 + nglobal + nbndry);
  int         *adjncy, *xadj, *perm, *mask, *xls, *invperm, *bsave;

  adjncy  = work();
  xadj    = adjncy  + tabSize + 1;
  perm    = xadj    + nsolve  + 1;
  mask    = perm    + nsolve;
  xls     = mask    + nsolve;
  invperm = xls     + nsolve;
  bsave   = invperm + nglobal;
  
  Veclib::copy (nbndry, bndmap(), 1, bsave, 1);
  for (i = nsolve; i < nglobal; i++) invperm[i] = i;

  fillAdjncy (adjncyList, adjncy, xadj, tabSize);
  delete   [] adjncyList;

  switch (optl) {
  case 1: {
    root = 1;
    Veclib::fill   (nsolve, 1, mask, 1);
    Femlib::fnroot (root, xadj, adjncy, mask, nlvl, xls, perm);
    Femlib::rcm    (root, xadj, adjncy, mask, perm, nlvl, xls);
    break;
  }
  case 2: {
    int  rtest, BWtest, BWmin = INT_MAX, best;

    if (verb) cout << ":";

    for (root = 1; root <= nsolve; root += 10) {
      rtest = root;
      Veclib::fill   (nsolve, 1, mask, 1);
      Femlib::fnroot (rtest, xadj, adjncy, mask, nlvl, xls, perm);
      Femlib::rcm    (rtest, xadj, adjncy, mask, perm, nlvl, xls);

      Veclib::sadd (nsolve, -1, perm, 1, perm, 1);
      for (i = 0; i < nsolve; i++) invperm[perm[i]] = i;
      Veclib::gathr (nbndry, invperm, bsave, bndmap());

      BWtest = globalBandwidth ();
      if (BWtest < BWmin) {
	BWmin = BWtest;
	best  = rtest;
	if (verb) cout << " " << BWmin;
      }
    }

    Veclib::fill (nsolve, 1, mask, 1);
    Femlib::rcm  (best, xadj, adjncy, mask, perm, nlvl, xls );

    break;
  }
  case 3: {
    int  BWtest, BWmin = INT_MAX, best;

    if (verb) cout << ":";

    for (root = 1; root <= nsolve; root++) {
      Veclib::fill (nsolve, 1, mask, 1);
      Femlib::rcm  (root, xadj, adjncy, mask, perm, nlvl, xls);

      Veclib::sadd (nsolve, -1, perm, 1, perm, 1);
      for (i = 0; i < nsolve; i++) invperm[perm[i]] = i;
      Veclib::gathr (nbndry, invperm, bsave, bndmap());

      BWtest = globalBandwidth ();
      if (BWtest < BWmin) {
	BWmin = BWtest;
	best  = root;
	if (verb) cout << " " << BWmin;
      }
    }

    Veclib::fill (nsolve, 1, mask, 1);
    Femlib::rcm  (best, xadj, adjncy, mask, perm, nlvl, xls );

    break;
  }
  default:
    break;
  }

  Veclib::sadd (nsolve, -1, perm, 1, perm, 1);
  for (i = 0; i < nsolve; i++) invperm[perm[i]] = i;
  Veclib::gathr (nbndry, invperm, bsave, bndmap());
  
  if (verb) cout << endl;
}


int Nsys::buildAdjncy (List<int>* adjncyList) const
// ---------------------------------------------------------------------------
// Traverse elements and build up a vector of linked lists that
// describe the global nodes adjacent to each global node.
//
// Return the total amount of storage required when information is packed
// into an integer vector, as required by genrcm.
// ---------------------------------------------------------------------------
{
  register int k, ntab;
  const int    next = nbndry / nel;

  for (k = 0, ntab = 0; k < nel; k++) {
    connectivSC (adjncyList, bndmap() + ntab, bndmsk() + ntab, next);
    ntab += next;
  }

  for (k = 0, ntab = 0; k < nsolve; k++)
    ntab += adjncyList[k].length();

  return ntab;
}


void Nsys::fillAdjncy (List<int>* adjncyList,
		       int*       adjncy    ,
		       int*       xadj      ,
		       const int  tabSize   ) const
// ---------------------------------------------------------------------------
// Load the information contained in adjncyList into the two vectors
// adjncy & xadj required by genrcm.
// ---------------------------------------------------------------------------
{
  char         routine[] = "Nsys::fillAdjncy";
  register int i, k;

  for (i = 0, k = 1; i < nsolve; i++) {
    xadj[i] = k;
    for (ListIterator<int> p(adjncyList[i]); p.more(); p.next(), k++)
      adjncy[k - 1] = p.current() + 1;
  }
  
  if (k != tabSize + 1)
    message (routine, "after traversing list, k != tabSize + 1", ERROR);

  adjncy[tabSize] = 0;
  xadj  [ nsolve] = k;
}


void Nsys::connectivSC (List<int>* adjList,
			const int* bmap   ,
			const int* mask   ,
			const int  next   ) const
// ---------------------------------------------------------------------------
// AdjList is an array of linked lists, each of which describes the global
// nodes that have connectivity with the the current node, i.e. which make
// a contribution to the weighted-residual integral for this node.
// This routine fills in the contribution from the current element.
//
// For general finite elements, all nodes of an element are interconnected,
// while for statically-condensed elements, only the boundary nodes are
// considered (since internal nodes are not global).
//
// Essential-BC nodes are ignored, since we're only interested in mimimizing
// bandwidths of global matrices.
// ---------------------------------------------------------------------------
{
  register int i, j, found, gidCurr, gidMate;
  
  for (i = 0; i < next; i++) {
    if (! mask[i]) {
      gidCurr = bmap[i];
      
      for (j = 0; j < next; j++) {
	if (i != j && ! mask[j]) {
	  ListIterator<int> a (adjList[gidCurr]);

	  for (gidMate = bmap[j], found = 0; !found && a.more(); a.next())
	    found = a.current() == gidMate;
	  if (!found) adjList[gidCurr].add (gidMate);
	}
      }
    }
  }
}


int Nsys::globalBandwidth () const
// --------------------------------------------------------------------------
// Precompute the bandwidth of assembled global matrix (including diagonal).
// --------------------------------------------------------------------------
{
  register int k, noff, nband = 0;
  const int    next = nbndry / nel;

  for (k = 0, noff = 0; k < nel; k++) {
    nband = max (bandwidthSC (bndmap() + noff, bndmsk() + noff, next), nband);
    noff += next;
  }

  ++nband; // -- Diagonal.

  return nband;
}



int Nsys::bandwidthSC (const int* bmap,
		       const int* mask,
		       const int  next) const
// ---------------------------------------------------------------------------
// Find the global equation bandwidth of this element, excluding diagonal.
// ---------------------------------------------------------------------------
{
  register int i;
  register int Min  = INT_MAX;
  register int Max  = INT_MIN;

  for (i = 0; i < next; i++) {
    if (!mask[i]) {
      Min = min (bmap[i], Min);
      Max = max (bmap[i], Max);
    }
  }

  return Max - Min;
}
