///////////////////////////////////////////////////////////////////////////////
// outline.C: generate input (initial loop) for qmesh.
//
// Usage: outline [-h] [file]
//
// Input to outline defines a number of lines or arcs, along which
// points are spaced according to rules.
//
// -- Sample input file:
//  0.0  0.0
// straight 10 uniform
// 12.0  0.0
// straight  8 geometric 1.1
// 12.0  5.0
// push
// straight  5 uniform
//  4.0  5.0
// pop
// straight  2 uniform
// 12.0  8.0
// straight 20 uniform
//  0.0  8.0
// straight  2 uniform
//  0.0  8.0
// arc       4 -2.0
//  0.0  4.0
// straight  4 uniform
// -- End.
//
// Input consists of:
// 1. points   (x,y pairs);
// 2. lines    (straight, arc);
// 3. keywords (push, pop);
// 4. options  (uniform, geometric, numeric values)
// and is free-format.
//
// Options act to modify line-generation commands.  The generation commands
// are always immediately followed by the number of new points to be
// generated, then by modifiers which are used to determine placement.
// 
// Line-generation commands are generally bracketed by points which define
// the start and end locations of the line to be generated.
// 
// Exceptions to this occur when "push" or "pop" follow the start
// point.  Push & pop should always occur in matching pairs.  After
// push is issued, generated points are pushed on a duplication stack
// as they are made.  When the matching pop occurs, the duplicate
// points are produced and the last one is used as the starting point
// for the next line. This allows the generation of internal lines of
// nodes.  Push and pop pairs cannot be nested (yet).
//
// On encountering EOF (or any input that can't be interpreted as a Point)
// the loop of Points is closed to the first point read and execution
// terminates after printing output.
// 
// Output consists of a list of the unique points generated, and a list
// of integer tags into this list of points that defines the generated loop.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <qmesh.h>
#include <Stack.h>


static char      prog[] = "outline";
static const int INS_MAX = 4096;

enum key { UNDEFINED, STRAIGHT, ARC, PUSH, POP };

static void getArgs       (int, char**, istream*& file);
static key  parse         (istream&);
static int  generateNodes (istream&, const key&, Node*&, Node*&,
			   const int&, List<Node*>&, List<int>&,
			   int&, const int&, Stack<int>&);
static void setSizes      (List<Node*>&, const List<int>&);
static void printup       (ostream&, const List<Node*>&, const List<int>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  istream*    file;
  char        err[StrMax];

  List<Node*> nodes;
  List<int>   loop;
  Stack<int>  save;

  int         id, nP, pushing = 0, moreInput = 1, id_max = 0;
  Point       P;
  Node        *startN, *homeN, *precN;
  key         Keyword;

  getArgs (argc, argv, file);

  *file >> P;
//startN = new Node (++id_max, P, 1.0, Node::BOUNDARY);
  startN = new Node (++id_max, P, 1.0, Node::DOMAIN_BOUNDARY_FIXED);
  homeN  = startN;
  nodes.add (homeN);
  loop.add  (homeN -> ID());

  do {
    switch (Keyword = parse (*file)) {
    case PUSH:
      if (pushing) error (prog, "pushes cannot be nested (yet)", ERROR);
      pushing = 1;
      precN   = startN;
      break;
    case STRAIGHT: case ARC:
      *file >> nP;
      if (nP < 1 || nP > INS_MAX) {
	sprintf (err, "Number of points (%1d) out of range", nP);
	error (prog, err, ERROR);
      }
      moreInput = generateNodes (*file, Keyword, homeN, startN, nP,
				 nodes, loop, id_max, pushing, save);
      break;
    case POP:
      if (!pushing) error (prog, "pop not matched to push", ERROR);
      pushing = 0;
      nP      = 0;
      while (!save.isEmpty()) {
	id = save.pop();
	if (nP) loop.add (id);
	nP++;
      }
      loop.add (precN -> ID());
      startN =  precN;
      break;
    default:
      error (prog, "never get here", ERROR);
      break;
    }
  } while (moreInput);

  setSizes (nodes, loop);
  printup  (cout, nodes, loop);

  return EXIT_SUCCESS;
}


static void getArgs (int       argc ,
		     char**    argv ,
		     istream*& input)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char buf[StrMax], c;
  char usage[] = "Usage: %s [options] [file]\n"
    "  [options]:\n"
    "  -h ... print this message\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog);
      cerr << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) {
      cerr << usage;
      sprintf (buf, "unable to open file: %s", *argv);
      message (prog, buf, ERROR);
    }
  } else input = &cin;
}


static key parse (istream& S)
// ---------------------------------------------------------------------------
// Look for (case-insensitive) string corresponding to keyword.
// ---------------------------------------------------------------------------
{
  char routine[] = "parse", org[StrMax], buf[StrMax];

  S >> buf;
  strcpy (org, buf);
  upperCase (buf);

  if (strstr (buf, "STRAIGHT"))
    return STRAIGHT;
  else if (strstr (buf, "ARC"))
    return ARC;
  else if (strstr (buf, "PUSH"))
    return PUSH;
  else if (strstr (buf, "POP"))
    return POP;
  else {
    sprintf (buf, "can't understand a keyword in %s", org);
    error   (routine, buf, ERROR);
  }
  
  return UNDEFINED;
}


static int generateNodes (istream&     infile ,
			  const key&   Keyword,
			  Node*&       homeN  ,
			  Node*&       startN ,
			  const int&   nP     ,
			  List<Node*>& nodes  ,
			  List<int>&   loop   ,
			  int&         id_max ,
			  const int&   push   ,
			  Stack<int>&  save   )
// ---------------------------------------------------------------------------
// Generate Nodes along a line, according to options.
//
// -- Example inputs:
// straight 10 uniform
// straight  8 geometric 1.1
// arc       4 -2.0
// 
// After the line-type command (which has already been read and is given
// by input Keyword) the following integer declares the number of Points
// to be generated (not including input starting Point startP).  This
// has also been read and is supplied as nP.
//
// For straight lines, points may be uniformly spaced, or computed as a
// geometric progression according to the supplied growth rate.
//
// For circular arcs, the supplied argument gives the radius of the arc.
// Negative values imply that the arc is concave with respect to the loop
// (i.e. reduces area compared to a straight line), and positive values
// imply the arc is convex.  Arc divisions are uniform.  The maximum angle
// that may be subtended by any arc segment is PI; larger arcs must be
// constructed as a series of smaller ones.
//
// Reset startP to be the final Point generated.
// ---------------------------------------------------------------------------
{
  char          routine[] = "generateNodes", err[StrMax], buf[StrMax];
  register int  i;
  int           goHome = 0;
  Point         dP, endP, startP = startN -> pos(), homeP = homeN -> pos();
  Node*         N;
  vector<Point> point(nP);

  // -- Compute points.

  switch (Keyword) {

  case STRAIGHT:
    {
      infile >> buf;
      strcpy (err, buf);
      upperCase (err);

      if (strstr (err, "UNIFORM")) {
	infile >> endP;
	if (!infile) { goHome = 1; endP = homeP; }
	dP = endP - startP;
	for (i = 0; i < nP; i++)
	  point[i] = (i + 1) / (real) nP * dP + startP;

      } else if (strstr (err, "GEOMETRIC")) {
	real initial, growth, len, dl, frac;
	
	infile >> growth;
	if (!infile) error (routine, "can't find geom. growth rate", ERROR);
	infile >> endP;
	if (!infile) { goHome = 1; endP = homeP; }
	dP  = endP - startP;
	len = endP.distance (startP);

	if (fabs (growth - 1.0) < EPSSP) { // -- Use uniform spacing.
	  for (i = 0; i < nP; i++)
	    point[i] = (i + 1) / (real) nP * dP + startP;
	} else {			         // -- Geometric progression.
	  initial = len * (1.0 - growth) / (1.0 - pow (growth, nP));
	  dl      = initial;
	  for (i = 0; i < nP; i++) {
	    frac     = dl / len;
	    point[i] = frac * dP + startP;
	    initial *= growth;
	    dl      += initial;
	  }
	}

      } else {
	sprintf (err, "invalid straight line modifier: %s", buf);
	error (routine, err, ERROR);
      }
    }
    break;

  case ARC:
    {
      real  radius, len, semiangle, angle;
      int   sign;
      Point midPoint, n, centre;

      infile >> radius;
      if (!infile) error (routine, "can't find arc's radius", ERROR);
      infile >> endP;
      if (!infile) { goHome = 1; endP = homeP; }
    
      len    = endP.distance (startP);
      sign   = (radius < 0.0) ? -1      :      1;
      radius = (sign   <   0) ? -radius : radius;
      
      if (2.0 * radius < len - EPSSP) {
	sprintf (err, "arc, radius %g can't span gap %g", radius, len);
	error (routine, err, ERROR);
      }

      dP = endP - startP;
      n  = unitNormal (startP, endP);
      midPoint  = 0.5 * (startP + endP);
      semiangle = asin (0.5 * len / radius);

      centre = -sign * cos (semiangle) * radius * n + midPoint;

      if (sign == 1)
	for (i = 0; i < nP; i++) {
	  angle = (i + 1) / (real) nP * 2.0 * semiangle;
	  point[i] = centre.relative (startP, radius, angle);
	}
      else 
	for (i = 0; i < nP; i++) {
	  angle = (nP - i - 1) / (real) nP * 2.0 * semiangle;
	  point[i] = centre.relative (endP, radius, angle);
	}
    }
    break;

  default:
    error (routine, "can't get here", ERROR);
    break;
  }

  // -- Load computed Points into appropriate storage.

  int nn = (goHome) ? nP - 1 : nP;

  if (push)
    for (i = 0; i < nn; i++) {
      N = new Node (++id_max, point[i], 1.0, Node::INTERIOR);
      nodes.add (N);
      loop.add  (N -> ID());
      save.push (N -> ID());
    }
  else
    for (i = 0; i < nn; i++) {
//    N = new Node (++id_max, point[i], 1.0, Node::BOUNDARY);
      N = new Node (++id_max, point[i], 1.0, Node::DOMAIN_BOUNDARY_FIXED);
      nodes.add (N);
      loop.add  (N -> ID());
    }

  startN = N;
  
  return !goHome;
}


static void setSizes (List<Node*>&     nodes,
		      const List<int>& loop )
// ---------------------------------------------------------------------------
// Set ideal sizes for Nodes as average of distance to contacting Nodes.
// ---------------------------------------------------------------------------
{
  const int           N = loop.length();
  vector<Node*>       finalLoop(N);
  register int        found, id, i = 0;
  Node                *prec, *curr, *next;
  real                len;
  ListIterator<Node*> n(nodes);
  ListIterator<int>   l(loop);

  for (; l.more(); l.next(), i++) {
    id    = l.current();
    found = 0;
    for (n.reset(); !found && n.more(); n.next()) {
      curr  = n.current();
      found = curr -> ID() == id;
    }
    finalLoop[i] = curr;
  }

  for (i = 0; i < N; i++) {
    prec = finalLoop[(i + N - 1) % N];
    curr = finalLoop[i];
    next = finalLoop[(i     + 1) % N];
    len  = curr -> pos().distance (prec -> pos());
    len += curr -> pos().distance (next -> pos());
    len *= 0.5;
    curr -> setSize (len);
  }
}


static void printup (ostream&           S    ,
		     const List<Node*>& nodes,
		     const List<int>&   loop )
// ---------------------------------------------------------------------------
// Print up input used by qmesh.
// ---------------------------------------------------------------------------
{
  S << nodes.length() << " BOUNDARY NODES" << endl;
  for (ListIterator<Node*> n(nodes); n.more(); n.next())
    S << *n.current() << endl;
  S << endl;

  S << loop.length() << " NODE LOOP" << endl;
  for (ListIterator<int> l(loop); l.more(); l.next())
    S << setw(5) << l.current();
  S << endl;
}
