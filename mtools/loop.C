///////////////////////////////////////////////////////////////////////////////
// loop.C
//
// The Loop class is central to the mesh generation algorithm.  All
// new Loops are created by dividing a parent Loop into left and right
// subLoops, hence it is natural to use a binary tree data structure
// as the underlying paradigm.
// 
// Loops always have an even number of Nodes.  Loops are recursively
// subdivided by a dividing-line strategy until they end up as either
// 4-Noded or 6-Noded.  Each 6-Noded Loop is then split into two
// subLoops of various sizes according to a set of rules enshrined in
// Loop::splitSix: these subLoops are passed back to the recursion.
// Finally all Loops are 4-Noded and the process terminates.
//
// Exceptions to this sequence occur for "offset" type boundary nodes,
// where 4-Noded subLoops are created before the splitting line
// strategy is employed, however the division process is the same as
// above in that the 4-Noded subLoops are split off the parent Loop
// one at a time until all offset Nodes are consumed.  This is a
// recursive process too.
// 
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <qmesh.h>


static inline int CW  (const int i, const int n) { return (i + n - 1) % n; }
static inline int CCW (const int i, const int n) { return (i     + 1) % n; }


Loop::Loop (const int& numnodes) :
 left  (0),
 right (0)
// ---------------------------------------------------------------------------
// Default constructor.  It is assumed that Loops are created sequentially
// and never destroyed.
// ---------------------------------------------------------------------------
{
  id = ++Global::loopIdMax;

  nodes    .setSize (numnodes);
  splitline.setSize (0);
}


Loop::Loop (vector<Node*>& bound,
	    vector<Node*>& split,
	    const int&     begin,
	    const int&     end  ,
	    const int&     dir  ) :
  left  (0),
  right (0)
// ---------------------------------------------------------------------------
// This is used to create a child Loop from a parent Loop (bound)
// and a splitline (split).  Split may be empty.
// 
// 1. Find indices of start and finish Nodes in bound.
// 2. Determine number of Nodes which come from bound and split.
// 3. Create vector of Nodes for this loop.
// 4. Start fill by traversing split, forwards if dir > 0, else backwards.
// 5. Complete fill by installing Nodes from bound.
// ---------------------------------------------------------------------------
{
  id = ++Global::loopIdMax;

  char routine[] = "Loop::Loop";

  register int i;
  int          nBound;
  const int    N = bound.getSize();
  const int    M = split.getSize();

  if (begin == end)
    if (M == 3 && dir == 0) {
      nodes.setSize (4);       // -- Special (last) case in Loop::splitSix.
      for (i = 0; i < M; i++) nodes[i] = split[M - i - 1];
      nodes[3] = nodes[begin];
      return;
    } else if (M != 3 && dir != 0) {
      message (routine, "start, finish Nodes same", ERROR);
    }

  nBound = (end > begin) ? end - begin + 1 : N - begin + end + 1;

  nodes.setSize (nBound + M);

  if   (dir > 0) for (i = 0; i < M; i++) nodes[i] = split[i];
  else           for (i = 0; i < M; i++) nodes[i] = split[M - i - 1];

  for (i = 0; i < nBound; i++) nodes[M + i] = bound[(begin + i) % N];
}


void Loop::split ()
// ---------------------------------------------------------------------------
// This is the heart of the generation algorithm.  Each loop with more
// than six nodes is recursively split into two subloops, as defined by
// a splitting line.
//
// The regular sequence of events is:
//
// 0. Return if this Loop has four Nodes.
// 1. For each node in Loop, create list of visible Nodes.
// 2. Determine the start and end Nodes of best splitting line.
// 3. Generate new Nodes on this splitting line.
// 4. Generate "left" and "right" subloops.
// 5. Split left and right subloops.
//
// However, if the Loop has six Nodes:
//
// 1. Classify into subgroups based on convexity and number of
//    180-degree interior angles.
// 2. If Loop is concave then immediately split in two four-Noded
//    Loops from concave corner.  A concave Loop has at least one interior
//     angle greater than 180 degrees.
// 3. Otherwise, classify into based on the number of 180-degree
//    interior angles, generate splitting to produce:
//      2 x 4-Noded subloops;
//      1 x 4-Noded subloop, 1 x 6-Noded concave subloop;
//      2 x 6-Noded concave subloops;
//      1 x 4-Noded subloop, 1 x 10-Noded subloop.
// 4. Attempt to split left and right subloops.
// 
// It is possible that the Loop has already been split when this routine
// is called at top level (by generation of offset Nodes), so we descend
// to bottom level before starting split.
// ---------------------------------------------------------------------------
{
  char routine[] = "Loop::split", err[StrMax];

  register int i;
  int          begin, end, insert;
  const int    N = nodes.getSize ();
  
  if (N & 1) {
    sprintf (err, "Loop %1d has odd number of Nodes (%1d)", id, N);
    message (routine, err, WARNING);
  }

  if (graphics) drawLoop (this);
  if (Global::verbose)
    cout << routine << ": loop " << id <<" ("<< N << " nodes)" << endl;
  if (Global::verbose > 1)
    cout << *this;
    
  if (N == 4)			// -- Can't go any smaller.
    return;

  if (!left && !right) {	// -- This may already have been split.
  
    if (N == 6)        // -- Split 6-Noded loop based on type.
      
      splitSix (begin, end);

    else {		        // -- Regular splitting on straight line.

      // -- Create list of visible nodes for each node in nodes.

      List<Node*>* visible = new List<Node*> [nodes.getSize ()];
      int          OK      = 0;

      for (i = 0; i < N; i++) {
	visibleNodes (visible[i], i);
	OK += visible[i].length ();
      }

      if (!OK) {
	sprintf (err, "Loop %1d cannot be split by straight line", id);
	message (routine, err, ERROR);
      }

      // -- Determine start and end nodes of best splitting line, create it.

      bestSplit   (visible, begin, end, insert);
      insertNodes (nodes[begin], nodes[end], insert);
    }
    
    // -- Create left and right subloops.

    left  = new Loop (nodes, splitline, end, begin, 1);
    right = new Loop (nodes, splitline, begin, end, 0);
  }

#ifdef PROMPT
  qpause ();
#endif

  if        (left  -> area() < EPSSP) {
    message (routine, "--> This Loop has zero area:", WARNING);
    cout << *left;
    sprintf (err, "LH Loop %1d has zero area, deleting", left -> ID());
    message (routine, err, WARNING);
    left = 0;
  } else if (right -> area() < EPSSP) {
    message (routine, "--> This Loop has zero area:", WARNING);
    cout << *right;
    sprintf (err, "RH Loop %1d has zero area, deleting", left -> ID());
    message (routine, err, WARNING);
    right = 0;
  }
    
  if (left)  left  -> split ();
  if (right) right -> split ();
}


void Loop::visibleNodes (List<Node*>& V,
			 const int&   i) const
// ---------------------------------------------------------------------------
// Into V place a list of all nodes in this loop which are visible
// from nodes[i].  Here we use  a different approach to that
// explained in Ref. [1], which is not robust.
//
// The Node angles are computed for nodes from the next CCW Node around
// to the next CW Node.  Use brute force approach of checking for 
// intersection of line segments, per Ref[3].  In addition, the unit normals
// of the two line segments either side of the point under examination
// for visibility have to point towards nodes[i].
// ---------------------------------------------------------------------------
{
  register int j, k, l, m, n;
  const int    N  = nodes.getSize ();
  const int    N1 = N - 1;
  int          visible;
  real         rangle, tangle;

  rangle = nodes[i] -> pos().angle (nodes[CCW(i, N)] -> pos(), 
				    nodes[ CW(i, N)] -> pos());
  rangle = (rangle < EPSSP) ? TWOPI - EPSSP : rangle - EPSSP;

  for (j = 2; j < N1; j++) {
    k = (i + j) % N;
    if (!cull (nodes[       i] -> pos(), nodes[        k] -> pos(),
	       nodes[CW(k, N)] -> pos(), nodes[CCW(k, N)] -> pos()) &&
	!cull (nodes[       k] -> pos(), nodes[        i] -> pos(),
	       nodes[CW(i, N)] -> pos(), nodes[CCW(i, N)] -> pos())) {

      tangle = nodes[i] -> pos().angle (nodes[CCW(i, N)] -> pos(),
					nodes[        k] -> pos());
      if (nodes[i] != nodes[k] && tangle > EPSSP && tangle < rangle) {
	visible = 1;
	for (l = 2; visible && l < N1; l++) {
	  m = (i + l) % N;
	  n = CCW (m, N);
	  if (nodes[m] == nodes[k] || nodes[n] == nodes[k])
	    continue;
	  if (fabs (nodes[m] -> pos().turn (nodes[i] -> pos(),
					    nodes[n] -> pos())) < EPSSP)
	    continue;
	  visible = visible && !cross (nodes[i] -> pos(), nodes[k] -> pos(), 
				       nodes[m] -> pos(), nodes[n] -> pos()); 
	}
	if (visible) V.xadd (nodes[k]);
      }
    }
  }

#ifdef DEBUG
  cout << "Node " << nodes[i] -> ID () << " (angle = " << rangle << "):";
  cout << " Nodes visible :";
  for (ListIterator<Node*> c (V); c.more (); c.next ())
    cout << "  " << c.current () -> ID ();
  cout << endl;
#endif
}


void Loop::bestSplit (List<Node*>* visible,
		      int&         start  ,
		      int&         end    ,
		      int&         NNodes ) const
// ---------------------------------------------------------------------------
// Return indices of start and end nodes of best splitting line in loop.
// For each node in loop, find best candidate line.  Return the best of these.
// 
// See the description in \S\,4.4 of Reference [1].
// ---------------------------------------------------------------------------
{
  char routine[] = "Loop::bestSplit";

  const int N = nodes.getSize();
  int       i, best   = 0;
  real      param     = 1.0e30;
  real      reflength = lengthScale();

  vector<real> SL     (N);	// -- Performance index for each split.
  vector<int>  End    (N);	// -- Matching end node index for start node.
  vector<int>  insert (N);	// -- Number of nodes to insert.

  SL = 1.0e30;

  for (i = 0; i < N; i++) {
    bestLine (visible[i], reflength, i, End[i], insert[i], SL[i]);
    if (SL[i] < param) {
      param = SL[i];
      best  = i;
    }
  }

  start  = best;
  end    = End[best];
  NNodes = insert[best];

  if (Global::verbose > 1) cout << routine << ": Loop " << id
    << ", Node " << nodes[start] -> ID() << " <--> Node " << nodes[end] -> ID()
      << ", " << NNodes << " insertions" << endl;
}


void Loop::bestLine (List<Node*>& visible ,
		     const real&  lenscale,
		     const int&   start   ,
		     int&         end     ,
		     int&         numnodes,
		     real&        pindex  ) const
// ---------------------------------------------------------------------------
// For the given starting Node, traverse the list of visible Nodes, estimate
// the number of Nodes to insert for the least size gradient error, compute
// the angle and length errors, combine to form an overall performance index.
// Return the optimal combination of end Node and number of Nodes to
// insert, as well as matching performance index.
//
// The loop which we are going to split has an even number of Nodes.
// Likewise the subloops which will be formed have to have an even number
// of Nodes.  This constrains how many Nodes can be inserted: e.g. if the
// proposed subloop has an odd number of Nodes, then an odd number of
// Nodes have to be placed along the splitting line to ensure subloops
// with even numbers of Nodes.  (At this stage, no new Nodes are actually
// generated.)
//
// Angle errors are assigned zero for cusps.
// ---------------------------------------------------------------------------
{
  const int M = visible.length ();

  if (!M) return;

  char routine[] = "Loop::bestLine", err[StrMax];

  int       i, j = 0, odd, nIns, cusp;
  const int N = nodes.getSize  ();

  real  theta, L, epsilon;
  real  da, dc, dx, dy, actualLength, epsOld, epsNew, SL, best = 1.0e30;
  real  theta1, theta2, theta3, theta4;
  Point pi, pj, pk, pa, pc, na, nc;

  vector<int> index (M);
  vector<int> ins   (M);

  // -- Generate indices in Loop vector "nodes" of visible Nodes.

  for (ListIterator<Node*> n (visible); n.more (); n.next (), j++)
    for (i = 0; i < N; i++)
      if (n.current () == nodes[i]) index[j] = i;

  // -- For each end Node, find:

  for (j = 0; j < M; j++) {
    end = index[j];
    
    // -- the relative length of the splitting line;

    dx = nodes[end] -> pos().x - nodes[start] -> pos().x;
    dy = nodes[end] -> pos().y - nodes[start] -> pos().y;
    actualLength = hypot (dx, dy);

    L = actualLength / lenscale;

    // -- the number of insertions that generates the smallest length error;

    odd    = (start & 1 && end & 1) || (!(start & 1) && !(end & 1));
    nIns   = (odd) ? 1 : 0;
    epsOld = spaceNodes (nodes[start], nodes[end], actualLength, nIns);

    while (nIns < InsMax) {	// -- Monotone search (inefficient?).
      nIns  += 2;
      epsNew = spaceNodes (nodes[start], nodes[end], actualLength, nIns);
      if (epsNew > epsOld) break;
      epsOld = epsNew;
    }
    
    if (nIns >= InsMax) {
      sprintf (err, "exceeded max insertions (%1d), Loop %1d, "
	       "Nodes %1d <--> %1d", InsMax, id,
	       nodes[start] -> ID(), nodes[end] -> ID());
      message (routine, err, ERROR);
    } else {
      nIns   -= 2;
      ins[j]  = nIns;
      epsilon = epsOld;
    }

    // -- the splitting angle deviation error;
    
    pi = nodes[CW  (start, N)] -> pos();
    pj = nodes[     start    ] -> pos();
    pk = nodes[CCW (start, N)] -> pos();

    da = pi.distance (pj);
    dc = pk.distance (pj);
    pa = 1.0 / da * (pi - pj);
    pc = 1.0 / dc * (pj - pk);
    na = unitNormal (pi,  pj);
    nc = unitNormal (pj,  pk);
    
    if (cusp = pj.turn (pi, pk) < EPSSP && na.dot (nc) < EPSSP) {
      theta1 = theta2 = M_PI_2;
    } else {
      theta1 = pj.angle (pk, nodes[end]->pos());
      theta2 = pj.angle (nodes[end]->pos(), pi);
    }
    
    pi = nodes[CW  (end, N)] -> pos();
    pj = nodes[     end    ] -> pos();
    pk = nodes[CCW (end, N)] -> pos();

    da = pi.distance (pj);
    dc = pk.distance (pj);
    pa = 1.0 / da * (pi - pj);
    pc = 1.0 / dc * (pj - pk);
    na = unitNormal (pi,  pj);
    nc = unitNormal (pj,  pk);
    
    if (cusp = pj.turn (pi, pk) < EPSSP && na.dot (nc) < EPSSP) {
      theta3 = theta4 = M_PI_2;
    } else {
      theta3 = pj.angle(pk, nodes[start]->pos());
      theta4 = pj.angle(nodes[start]->pos(), pi);
    }

    theta1 = fabs (theta1 - M_PI_2);
    theta2 = fabs (theta2 - M_PI_2);
    theta3 = fabs (theta3 - M_PI_2);
    theta4 = fabs (theta4 - M_PI_2);

    theta  = (theta1 + theta2 + theta3 + theta4) / TWOPI;

    // -- the combined performance index, eq. (1).  Keep index of least SL.

    SL = Global::C1*theta + Global::C2*L + Global::C3*epsilon;
    
    if (SL < best) { best = SL; i = j; }
  }

  pindex   = best;  
  end      = index[i];
  numnodes = ins  [i];
}


istream& operator >> (istream& s,
		      Loop&    l)
// ---------------------------------------------------------------------------
// Read in Loop from Node list in stream.
//
// 1 2 3 4 5 ... NN ...     { Loop is assumed closed by return to start tag.
//
// MM specifiers must be on lines of their own, but input is
// otherwise free-format.
//
// The first tag does not need to be re-specified at the end of the node loop
// tag list: MM is the number of nodes in the loop not including the return
// to the first node (e.g. a quad loop would have MM = 4, not 5).
// ---------------------------------------------------------------------------
{
  char  routine[] = "operator >> (istream&, Loop&)";
  char                err[StrMax];
  int                 i, id, npts, found;
  register Node*      N;
  ListIterator<Node*> n (Global::nodeList);

  npts = l.nodes.getSize ();

  for (i = 0; i < npts; i++) {
    if (!(s >> id)) {
      sprintf (err, "problem reading tag for loop node %1d", i);
      message (routine, err, ERROR);
    }
    
    for (found = 0, n.reset(); !found && n.more(); n.next()) {
      N = n.current();
      found = N -> ID() == id;
    }

    if (!found) {
      sprintf (err, "couldn't find node tag %1d in input list", id);
      message (routine, err, ERROR);
    } else
      l.nodes[i] = N;
  }

  return s;
}


ostream& operator << (ostream& s,
		      Loop&    l)
// ---------------------------------------------------------------------------
// Recursively print up Loop information on s.
// ---------------------------------------------------------------------------
{
  if (l.left) s << *(l.left);

  s << "-- Loop " << l.id << " --" << endl;
 
  int       i;
  const int N = l.nodes.getSize ();

  for (i = 0; i < N; i++) s << *l.nodes[i] << endl;
  
  if (l.right) s << *(l.right);
  
  return s;
}


int Loop::points (vector<real>& x,
		  vector<real>& y) const
// ---------------------------------------------------------------------------
// Load locations of Nodes into x & y.
// ---------------------------------------------------------------------------
{
  register int   i;
  register Point p;
  const int      N = nodes.getSize ();
  
  x.setSize (N);
  y.setSize (N);

  for (i = 0; i < N; i++) {
    p = nodes[i] -> pos();
    x[i] = p.x;
    y[i] = p.y;
  }

  return N;
}


int Loop::line (vector<real>& x,
		vector<real>& y) const
// ---------------------------------------------------------------------------
// Load locations of splitline Nodes into x & y.
// ---------------------------------------------------------------------------
{
  register int   i;
  register Point p;
  const int      N = splitline.getSize ();
  
  x.setSize (N);
  y.setSize (N);

  for (i = 0; i < N; i++) {
    p = splitline[i] -> pos();
    x[i] = p.x;
    y[i] = p.y;
  }

  return N;
}


real Loop::spaceNodes (const Node*    begin   ,
		       const Node*    end     ,
		       const real&    distance,
		       const int&     insertN ) const
// ---------------------------------------------------------------------------
// This routine determines mesh refinement along a splitting line.
//
// Given a beginning node, an end node, and a number of nodes to insert,
// compute and return the error index epsilon in eq. (1), Reference [1].
//
// See \S\,4.5 of Reference [1].
//
// Refinement coefficient defines approximate proportion of splitting line
// that is to be occupied by uniform subdivisions at smallest size.
// ---------------------------------------------------------------------------
{
  int  i;
  real scale;

  if (insertN == 0) {
    scale = distance / (end -> prefSize() + begin -> prefSize());
    return fabs (1.0 - scale);
  }

  const Node   *N1, *N2;
  vector<real> spacing (insertN + 1);

  if   (begin -> prefSize() <= end -> prefSize()) { N1 = begin; N2 = end;   } 
  else                                            { N1 = end;   N2 = begin; }

  const real b = N1 -> prefSize ();
  const int  M = (int) (Global::refCoeff * (real) insertN);

  real idealLength   = 0;
  real idealGradient = (N2 -> prefSize() - N1 -> prefSize()) / (insertN - M);

  for (i = 0; i < M; i++) {
    spacing[i]   = N1 -> prefSize();
    idealLength += spacing[i];
  }
  for (i = M; i <= insertN; i++) {
    spacing[i]   = (i - M) * idealGradient + b;
    idealLength += spacing[i];
  }
  scale = distance / idealLength;

  return fabs (1.0 - scale);
}


void Loop::insertNodes (const Node* begin  ,
		        const Node* end    ,
		        const int&  insertN)
// ---------------------------------------------------------------------------
// This routine inserts Nodes along a splitting line.
//
// See \S\,4.5 of Reference [1].
// ---------------------------------------------------------------------------
{
  if (insertN < 1) return;

  splitline.setSize (insertN);

  int          i;
  real         propn, idealGradient, idealLength = 0.0, sumlen = 0.0;
  vector<real> spacing (insertN + 1);

  const Node   *N1, *N2;
  Node         *newNode, *oldNode;

  if   (begin -> prefSize() <= end -> prefSize()) { N1 = begin; N2 = end;   } 
  else                                            { N1 = end;   N2 = begin; }

  const real b = N1 -> prefSize();
  const int  M = (int) (Global::refCoeff * (real) insertN);
  Point      P, D = N2 -> pos() - N1 -> pos();

  idealGradient = (N2 -> prefSize() - N1 -> prefSize()) / (insertN - M);

  for (i = 0; i < M; i++) {
    spacing[i]   = N1 -> prefSize();
    idealLength += spacing[i];
  }
  for (i = M; i <= insertN; i++) {
    spacing[i]   = (i - M) * idealGradient + b;
    idealLength += spacing[i];
  }

  for (i = 0; i < insertN; i++) {
    sumlen += spacing[i];
    propn   = sumlen / idealLength ;
    P = N1 -> pos() + propn * D;
    newNode = new Node (++Global::nodeIdMax, P, 0.5*(spacing[i]+spacing[i+1]));
    if (oldNode = Global::exist (newNode)) {
      delete newNode;
      Global::nodeIdMax--;
      newNode = oldNode;
    }
    if (N1 == begin)
      splitline[i] = newNode;
    else
      splitline[insertN - i - 1] = newNode;
    Global::nodeList.add (newNode);
  }
}


Point Loop::centroid () const
// ---------------------------------------------------------------------------
// Return centroid of this Loop's Nodal positions.
// ---------------------------------------------------------------------------
{
  register int   i;
  register Point P;
  const int      N = nodes.getSize();

  P = 0.0;
  for (i = 0; i < N; i++) P += nodes[i] -> pos();
  P *= 1.0 / N;

  return P;
}


void Loop::splitSix (int& begin, int& end)
// ---------------------------------------------------------------------------
// Find what type of six-Noded Loop *this is.  Set up splitline to
// split *this into two subloops.  Return indices of begin & end Nodes.
//
// For all cases where new Nodes are created, their ideal mesh size is
// set as the size of the Loop.
//
// The code deals with a bunch of special cases, so it's fairly hideous.
// ---------------------------------------------------------------------------
{
  char         routine[] = "Loop::splitSix";

  register int i1, i2, i3, i4, i5, i6;
  Point        P, P1, P2, D, D1, D2, D3;
  real         d1, d2, d3, dir;
  const real   EPS = (sizeof (real) == sizeof (float)) ? EPSSP : EPSDP;
  int          n180 = 0;
  vector<int>  index180 (6);

  // -- Find number of 180-degree interior angles, and their indices.
  //    If any concave corners are found, deal with immediately & return.

  for (i1 = 0; i1 < 6; i1++) {
    i2 = (i1     + 1) % 6;
    i3 = (i1 + 6 - 1) % 6;
    dir = nodes[i1] -> pos().turn (nodes[i3] -> pos(), nodes[i2] -> pos());
    if (fabs (dir) < EPS)
      index180[n180++] = i1;
    else if (dir <= -EPS) {	// -- Bingo.  A concave loop.
      if (Global::verbose)
	cout << routine << ": Loop " << id << " is concave" << endl;
      splitline.setSize (0);
      begin =  i1;
      end   = (i1 + 3) % 6;
      return;
    }
  }

  if (Global::verbose) cout << routine << ": Loop " << id << " has " << n180
      << " 180-degree interior angles" << endl;

  // -- All inputs guaranteed convex from now on.

  Global::limits (P1, P2);
  const real psize = 0.5 * P2.distance (P1);
  Node       *oldNode, *newNode;
  int        nodd, nsum;

  switch (n180) {

  case 0:
    // -- Split into 2 x 4 on shortest line.
    //
    d2 = 1.0e30;

    for (i1 = 0; i1 < 3; i1++) {
      D   = nodes[(i1 + 3) % 6] -> pos() - nodes[i1] -> pos();
      d1 = D.magnitude ();
      if (d1 < d2) {
	d2 = d1;
	i2 = i1;
      }
    }
    splitline.setSize (0);
    begin =  i2;
    end   = (i2 + 3) % 6;
    return;

  case 1:
    // -- Split into 2 x 4 across from 180-degree.
    //
    splitline.setSize (0);
    begin = index180[0];
    end   = (begin + 3) % 6;
    return;
      
  case 2:
    // -- 3 options, choose on number of nodes between.
    //
    i5 = index180[0],
    i6 = index180[1];

    switch (i6 - i5) {	// -- No. of Nodes (0 < i < 6) to next 180-deg.
    case 1:
      // -- Insert two new nodes, spaced along line subdividing element;

      i3 = (i5 + 6 - 1) % 6;
      i4 = (i5 + 6 - 2) % 6;
      P1 = 0.5 * (nodes[i3] -> pos() + nodes[i4] -> pos());

      i3 = (i6 + 1) % 6;
      i4 = (i6 + 2) % 6;
      P2 = 0.5 * (nodes[i3] -> pos() + nodes[i4] -> pos());

      D1 = P2 - P1;
      d1 = D1.magnitude ();

      i3 = (i5 + 6 - 1) % 6;
      i4 = (i6     + 1) % 6;

      D2 = nodes[i4] -> pos() - nodes[i3] -> pos();
      d2 = D2.magnitude ();

      splitline.setSize (2);

      D3 = nodes[i5] -> pos() - nodes[i3] -> pos();
      d3 = D3.magnitude ();

      P2 = P1 + d3 / d2 * D1;

      newNode = new Node (++Global::nodeIdMax, P2, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);	
      
      splitline[0] = newNode;
      
      D3 = nodes[i6] -> pos() - nodes[i3] -> pos();
      d3 = D3.magnitude ();

      P2 = P1 + d3 / d2 * D1;

      newNode = new Node (++Global::nodeIdMax, P2, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);
      splitline[1] = newNode;
      
      begin = i5;
      end   = (i6 + 2) % 6;
      return;
    case 2:
      // -- Insert 1 new node at centroid of remaining Nodes.
      P = 0.0;
      for (i2 = 0; i2 < 6; i2++)
	if (i2 != i5 && i2 != i6)
	  P += nodes[i2] -> pos();
      P *= 0.25;

      splitline.setSize (1);
      newNode = new Node (++Global::nodeIdMax, P, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);
      splitline[0] = newNode;

      begin = i6;
      end   = i5;
      return;
    case 3:
      // -- Break across middle.
      splitline.setSize (0);
      begin = i5;
      end   = i6;
      return;
    case 4:
      // -- Same as case 2, but other order.
      P = 0.0;
      for (i2 = 0; i2 < 6; i2++)
	if (i2 != i5 && i2 != i6)
	  P += nodes[i2] -> pos();
      P *= 0.25;

      splitline.setSize (1);
      newNode = new Node (++Global::nodeIdMax, P, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);
      splitline[0] = newNode;

      begin = i5;
      end   = i6;
      return;
    case 5:
      // -- Same as case 1, but other order.
      i3 = (i5 + 1) % 6;
      i4 = (i5 + 2) % 6;
      P1 = 0.5 * (nodes[i3] -> pos() + nodes[i4] -> pos());

      i3 = (i6 + 6 - 1) % 6;
      i4 = (i6 + 6 - 2) % 6;
      P2 = 0.5 * (nodes[i3] -> pos() + nodes[i4] -> pos());

      D1 = P2 - P1;
      d1 = D1.magnitude ();

      i3 = (i5  + 1) % 6;
      i4 = (i6 + 6 - 1) % 6;

      D2 = nodes[i4] -> pos() - nodes[i3] -> pos();
      d2 = D2.magnitude ();

      splitline.setSize (2);

      D3 = nodes[i5] -> pos() - nodes[i3] -> pos();
      d3 = D3.magnitude ();

      P2 = P1 + d3 / d2 * D1;

      newNode = new Node (++Global::nodeIdMax, P2, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);	
    
      splitline[0] = newNode;

      D3 = nodes[i6] -> pos() - nodes[i3] -> pos();
      d3 = D3.magnitude ();

      P2 = P1 + d3 / d2 * D1;

      newNode = new Node (++Global::nodeIdMax, P2, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);	
  
      splitline[1] = newNode;

      begin = i5;
      end   = (i6 + 6 - 2) % 6;
      return;

    default:
      message (routine, "impossible! two 180-deg. angles at same Node", ERROR);
      break;
    }
    break;

  case 3:

    // -- We're splitting a triangular shape.  Again, there are three options.
    //
    //    These we classify according to the maximum number of 
    //    consecutive 180-degree interior angles (1, 2, or 3).
    //    We can pick between these on basis of index180.
    
    nsum = 0;
    nodd = 0;
    for (i1 = 0; i1 < 3; i1++) {
      nodd += (index180[i1] & 1) ? 1 : 0;
      nsum +=  index180[i1];
    }
    
    if (nodd == 0 || nodd == 3) {
      // -- Rotational symmetry case.  New Node at centroid of other Nodes.
      
      P = 0.0;
      if (index180[0] == 0) {
	P += nodes[1] -> pos(); P += nodes[3] -> pos(); P += nodes[5] -> pos();
      } else {
	P += nodes[0] -> pos(); P += nodes[2] -> pos(); P += nodes[4] -> pos();
      }
      P *= 0.3333333333333;
        
      splitline.setSize (1);
      newNode = new Node (++Global::nodeIdMax, P, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);
      splitline[0] = newNode;

      begin = index180[0];
      end   = index180[1];

    } else if (nsum % 3) {
      // -- Here, there are 2 consecutive 180-degree nodes on a side.  Ugh.

      splitline.setSize (0);
      
      switch (100 * index180[0] + 10 * index180[1] + index180[2]) {
      case 13: case 23: case 34: case 35: 
	begin = 0; end = 3;
	break;
      case 14: case 124: case 134: case 145:
	begin = 1; end = 4;
	break;
      case 25: case 125: case 235: case 245:
	begin = 2; end = 5;
	break;
      default:
	message (routine, "unknown case splitting triangular loop", ERROR);
	break;
      }

    } else {
      // -- Three 180-degree Nodes on a side.  Make 1 x 4-Node, 1 x 10-Node.

      splitline.setSize (3);

      switch (100 * index180[0] + 10 * index180[1] + index180[2]) {
      case 123:
	i1 = 5; i2 = 0; i3 = 1; i4 = 2; i5 = 3; i6 = 4;
	break;
      case 12:
	i1 = 4; i2 = 5; i3 = 0; i4 = 1; i5 = 2; i6 = 3;
	break;
      case 15:
	i1 = 3; i2 = 4; i3 = 5; i4 = 0; i5 = 1; i6 = 2;
	break;
      case 45:
	i1 = 2; i2 = 3; i3 = 4; i4 = 5; i5 = 0; i6 = 1;
	break;
      case 345:
	i1 = 1; i2 = 2; i3 = 3; i4 = 4; i5 = 5; i6 = 0;
	break;
      case 234:
	i1 = 0; i2 = 1; i3 = 2; i4 = 3; i5 = 4; i6 = 5;
	break;
      default:
	message (routine, "unknown case splitting triangular loop", ERROR);
	break;
      }  
      P  = nodes[i1] -> pos();
      P1 = nodes[i5] -> pos();
      D1 = nodes[i6] -> pos();
      D  = P1 - D1;
      D1 = P  + D;
      P2 = nodes[i3] -> pos();
      D2 = nodes[i2] -> pos();
      D  = P2 - D2;
      D2 = P  + D;
      D3 = intersect (P1, D1, P2, D2);
      D1 = P1 + 0.66666666666 * (D3 - P1);
      D2 = P2 + 0.66666666666 * (D3 - P2);
      D3 =      0.33333333333 * (D1 + D2 + nodes[i4] -> pos());

      newNode = new Node (++Global::nodeIdMax, D1, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);	
      splitline[0] = newNode;

      newNode = new Node (++Global::nodeIdMax, D3, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);	
      splitline[1] = newNode;

      newNode = new Node (++Global::nodeIdMax, D2, psize);
      if (oldNode = Global::exist (newNode)) {
	delete newNode;
	Global::nodeIdMax--;
	newNode = oldNode;
      } else 
	Global::nodeList.add (newNode);
      splitline[2] = newNode;

      begin = i1;
      end   = i1;
    }
    return;

  default:
    message (routine, "impossible! More than three 180-degree angles", ERROR);
    break;
  }
}


void Loop::connect ()
// ---------------------------------------------------------------------------
// Recursively visit all Loops and for each Node, add information about
// the Nodes it is connected to.  We only do this at the level for which
// there are four Nodes in each Loop.
// ---------------------------------------------------------------------------
{
  if (nodes.getSize() == 4) {

    register int        i, found1, found2;
    register Node*      N;
    register Node*      N1;
    register Node*      N2;
    ListIterator<Node*> n (Global::nodeList);
    
    for (i = 0; i < 4; i++) {
      N1 = nodes[i];
      N2 = nodes[(i + 1) % 4];
      
      found1 = found2 = 0;
      for (n.reset (); n.more (); n.next ()) {
	N = n.current ();
	if (!found1)
	  if (found1 = (N == N1)) N -> xadd (N2);
	if (!found2)
	  if (found2 = (N == N2)) N -> xadd (N1);
	if ( found1 && found2 ) break;
      }
    }
    return;
  }

  if (left)  left  -> connect ();
  if (right) right -> connect ();
}


void Loop::smooth ()
// ---------------------------------------------------------------------------
// Laplacian smoothing.  Visit each Node, move non-boundary nodes to
// centroid of connected Nodes.
// ---------------------------------------------------------------------------
{
  ListIterator<Node*> n (Global::nodeList);
  register Node*      N;
  Point               cen;

  for (; n.more (); n.next ()) {
    N   = n.current ();
    cen = N -> centroid ();
    N -> setPos (cen);
  }
}


void Loop::drawQuad (const int& numbers) const
// ---------------------------------------------------------------------------
// This is for completed Loops (only do 4-Noded Loops).
// ---------------------------------------------------------------------------
{
  const int N = nodes.getSize ();

  if (N == 4) drawLoop (this, numbers);

  if (left)  left  -> drawQuad (numbers);
  if (right) right -> drawQuad (numbers);
}


void Loop::quads (List<Quad*>& elements) const
// ---------------------------------------------------------------------------
// Load 4-noded loops into elements.  We could just point to the nodes
// vector for each 4-noded loop, but this way, elements is subsequently
// completely independent of the Loop tree.
// ---------------------------------------------------------------------------
{
  int       i;
  Quad*     E;
  const int N = nodes.getSize ();

  if (N == 4) {
    E = new Quad;
    for (i = 0; i < 4; i++) E -> vertex[i] = nodes[i];
    elements.add (E);
  }

  if (left)  left  -> quads (elements);
  if (right) right -> quads (elements);
}


real Loop::area () const
// ---------------------------------------------------------------------------
// Compute loop area using cross-products.  See Ref. [2] \S\,7.8.5.
// ---------------------------------------------------------------------------
{
  register int  i;
  register real SA = 0.0;
  const int     N1 = nodes.getSize() - 1;
  const Point   P  = nodes[0] -> pos();

  for (i = 1; i < N1; i++)
    SA += P.turn (nodes[i + 1] -> pos(), nodes[i] -> pos());
  
  return 0.5 * SA;
}


void Loop::offset ()
// ---------------------------------------------------------------------------
// Traverse loop CCW and generate offset Nodes if specified.
// See \S\,4.8, and Table~I, Ref[1].
// 
// Offset nodes are produced in pairs, and the present Loop is split
// recursively into a quad (left) and the remainder (right).
//
// A pair of offset Nodes is needed in order to create a 4-Noded Loop.
//
// What happens is determined by the nature of the interior angles at
// each of these Nodes (classification and associated rules per Table I):
//   Case 1: 120 degrees <= angle <= 240 degrees.
//   Case 2:   0 degrees <= angle <  120 degrees.
//   Case 3: 240 degrees <  angle <  360 degrees.
//
// Associated rules for placement of new Nodes are:
//   Case 1: along bisecting line of interior angle, at current ideal length.
//   Case 2: resultant of vectors which point to the two adjoining Nodes.
//   Case 3: same, except use negative of resultant.
//
// What occurs is summarized in the following table.  In cases when only a
// single new Node is produced, an existing Node is used to in making up
// the four that for the new left subLoop.
//   +------+------+---------------------------------------------------------+
//   | curr | CCW  |                        Action                           |
//   | Case | Case |                                                         |
//   +------+------+---------------------------------------------------------+
//   |   1  |   1  |  Use Rule 1 twice to make 2 new Nodes.                  |
//   |   2  |   2  |  SubLoop with no new Nodes (forwards & back).
//   |   3  |   3  |  Use Rule 3 twice to make 2 new Nodes.                  |
//   |   1  |   2  |  Use Rule 2 from 2nd Node to make 1 new Node.           |
//   |   2  |   1  |    Reverse.                                             |
//   |   1  |   3  |  Use Rule 1 and Rule 3, in order, to make 2 new Nodes.  |
//   |   3  |   1  |    Reverse.                                             |
//   |   2  |   3  |  Use Rule 3 to create 1 new Node.                       |
//   |   3  |   2  |    Reverse.                                             |
//   +------+------+---------------------------------------------------------+
// ---------------------------------------------------------------------------
{
  const int N = nodes.getSize ();

  if (N < 6) return;
  
  char routine[] = "Loop::offset";

  register int i;

  int          caseA, caseB;
  real         anglA, anglB;
  Node         *Ni, *Nj, *Np, *Nq, *newNode, *oldNode;
  Point        P, P1, P2;

  // -- Search for first pair of offset Nodes in *this.

  for (i = 0; i < N; i++) {
    Ni = nodes[i];
    Np = nodes[CCW (i, N)];
    if (Ni -> offset() && Np -> offset()) {
      Nj    = nodes [CW (i, N)];
      Nq    = nodes[CCW(CCW(i,N),N)];
      anglA = Ni -> pos().angle (Np -> pos(), Nj -> pos());
      anglB = Np -> pos().angle (Nq -> pos(), Ni -> pos());

      if      (anglA <       TWOPI / 3.0) caseA = 2;
      else if (anglA > 2.0 * TWOPI / 3.0) caseA = 3;
      else                                caseA = 1;

      if      (anglB <       TWOPI / 3.0) caseB = 2;
      else if (anglB > 2.0 * TWOPI / 3.0) caseB = 3;
      else                                caseB = 1;
      
      if (Nj -> interior () && Nq -> interior () || caseA == 2 && caseB == 2) {
	splitline.setSize (0);
	left  = new Loop (nodes, splitline, CW (i, N), CCW (CCW (i, N), N), 1);
	right = new Loop (nodes, splitline, CCW (CCW (i, N), N), CW (i, N), 0);
	right -> offset  ();
	break;
      
      } else if (caseA == 1 && caseB == 1) {
	if (Nj -> interior ()) {
	  splitline.setSize (1);
	  P  = Np -> pos().relative (Nq -> pos(), Np -> prefSize(), 0.5*anglB);
	  newNode = new Node (++Global::nodeIdMax, P, Np -> prefSize ());
	  if (oldNode = Global::exist (newNode)) {
	    delete newNode;
	    Global::nodeIdMax--;
	    newNode = oldNode;
	  } else 
	    Global::nodeList.add (newNode);
	  splitline[0] = newNode;

	  left  = new Loop (nodes, splitline, CW (i, N), CCW (i, N), 1);
	  right = new Loop (nodes, splitline, CCW (i, N), CW (i, N), 0);
	} else {
	  splitline.setSize (2);
	  P1 = Ni -> pos().relative (Np -> pos(), Ni -> prefSize(), 0.5*anglA);
	  P2 = Np -> pos().relative (Nq -> pos(), Np -> prefSize(), 0.5*anglB);

	  newNode = new Node (++Global::nodeIdMax, P2, Np -> prefSize ());
	  if (oldNode = Global::exist (newNode)) {
	    delete newNode;
	    Global::nodeIdMax--;
	    newNode = oldNode;
	  } else 
	    Global::nodeList.add (newNode);
	  splitline[0] = newNode;

	  newNode = new Node (++Global::nodeIdMax, P1, Ni -> prefSize ());
	  if (oldNode = Global::exist (newNode)) {
	    delete newNode;
	    Global::nodeIdMax--;
	    newNode = oldNode;
	  } else 
	    Global::nodeList.add (newNode);
	  splitline[1] = newNode;

	  left  = new Loop (nodes, splitline, i, CCW (i, N), 1);
	  right = new Loop (nodes, splitline, CCW (i, N), i, 0);
	}
	right -> offset  ();
	break;

      } else if (caseA == 3 && caseB == 3) {
	splitline.setSize (2);
	P1 = Np -> pos() - Nq -> pos();
	P2 = Np -> pos() - Ni -> pos();
	P  = Np -> pos() + (P1 + P2);
	newNode = new Node (++Global::nodeIdMax, P, Np -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);
	splitline[0] = newNode;

	P1 = Ni -> pos() - Np -> pos();
	P2 = Ni -> pos() - Nj -> pos();
	P  = Ni -> pos() + (P1 + P2);

	newNode = new Node (++Global::nodeIdMax, P, Ni -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);
	splitline[1] = newNode;

	left  = new Loop (nodes, splitline, i, CCW (i, N), 1);
	right = new Loop (nodes, splitline, CCW (i, N), i, 0);
	right -> offset  ();
	break;

      } else if (caseA == 1 && caseB == 2) {
	splitline.setSize (1);
	P1 = Ni -> pos() - Np -> pos();
	P2 = Nq -> pos() - Np -> pos();
	P  = Np -> pos() + (P1 + P2);

	newNode = new Node (++Global::nodeIdMax, P, Np -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);
	splitline[0] = newNode;

	left  = new Loop (nodes, splitline, i, CCW (CCW (i, N), N), 1);
	right = new Loop (nodes, splitline, CCW (CCW (i, N), N), i, 0);
	right -> offset  ();
	break;

      } else if (caseA == 2 && caseB == 1) {
	splitline.setSize (1);

	if (Nj -> interior()) {
	  P = Np -> pos().relative (Nq -> pos(), Np -> prefSize(), 0.5*anglB);
	} else {
	  P1 = Nj -> pos() - Ni -> pos();
	  P2 = Np -> pos() - Ni -> pos();
	  P  = Ni -> pos() + (P1 + P2);
	}

	newNode = new Node (++Global::nodeIdMax, P, Ni -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);
	splitline[0] = newNode;

	left  = new Loop (nodes, splitline, CW (i, N), CCW (i, N), 1);
	right = new Loop (nodes, splitline, CCW (i, N), CW (i, N), 0);
	right -> offset  ();

	break;

      } else if (caseA == 1 && caseB == 3) {
	splitline.setSize (2);
	P1 = Np -> pos() - Nq -> pos();
	P2 = Np -> pos() - Ni -> pos();
	P  = Np -> pos() + (P1 + P2);
	
	newNode = new Node (++Global::nodeIdMax, P, Np -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);

	splitline[0] = newNode;

	P  = Ni -> pos().relative (Np -> pos(), Ni -> prefSize(), 0.5 * anglA);

	newNode = new Node (++Global::nodeIdMax, P, Ni -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);

	splitline[1] = newNode;

	left  = new Loop (nodes, splitline, i, CCW (i, N), 1);
	right = new Loop (nodes, splitline, CCW (i, N), i, 0);
	right -> offset  ();
	break;

      } else if (caseA == 3 && caseB == 1) {
	splitline.setSize (2);
	P  = Np -> pos().relative (Nq -> pos(), Nq -> prefSize(), 0.5 * anglB);

	newNode = new Node (++Global::nodeIdMax, P, Ni -> prefSize());
 	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);
	splitline[0] = newNode;

	P1 = Ni -> pos() - Nj -> pos();
	P2 = Ni -> pos() - Np -> pos();
	P  = Ni -> pos() + (P1 + P2);

	newNode = new Node (++Global::nodeIdMax, P, Ni -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);
	splitline[1] = newNode;

	left  = new Loop (nodes, splitline, i, CCW (i, N), 1);
	right = new Loop (nodes, splitline, CCW (i, N), i, 0);
	right -> offset  ();
	break;

      } else if (caseA == 2 && caseB == 3) {
	splitline.setSize (1);
	P1 = Np -> pos() - Nq -> pos();
	P2 = Np -> pos() - Ni -> pos();
	P  = Np -> pos() + (P1 + P2);

	newNode = new Node (++Global::nodeIdMax, P, Np -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode);
	splitline[0] = newNode;

	left  = new Loop (nodes, splitline, CW (i, N), CCW (i, N), 1);
	right = new Loop (nodes, splitline, CCW (i, N), CW (i, N), 0);
	right -> offset  ();
	break;

      } else if (caseA == 3 && caseB == 2) {
	splitline.setSize (1);
	P1 = Ni -> pos() - Nj -> pos();
	P2 = Ni -> pos() - Np -> pos();
	P  = Ni -> pos() + (P1 + P2);

	newNode = new Node (++Global::nodeIdMax, P, Ni -> prefSize());
	if (oldNode = Global::exist (newNode)) {
	  delete newNode;
	  Global::nodeIdMax--;
	  newNode = oldNode;
	} else 
	  Global::nodeList.add (newNode); 
	splitline[0] = newNode;

	left  = new Loop (nodes, splitline, i, CCW (CCW (i, N), N), 1);
	right = new Loop (nodes, splitline, CCW (CCW (i, N), N), i, 0);
	right -> offset  ();
	break;

      } else
	message (routine, "classification error", ERROR);
    }
  }
}


void Loop::printMesh (ostream& s) const
// ---------------------------------------------------------------------------
// Print up list of Nodes and 4-noded Loops (elements) on s.
// ---------------------------------------------------------------------------
{
  int                 num;
  Node*               start;
  ListIterator<Node*> n (Global::nodeList);

  s << "Mesh {" << endl;
  s << Global::nodeList.length() << "  Vertices" << endl;
  for (n .reset(); n.more(); n.next())
    s << *(n.current()) << setw (10) << 0.0 << endl;

  num = 0;
  countElements (num);
  s << endl << num << "  Elements" << endl;
  
  num = 0;
  printQuad (s, num);

  // -- Print out the non-interior Nodes in order given, assumed
  //    to provide a continuous CCW boundary to the region.

  s << endl << "1  Boundary" << endl << "1  1  ";
  num = 0;
  for (n.reset(); n.more(); n.next()) {
    num += (n.current() -> interior()) ? 0 : 1;
    if (num == 1) start = n.current();
  }
  s << num + 1;

  for (n.reset(); n.more(); n.next())
    if (!n.current() -> interior()) s << setw (5) << n.current() -> ID();
  s << setw (5) << start -> ID() << endl;

  s << endl << "0  Curves" << endl;
  s << "}" << endl;
}


void Loop::printQuad (ostream& s, int& n) const
// ---------------------------------------------------------------------------
// Recursively print up 4-noded Loops on s.
// ---------------------------------------------------------------------------
{
  const int N = nodes.getSize();
  
  if (N == 4) {
    s << setw (5) << ++n << "\t4";
    for (register int i = 0; i < N; i++)
      s << setw (5) << nodes[i] -> ID();
    s << endl;
  }

  if (left)  left  -> printQuad (s, n);
  if (right) right -> printQuad (s, n);
}


void Loop::countElements (int& sum) const
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  if (nodes.getSize() == 4) ++sum;
  
  if (left)  left  -> countElements (sum);
  if (right) right -> countElements (sum);
}


void Loop::nodeLabel (const int& num  ,
		      char*      label) const
// ---------------------------------------------------------------------------
// Print ID for nodes[i] in label.
// ---------------------------------------------------------------------------
{
  char routine[] = "Loop::nodeLabel";

  if (num >= nodes.getSize ())
    error (routine, "input exceeds number of nodes in Loop", ERROR);

  sprintf (label, "%1d", nodes[num] -> ID());
}


real Loop::lengthScale () const
// ---------------------------------------------------------------------------
// Return a descriptive length scale for loop.
// As a first hack, return hypotenuse of smallest spanning x-y rectangle.
// ---------------------------------------------------------------------------
{
  Point Pmin, Pmax;

  limits (Pmin, Pmax);

  return hypot (Pmax.x - Pmin.x, Pmax.y - Pmin.y);
}


void Loop::limits (Point& Pmin,
		   Point& Pmax) const
// ---------------------------------------------------------------------------
// Get x & y limits for this loop.
// ---------------------------------------------------------------------------
{
  register int i;
  const int    N = nodes.getSize (); 

  real X, Y, xmin, ymin, xmax, ymax;

  xmin = ymin =  1.0e30;
  xmax = ymax = -1.0e30;

  for (i = 0; i < N; i++) {
    X = nodes[i] -> pos () . x;
    Y = nodes[i] -> pos () . y;
    if      (X < xmin) xmin = X;
    else if (X > xmax) xmax = X;
    if      (Y < ymin) ymin = Y;
    else if (Y > ymax) ymax = Y;
  }

  Pmin.x = xmin;  Pmin.y = ymin;
  Pmax.x = xmax;  Pmax.y = ymax;
}
