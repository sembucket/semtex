///////////////////////////////////////////////////////////////////////////////
// mesh.C: read information from a FEML stream, provide
// facilities for generation of mesh knots and initial connectivity.
//
// Example/required parts of a FEML file:
//
// <NODES NUMBER=9>
// #	tag	x	y	z
// 	1	0.0	0.0	0.0
// 	2	2.0	0.0	0.0
// 	3	4.0	0.0	0.0
// 	4	0.0	0.5	0.0
// 	5	2.0	0.5	0.0
// 	6	4.0	0.5	0.0
// 	7	0.0	1.0	0.0
// 	8	2.0	1.0	0.0
// 	9	4.0	1.0	0.0
// </NODES>
// 
// <ELEMENTS NUMBER=4>
// #	tag	type    nodes
// 	1	<Q>	1 2 5 4	 </Q>
// 	2	<Q>	2 3 6 5  </Q>
// 	3	<Q>	4 5 8 7  </Q>
// 	4	<Q>	5 6 9 8  </Q>
// </ELEMENTS>
// 
// <SURFACES NUMBER=8>
// #	tag	elmt	face	type
//	1	1	1	<B>	w	</B>
//	2	2	1	<B>	w	</B>
//	3	2	2	<B>	o	</B>
//	4	4	2	<B>	o	</B>
//	5	4	3	<B>	w	</B>
//	6	3	3	<B>	w	</B>
//	7	3	4	<B>	v	</B>
//	8	1	4	<B>	v	</B>
// </SURFACES>
//
// Optional:
//
// <CURVES NUMBER=1>
// #    tag     elmt    face    specification
//      1       4       3       <arc>      -1.0       </arc>
// </CURVES>
//
// <GROUPS NUMBER=3>
// #	tag	name	descriptor
// 	1	v	value
// 	2	w	wall
// 	3	o	outflow
// </GROUPS>
// 
// <BCS NUMBER=3>
// #	tag	group	number, followed by BCs.
// 	1	v	4
// 			<D>	u = 1.0-4.0*(y-0.5)^2.0	</D>
// 			<D>	v = 0.0			</D>
// 			<D>	w = 0.0			</D>
// 			<H>	p			</H>
// 	2	w	4
// 			<D>	u = 0.0			</D>
// 			<D>	v = 0.0			</D>
// 			<D>	w = 0.0			</D>
// 			<H>	p 			</H>
// 	3	o	4
// 			<N>	u = 0.0			</N>
// 			<N>	v = 0.0			</N>
// 			<N>	w = 0.0			</N>
// 			<D>	p = 0.0			</D>
// </BCS>
// 
// NB: Node, Element and Side IDs are internally held as one less than
// input value, i.e. commence at 0 instead of 1.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include <iostream.h>
#include <fstream.h>
#include <strstream.h>
#include <iomanip.h>

#include <Mesh.h>
#include <Utility.h>
#include <Femlib.h>

#ifdef __DECCXX
  #pragma define_template vector<int>
  #pragma define_template vector<Mesh::Node*>
  #pragma define_template vector<Mesh::Side*>
  #pragma define_template vector<Curve*>
#endif


static inline int rma (int i, int j, int n)
// -- Row-major offsetting for 2D arrays with 0-based indexing.
{ return j + i * n; }


Mesh::Mesh (FEML& f) : feml (f)
// ---------------------------------------------------------------------------
// Create a Mesh using information available in feml.
// ---------------------------------------------------------------------------
{
  char      routine[] = "Mesh::Mesh", err[StrMax], tag[StrMax], nextc;
  int       i, j, k, K, Nn;
  const int verb = (int) Femlib::value ("VERBOSE");
  Node*     N;
  Elmt*     E;

  nodeTable .setSize (0);
  elmtTable .setSize (0);
  curveTable.setSize (0);

  // -- Input Nodes.

  nodeTable.setSize (Nn = feml.attribute ("NODES", "NUMBER"));

  if (Nn < 5) {
    sprintf (err, "At least 4 Nodes are needed, found %1d declared", Nn);
    message (routine, err, ERROR);
  }

  if (verb) cout << "Reading vertices ...";

  for (i = 0; i < Nn; i++) {

    while ((nextc = feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');

    N = new Mesh::Node;
    feml.stream() >> N -> ID >> N -> loc.x >> N -> loc.y >> N -> loc.z;
    N -> ID--;
    N -> gID = UNSET;

    if (N -> ID >= Nn) {
      sprintf (err, "Node ID %1d exceeds attribution (%1d)", N -> ID + 1, Nn);
      message (routine, err, ERROR);
    } else 
      nodeTable[N -> ID] = N;
  }

  if (verb) cout << " done" << endl;

  // -- Input Elmt corner vertex nodes.
  //    Presently, only quad (<Q>) elements are allowed.
  
  elmtTable.setSize (K = feml.attribute ("ELEMENTS", "NUMBER"));

  if (K < 1) {
    sprintf (err, "at least 1 element needed, %1d attributed", K);
    message (routine, err, ERROR);
  }

  if (verb) cout << "Reading elements ...";

  for (i = 0; i < K; i++) {

    while ((nextc = feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');

    E = new Mesh::Elmt;

    feml.stream() >> E -> ID >> tag;
    E -> ID--;

    if (strcmp (tag, "<Q>") == 0) {
      E -> node.setSize (4);
      E -> side.setSize (4);

      for (j = 0; j < 4; j++) {
	feml.stream() >> k;
	k--;
	if (k >= Nn) {
	  sprintf (err, "in element %1d, node tag %1d exceeds maximum (%1d)",
		   E -> ID + 1, k + 1, Nn);
	  message (routine, err, ERROR);
	} else 
	  E -> node[j] = nodeTable[k];
      }
    } else {
      sprintf (err, "unrecognized element tag: %s", tag);
      message (routine, err, ERROR);
    }

    feml.stream() >> tag;
    if (strcmp (tag, "</Q>") != 0) {
      sprintf (err, "closing tag </Q> missing for element %1d", E -> ID + 1);
      message (routine, err, ERROR);
    }

    if (E -> ID >= K) {
      sprintf (err, "element ID (%1d) exceeds attribution (%1d)", E -> ID+1,K);
      message (routine, err, ERROR);
    } else
      elmtTable[E -> ID] = E;
  }

  if (verb) cout << " done" << endl;

  if (verb) cout << "Setting up mesh internal connectivity ...";
  assemble ();
  if (verb) cout << " done" << endl;
  
  if (verb) cout << "Installing mesh external surface data ...";
  surfaces ();
  if (verb) cout << " done" << endl;

  if (verb) cout << "Checking mesh connectivity ...";
  checkAssembly ();
  if (verb) cout << " done" << endl;

  if (verb) cout << "Installing mesh curved sides ...";
  curves ();
  if (verb) cout << " done" << endl;
}


void Mesh::assemble ()
// ---------------------------------------------------------------------------
// Traverse Elmts and fill in Side-based connectivity information.
//
// On exit, Sides that don't have mateElmt set should correspond to Surfaces.
// ---------------------------------------------------------------------------
{
  register int  i, j, r, s, found;
  const    int  Ne = nEl();
  register Elmt *E, *ME;	               // -- M <==> "mate".
  register Side *S, *MS;

  // -- First, build Elmt Sides.

  for (i = 0; i < Ne; i++) {
    E = elmtTable (i);
    const int Nn = E -> nNodes();
    
    for (j = 0; j < Nn; j++) {
      S = new Side (j, E -> node (j), E -> ccwNode (j));
      S -> thisElmt = E;
      E -> side (j) = S;
    }
  }

  // -- Now traverse Elmts and build Side--Side connections.

  for (i = 0; i < Ne; i++) {
    E = elmtTable (i);
    const int Nn = E -> nNodes();

    for (j = 0; j < Nn; j++) {
      S = E -> side (j);
      found = 0;

      for (r = 0; !found && r < Ne; r++) {
	ME = elmtTable (r);
	const int Nm = ME -> nNodes();
	
	for (s = 0; !found && s < Nm; s++) {
	  MS = ME -> side (s);

	  if (found = ( S -> startNode == MS -> endNode &&
		       MS -> startNode ==  S -> endNode )) {
	    S -> mateElmt = ME;
	    S -> mateSide = MS;
	  }
	}
      }
    }
  }
}


void Mesh::surfaces ()
// ---------------------------------------------------------------------------
// Surface information can either declare
//   a boundary group name <B> group </B>
// or set up
//   a periodic boundary   <P> elmt side </P>.
//
// Periodic boundaries can be set only once.
// ---------------------------------------------------------------------------
{
  char      routine[] = "Mesh::surfaces", err[StrMax], tag[StrMax], nextc;
  int       i, e, s, t;
  const int K = feml.attribute ("SURFACES", "NUMBER");

  for (i = 0; i < K; i++) {

    while ((nextc = feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');

    // -- Get element and side number information.

    feml.stream() >> t >> e >> s;
    if (t > K) {
      sprintf (err, "Surface tag no. %1d exceeds attribution (%1d)",
	       t, K);
      message (routine, err, ERROR);
    } else if (e > nEl()) {
      sprintf (err, "Surface %1d element no. %1d too large (%1d)",
	       t, e, nEl());
      message (routine, err, ERROR);
    } else if (s > elmtTable (e - 1) -> nNodes()) {
      sprintf (err, "Surface %1d elmt %1d side no. %1d too large (%1d)",
	       t, e, s, elmtTable (e - 1) -> nNodes());
      message (routine, err, ERROR);
    } else if (elmtTable (e - 1) -> side (s - 1) -> mateElmt) {
      Mesh::Elmt* ME = elmtTable (e - 1) -> side (s - 1) -> mateElmt;
      Mesh::Side* MS = elmtTable (e - 1) -> side (s - 1) -> mateSide;
      sprintf (err, "Surface %1d elmt %1d side %1d already set to mate "
	       "elmt %1d side %1d", t, e, s, ME -> ID + 1, MS -> ID + 1);
      message (routine, err, ERROR);
    } 
    
    // -- Set up either a boundary group or a periodic boundary.

    feml.stream() >> tag;
    e--; s--;

    if (strcmp (tag, "<B>") == 0) {
      
      // -- Boundary group.
      //    Group information for this side should not be set already.

      if (elmtTable (e) -> side (s) -> group) {
	sprintf (err, "Surface %1d: group already set (%c)",
		 t, elmtTable (e) -> side (s) -> group);
	message (routine, err, ERROR);
      }
      
      feml.stream() >> elmtTable (e) -> side (s) -> group;

      // -- Clean up.

      feml.stream() >> tag;
      if (strcmp (tag, "</B>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <B> with %s", t, tag);
	message (routine, err, ERROR);
      }
      
    } else if (strcmp (tag, "<P>") == 0) {
      
      // -- Periodic.
      //    These are set two at a time (this and indicated mate).
      //    Mate information for both this side and its indicated mate
      //    should not be previously set.

      int me, ms;
      feml.stream() >> me >> ms;
      
      if (me < 1 || me > nEl()) {
	sprintf (err, "Surface %1d, mating elmt no. %1d out of range (1--%1d)",
		 t, me, nEl());
	message (routine, err, ERROR);
      } else if (ms < 1 || ms > elmtTable (e) -> nNodes()) {
	sprintf (err, "Surface %1d, mating side no. %1d out of range (1--%1d)",
		 t, me, elmtTable (e) -> nNodes());
	message (routine, err, ERROR);
      } else if (elmtTable (me - 1) -> side (ms - 1) -> mateElmt ||
		 elmtTable (me - 1) -> side (ms - 1) -> group    ) {
	sprintf (err, "Surface %1d, mating elmt %1d, side %1d already set",
		 t, me, ms);
	message (routine, err, ERROR);
      }

      me--; ms--;

      elmtTable (e)  -> side (s)  -> mateElmt = elmtTable (me);
      elmtTable (e)  -> side (s)  -> mateSide = elmtTable (me) -> side (ms);
      elmtTable (me) -> side (ms) -> mateElmt = elmtTable (e);
      elmtTable (me) -> side (ms) -> mateSide = elmtTable (e)  -> side (s);
      
      // -- Clean up.
      
      feml.stream() >> tag;
      if (strcmp (tag, "</P>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <P> with %s", t, tag);
	message (routine, err, ERROR);
      }
      
    } else {
      sprintf (err, "couldn't recognize Surface tag %s", tag);
      message (routine, err, ERROR);
    }
  }
}


void Mesh::showAssembly (Mesh& m)
// ---------------------------------------------------------------------------
// Static debugging function: Print edge--edge connectivity information.
// ---------------------------------------------------------------------------
{
  int   i, j, Ne, Nn;
  Elmt* E;
  Side* S;

  Ne = m.nEl();
  cout << "# " << Ne << " Elmts" << endl;

  for (i = 0; i < Ne; i++) {
    E  = m.elmtTable (i);
    Nn = E -> nNodes ();

    cout << "# Elmt: " << E -> ID + 1 << ", Vertices:";
    for (j = 0; j < Nn; j++)
      cout << " " << E -> node (j) -> ID + 1;

    cout << ", Mating:";
    for (j = 0; j < Nn; j++) {
      S = E -> side (j);
      cout << '\t' << S -> ID + 1;
      cout << "->";
      if   (!S -> mateElmt) cout << S -> group;
      else                  cout << S -> mateElmt -> ID + 1 << "."
				 << S -> mateSide -> ID + 1;
    }
    cout << endl;
  }
}


void Mesh::checkAssembly ()
// ---------------------------------------------------------------------------
// All element sides have to either mate an adjoining element or fall on
// a boundary.  Check it out.  But surfaces() must have been called first.
// ---------------------------------------------------------------------------
{
  char     routine[] = "Mesh::checkAssembly", err[StrMax];
  Elmt*    E;
  Side*    S;
  register int i, j;
  const    int Ne = nEl();

  for (i = 0; i < Ne; i++) {
    E = elmtTable (i);
    const int Ns = E -> nNodes();

    for (j = 0; j < Ns; j++) {
      S = E -> side (j);
      if (S -> mateSide == 0) {
	sprintf (err, "Elmt %1d Side %1d not set",
		 S -> thisElmt -> ID + 1, S -> ID + 1);
	message (routine, err, ERROR);
      }
    }
  }
  
  if (Femlib::value ("VERBOSE") > 1) {
    cout << "Mesh connectivity summary:" << endl;
    showAssembly (*this);
  }
}


void Mesh::curves ()
// ---------------------------------------------------------------------------
// Read in curved edge information and store in Mesh curveTable.
//
// Curved edges are specified by lines like:
//   curveID  elementID  sideID <C> ... </C>
//
// Tags <C> and </C> delimit the kind of curve to be defined; presently
// the only defined kind is <ARC>.  Between the delimiters the number
// of parameters is user-defined.
// ---------------------------------------------------------------------------
{
  if (!feml.seek ("CURVES")) return;
  
  char   routine[] = "Mesh::curves";
  char   err[StrMax], buf[StrMax];
  int    i, K, id, elmt, side, ns;
  Curve* C;
  Elmt*  E;
  Side*  S;

  curveTable.setSize (K = feml.attribute ("CURVES", "NUMBER"));

  for (i = 0; i < K; i++) {
    feml.stream() >> id >> elmt >> side;

    if (id > K) {
      sprintf (err, "Curve ID %1d exceeds attribution (%1d)", id, K);
      message (routine, err, ERROR);
    } else if (elmt > nEl()) {
      sprintf (err, "Curve ID %1d, Elmt no. %1d too large (%1d)", id, elmt, K);
      message (routine, err, ERROR);
    } else if (side > (ns = elmtTable (elmt - 1) -> nNodes())) {
      sprintf (err, "Curve ID %1d, Side no. %1d too large (%1d)", id,side, ns);
      message (routine, err, ERROR);
    }
    
    S = elmtTable (elmt - 1) -> side (side - 1);

    feml.stream() >> buf;

    if (strcmp (buf, "<ARC>") == 0) {
      real radius;
      feml.stream() >> radius;

      C = new CircularArc (id, S, radius);

      feml.stream() >> buf;
      if (strcmp (buf, "</ARC>") != 0) {
	sprintf (err, "Curve ID %1d, couldn't close <ARC> with </ARC>", id);
	message (routine, err, ERROR);
      }
    } else {
      sprintf (err, "Curve %1d, unknown curve kind %s", buf);
      message (routine, err, ERROR);
    }

    curveTable (i) = C;
  }
}


CircularArc::CircularArc (const int   id,
			  Mesh::Side* S ,
			  const real  R )
// ---------------------------------------------------------------------------
// Constructor for CircularArc.  R is the radius of arc, and its sign
// specifies the convexity of the element edge.
//
// R +ve ==> arc increases area enclosed by element (cf straight line),
// R -ve ==> arc decreases area enclosed by element.
// ---------------------------------------------------------------------------
{
  char err[StrMax], routine[] = "CircularArc::CircularArc";

  curveSide = S;

  convexity = (R < 0.0) ? -1 : 1;
  radius    = fabs (R);
  Point P1  = curveSide -> startNode -> loc;
  Point P2  = curveSide -> endNode   -> loc;
  Point unitNormal, link, midpoint, centroid = {0.0, 0.0, 0.0};
  real  dx, dy, l, sign = 0.0;

  midpoint.x   = 0.5 * (P2.x + P1.x);
  midpoint.y   = 0.5 * (P2.y + P1.y);
  dx           =        P2.x - P1.x;
  dy           =        P2.y - P1.y;
  l            = hypot (dx, dy);
  unitNormal.x = -dy / l;
  unitNormal.y =  dx / l;

  if (2.0 * radius < l) {
    sprintf (err, "curve %1d:\narc, radius %f, can't span nodes %1d & %1d",
	     id, radius, 
	     curveSide -> startNode -> ID + 1, curveSide -> endNode -> ID + 1);
    message (routine, err, ERROR);
  } else
    semiangle = asin (0.5*l / radius);

  centroid = curveSide -> thisElmt -> centroid ();
  
  link.x = centroid.x - midpoint.x;
  link.y = centroid.y - midpoint.y;

  // -- Sign +1 if centre lies in direction of centroid from side midpoint.

  sign = link.x * unitNormal.x + link.y * unitNormal.y;
  sign = convexity * sign / fabs (sign);

  centre.x = midpoint.x + sign * cos (semiangle) * radius * unitNormal.x;
  centre.y = midpoint.y + sign * cos (semiangle) * radius * unitNormal.y;
}


void CircularArc::compute (const int   np     ,
			   const real* spacing,
			   Point*      knot   ) const
// ---------------------------------------------------------------------------
// Distribute np knots along arc according to spacing on -1, 1.
// ---------------------------------------------------------------------------
{
  Point     P1 = curveSide -> startNode -> loc;
  Point     P2 = curveSide -> endNode   -> loc;
  real      theta1, theta2, dtheta, phi;
  const int nm = np - 1;

  theta1 = atan2 (P1.y - centre.y, P1.x - centre.x);
  theta2 = atan2 (P2.y - centre.y, P2.x - centre.x);
  dtheta = theta2 - theta1;

  if (fabs (dtheta) > 2.0*semiangle + EPSSP)
    dtheta += (dtheta < 0.0) ? TWOPI : -TWOPI;

  knot[ 0].x = P1.x;  knot[ 0].y  = P1.y;
  knot[nm].x = P2.x;  knot[nm].y = P2.y;

  for (int i(1); i < nm; i++) {
    phi = theta1 + dtheta * 0.5 * (spacing[i] + 1.0);
    knot[i].x = centre.x + radius * cos (phi);
    knot[i].y = centre.y + radius * sin (phi);
  }
}


Point Mesh::Elmt::centroid () const
// ---------------------------------------------------------------------------
// Return point that is centroid of element Node points.
// ---------------------------------------------------------------------------
{
  register int i;
  const    int K = nNodes();
  Point        C = {0.0, 0.0, 0.0};

  for (i = 0; i < K; i++) {
    Point P = node (i) -> loc;
    C.x += P.x;
    C.y += P.y;
  }

  C.x /= K;
  C.y /= K;

  return C;
}


void Mesh::meshSide (int         np     ,
		     int         elmt   ,
		     int         side   ,
		     const real* spacing,
		     Point*      knot   ) const
// ---------------------------------------------------------------------------
// If a curved side can be identified for the nominated element and side,
// compute the points using appropriate routine.  Otherwise compute points
// along a straight side.
//
// Spacing gives location of knots in master coordinates [-1, 1].
// ---------------------------------------------------------------------------
{
  char      routine[] = "Mesh::meshSide";
  char      i, j;
  const int Nc = curveTable.getSize();
  const int Ne = elmtTable .getSize();

  if (np < 2) message (routine, "must have at least two points", ERROR);

  for (i = 0; i < Nc; i++) {
    Curve* C = curveTable (i);
    if (C -> ismatch (elmt, side)) {
      C -> compute (np, spacing, knot);
      return;
    }
  }

  // -- Fall though default: straight line.

  const Side* S  = elmtTable (elmt) -> side (side);
  const Point P1 = S -> startNode -> loc;
  const Point P2 = S -> endNode   -> loc;
  const real  dx = P2.x - P1.x;
  const real  dy = P2.y - P1.y;

  for (i = 0; i < np; i++) {
    knot[i].x = P1.x + dx * 0.5 * (spacing[i] + 1.0);
    knot[i].y = P1.y + dy * 0.5 * (spacing[i] + 1.0);
  }
  return;
}


void Mesh::meshElmt (const int   ID,
		     const int   np,
		     const real* z ,
		     real*       x ,
		     real*       y ) const
// ---------------------------------------------------------------------------
// Generate mesh points for Elmt No ID (IDs begin at 1).
// Generate element-edge points, then internal points using a Coons patch.
// Input z contains the spacing of edge knot points along interval [-1, 1].
//
// For a quad mesh, equal-order on each side, x & y have row-major ordering.
// ---------------------------------------------------------------------------
{
  char     routine[] = "Mesh::meshElmt";
  register int  i, j;
  const    int  nm = np - 1;
  const    int  ns = elmtTable (ID) -> nNodes();
  vector<Point> P (np);

  // -- Compute and load peripheral points.

  for (j = 0; j < ns; j++) {
    meshSide (np, ID, j, z, P());
    for (i = 0; i < nm; i++) {
      switch (j) {
      case 0:
	x[rma (     0,      i, np)] = P[i].x;
	y[rma (     0,      i, np)] = P[i].y;
	break;
      case 1:
	x[rma (     i,     nm, np)] = P[i].x;
	y[rma (     i,     nm, np)] = P[i].y;
	break;
      case 2:
	x[rma (    nm, nm - i, np)] = P[i].x;
	y[rma (    nm, nm - i, np)] = P[i].y;
	break;
      case 3:
	x[rma (nm - i,      0, np)] = P[i].x;
	y[rma (nm - i,      0, np)] = P[i].y;
	break;
      default:
	message (routine, "never happen", ERROR);
	break;
      }
    }
  }

  // -- Coons patch on (-1, 1) X (-1, 1) to make interior points.

  for (i = 1; i < nm; i++)
    for (j = 1; j < nm; j++) {

      x[rma (i, j, np)] = 0.50 * ( (1.0 - z[j]) * x[rma ( i,  0, np)] +
				   (1.0 + z[j]) * x[rma ( i, nm, np)] +
				   (1.0 - z[i]) * x[rma ( 0,  j, np)] +
				   (1.0 + z[i]) * x[rma (nm,  j, np)] )
	     
	 - 0.25 * ( (1.0 - z[j]) * (1.0 - z[i]) * x[rma ( 0,  0, np)] +
	            (1.0 - z[j]) * (1.0 + z[i]) * x[rma (nm,  0, np)] +
		    (1.0 - z[i]) * (1.0 + z[j]) * x[rma ( 0, nm, np)] +
		    (1.0 + z[i]) * (1.0 + z[j]) * x[rma (nm, nm, np)] );

      y[rma (i, j, np)] = 0.50 * ( (1.0 - z[j]) * y[rma ( i,  0, np)] +
				   (1.0 + z[j]) * y[rma ( i, nm, np)] +
				   (1.0 - z[i]) * y[rma ( 0,  j, np)] +
				   (1.0 + z[i]) * y[rma (nm,  j, np)] )
	     
	 - 0.25 * ( (1.0 - z[j]) * (1.0 - z[i]) * y[rma ( 0,  0, np)] +
                    (1.0 - z[j]) * (1.0 + z[i]) * y[rma (nm,  0, np)] +
		    (1.0 - z[i]) * (1.0 + z[j]) * y[rma ( 0, nm, np)] +
		    (1.0 + z[i]) * (1.0 + z[j]) * y[rma (nm, nm, np)] );
    }
}


int Mesh::globalID (const int np  ,
		    int*      btog)
// ---------------------------------------------------------------------------
// Generate connectivity (i.e. global knot numbers) for a mesh with np
// knot points (i.e. Lagrange knots) along each element side, ignoring
// internal points (i.e. generate connectivity for static-condensation form).
//
// Fill btog (element-by-element storage of these global numbers) for whole
// mesh: for a mesh of quad elements, btog must hold 4*(np-1)*nEl integers.
// Return the number of global knots (maximum global knot number + 1). 
//
// NB: np >= 2, also global numbers generated here start at 0.
// NB: this connectivity information is generated without reference to BCs.
// ---------------------------------------------------------------------------
{
  char routine[] = "Mesh::globalID";

  if (np < 2) message (routine, "need at least 2 knots", ERROR);

  // -- Create element-side based gID storage, if required, & initialize gIDs.
  
  register int i, j, k, ns;
  const    int nel = nEl(), ni = np - 2;
  int          nGid = 0, nb = 0;
  Elmt*        E;
  Side*        S;

  // -- Allocate space, unset all knot numbers.

  for (i = 0; i < nel; i++) {
    E  = elmtTable (i);
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side (j);
      S -> gID.setSize (ni);      
      S -> startNode -> gID = UNSET;
      S -> endNode   -> gID = UNSET;
      if (ni)      S -> gID = UNSET;
    }
  }

  // -- Generate connectivity information.

  for (i = 0; i < nel; i++) {
    E  = elmtTable (i);
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side (j);
      S -> connect (ni, nGid);
    }
  }

  // -- Fill btog.

  for (i = 0; i < nel; i++) {
    E  = elmtTable (i);
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side (j);
      btog[nb++] = S -> startNode -> gID;
      for (k = 0; k < ni; k++)
	btog[nb++] = S -> gID (k);
    }
  }

  // -- Deallocate internal knot number storage.

  for (i = 0; i < nel; i++) {
    E  = elmtTable (i);
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side (j);
      S -> gID.setSize (0);      
      S -> startNode -> gID = UNSET;
      S -> endNode   -> gID = UNSET;
    }
  }

  return nGid;
}


void Mesh::Side::connect (const int ni ,
			  int&      gid)
// ---------------------------------------------------------------------------
// Fill in connectivity for this element side, updating global number gid.
// ---------------------------------------------------------------------------
{
  register int   i, k;
  register Side* otherSide;

  if (startNode -> gID == UNSET) startNode -> gID = gid++;

  if (ni) {			// -- Do side-internal gids.
    if (mateElmt) {
      otherSide = mateSide;
      if (otherSide -> gID[0] == UNSET)
	for (i = 0; i < ni; i++)
	  gID[i] = gid++;
      else
	for (i = 0, k = ni - 1; i < ni; i++, k--)
	  gID[i] = otherSide -> gID[k];
    } else
      for (i = 0; i < ni; i++)
	gID[i] = gid++;
  }

  if (endNode -> gID == UNSET) endNode -> gID = gid++;
}


void Mesh::showGlobalID (Mesh& m)
// ---------------------------------------------------------------------------
// Print knot connectivity information (global node numbers).
// As things are set up now, this function is of no use since gIDs are
// always UNSET except during calls to Mesh::globalID.
// ---------------------------------------------------------------------------
{
#if 0
  register int i, j, k;
  const    int nel = m.nEl();
  int      ni, ns;
  Elmt*    E;
  Side*    S;

  cout << "# " << nel << " NEL" << endl;
  cout << "# 0 BANDWIDTH" << endl;

  for (i = 0; i < nel; i++) {
    E  = m.elmtTable (i);
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side (j);

      cout << setw (5) << S -> thisElmt -> ID + 1;
      cout << setw (5) << S -> ID + 1;
      cout << setw (5) << 1;
      cout << setw (5) << S -> startNode -> gID;
      cout << endl;

      ni = S -> gID.getSize ();

      if (ni)
	for (k = 0; k < ni; k++) {
	  cout << setw (5) << S -> thisElmt -> ID + 1;
	  cout << setw (5) << S -> ID + 1;
	  cout << setw (5) << k + 2;
	  cout << setw (5) << S -> gID[k];
	  cout << endl;
	}

      cout << setw (5) << S -> thisElmt -> ID + 1;
      cout << setw (5) << S -> ID + 1;
      cout << setw (5) << ni + 2;
      cout << setw (5) << S -> endNode -> gID;
      cout << endl;
    }
  }
#endif
}


void CircularArc::printNek () const
// ---------------------------------------------------------------------------
// Print out information in NEKTON format.
// ---------------------------------------------------------------------------
{
  cout << setw (2)  << curveSide -> ID + 1
       << setw (5)  << curveSide -> thisElmt -> ID + 1
       << setw (14) << 1.0*convexity*radius
       << setw (14) << 0.0
       << setw (14) << 0.0
       << setw (14) << 0.0
       << setw (14) << 0.0
       << " C"
       << endl; 
}


void Mesh::printNek () const
// ---------------------------------------------------------------------------
// Print out mesh information in NEKTON format.
// ---------------------------------------------------------------------------
{
  char       routine[] = "Mesh::printNek";
  char       err [StrMax], buf[StrMax];
  ostrstream os  (err, StrMax);

  int        i, j, ns, nel = nEl();
  Elmt       *E, *ME;
  Side       *S;

  cout.precision (5);
  cout.setf      (ios::scientific,ios::floatfield);

  // -- Elements.

  cout << setw(14) << 10.0
       << setw(14) << 10.0
       << setw(14) << 0.0 
       << setw(14) << 0.0
       << " XFAC,YFAC,XZERO,YZERO"
       << endl;

  cout << "**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y."
       << endl;

  cout << setw(10) << nel
       << setw(10) << 2 
       << setw(10) << nel
       << " NEL,NDIM,NELV"
       << endl;

  for (i = 0; i < nel; i++) {
    E = elmtTable (i);

    cout << "ELEMENT   "
         << setw(10) << E -> ID + 1
         << " [  1A]  GROUP 0"
         << endl;

    ns = E -> nNodes();

    for (j = 0; j < ns; j++)
      cout << setw(14) << E -> node (j) -> loc.x;
    cout << endl;

    for (j = 0; j < ns; j++)
      cout << setw(14) << E -> node (j) -> loc.y;
    cout << endl;
  }

  // -- Curved sides.

  ns = curveTable.getSize();
  cout << "***** CURVED SIDE DATA *****" << endl;
  cout << setw(5) << ns
       << " Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE"
       << endl;
  for (i = 0; i < ns; i++)
    curveTable (i) -> printNek ();

  // -- Boundary conditions.

  cout << "***** BOUNDARY CONDITIONS *****" << endl;
  cout << "***** FLUID BOUNDARY CONDITIONS *****" << endl;

  for (i = 0; i < nel; i++) {
    E  = elmtTable (i);
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S  = E -> side (j);
      ME = S -> mateElmt;
      if (ME) {
	cout << "E  "
	     << setw (5)  << E -> ID + 1
	     << setw (3)  << S -> ID + 1
	     << setw (14) << 1.0*ME -> ID + 1
	     << setw (14) << 1.0*S  -> mateSide -> ID + 1
	     << setw (14) << 1.0
	     << endl;
      } else {
	describeGrp (S -> group, buf);
	if (strstr (buf, "value")) {
	  describeBC (S -> group, 'u', buf);
	  describeBC (S -> group, 'v', err);
	  cout << "V  "
	       << setw (5) << E -> ID + 1
	       << setw (3) << S -> ID + 1
	       << " " << buf
	       << " " << err
	       << setw (14) << 0.0
	       << endl;
	} else if  (strstr (buf, "wall")) {
	  cout << "W  "
	       << setw (5)  << E -> ID + 1
	       << setw (3)  << S -> ID + 1
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << endl;
	} else if (strstr (buf, "outflow")) {
	  cout << "O  "
	       << setw (5)  << E -> ID + 1
	       << setw (3)  << S -> ID + 1
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << endl;
	} else {
	  os << "Elmt " << E -> ID + 1 << " side " << S -> ID + 1
	     << " --- B.C. type "  << buf
	     << " not implemented" << ends;
	  message (routine, err, ERROR);
	}
      }
    }
  }
}


void Mesh::describeGrp (char  G,
			char* S) const
// ---------------------------------------------------------------------------
// Search feml file info for string descriptor matching G, load into S.
// ---------------------------------------------------------------------------
{
  char      routine[] = "Mesh::describeGrp";
  char      nextc, groupc, err[StrMax], buf[StrMax];
  int       i, id, found = 0;
  const int N = feml.attribute ("GROUPS", "NUMBER");
  
  for (i = 0; !found && i < N; i++) {
    while ((nextc = feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');
    feml.stream() >> id >> groupc >> buf;
    if (found = (groupc == G)) strcpy (S, buf);
  }

  if (!found) {
    sprintf (err, "no group found to match '%c'", G);
    message (routine, err, ERROR);
  }
}


void Mesh::describeBC (char  G,
		       char  F, 
		       char* S) const
// ---------------------------------------------------------------------------
// Find the BC description string matching Group G for Field F, load into S.
// ---------------------------------------------------------------------------
{
  char      routine[] = "Mesh::describeBC";
  char      eql, groupc, fieldc, nextc, err[StrMax], buf[StrMax];
  int       i, j, id, nbcs, found = 0;
  const int N = feml.attribute ("BCS", "NUMBER");

  for (i = 0; !found && i < N; i++) {

    while ((nextc = feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');

    feml.stream() >> id >> groupc >> nbcs;

    for (j = 0; !found && j < nbcs; j++) {

      // -- Trash tag. (i.e. ignore type of BC. This is an error!).

      feml.stream() >> buf;

      // -- Check for match and take appropriate action if true.

      feml.stream() >> fieldc;

      if (found = (groupc == G) && (fieldc == F)) {
	feml.stream() >> eql;
	if (eql == '=') {
	  S[0] = F;
	  S[1] = '\0';
	  strcat (S, " = ");
	  feml.stream() >> buf;
	  strcat (S, buf);
	} else {
	  sprintf (err, "Group '%c', Field '%c', expected '=', got '%c",
		   G, F, eql);
	  message (routine, err, ERROR);
	}
      } else {
	feml.stream().ignore (StrMax, '\n');
      }
    }
  }
      
  if (!found) {
    sprintf (err, "couldn't find BC to match Group '%c', Field '%c'", G, F);
    message (routine, err, ERROR);
  }
}
