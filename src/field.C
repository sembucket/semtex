/*****************************************************************************
 * field.C:  routines to deal with Fields.
 *****************************************************************************/

static char RCSid[] = "$Id$";

#include "Fem.h"


static int  nOrder (const void *a, const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare first element (global node number) of two arrays.
// ---------------------------------------------------------------------------
{ return ((int *)a)[0] - ((int *)b)[0]; }


static int  sOrder (const void *a, const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare second element (solve mask) of two arrays.
// ---------------------------------------------------------------------------
{ return ((int *)b)[1] - ((int *)a)[1]; }


Field::Field () :
// ---------------------------------------------------------------------------
// Default constructor.
// ---------------------------------------------------------------------------
field_name      (0),
n_data          (0),
n_elmt_bnodes   (0),
n_gid           (0),
n_solve         (0),
data            (0),
elmt_bndry_gid  (0),
elmt_bndry_mask (0),
n_cons          (0),
n_pack          (0),
n_band          (0),
Hp              (0),
Hc              (0),
geometry_economized (0),
matrices_economized (0)     
{ }


Field::Field (const Mesh& M, int np)
// ---------------------------------------------------------------------------
// Use M to allocate storage for Field and its Elements.  Name Field 'u'.
//
// Compute Mesh knot locations (xmesh & ymesh), but not geometric partials.
// ---------------------------------------------------------------------------
{
  Element* E;

  n_data = n_mesh = n_elmt_bnodes = 0;

  for (ElmtsOfMesh e(M); e.more(); e.next()) {
    E = new Element  (e.ID(), np, e.nSides());
    element_list.add (E);
    n_data        += E -> nTot ();
    n_mesh        += E -> nMsh ();
    n_elmt_bnodes += E -> nExt ();
  }

  data            = rvector (n_data);
  mesh            = rvector (n_mesh);
  elmt_bndry_gid  = ivector (n_elmt_bnodes);
  elmt_bndry_mask = ivector (n_elmt_bnodes);

  int    offset = 0;
  real*  m      = mesh;
  int*   g      = elmt_bndry_gid;
  int*   s      = elmt_bndry_mask;

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E = k.current ();
    E -> install  (offset, m, g, s);
    E -> mesh     (M);
    offset += E -> nTot ();
    m      += E -> nMsh ();
    g      += E -> nExt ();
    s      += E -> nExt ();
  }

  field_name = 'u';
}

 
Field::Field (const Field& f, const Mesh& M, char tag)
// ---------------------------------------------------------------------------
// Copy constructor.
//
// New Field has all same internal data and pointers of old Field, *except*
//
// 1.  tag != 0.  Make a new, independent, Field. New storage is allocated
//     for data, map & mask areas, new elements and boundaries are linked.
// 
// 2.  tag == 0.  Auxillary.  New storage is allocated for data area only,
//     and a new element list is also built, with pointers into new storage.
//     NB: in this case, boundary list is empty: Field has geometry but no
//     boundary conditions.  This kind of Field is intended as auxillary
//     storage, for forcing terms.
//
// In both cases, new storage areas are intialized from input Field. 
// ---------------------------------------------------------------------------
{
  Element*  E;
  real*     d;
  int*      g;
  int*      s;

  memcpy (this, &f, sizeof (Field));

  element_list.clear();
  boundary_list.clear();

  field_name = tag;
  d = data   = rvector (n_data);
  memcpy (data, f.data, n_data * sizeof (real));

  if (tag) {
    elmt_bndry_gid  = ivector (n_elmt_bnodes);
    elmt_bndry_mask = ivector (n_elmt_bnodes);
    memcpy (elmt_bndry_gid,  f.elmt_bndry_gid,  n_elmt_bnodes * sizeof (int)); 
    memcpy (elmt_bndry_mask, f.elmt_bndry_mask, n_elmt_bnodes * sizeof (int));

    g = elmt_bndry_gid;
    s = elmt_bndry_mask;
    
    int offset = 0;
    for (ListIterator<Element*> k(f.element_list); k.more(); k.next()) {
      E = new Element (*k.current ());
      E -> install (offset, 0, g, s);
      element_list.add (E);
      offset += E -> nTot();
      g      += E -> nExt();
      s      += E -> nExt();
    }
    
    if   (tag == PRESSURE) buildBoundaries (f);
    else                   buildBoundaries (M);

  } else {
    int offset = 0;
    for (ListIterator<Element*> k(f.element_list); k.more(); k.next()) {
      E = new Element (*k.current ());
      E -> install (offset, 0, 0, 0);
      element_list.add (E);
      offset += E -> nTot();
    }
  }
}


Field&  Field::operator = (real val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (n_data,      data, 1);
  else              Veclib::fill (n_data, val, data, 1);

  return *this;
}


Field&  Field::operator += (real val)
// ---------------------------------------------------------------------------
// Add val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (n_data, val, data, 1, data, 1);

  return *this;
}


Field&  Field::operator -= (real val)
// ---------------------------------------------------------------------------
// Add -val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::ssub (n_data, val, data, 1, data, 1);

  return *this;
}


Field&  Field::operator *= (real val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  if (val == 0.0) Veclib::zero (n_data, data, 1);
  else            Blas  ::scal (n_data, val, data, 1);

  return *this;
}


Field&  Field::operator /= (real val)
// ---------------------------------------------------------------------------
// Divide field storage area by val.
// ---------------------------------------------------------------------------
{
  if (val == 0.0) message ("Field::op /= (real)", "divide by zero", ERROR);
  else            Blas::scal (n_data, 1.0/val, data, 1);

  return *this;
}


Field&  Field::operator = (const Field& f)
// ---------------------------------------------------------------------------
// This presumes the two fields conform, and copies f's value storage to LHS.
// ---------------------------------------------------------------------------
{
  if (f.n_data != n_data)
    message ("Field::operator = (const Field&)", "size mismatch", ERROR);
  else
    Veclib::copy (n_data, f.data, 1, data, 1);
  
  return *this;
}


Field& Field::operator += (const Field& f)
// ---------------------------------------------------------------------------
// Add f's value to this Field's.
// ---------------------------------------------------------------------------
{
  if (f.n_data != n_data)
    message ("Field::operator += (const Field&)", "size mismatch", ERROR);
  else
    Veclib::vadd (n_data, data, 1, f.data, 1, data, 1);

  return *this;
}


Field&  Field::operator -= (const Field& f)
// ---------------------------------------------------------------------------
// Subtract f's value from this Field's.
// ---------------------------------------------------------------------------
{
  if (f.n_data != n_data)
    message ("Field::operator -= (const Field*)", "size mismatch", ERROR);
  else
    Veclib::vsub (n_data, data, 1, f.data, 1, data, 1);

  return *this;
}


Field&  Field::operator *= (const Field& f)
// ---------------------------------------------------------------------------
// Multiply this Field's value by f's value.
// ---------------------------------------------------------------------------
{
  if (f.n_data != n_data)
    message ("Field::operator *= (const Field&)", "size mismatch", ERROR);
  else
    Veclib::vmul (n_data, data, 1, f.data, 1, data, 1);

  return *this;
}


Field&  Field::operator /= (const Field& f)
// ---------------------------------------------------------------------------
// Divide this Field's value by f's value.
// ---------------------------------------------------------------------------
{
  if (f.n_data != n_data)
    message ("Field::operator /= (const Field&)", "size mismatch", ERROR);
  else
    Veclib::vdiv (n_data, data, 1, f.data, 1, data, 1);

  return *this;
}


Field&  Field:: operator = (const char* function)
// ---------------------------------------------------------------------------
// Set Field's value to temporo-spatially varying function.
// ---------------------------------------------------------------------------
{
  Element*  E;
  int       offset = 0;

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E      = k.current ();
    offset = E -> nOff ();
    E -> evaluate (function, data + offset);

  }

  return *this;
}


Field&  Field::prod (const Field& a, const Field& b)
// ---------------------------------------------------------------------------
// Set this Field's value as the product of a & b.
// ---------------------------------------------------------------------------
{
  Veclib::vmul (n_data, a.data, 1, b.data, 1, data, 1);

  return *this;
}


Field&  Field::sum (const Field& a, const Field& b)
// ---------------------------------------------------------------------------
// Set this Field's value as the sum of a & b.
// ---------------------------------------------------------------------------
{
  Veclib::vadd (n_data, a.data, 1, b.data, 1, data, 1);

  return *this;
}


Field&  Field::addprod (const Field& a, const Field& b)
// ---------------------------------------------------------------------------
// Add the product of a & b to this Field.
// ---------------------------------------------------------------------------
{
  Veclib::vvtvp (n_data, a.data, 1, b.data, 1, data, 1, data, 1);

  return *this;
}


Field&  Field::axpy (real alpha, const Field& x)
// ---------------------------------------------------------------------------
// Add alpha * x to this Field.
// ---------------------------------------------------------------------------
{
  Blas::axpy (n_data, alpha, x.data, 1, data, 1);

  return *this;
}


void  Field::printMesh (Field* F)
// ---------------------------------------------------------------------------
// Static member function.
//
// Mesh location information is written out element-by-element.
// ---------------------------------------------------------------------------
{
  cout << F -> element_list.first() -> nKnot() << " "
       << F -> element_list.first() -> nKnot() << " 1 " 
       << F -> nEl() << " NR NS NZ NEL" << endl;

  for (ListIterator<Element*> k(F -> element_list); k.more(); k.next())
    k.current() -> printMesh ();
}


void  Field::mapElements ()
// ---------------------------------------------------------------------------
// Carry out computation of mesh geometric factors for all elements.
// ---------------------------------------------------------------------------
{
  Element* E;

  for (ListIterator<Element*> i(element_list); i.more(); i.next()) {
    E = i.current ();
    E -> map ();
    if (FamilyMgr::active) E -> economizeGeo ();
  }
  
  geometry_economized = FamilyMgr::active;
}


void  Field::buildBoundaries (const Field& f)
// ---------------------------------------------------------------------------
// Construct a new list of Boundary edges, copying BCs from f.
//
// The list is made by traversing the boundary_list in f, and copying 
// each Boundary, while making pointers into the current element_list.
// ---------------------------------------------------------------------------
{
  Boundary*  N;

  for (ListIterator<Boundary*> b(f.boundary_list); b.more(); b.next()) {
    N = new Boundary (*b.current (), element_list);
    boundary_list.add (N);
  }
}


void  Field::buildBoundaries (const Mesh& M)
// ---------------------------------------------------------------------------
// Construct a new list of Boundary edges, using BCmanager and Mesh M.
//
// NB: All boundary_lists made from the same mesh traverse elements/edges
// in the same order.
// ---------------------------------------------------------------------------
{
  Element*   E;
  Boundary*  B;
  BC*        bc;
  int        tag;
  int        id = 0;

  BCmanager::incVar ();

  for (ElmtsOfMesh e(M); e.more(); e.next()) {
    for (Mesh::SidesOfElmt s(e.current()); s.more(); s.next()) {
      if (s.current().isBoundary (tag)) {
	bc = BCmanager::getBC (tag);
	for (ListIterator<Element*> i(element_list); i.more(); i.next()) {
	  E = i.current ();
	  if (E -> ID() == e.current().ID) {
	    B = new Boundary (++id, E, s.current().ID, bc);
	    boundary_list.add (B);
	    break;
	  }
	}
      }
    }
  }
}


void  Field::printBoundaries (Field* F)
// ---------------------------------------------------------------------------
// Static member function.
//
// (Debugging) Utility to print information contained in a Boundary list.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Field::printBoundaries";
  char  s[StrMax];

  if (!F->boundary_list.length()) {
    message (routine, "empty Boundary list", WARNING); return;
  }

  sprintf (s, "# -- Field '%c' BOUNDARY LIST INFORMATION:", F -> field_name);
  message (routine, s, REMARK);

  for (ListIterator<Boundary*> k(F->boundary_list); k.more(); k.next())
    k.current () -> print ();
}


void  Field::evaluateBoundaries (int step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
// ---------------------------------------------------------------------------
{
  for (ListIterator<Boundary*> i(boundary_list); i.more(); i.next())
    i.current () -> evaluate (step);
}


void  Field::connect (Mesh& M, int np)
// ---------------------------------------------------------------------------
// Generate global node numbers for element boundaries.
//
// Method:
// 1)  Generate an initial set of numbers using Mesh functions;
// 2)  Sort the numbers so that Essential boundary conditions come last;
// 3)  Scatter the rearranged list into the bmap storage of the elements.
//
// Then carry out RCM renumbering and build element-boundary smoothing
// weights.
// ---------------------------------------------------------------------------
{
  char      routine[] = "Field::connect";
  Element*  E;
  int*      b;

  n_gid = M.connectSC (np);	// -- Set up M's internal storage.

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E = k.current();
    b = ivector (E -> nExt());
    for (ElmtsOfMesh e(M); e.more(); e.next()) {
      Mesh::Elmt& ME = e.current();
      if (E -> ID() == ME.ID) {
	ME.getGid (b);
	E -> gidInsert (b);
      }
    }
    freeVector (b);
  }

  M.connectSC   (2);		// -- Minimize M's internal storage.

  setMask       ();		// -- Mask ESSENTIAL BC nodes.
  sortGid       ();		// -- Sort ESSENTIAL BCs last.
  renumber      ();		// -- RCM renumbering.
  buildSmoother ();		// -- Global mass smoothing weights.

  if (option ("VERBOSE") > 1) {
    char  s[StrMax];
    sprintf (s, "-- GLOBAL NUMBERING FOR FIELD %c:", field_name);
    message (routine, s, REMARK);
    Field::printConnect (this);
  }
}


void  Field::sortGid ()
// ---------------------------------------------------------------------------
// Global node numbers get sorted to place essential-BC nodes last; this
// simplifies the later partition of global matrices---see the header for
// Field::buildSys.  The non-essential type node numbers can be sorted to
// optimize global matrix bandwidths, but this is not done here.
// ---------------------------------------------------------------------------
{
  int** gOrder = imatrix (n_gid, 2);        // For sorting global numbers.
  int*  bOld   = ivector (n_elmt_bnodes);   // Copy of original boundary map.
  int*  ramp   = ivector (n_gid);           // Ordered global node numbers.
  int*  mask   = ivector (n_gid);           // Corresponding solve mask.
  int*  shufl  = ivector (n_gid);           // Gids sorted by mask.

  Veclib::copy  (n_elmt_bnodes, elmt_bndry_gid, 1, bOld, 1);

  Veclib::ramp  (n_gid, 0,    1,  ramp, 1);
  Veclib::scatr (n_elmt_bnodes, elmt_bndry_mask, elmt_bndry_gid, mask);
  Veclib::copy  (n_gid, ramp, 1, *gOrder,   2);
  Veclib::copy  (n_gid, mask, 1, *gOrder+1, 2);

  n_solve = Veclib::count (n_gid, mask, 1);
  
  if (n_gid != n_solve) {
    /* -- Partition into "unknown" nodes & "essential BC" nodes,
          Then sort each part into ascending node number order.
	  Finally replace the old gids with the new, sorted ones. */

    qsort (*gOrder,           n_gid,           2*sizeof (int), sOrder);
    qsort (*gOrder,           n_solve,         2*sizeof (int), nOrder);
    qsort (*gOrder+2*n_solve, n_gid - n_solve, 2*sizeof (int), nOrder);

    Veclib::copy  (n_gid, *gOrder, 2, shufl, 1);
    Veclib::scatr (n_gid, ramp, shufl, mask);   // -- Mask <- new node numbers.
    Veclib::gathr (n_elmt_bnodes, mask, bOld, elmt_bndry_gid); // Replace.
  }

  freeMatrix (gOrder);
  freeVector (bOld);
  freeVector (ramp);
  freeVector (mask);
  freeVector (shufl);
}


void  Field::printConnect (Field* F)
// ---------------------------------------------------------------------------
// Static member function.
//
// (Debugging) utility to print up mesh connectivity information.
// ---------------------------------------------------------------------------
{
  Element  *E;
  int       offset;

  for (ListIterator<Element*> k (F -> element_list); k.more (); k.next ()) {
    E      = k.current ();
    offset = E -> nOff ();
    E -> printBndry (F -> data + offset);
  }
} 


void  Field::setMask ()
// ---------------------------------------------------------------------------
// Set up element solve masks.  The solve mask is really a shorthand for an
// essential-type boundary condition: if the solve-mask is zero, the BC is
// an essential one.  All other element edges have a solve mask of 1.
// ---------------------------------------------------------------------------
{
  int  ntot    = nBndry();
  int  gid_max = elmt_bndry_gid[Veclib::imax (ntot, elmt_bndry_gid, 1)] + 1;
  int* gmask   = ivector (gid_max);
  Veclib::fill (gid_max, 1, gmask, 1);

  Boundary*  B;

  for (ListIterator<Boundary*> b(boundary_list); b.more(); b.next()) {
    B = b.current();
    if (B -> isEssential()) B -> mask (gmask);
  }

  Veclib::gathr (ntot, gmask, elmt_bndry_gid, elmt_bndry_mask);

  freeVector (gmask);
}


void  Field::renumber ()
// ---------------------------------------------------------------------------
// From the initial ordering in specified by elmt_bndry_gid and
// elmt_bndry_mask, use RCM to generate a reduced-bandwidth numbering scheme.
// Reload into elmt_bndry_gid.
// ---------------------------------------------------------------------------
{
  // -- Build node adjacency tables.
  
  List<int>* adjncyList = new List<int> [n_solve];
  int        tabSize    = buildAdjncy (adjncyList);

  // -- Allocate memory.

  int* adjncy = ivector (tabSize+1);
  int* xadj   = ivector (n_solve+1);
  int* perm   = ivector (n_solve);
  int* mask   = ivector (n_solve);
  int* xls    = ivector (n_solve);

  int* invperm = ivector (n_gid);
  int* oldmap  = ivector (n_elmt_bnodes);

  // -- Pack adjacency tables into form required by genrcm routine, run it.

  fillAdjncy (adjncyList, adjncy, xadj, tabSize);

  genrcm (n_solve, xadj, adjncy, perm, mask, xls);
  Veclib::sadd (n_solve, -1, perm, 1, perm, 1);

  // -- Invert rcm's perm vector (identity for nodes that won't be changed).

  register int  i;
  for (i = 0;       i < n_solve; i++) invperm[perm[i]] = i;
  for (i = n_solve; i < n_gid;   i++) invperm[i]       = i;

  // -- Gather permuted node numbers into elmt_bndry_gid.

  Veclib::copy  (n_elmt_bnodes, elmt_bndry_gid, 1, oldmap, 1);
  Veclib::gathr (n_elmt_bnodes, invperm, oldmap, elmt_bndry_gid);

  // -- Clean up.

  delete [] adjncyList;

  freeVector (adjncy);
  freeVector (xadj);
  freeVector (perm);
  freeVector (mask);
  freeVector (xls);
  freeVector (invperm);
  freeVector (oldmap);

}


int  Field::buildAdjncy (List<int>* adjncyList) const
// ---------------------------------------------------------------------------
// Traverse list of elements and build up a vector of linked lists that
// describe the global nodes adjacent to each global node.
//
// Return the total amount of storage required when information is packed
// into an integer vector, as required by genrcm.
// ---------------------------------------------------------------------------
{
  for (ListIterator<Element*> k(element_list); k.more(); k.next())
    k.current() -> connectivSC (adjncyList);

  register int  i, ntab;
  for (i = 0, ntab = 0; i < n_solve; i++) ntab += adjncyList[i].length();

  return ntab;
}


void  Field::fillAdjncy (List<int>* adjncyList,
			 int*       adjncy    ,
			 int*       xadj      ,
			 int        tabSize   ) const
// ---------------------------------------------------------------------------
// Load the information contained in adjncyList into the two vectors
// adjncy & xadj required by genrcm.
// ---------------------------------------------------------------------------
{
  char          routine[] = "Field::fillAdjncy";
  register int  i, k;

  k = 0;
  for (i = 0; i < n_solve; i++) 
    for (ListIterator<int> p(adjncyList[i]); p.more(); p.next()) {
      adjncy[k] = p.current() + 1;
      xadj[i]   = k + 1;
      k++;
    }
  
  if (k != tabSize)
    message (routine, "after traversing list, k != tabSize", ERROR);

  adjncy[tabSize] = 0;
  xadj  [n_solve] = k + 1;
}


void  Field::buildSmoother ()
// ---------------------------------------------------------------------------
// Make element-boundary mass-averaged smoothing denominator in inv_mass
// storage area.  Call this procedure after creating a new field with its
// own global numbering system (after calling connect, so that n_gid is set).
// ---------------------------------------------------------------------------
{
  int       ntot;
  real*     unity;
  Element*  E;

  inv_mass = rvector (n_gid);
  Veclib::zero (n_gid, inv_mass, 1);

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E     = k.current ();
    ntot  = E -> nTot ();
    unity = rvector (ntot);
    Veclib::fill (ntot, 1.0, unity, 1);
    E -> bndryDsSum (unity, inv_mass);
    freeVector (unity);
  }

  Veclib::vrecp (n_gid, inv_mass, 1, inv_mass, 1);
}


Field&  Field::smooth ()
// ---------------------------------------------------------------------------
// Smooth values along element boundaries by direct stiffness summation.
//
// Mass-average smoothing.
// ---------------------------------------------------------------------------
{
  register int  offset;
  Element*      E;

  real* dssum = rvector (n_gid);
  Veclib::zero (n_gid, dssum, 1);

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E      = k.current ();
    offset = E -> nOff ();
    E -> bndryDsSum  (data + offset, dssum);
  }

  Veclib::vmul (n_gid, dssum, 1, inv_mass, 1, dssum, 1);

  for (k.reset(); k.more(); k.next()) {
    E      = k.current ();
    offset = E -> nOff ();
    E -> bndryInsert (dssum, data + offset);
  }

  freeVector (dssum);

  return *this;
}


Field&  Field::grad (int flag1, int flag2)
// ---------------------------------------------------------------------------
// Operate on Field(s) to produce the nominated index of the gradient.
// ---------------------------------------------------------------------------
{
  char      routine[] = "Field::grad";
  Element*  E;
  int       offset;
  ListIterator<Element*> i(element_list);

  if (flag1 && flag2) 
    message (routine, "can't make both components simultaneously", ERROR);
  else if (flag1)
    for (; i.more(); i.next()) {
      E      = i.current ();
      offset = E -> nOff ();
      E -> grad (data + offset, 0);
    }
  else if (flag2)
    for (; i.more(); i.next()) {
      E      = i.current ();
      offset = E -> nOff ();
      E -> grad (0, data + offset);
    }

  return *this;
}


void  Field::buildSys (real lambda2)
// ---------------------------------------------------------------------------
// Build a direct-solve discrete Helmholtz matrix system for this Field.
//
// Problem for solution is the statically condensed discrete version of the
// weak form of
//                                          2
//                       div grad u - lambda  u = f.
//
// This routine constructs the LHS matrices for this problem.
//
//
// CONVENTIONS:
// -----------
//   global    Helmholtz matrices: uppercase H,
//   elemental Helmholtz matrices: lowercase h.
//   SC <==> "Statically-Condensed".
//
// BACKGROUND:
// ----------
// A global Helmholtz matrix problem looks like this before partitioning:
//
//   +--------+-------------+ /  \     /  \
//   |        |             | |  |     |  |
//   |   Hp   |     Hc      | |u |     |f |   u : nodal values for solution.
//   |        |(constraint) | | s|     | s|    s
//   |        |             | |  |     |  |       (n_solve values)
//   +--------+-------------+ +--+     +--+
//   |        | Hess: this  | |  |  =  |  |
//   |        | partition   | |  |     |  |
//   |    T   | relates to  | |u |     |f |   u : are given essential BCs.
//   |  Hc    | essential   | | g|     | g|    g
//   |        | BC nodes &  | |  |     |  |       (n_gid - n_solve values)
//   |        | is not used | |  |     |  |
//   |        | or computed | |  |     |  |
//   +--------+-------------+ \  /     \  /
//
// Partition out the parts of the matrix corresponding to the known nodal
// values (ESSENTIAL BCs), and solve instead
//
//   +--------+ /  \     /  \     +-------------+ /  \
//   |        | |  |     |  |     |             | |  |
//   |   Hp   | |u |     |f |     |     Hc      | |  |
//   |        | | s|  =  | s|  -  |             | |u |
//   |        | |  |     |  |     |             | | g|
//   +--------+ \  /     \  /     +-------------+ |  |
//                                                |  |
//                                                |  |
//                                                \  /
//
// Here n_gid is the number of nodes that receive global node numbers,
// typically those on the mesh edges.  N_solve is the number of these
// nodes that have values that must be solved for, i.e. n_gid minus the
// number of global nodes situated on essential-type boundaries.
//
// If there are no given essential BCs and the Helmholtz constant lambda2 is
// zero, the last row and column of Hp are never constructed (n_solve=n_gid-1)
// since the full matrix is singular: the last value is arbitrarily assigned
// zero in order to pin the solution and Hc, nominally n_gid x 1, is not used.
//
// CONSTRUCTION STRATEGY:
// ---------------------
// 1.  Form elemental Helmholtz matrices.
//
// 2.  If possible, construct the statically condensed and partitioned forms
//     of the elemental Helmholtz matrices, keeping the elemental
//     edge/interior coupling and factored interior resolution matrices in a
//     list.  Can't do this step for bilinear elements (no interior nodes).
//
// 3.  Post the (condensed) elemental Helmholtz matrices to the global
//     Helmholtz matrices Hp & Hc, accounting for any special packing scheme
//     (band/symmetric) in use.
//
// 4.  Factor the global Helmholtz matrix.
//
// CONSTRAINT PARTITION:
// --------------------
// As far as the constraint partitioning is concerned, there are 4 cases:
// case 0:  n_solve == n_gid:  For all-natural boundaries, OK if lambda2 != 0,
//          since global Helmholtz matrix will be non-singular.
//          In this case there is no constraint matrix partition.
// case 1:  n_solve == n_gid - 1: This is the case imposed for all-natural
//          boundaries when lambda2 == 0.  Constraint partition nominally
//          n_solve x 1, but is not needed (or constructed) since value
//          of variable for last gid is assumed to be zero.
// case 2:  n_solve < n_gid - 1:  Constraint matrix partition is of size
//          n_solve x (n_solve - n_gid).
// case 3:  n_solve == 0:  This for a single-element problem when all
//          boundaries are of essential type.  Here there also there
//          will be no constraint matrix since the internal SC matrix
//          is used to obtain solution for internal nodes (if any).
// ---------------------------------------------------------------------------
{
  char    routine[] = "Field::buildSys";
  char    s[StrMax];
  int     info, verbose = option ("VERBOSE");
  real**  hbb;
  real**  rmat;
  real*   rwrk;
  int*    iwrk;
  real    condition;

  // -- Set up. n_band additionally flags use of packed or band storage.

  n_band = n_pack = n_cons = 0;

  if (n_solve) {
    if (n_gid == n_solve && lambda2 == 0.0) n_solve--;     // -- Pin solution.
    n_band = globalBandwidth ();
    n_band = (n_band < (n_solve + 1) >> 1) ? n_band : 0;
    n_pack = (n_band) ? n_solve * n_band : ((n_solve + 1) * n_solve) >> 1;
    n_cons = n_gid - n_solve;
  
    Hp = rvector (n_pack);
    Veclib::zero (n_pack, Hp, 1);

    if (n_cons > 1) {
      Hc = rmatrix (n_solve, n_cons);
      Veclib::zero (n_solve*n_cons, *Hc, 1);
    }

    if (verbose > 0) {
      if (n_band)
	sprintf (s, ": Banded system matrix:     %1dx%1d\t(%1d words)",
		 n_solve, n_band, n_pack);
      else
	sprintf (s, ": Triangular system matrix: %1dx%1d\t(%1d words)",
		 n_solve, n_solve, n_pack);
      message (routine, s, REMARK);
      sprintf (s, "Constraint matrix:        %1dx%1d\t(%1d words)",
	       n_solve, n_cons, n_solve * n_cons);
      message ("                 ", s, REMARK);
    }
  }

  // -- Loop over elements, creating & posting elemental Helmholtz matrices.

  Element*  E;

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E = k.current();

    int nExt  = E -> nExt();
    int nInt  = E -> nInt();
    int nTot  = E -> nTot();
    int nKnot = E -> nKnot();
		  
    hbb  = rmatrix (nExt, nExt);
    rmat = rmatrix (nKnot, nKnot);
    rwrk = rvector (nExt*nTot);
 
    E -> HelmholtzSC (lambda2, hbb, rmat, rwrk);
    E -> post        (hbb, Hp, Hc, n_solve, n_cons, n_band);

    freeMatrix (hbb );
    freeMatrix (rmat);
    freeVector (rwrk);
 
    if (FamilyMgr::active) E -> economizeMat ();
  }

  matrices_economized = FamilyMgr::active;

  // -- Factor global Helmholtz matrix.

  if (n_solve) {
    if   (n_band) Lapack::pbtrf ("U", n_solve, n_band - 1, Hp, n_band, info);
    else          Lapack::pptrf ("U", n_solve,             Hp,         info);

    if (info) message (routine, "failed to factor matrix", ERROR);

    if (verbose > 0) {
      rwrk = rvector (3 * n_solve);
      iwrk = ivector (n_solve);
      if   (n_band) Lapack::pbcon ("U", n_solve, n_band - 1, Hp, n_band, 1.0,
				   condition, rwrk, iwrk, info);
      else          Lapack::ppcon ("U", n_solve,             Hp,         1.0,
				   condition, rwrk, iwrk, info);
      sprintf (s, "System condition number:  %g", condition);
      message ("                 ", s, REMARK);
      freeVector (rwrk);
      freeVector (iwrk);
    }
  }
}


void  Field::solveSys (Field* F)
// ---------------------------------------------------------------------------
// Carry out direct solution of this Field using F as forcing.
//
// Problem for solution is the statically condensed discrete version of the
// weak form of
//                                          2
//                       div grad u - lambda  u = F.
//
// This routine creates the RHS vector from the input forcing field F
// and the Field's boundary conditions.  Then solution is carried out as
// indicated in the header to buildSys.
//
// The RHS vector is constructed with length of the number of element-edge
// nodes in the problem (n_gid).  The first n_solve values contain forcing
// terms for the free (not ESSENTIAL BC) nodes in the problem, derived from
// the forcing field and the NATURAL BCs, while the remaining n_cons values
// get loaded from ESSENTIAL BC values.
//
// Forcing field F is overwritten/destroyed during processing.
// ---------------------------------------------------------------------------
{
  int  info;

  // -- Allocate, fill values of RHS vector.

  real* RHS  = rvector (n_gid);
  Veclib::zero (n_gid, RHS, 1);
  buildRHS     (F, RHS);

#ifdef DEBUG
  char  routine[] = "Field::solveSys";
  char  s[StrMax];
  if (option ("VERBOSE") > 2) {
    sprintf (s, ": -- Solution for Field: %'%c', n_gid: %1d, n_solve: %1d",
	     field_name, n_gid, n_solve);
    message (routine, s, REMARK);
    message ("", "-- (unconstrained) RHS:", REMARK);
    for (int i = 0; i < n_gid; i++) {
      sprintf (s, "  RHS %5d: %g", i, RHS[i]);
      message ("", s, REMARK);
    }
  }
#endif

  // -- Solve for unknown global-node values (if any), applying Hc, then Hp.

  if (n_solve) {
    if (n_cons > 1)   // -- Apply ESSENTIAL BC constraint matrix to RHS.
      Blas::gemv ("T", n_cons, n_solve, -1.0, *Hc, n_cons,
		  RHS + n_solve, 1, 1.0, RHS, 1);

#ifdef DEBUG
    if (option ("VERBOSE") > 2) {
      message ("", "-- (constrained) RHS:", REMARK);
      for (int i = 0; i < n_gid; i++) {
	sprintf (s, "  RHS %5d: %g", i, RHS[i]);
	message ("", s, REMARK);
      }
    }
#endif

    // -- Solve unknown global-node values using factored matrix Hp.
    if (n_band) 
      Lapack::pbtrs ("U", n_solve, n_band-1, 1, Hp, n_band, RHS, n_gid, info);
    else     
      Lapack::pptrs ("U", n_solve,           1, Hp,         RHS, n_gid, info);
  }

#ifdef DEBUG
  if (option ("VERBOSE") > 2) {
    message ("", "-- Solution vector:", REMARK);
    for (int i = 0; i < n_gid; i++) {
      sprintf (s, "  RHS %5d: %g", i, RHS[i]);
      message ("", s, REMARK);
    }
  }
#endif

  // -- Resolve element external (and internal, if S-C) nodes.

  Element*  U;
  int       offset;

  for (ListIterator<Element*> u(element_list); u.more(); u.next()) {
    U      = u.current ();
    offset = U -> nOff ();
    U -> resolveSC (RHS, F -> data + offset, data + offset);
  }

  freeVector (RHS);
}


void Field::buildRHS (Field* F, real* RHS) const
// ---------------------------------------------------------------------------
// RHS is true RHS for solvable global nodes, LHS for essential BC ones.
//
// Both forcing and NATURAL BCs make a contribution to the true RHS
// (stored in the first n_solve locations of RHS), while ESSENTIAL BCs
// determine LHS values (stored in the remaining locations of RHS).
//
// Various RHS terms are simplified by use of diagonal mass matrix.
//
// On input, F contains the forcing field values.
// NB: Side-effect: F is changed in dsForcingSC, using internal matrices
// for the current Field.
// ---------------------------------------------------------------------------
{
  Element*  E;
  int       offset;

  for (ListIterator<Element*> u (element_list); u.more (); u.next ()) {
    E      = u.current ();
    offset = E -> nOff ();
    E -> dsForcingSC (F -> data + offset, RHS);
  }

  Boundary *B;
  for (ListIterator<Boundary*> b(boundary_list); b.more(); b.next()) {
    B = b.current();
    if   (B -> isEssential ()) B -> enforce (RHS);
    else                       B -> dsSum   (RHS);
  }

  if (n_solve == n_gid - 1) RHS[n_gid - 1] = 0.0;  // Pin last value to earth.
}


int Field::globalBandwidth () const
// --------------------------------------------------------------------------
// Find the bandwidth of the assembled global matrix (including diagonal).
// --------------------------------------------------------------------------
{
  int  nband = 0;

  for (ListIterator<Element*> k(element_list); k.more(); k.next())
    nband = max (k.current() -> bandwidthSC (), nband);

  ++nband; // -- Diagonal.

  return nband;
}


int  Field::switchPressureBCs (const BC* hopbc, const BC* zero)
// ---------------------------------------------------------------------------
// This is called by PBCmanager.  Traverse the Field boundaries, and switch
// OUTFLOW boundaries to be ESSENTIAL (value 0.0), all other types to be
// HOPBC.  Return the total amount of edge-based storage on traverse.
// ---------------------------------------------------------------------------
{
  Boundary* B;
  int       ntot(0);

  for (ListIterator<Boundary*> j(boundary_list); j.more(); j.next()) {
    B = j.current ();
    B -> resetKind (hopbc, zero);
    ntot += B -> nKnot ();
  }

  return ntot;
}


void  Field::errors (const char* function)
// ---------------------------------------------------------------------------
// Compare F with function, print the infinity-norm Li, the 2-norm L2
// and the Sobolev 1-norm H1.
//
// The norms are found element-by-element, using projection onto higher-order
// elements and high-order quadrature.
// ---------------------------------------------------------------------------
{
  const int NQUAD = 15;

  Element* E;
  Element* P;
  real     area = 0.0;
  real     Li   = 0.0;
  real     L2   = 0.0;
  real     H1   = 0.0;
  real*    sol;
  real*    v;
  real*    x;
  int      ntot, nmsh;

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E = k.current ();
    P = new Element (*E, NQUAD);

    ntot = P -> nTot ();
    nmsh = P -> nMsh ();

    v   = rvector (ntot);
    x   = rvector (nmsh);
    sol = rvector (ntot);

    P -> install (0, x, 0, 0);
    P -> project (*E, data + E -> nOff (), v);
    P -> map     ();

    P -> evaluate (function, sol);
    Veclib::vsub  (ntot, v, 1, sol, 1, v, 1);

    area += P -> area ();
    Li    = max (Li, P -> norm_inf (v));
    L2   += P -> norm_L2 (v);
    H1   += P -> norm_H1 (v);

    freeVector (v);
    freeVector (x);
    freeVector (sol);
  }
  
  L2 /= area;
  H1 /= area;

  char  s[StrMax];
  sprintf (s, "Field '%c' error norms (inf, L2, H1): %.3e %.3e %.3e",
	   field_name, Li, L2, H1);
  message ("", s, REMARK);
}


Field&  Field::evaluate (const char* function)
// ---------------------------------------------------------------------------
// Evaluate function over each element.
// ---------------------------------------------------------------------------
{
  vecInit   ("x y", function);
  vecInterp (n_data, mesh, mesh + n_data, data);

  return *this;
}


Vector  Field::normalTraction (const Field& pre, const Field& vel)
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute normal tractive forces on all WALL boundaries, taking pre to be
// the pressure field and getting WALL BCs from vel (a velocity Field).
//
// Rely on the fact that all boundary lists traverse mesh in same order
// (as made by field::buildBoundaries).
// ---------------------------------------------------------------------------
{
  ListIterator<Boundary*>  p(pre.boundary_list);
  ListIterator<Boundary*>  u(vel.boundary_list);
  Boundary*                P;
  Boundary*                U;
  Vector                   secF, F = {0.0, 0.0, 0.0};
  
  for (; u.more(); u.next(), p.next()) {
    P = p.current ();
    U = u.current ();
    if (U -> isWall ()) {
      secF = P -> normalTraction (pre.data + P -> nOff (), P -> nSkip ());
      F.x += secF.x;
      F.y += secF.y;
    }
  }

  return F;
}


Vector  Field::tangentTraction (const Field& U , const Field& V ,
			              Field& Dx,       Field& Dy)
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute (2D) tangential viscous tractive forces on all WALL boundaries,
// treating U & V as first and second velocity components, respectively.
// Du & Dv are workspace Field storage used for computation of velocity
// gradients.
//
// Compute viscous tractive forces on wall from
//
//  t_i  = - T_ij * n_j       (minus sign for force exerted BY fluid ON wall),
//
// where
//
//  T_ij = viscous stress tensor
//                          dU_i    dU_j
//       = RHO * KINVIS * ( ----  + ---- ) .
//                          dx_j    dx_i
//
// ---------------------------------------------------------------------------
{
  real                     mu = dparam ("KINVIS") * dparam ("RHO");
  ListIterator<Boundary*>  u(U.boundary_list);
  ListIterator<Boundary*>  v(V.boundary_list);
  Boundary*                Bu;
  Boundary*                Bv;
  Vector                   secF, F = {0.0, 0.0, 0.0};

  // -- Terms from du_1/dx_i.

  Dx = Dy = U;
  
  Dx.grad(1, 0).smooth();
  Dy.grad(0, 1).smooth();
 
  for (; u.more(); u.next(), v.next ()) {
    Bu = u.current ();
    Bv = v.current ();
    if (Bu -> isWall ()) {
      secF = Bu -> tangentTraction (Dx.data + Bu -> nOff (),
				    Dy.data + Bv -> nOff (),
				    Bv -> nSkip (),  mu, 1);
//      printf ("sec: %1d, Fx: %f\n", Bu -> ID(), secF.x);
      F.x += secF.x;
      F.y += secF.y;
    }
  }

  // -- Terms form du_2/dx_i.

  Dx = Dy = V;
  
  Dx.grad(1, 0).smooth();
  Dy.grad(0, 1).smooth();

  for (u.reset(), v.reset(); u.more(); u.next(), v.next ()) {
    Bu = u.current ();
    Bv = v.current ();
    if (Bu -> isWall ()) {
      secF = Bu -> tangentTraction (Dx.data + Bu -> nOff (),
				    Dy.data + Bv -> nOff (), 
				    Bv -> nSkip (),  mu, 2);
//      printf ("sec: %1d, Fx: %f\n", Bu -> ID(), secF.x);
      F.x += secF.x;
      F.y += secF.y;
    }
  }

  return F;
}
