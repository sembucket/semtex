/*****************************************************************************
 * field.C:  routines to deal with Fields.
 *****************************************************************************/

// $Id$

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
Hc              (0)
{ }




 
Field::Field (const Field& f, char tag)
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
  Element* E;
  real*    d;
  int*     g;
  int*     s;

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
    
    for (ListIterator<Element*> k(f.element_list); k.more(); k.next()) {
      E = new Element (*k.current());
      E -> install (d, 0, g, s);
      element_list.add (E);
      d += E -> nTot();
      g += E -> nExt();
      s += E -> nExt();
    }
    
    if (tag == 'p')
      buildBoundaries (f);
    else
      buildBoundaries ();

  } else {
    for (ListIterator<Element*> k(f.element_list); k.more(); k.next()) {
      E = new Element (*k.current());
      E -> install (d, 0, 0, 0);
      element_list.add (E);
      d += E -> nTot();
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
  for (ListIterator<Element*> k(element_list); k.more(); k.next())
    k.current () -> evaluate (function);

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





void  Field::readMesh (istream& strm, int np)
// ---------------------------------------------------------------------------
// Read mesh description from mesh stream strm.
//
// Mesh node locations (xmesh & ymesh) are computed, but geometric partials
// are not.
//
// Mesh file description:
// ---------------------
// The first line of the file contains the number of element corner
// vertices: <num> VERTICES.  Following this is a list of corner vertices
// given in Cartesian form, one vertex per line.  E.g.:
// 4 VERTICES
// 0.0  0.0
// 1.0  0.0
// 1.0  1.0
// 0.0  1.0
//
// The list of vertices is followed by a blank line, then a section which
// describes the elements, starting with the line <num> ELEMENTS.
// This is also followed by a blank line and domain information is then
// given on an element-wise basis.  No global node-numbering is needed; this
// is generated automatically (& dealt with elsewhere).
//
// Corner vertex numbering proceeds CCW round the element.  Side 1 lies
// between corner vertices 1 & 2 and side numbering is also CCW-sequential
// around the element.  Sides can either lie on the domain boundary (BC) or
// mate with another element (EL).  Internal numbering of elements and sides
// starts at zero.
//
// Information for each element is followed by a blank line, except that the
// last element may be followed by EOF.
//
// Element data for a biquadratic element:
// ELEMENT number ORDER 2
// 1 2 3 4                             (these are indices in the vertex list)
// SIDE 1  BC  tag-num
// side 2  EL  mate-el SIDE mate-side
// SI   3  EL  mate-el si   mate-side
// SIDE 4  EL  mate-el SIDE mate-side
// ---------------------------------------------------------------------------
{
  char  routine[] = "Field::readMesh";
  char  s1[StrMax], s2[StrMax];
  int   nvert;

  // -- Input vertex information.

  while (strm.getline (s1, StrMax)) {
    upperCase (s1);
    if ( (strstr (s1, "VERTICES")) && (sscanf (s1, "%d", &nvert))) {
      strm.getline (s1, StrMax);
      break;
    } else {
      sprintf (s2, "couldn't get number of vertices: %s", s1);
      message (routine, s2, ERROR);
    }
  }

  Point* vertexList = new Point[nvert];

  int  i(0);
  while (strm.getline (s1, StrMax) && i < nvert) {
    if (sscanf (s1, "%*s%lf%lf", &vertexList[i].x, &vertexList[i].y) != 2) {
      sprintf (s2, "expected info for vertex %1d, got: %s", i+1, s1);
      message (routine, s2, ERROR);
    }
    ++i;
  }
  
  // -- Check that vertex input completed OK.

  if (i < nvert) {
    sprintf (s1, "expected another vertex, got newline at vertex #%1d", i);
    message (routine, s1, ERROR);
  } else if (s1[0]) {
    sprintf (s2, "read all vertices, but next line not blank: %s", s1);
    message (routine, s2, ERROR);
  }

  // -- Input element information.

  int nel;

  while (strm.getline (s1, StrMax)) {
    upperCase (s1);
    if ( (strstr (s1, "ELEMENTS")) && (sscanf (s1, "%d", &nel)))
      break;
    else{
      sprintf (s2, "couldn't get number of elements: %s", s1);
      message (routine, s2, ERROR);
    }
  }

  int       id;
  int**     vertexTable = new int* [nel];
  Element*  E;

  for (i = 0; i < nel; i++) {
    strm.getline(s1, StrMax).getline(s1, StrMax);
    upperCase (s1);

    if (strstr (s1, "ELEMENT")) {

      element_list.add (E = new Element);
      
      sscanf (s1, "%*s %d", &id);
      if (id > nel) {
	sprintf (s2, "element id: %1d exceeds no. of elements: %1d", id, nel);
	message (routine, s2, ERROR);
      }

      E -> setState (--id, np, 4);
      E -> read     (strm, vertexTable[id]);

    } else {
      sprintf (s2, "element descriptor for #%1d? : %s", i+1, s1);
      message (routine, s2, ERROR);
    }
  }

  // -- Check element input and nominated element connectivity.

  int* check = ivector (nel);
  Veclib::fill (nel, 1, check, 1);
  for (ListIterator<Element*> k(element_list); k.more(); k.next())
    check[k.current() -> ID ()] = 0;
  if (Veclib::any (nel, check, 1)) {
    sprintf (s1, "no input for element %1d", Veclib::first (nel, check, 1)+1);
    message (routine, s1, ERROR);
  }
  freeVector (check);

  for (k.reset(); k.more(); k.next()) k.current() -> checkMate (element_list);

  // -- All the mesh input file has been read; allocate & install storage.

  n_data = n_mesh = n_elmt_bnodes = 0;

  for (k.reset(); k.more(); k.next()) {
    E = k.current();
    n_data        += E -> nTot ();
    n_mesh        += E -> nMsh ();
    n_elmt_bnodes += E -> nExt ();
  }

  data            = rvector (n_data);
  mesh            = rvector (n_mesh);
  elmt_bndry_gid  = ivector (n_elmt_bnodes);
  elmt_bndry_mask = ivector (n_elmt_bnodes);

  double* d = data;
  double* m = mesh;
  int*    g = elmt_bndry_gid;
  int*    s = elmt_bndry_mask;

  for (k.reset(); k.more(); k.next(), i++) {
    E = k.current();
    E -> install (d, m, g, s);
    d += E -> nTot ();
    m += E -> nMsh ();
    g += E -> nExt ();
    s += E -> nExt ();
  }

  // -- Fill mesh internal node locations.

  for (k.reset(), i = 0; k.more(); k.next(), i++)
    k.current() -> mesh (vertexTable[i], vertexList);

  for (i = 0; i < nel; i++) delete [] vertexTable[i];
  delete [] vertexTable;
  delete [] vertexList;

  field_name = 'u';
}





void  Field::printMesh (Field* F)
// ---------------------------------------------------------------------------
// Mesh location information is written out element-by-element.
// ---------------------------------------------------------------------------
{
  cout
    << F -> element_list.first() -> nKnot() << " "
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
  for (ListIterator<Element*> i(element_list); i.more(); i.next())
    i.current() -> map();
}





void  Field::buildBoundaries (const Field& f)
// ---------------------------------------------------------------------------
// Construct a new list of Boundary edges, copying BCs from f.
//
// The list is made by traversing the boundary_list in f, and copying 
// each Boundary, while making pointers into the current element_list.
// ---------------------------------------------------------------------------
{
  Boundary* N;

  for (ListIterator<Boundary*> b(f.boundary_list); b.more(); b.next()) {
    N = new Boundary (*b.current (), element_list);
    boundary_list.add (N);
  }
}





void  Field::buildBoundaries ()
// ---------------------------------------------------------------------------
// Construct a new list of Boundary edges, using BCmanager.
//
// The list is made by searching the list of elements: each time an element
// edge is identified as lying on the domain boundary, a new Boundary structure
// is made and added to the list.  Boundary condition information is found
// by searching the list of boundary condition specifiers for a matching tag,
// after which the a pointer to its specifier is placed in Bedge structure.
//
// NB: estart values assume BLAS-conformant behaviour with negative skips.
// ---------------------------------------------------------------------------
{
  int       side, id  = 0;
  Element*  E;
  Boundary* B;
  BC*       bc;

  BCmanager::incVar ();

  for (ListIterator<Element*> i(element_list); i.more(); i.next()) {
    E = i.current();
    for (side = 0; side < E -> nSide(); side++)
      if (E -> sideKind (side) == DOMAIN_BOUNDARY) {
	bc = BCmanager::getBC (E -> sideTag (side));
	B  = new Boundary (id++, E, side, bc);
	boundary_list.add(B);
      }
  }
}





void  Field::printBoundaries (Field* F)
// ---------------------------------------------------------------------------
// (Debugging) Utility to print information contained in a Boundary list.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Field::printBoundaries";

  if (!F->boundary_list.length()) {
    message (routine, "empty Boundary list", WARNING); return;
  }

  cout << "# -- Field '" << F->field_name <<
    "' BOUNDARY LIST INFORMATION:" << endl;

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





void  Field::readConnect (const char *session)
// ---------------------------------------------------------------------------
// Read in global mesh-numbering information generated by CONNECT, and
// re-arrange so that Essential BCs have highest numbers.
//
// Input "name" is the root of session, for connectivity file "name.con".
//
// The data then contained in element bmaps are the global equation numbers
// (if solve mask is non-zero).
//
// Method:
// 1)  Read in a list of global node numbers generated by CONNECT.
// 2)  Sort the numbers so that Essential boundary conditions come last.
// 3)  Scatter the rearranged list into the bmap storage of the elements.
// ---------------------------------------------------------------------------
{
  char       routine[] = "Field::readConnect";
  char       s[StrMax];

  int                     elmt, nLines = 0, nExpect = 0;
  ListIterator<Element*>  e(element_list);
  Element*                Esave   = 0;
  int                     nel     = nEl();

  ifstream file (strcat (strcpy (s, session), ".con"));
  if (!file) message (routine, "can't open connectivity file", ERROR);

  n_gid = 0;
  while (file.getline (s, StrMax)) {
    if (s[0] == '#') continue;

    if (sscanf (s, "%d", &elmt) != 1)
      message (routine, "failed scanning connectivity data",  ERROR);
    if (elmt > nel)
      message (routine, "element number exceeds declaration", ERROR);
    --elmt;
    ++nLines;
 
    if (Esave && elmt == Esave -> ID ())
      n_gid = max (n_gid, Esave -> gidInsert (s));
    else {
      Esave = 0;
      for (e.reset (); !Esave && e.more (); e.next ())
	if (elmt  == e.current () -> ID ()) {
	  Esave    = e.current ();
	  nExpect += Esave -> nKnot () * Esave -> nSide ();
	}
      if (!Esave) {
	sprintf (s, "couldn't locate element %1d in list", elmt + 1);
	message (routine, s, ERROR);
      } else
	n_gid = max (n_gid, Esave -> gidInsert (s));
    }
  }

  file.close();

  if (nLines != nExpect)
    message (routine, "number of connectivity data mismatch", ERROR);

  setMask ();

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
// (Debugging) Utility to print up mesh connectivity information.
// ---------------------------------------------------------------------------
{
  cout << "# -- ELEMENT-BOUNDARY CONNECTIVITY & VALUE INFORMATION --" << endl;

  for (ListIterator<Element*> k(F->element_list); k.more(); k.next())
    k.current() -> printBndry ();
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

  int* adjncy  = ivector (tabSize);
  int* xadj    = ivector (n_solve+1);
  int* perm    = ivector (n_solve);
  int* mask    = ivector (n_solve);
  int* xls     = ivector (n_solve);

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
  char  routine[] = "Field::fillAdjncy";
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
  xadj[n_solve]   = k + 1;
}




Field&  Field::dsSmooth ()
// ---------------------------------------------------------------------------
// Smooth values along element boundaries by direct stiffness summation.
//
// Simple average smoothing.
// ---------------------------------------------------------------------------
{
  register int            i;
  ListIterator<Element*>  k(element_list);

  real* dssum = rvector (n_gid);
  real* denom = rvector (n_gid);

  Veclib::zero (n_gid, dssum, 1);
  Veclib::zero (n_gid, denom, 1);

  for (i = 0; i < n_elmt_bnodes; i++) denom[elmt_bndry_gid[i]] += 1.0;

  for (k.reset(); k.more(); k.next()) k.current() -> bndryDsSum  (dssum);

  Veclib::vdiv (n_gid, dssum, 1, denom, 1, dssum, 1);

  for (k.reset(); k.more(); k.next()) k.current() -> bndryInsert (dssum);

  freeVector (dssum);
  freeVector (denom);

  return *this;
}




Field&  Field::grad (int index)
// ---------------------------------------------------------------------------
// Operate on Field to produce the nominated index of the gradient.
// ---------------------------------------------------------------------------
{
  ListIterator<Element*> i(element_list);

  switch (index) {
  case 0:  for     (; i.more(); i.next())  i.current() -> d_dx (); break;
  case 1:  for     (; i.more(); i.next())  i.current() -> d_dy (); break;
  default: message ("Field::grad(int)", "illegal index", ERROR);   break;
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
  char  routine[] = "Field::buildSys";
  int   info;

  // -- Set up. n_band additionally flags use of packed or band storage.

  n_band = n_pack = n_cons = 0;

  if (n_solve) {
    if (n_gid == n_solve && lambda2 == 0.0) n_solve--; // Pin solution.
    n_band = globalBandwidth ();
    n_band = (n_band < (n_solve+1)>>1) ? n_band : 0;
    n_pack = (n_band) ? n_solve * n_band : ((n_solve + 1) * n_solve)>>1;
    n_cons = n_gid - n_solve;
  
    Hp = rvector (n_pack);
    Veclib::zero (n_pack, Hp, 1);

    if (n_cons > 1) {
      Hc = rmatrix (n_solve, n_cons);
      Veclib::zero (n_solve*n_cons, *Hc, 1);
    }
  }

  // -- Loop over elements, creating & posting elemental Helmholtz matrices.

  Element* E;

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E = k.current();

    int nExt  = E -> nExt();
    int nInt  = E -> nInt();
    int nTot  = E -> nTot();
    int nKnot = E -> nKnot();
		  
    real** hbb  = rmatrix (nExt, nExt);
    real** dmat = rmatrix (nKnot, nKnot);
    real*  dwrk = rvector (nExt*nTot);
 
    E -> HelmholtzSC (lambda2, hbb, dmat, dwrk);
    E -> post        (hbb, Hp, Hc, n_solve, n_cons, n_band);

    freeMatrix (hbb );
    freeMatrix (dmat);
    freeVector (dwrk);
  }

  // -- Factor global Helmholtz matrix.

  if (n_solve)
    if (n_band) {
      Lapack::pbtrf ("U", n_solve, n_band - 1, Hp, n_band, info);
      if (info) message (routine, "pbtrf failed to factor matrix", ERROR);
    } else {
      Lapack::pptrf ("U", n_solve, Hp, info);
      if (info) message (routine, "pptrf failed to factor matrix", ERROR);
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
  buildRHS (F, RHS);

  // -- Solve for unknown global-node values (if any), applying Hc, then Hp.

  if (n_solve) {
    if (n_cons > 1)   // -- Apply ESSENTIAL BC constraint matrix to RHS.
      Blas::gemv ("T", n_cons, n_solve, -1.0, *Hc, n_cons,
		  RHS + n_solve, 1, 1.0, RHS, 1);

    // -- Solve unknown global-node values using factored matrix Hp.
    if (n_band) 
      Lapack::pbtrs ("U", n_solve, n_band-1,1, Hp, n_band, RHS, n_solve, info);
    else     
      Lapack::pptrs ("U", n_solve,          1, Hp,         RHS, n_solve, info);
  }

  // -- Resolve element external (and internal, if S-C) nodes.

  ListIterator<Element*> f(F->element_list);
  for (ListIterator<Element*> u(element_list);
       u.more(), f.more();
       u.next(), f.next())
    u.current () -> resolveSC (RHS, f.current ());

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
  ListIterator<Element*>  u(element_list);
  ListIterator<Element*>  f(F->element_list);

  for (; u.more(); u.next(), f.next())
    f.current() -> dsForcingSC (u.current(), RHS);

  if (n_solve == n_gid - 1) RHS [n_gid - 1] = 0.0;  // Fix last value.

  Boundary *B;
  for (ListIterator<Boundary*> b(boundary_list); b.more(); b.next()) {
    B = b.current();
    if   (B -> isEssential()) B -> enforce (RHS);
    else                      B -> dsSum   (RHS);
  }
}





int Field::globalBandwidth () const
// --------------------------------------------------------------------------
// Find the bandwidth of the assembled global matrix (including diagonal).
// --------------------------------------------------------------------------
{
  int  nband = 0;

  for (ListIterator<Element*> k(element_list); k.more(); k.next())
    nband = max (k.current() -> bandwidthSC (), nband);

  ++nband; // Diagonal.

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
    B = j.current();
    B -> resetKind (hopbc, zero);
    ntot += B -> nKnot();
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
  const int Nquad = 15;

  Element* E;
  Element* P;
  real     area = 0.0;
  real     Li   = 0.0;
  real     L2   = 0.0;
  real     H1   = 0.0;
  real*    tmp;
  real*    sol;
  real*    v;
  real*    x;
  int      ntot, nmsh;

  for (ListIterator<Element*> k(element_list); k.more(); k.next()) {
    E = k.current ();
    
    P = new Element (*E, Nquad);
    ntot = P -> nTot ();
    nmsh = P -> nMsh ();

    v   = rvector (ntot);
    x   = rvector (nmsh);
    tmp = rvector (ntot);
    sol = rvector (ntot);

    P -> install (v, x, 0, 0);
    P -> project (*E);
    P -> map     ();

    P -> extract  (tmp);
    P -> evaluate (function, sol);
    Veclib::vsub  (ntot, tmp, 1, sol, 1, tmp, 1);
    P -> insert   (tmp);

    area += P -> area ();
    Li    = max (Li, P -> norm_inf ());
    L2   += P -> norm_L2 ();
    H1   += P -> norm_H1 ();

    freeVector (v);
    freeVector (x);
    freeVector (tmp);
    freeVector (sol);
  }
  
  L2 /= area;
  H1 /= area;

  cout << "-- Error norms for Field " << field_name << " (inf, L2, H1):";
  cout << Li << "  " << L2 << "  " << H1 << endl;

}





Field&  Field::evaluate (const char* function)
// ---------------------------------------------------------------------------
// Evaluate function over each element.
// ---------------------------------------------------------------------------
{
  for (ListIterator<Element*> k(element_list); k.more(); k.next())
    k.current () -> evaluate (function);

  return *this;
}
