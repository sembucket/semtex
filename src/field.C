///////////////////////////////////////////////////////////////////////////////
// field.C:  routines for Field class.
///////////////////////////////////////////////////////////////////////////////

static char 
RCSid[] = "$Id$";


#include <Fem.h>


#ifdef __DECCXX
  #pragma define_template max<int>
#endif


Field::Field (Mesh&      M ,
	      const int& np)
// ---------------------------------------------------------------------------
// Use M to allocate storage for Field and its Elements.
//
// Compute Mesh knot locations (xmesh & ymesh), but not geometric partials.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      k;
  const int         N = M.nEl ();

  Elmt.setSize (N);

  n_data = n_mesh = n_elmt_bnodes = 0;

  // -- Create all the Elements using their Mesh description.

  k = 0;
  for (ElmtsOfMesh e (M); e.more (); e.next (), k++) {
    E = Elmt[k] = new Element (e.ID (), np, e.nSides ());

    n_data        += E -> nTot ();
    n_mesh        += E -> nMsh ();
    n_elmt_bnodes += E -> nExt ();

    elmt_np_max = max (elmt_np_max, E -> nKnot ());
    elmt_nt_max = max (elmt_nt_max, E -> nTot  ());
    elmt_ne_max = max (elmt_ne_max, E -> nExt  ());
    elmt_ni_max = max (elmt_ne_max, E -> nInt  ());
  }

  data = rvector (n_data);
  mesh = rvector (n_mesh);
  n_elmt_inodes = n_data - n_elmt_bnodes;

  // -- Install mesh storage and compute mesh node locations in each element.

  register real* m      = mesh;
  register int   offset = 0;

  for (k = 0; k < N; k++) {
    E = Elmt[k];

    E -> install      (offset, m, 0, 0);
    E -> mesh         (M);
    E -> map          ();
    E -> economizeGeo ();

    offset += E -> nTot ();
    m      += E -> nMsh ();
  }
}


Field::Field (const Field& master) :
  Elmt          (0),
  data          (0),
  n_data        (master.n_data),
  n_mesh        (master.n_mesh),
  n_elmt_bnodes (master.n_elmt_bnodes),
  n_elmt_inodes (master.n_elmt_inodes),
  elmt_np_max   (master.elmt_np_max),
  elmt_nt_max   (master.elmt_nt_max),
  elmt_ne_max   (master.elmt_ne_max),
  elmt_ni_max   (master.elmt_ni_max),
  mesh          (master.mesh)
// ---------------------------------------------------------------------------
// Copy constructor: new Field has its own storage area, Elements,
// shares master's mesh geometry.
// ---------------------------------------------------------------------------
{
  Elmt.setSize   (master.nEl ());
  data = rvector (n_data);
  setName ();

  register int k;
  const int    N = master.nEl ();

  for (k = 0; k < N; k++)
    Elmt[k] = (new Element (*master.Elmt[k]));
}


Field& Field::operator = (const real& val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (n_data,      data, 1);
  else              Veclib::fill (n_data, val, data, 1);

  return *this;
}


Field& Field::operator += (const real& val)
// ---------------------------------------------------------------------------
// Add val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (n_data, val, data, 1, data, 1);

  return *this;
}


Field& Field::operator -= (const real& val)
// ---------------------------------------------------------------------------
// Add -val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (n_data, -val, data, 1, data, 1);

  return *this;
}


Field& Field::operator *= (const real& val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (n_data, data, 1);
  else              Blas::scal   (n_data, val, data, 1);

  return *this;
}


Field& Field::operator /= (const real& val)
// ---------------------------------------------------------------------------
// Divide field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) message ("Field::op /= (real)", "divide by zero", ERROR);
  else              Blas::scal (n_data, 1.0/val, data, 1);

  return *this;
}


Field& Field::operator = (const Field& f)
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


Field& Field::operator -= (const Field& f)
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


Field& Field::operator *= (const Field& f)
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


Field& Field::operator /= (const Field& f)
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


Field& Field::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set Field's value to temporo-spatially varying function.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      k, offset = 0;
  const int         N = nEl ();

  for (k = 0; k < N; k++) {
    E      = Elmt[k];
    offset = E -> nOff ();
    E -> evaluate (function, data + offset);
  }

  return *this;
}


Field& Field::prod (const Field& a,
		    const Field& b)
// ---------------------------------------------------------------------------
// Set this Field's value as the product of a & b.
// ---------------------------------------------------------------------------
{
  if (a.n_data != n_data || b.n_data != n_data)
    message ("Field::prod", "size mismatch", ERROR);

  Veclib::vmul (n_data, a.data, 1, b.data, 1, data, 1);

  return *this;
}


Field& Field::sum (const Field& a,
		   const Field& b)
// ---------------------------------------------------------------------------
// Set this Field's value as the sum of a & b.
// ---------------------------------------------------------------------------
{
  if (a.n_data != n_data || b.n_data != n_data)
    message ("Field::sum", "size mismatch", ERROR);

  Veclib::vadd (n_data, a.data, 1, b.data, 1, data, 1);

  return *this;
}


Field& Field::addprod (const Field& a,
		       const Field& b)
// ---------------------------------------------------------------------------
// Add the product of a & b to this Field.
// ---------------------------------------------------------------------------
{
  if (a.n_data != n_data || b.n_data != n_data)
    message ("Field::addprod", "size mismatch", ERROR);

  Veclib::vvtvp (n_data, a.data, 1, b.data, 1, data, 1, data, 1);

  return *this;
}


Field& Field::axpy (const real&  alpha,
		    const Field& x    )
// ---------------------------------------------------------------------------
// Add alpha * x to this Field.
// ---------------------------------------------------------------------------
{
  if (x.n_data != n_data)
    message ("Field::axpy", "size mismatch", ERROR);

  Blas::axpy (n_data, alpha, x.data, 1, data, 1);

  return *this;
}


Field& Field::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (n_data, data, 1, data, 1);

  return *this;
}


void Field::swapData (Field* x,
		      Field* y)
// ---------------------------------------------------------------------------
// Swap data areas of two fields.
// ---------------------------------------------------------------------------
{
  if (x -> n_data != y -> n_data)
    message ("Field::swapData", "size mismatch", ERROR);

  real* tmp = x -> data;
  x -> data = y -> data;
  y -> data = tmp;
}


void Field::printMesh (const Field* F)
// ---------------------------------------------------------------------------
// Static member function.
//
// Mesh location information is written out element-by-element.
// ---------------------------------------------------------------------------
{
  cout << F -> Elmt[0] -> nKnot() << " "
       << F -> Elmt[0] -> nKnot() << " 1 " 
       << F -> nEl() << " NR NS NZ NEL" << endl;

  int k, N = F -> nEl ();

  for (k = 0; k < N; k++)
    F -> Elmt[k] -> printMesh ();
}


Field& Field::grad (const int& flag1, 
		    const int& flag2)
// ---------------------------------------------------------------------------
// Operate on Field(s) to produce the nominated index of the gradient.
// ---------------------------------------------------------------------------
{
  char routine[] = "Field::grad";

  register Element* E;
  register int      k, offset;
  const int         N = nEl ();

  if (flag1 && flag2) 
    message (routine, "can't make both components simultaneously", ERROR);
  else if (flag1)
    for (k = 0; k < N; k++) {
      E      = Elmt[k];
      offset = E -> nOff ();
      E -> grad (data + offset, 0);
    }
  else if (flag2)
    for (k = 0; k < N; k++) {
      E      = Elmt[k];
      offset = E -> nOff ();
      E -> grad (0, data + offset);
    }

  return *this;
}


Field& Field::errors (const char* function)
// ---------------------------------------------------------------------------
// Compare F with function, print the infinity-norm Li, the 2-norm L2
// and the Sobolev 1-norm H1.
//
// The norms are found element-by-element, using projection onto higher-order
// elements and high-order quadrature.
// ---------------------------------------------------------------------------
{
  const int NQUAD = 15;

  Element*  E;
  Element*  P;
  real      area = 0.0;
  real      Li   = 0.0;
  real      L2   = 0.0;
  real      H1   = 0.0;
  real*     sol;
  real*     v;
  real*     x;
  int       k, ntot, nmsh;
  const int N = nEl ();

  for (k = 0; k < N; k++) {
    E = Elmt[k];
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
  ostrstream (s, StrMax) << "Field '"
                         << getName ()
                         << "' error norms (inf, L2, H1): "
			 << Li << "  " << L2 << "  " << H1 << ends;
  message ("", s, REMARK);

  return *this;
}


ostream& operator << (ostream& strm,
		      Field&   F   )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  strm.write ((char*) F.data, F.n_data * sizeof (real));

  return strm;
}


istream& operator >> (istream& strm,
		      Field&   F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.
// ---------------------------------------------------------------------------
{
  strm.read ((char*) F.data, F.n_data * sizeof (real));

  return strm;
}


void Field::describe (char* s)  const
// ---------------------------------------------------------------------------
// Load s with a (prism-compatible) description of field geometry:
// NR NS NZ NEL.
// ---------------------------------------------------------------------------
{
  ostrstream (s, StrMax) << elmt_np_max << " "
                         << elmt_np_max << " "
                         << 1           << " "
                         << nEl ()      << ends;
}
