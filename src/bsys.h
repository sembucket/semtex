#ifndef BSYS_H
#define BSYS_H

class BoundarySys
// ===========================================================================
// This class automates the retrieval of the boundary condition
// applicators (Boundary objects), global numbering schemes
// (NumberSys) and inverse mass matrix for a given Field and Fourier
// mode.
// ===========================================================================
{
public:
  BoundarySys (BCmgr*, const vector<Element*>&, const char);
  ~BoundarySys () { };

  char                     field () const { return field_name; }
  integer                  nSurf () const { return nbound; }
  integer                  mixBC () const { return mixed; }
  const vector<Boundary*>& BCs   (const integer) const;
  const NumberSys*         Nsys  (const integer) const;
  const real*              Imass (const integer) const;

private:
  char               field_name;
  integer            nbound    ;  // Number of element edges with BCs.
  integer            mixed     ;  // Flags presence of mixed BC type.
  vector<Boundary*>* boundary  ;  // Boundary*'s           for modes 0, 1, 2.
  NumberSys**        number    ;  // NumberSys*'s          for modes 0, 1, 2.

  void buildbcs (const BCmgr*, const vector<Element*>&);
};

#endif
