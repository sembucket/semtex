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

  char                     field () const { return _field_name; }
  int_t                    nSurf () const { return _nbound; }
  int_t                    sizecontrolbc() { return _uc; }
  int_t                    mixBC () const { return _mixed; }
  int_t                    ToutflowBC () const { return _Toutflow; }
  const vector<Boundary*>& BCs   (const int_t) const;
  const NumberSys*         Nsys  (const int_t) const;
  const real*              Imass (const int_t) const;

private:
  char               _field_name;
  int_t              _nbound    ;  // Number of element edges with BCs.
  int_t              _uc        ;  // Number of control boundary nodes in one plane.
  bool               _mixed     ;  // Flags presence of mixed BC type.
  bool               _Toutflow  ;  // Flags presence of Toutflow BC type.
  vector<Boundary*>* _boundary  ;  // Boundary*'s  for modes 0, 1, 2. and adjoint 0, 1, 2 
  NumberSys**        _number    ;  // NumberSys*'s for modes 0, 1, 2. and adjoint 0, 1, 2

  void buildbcs (const BCmgr*, const vector<Element*>&);
};

#endif
