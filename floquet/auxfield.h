#ifndef AUXFIELD_H
#define AUXFIELD_H


class AuxField
// ===========================================================================
// Physical field storage class, no BCs.
//
// An AuxField has a storage area for data, physical space mesh, and a
// matching list of elements which can access this storage area on an
// element-wise basis if required.  Functions to operate on fields as
// a whole are included.  AuxFields have no boundary conditions, no
// global numbering system, and no solution routines.
// ===========================================================================
{
friend istream& operator >> (istream&, AuxField&);
friend ostream& operator << (ostream&, AuxField&);
friend class     Field;
friend class     PBCmgr;

public:
  AuxField (real*, const integer, vector<Element*>&, const char = 0);

  char name     ()      const { return _name; }
  void describe (char*) const;

  AuxField& operator  = (const real);
  AuxField& operator += (const real);
  AuxField& operator -= (const real);
  AuxField& operator *= (const real);
  AuxField& operator /= (const real);

  AuxField& operator  = (const AuxField&);
  AuxField& operator += (const AuxField&);
  AuxField& operator -= (const AuxField&);
  AuxField& operator *= (const AuxField&);

  AuxField& operator  = (const char*);
  AuxField& axpy        (const real, const AuxField&);

  AuxField& times       (const AuxField&, const AuxField&);
  AuxField& timesPlus   (const AuxField&, const AuxField&);
  AuxField& timesMinus  (const AuxField&, const AuxField&);

  AuxField& addToPlane  (const integer, const real);
  AuxField& getPlane    (const integer, real*);
  AuxField& setPlane    (const integer, const real*);
  AuxField& setPlane    (const integer, const real);

  AuxField& update   (const integer, const real*, const real, const real); 

  AuxField& gradient (const integer);
  AuxField& mulY     ();
  AuxField& divY     ();

  real mode_L2     (const integer) const;
  real integral    ()          const;
  real CFL         (const integer) const;

  real probe (const Element*, const real, const real, const integer)  const;
  real probe (const Element*, const real, const real, const real) const;

  AuxField& reverse     ();
  static void swapData  (AuxField*, AuxField*);
  static void couple    (AuxField*, AuxField*, const integer);

protected:
  char              _name ;	// Identification tag.  '\0' by default.
  vector<Element*>& _elmt ;	// Quadrilateral elements.
  integer           _nz   ;	// number of data planes (per process).
  integer           _size ;	// _nz * Geometry::planeSize().
  real*             _data ;	// 2/3D data area, element x element x plane.
  real**            _plane;	// Pointer into data for each 2D frame.

private:
  AuxField& operator /= (const AuxField&) { return *this; }

};

#endif