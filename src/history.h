#ifndef HISTORY_H
#define HISTORY_H


class HistoryPoint
// ===========================================================================
// Class used to provide x,y,z history information.
// ===========================================================================
{
public:
  HistoryPoint (const integer id, const Element* e, const real r, 
		const real s, const real x, const real y, const real z):
    _id (id), _E (e), _r (r), _s (s), _x (x), _y (y), _z (z) { }

  integer               ID      () const { return _id; } 
  void                  extract (vector<AuxField*>&, real*) const;
  static const Element* locate  (const real, const real,
				 vector<Element*>&, real&, real&);

private:
  const integer  _id;		// Numeric identifier.
  const Element* _E ;		// Pointer to element.
  const real     _r ;		// Canonical-space r-location.
  const real     _s ;		// Canonical-space s-location.
  const real     _x ;		// x location.
  const real     _y ;		// y location.
  const real     _z ;		// Location in homogeneous direction.
};

#endif
