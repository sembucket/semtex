#ifndef EDGE_H
#define EDGE_H

class Edge
// ===========================================================================
// Element-edge data and operators.
// ===========================================================================
{
public:
  Edge (const char*, const Element*, const integer); 

  integer bOff  () const { return _elmt -> ID() * Geometry::nExtElmt(); }
  integer dOff  () const { return _doffset; }
  integer dSkip () const { return _dskip; }

  void   get       (const real*, real*)                                 const;
  void   geometry  (real*, real*, real* = 0, real* = 0, real* = 0)      const;
  void   curlCurl  (const integer,
		    const real*, const real*, const real*,
		    const real*, const real*, const real*,
		    real*, real*, real*, real*, real*)                  const;
  void   mulY      (real*)                                              const;
  void   divY      (real*)                                              const;

  void   addForGroup (const char*, const real,  real*)                  const;
  void   setForGroup (const char*, const real,  real*)                  const;

  real   normalFlux      (const char*, const real*, const real*, real*) const;
  real   gradientFlux    (const char*, const real*, real*)              const;
  Vector normalTraction  (const char*, const real*, real*)              const;
  Vector tangentTraction (const char*, const real*, const real*, real*) const;

protected:
  integer        _np     ;	// Matches Geometry::nP().
  char*          _group  ;	// Group string.

  const Element* _elmt   ;	// Corresponding element.
  integer        _side   ;	// Corresponding side.

  integer        _doffset;	// Offset in Field data plane (matches BLAS).
  integer        _dskip  ;	// Skip   in Field data plane (matches BLAS).

  real*          _x      ;	// Node locations.
  real*          _y      ;	//
  real*          _nx     ;	// Unit outward normal components at nodes.
  real*          _ny     ;	// 
  real*          _area   ;	// Weighted multiplier for parametric mapping.
};

#endif
