#ifndef BOUNDARY_H
#define BOUNDARY_H


class Boundary : public Edge
// ===========================================================================
// Physical field element-wise boundary class.
// ===========================================================================
{
public:
  Boundary (const integer id, const char* group, const Condition* bcond,
	    const Element* elmt, const integer side):
    Edge (group, elmt, side), _id (id), _bcond (bcond) { }

  integer ID       () const { return _id; }
  void    print    () const;

  void    evaluate (const integer, const integer,  real*)               const;

  // -- Impose essential BCs:
  void   set       (const real*, const integer*, real*)                 const;
  // -- Apply natural BCs:
  void   sum       (const real*, const integer*, real*, real*)          const;
  // -- Apply mixed BCs:
  void   augmentSC (const integer, const integer, const integer*,
		    real*, real*)                                       const;
  void   augmentOp (const integer*, const real*, real*)                 const;
  void   augmentDg (const integer*, real*)                              const;

private:
  integer          _id   ;	// Ident number.
  const Condition* _bcond;	// Boundary condition.
};

#endif
