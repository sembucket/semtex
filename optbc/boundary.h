#ifndef BOUNDARY_H
#define BOUNDARY_H

class Boundary : public Edge
// ===========================================================================
// Physical field element-wise boundary class.
// ===========================================================================
{
public:
  Boundary (const int_t, const char*, const Condition*,
	    const Element*, const int_t);

  int_t ID        () const { return _id; }
  void  print     () const;

  void  evaluate  (const int_t,const int_t,real_t*)                      const;

  // -- Impose essential BCs:
  void  set       (const real_t*,const int_t*,real_t*)                   const;

  // -- Retrieve essential BCs from globally numbered vector.
  void  get       (const real_t*,const int_t*,real_t*)                   const;

  // -- Apply natural BCs:
  void  sum       (const real_t*,const int_t*,real_t*,real_t*)           const;

  // -- Apply mixed BCs:
  void  augmentSC (const int_t,const int_t,const int_t*,real_t*,real_t*) const;
  void  augmentOp (const int_t*,const real_t*,real_t*)                   const;
  void  augmentDg (const int_t*,real_t*)                                 const;

  const Condition* bcond() const {return _bcond;}      

  real_t  controlnorm (real_t*)                                      const;
  real_t  controlnorm_mixed (real_t*, real_t*)                       const;
  real_t  controllength ()                                           const;
  void    switchK   (const real_t*,const real_t*, const bool)        const; 
  void    takeC     (const real_t* ) const;
  
  real_t* uc;

private:
  int_t            _id   ;	// Ident number.
  const Condition* _bcond;	// Boundary condition.
};

#endif
