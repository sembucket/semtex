#ifndef CONDITION_H
#define CONDITION_H

class Condition
// ===========================================================================
// Virtual base class for boundary condition application.
//
// Each concrete class is derived from the virtual base class Condition:
// 1. Essential         : essential BC with constant, supplied, value.
// 2. EssentialFunction : essential BC, value obtained by parsing a function.
// 3. Natural           : natural BC with constant, supplied, value.
// 4. NaturalFunction   : natural BC, value obtained by parsing a function.
// 5. NaturalHOPBC      : "high-order" pressure BC, natural, computed value.
// 6. Mixed             : transfer coefficient type, 2 supplied values.
//
// Note that for supplied-value BC types, the value is a physical-space
// description of the BC, and is now set once at run-time (cannot be reset).
// ===========================================================================
{
public:
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)       const=0;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                      const=0;
  virtual void sum       (const int_t, const int_t*,
		          const real_t*,const real_t*,real_t*,real_t*) const=0;
  virtual void augmentSC (const int_t,  const int_t, const int_t,
			  const int_t*,const real_t*,real_t*,real_t*)  const=0;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)       const=0;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                      const=0;
  virtual void describe  (char* tgt)                                   const=0;

  virtual ~Condition()   { }
};


class Essential : public Condition
// ===========================================================================
// Essential BC applicator.  This one is for plain (constant value)
// Dirichlet/essential BCs.
// ===========================================================================
{
public:
  Essential              (const char* v) : _value (strtod (v, 0)) { }
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const;
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const
    { };
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  real_t _value;
};


class EssentialFunction : public Condition
// ===========================================================================
// Essential BC applicator for specified function Dirichlet/essential BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  EssentialFunction      (const char*);
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const;
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const
    { };
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  char* _function;
};


class Natural : public Condition
// ===========================================================================
// Natural BC applicator.  This one is for plain (constant value)
// Neumann/natural BCs.
// ===========================================================================
{
public:
  Natural                (const char* v) : _value (strtod (v, 0)) { }
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  real_t _value;
};


class NaturalFunction : public Condition
// ===========================================================================
// Natural BC applicator for specified function Neumann/natural BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  NaturalFunction        (const char*);
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  char* _function;
};


class NaturalHOPBC : public Condition
// ===========================================================================
// High-order pressure BC.  This is a natural BC, which has its value
// set at each timestep by an extrapolative process.  The value of dP/dn
// is obtained from the N-S equations.
// ===========================================================================
{
public:
  NaturalHOPBC           () { }
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void describe  (char*)                                         const;
};


class Mixed : public Condition
// ===========================================================================
// Boundary condition class for mixed type BCs of form dc/dn + K(c - C) = 0.
// Mixed BCs affect problem Helmholtz matrices, but only on the diagonal, 
// and element-boundary, terms.
// ===========================================================================
{
public:
  Mixed                  (const char*);
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*, const real_t*,
			  const real_t*, real_t*, real_t*)               const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void describe  (char*)                                         const;
private:
  real_t _K_;		// -- This is "K" above.
  real_t _C_;		// -- This is "C" above.
};

#endif
