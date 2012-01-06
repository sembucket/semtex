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
// 7. Convective        : convective-type mixed BC with supplied velocity.
//
// Note that for supplied-value BC types, the value is a physical-space
// description of the BC, and is now set once at run-time (cannot be reset).
// This is not true for those obtained by parsing a function: this can
// be re-parsed every timestep (if so flagged by command-line option
// -t or -tt).
//
// Also, each condition class derived from the base has to define all the
// pure virtual functions listed below (except the destructor) but
// some of these will just be stubs that do nothing for any particular
// type.  Those stubs are indicated in the present header
// file with the function body "{ };".
//
// For essential/Dirichlet BCs, the method "set" must be defined;
// For natural/Neumann BCs, the method "sum" must be defined, while
// For mixed/Robin BCs, the three methods "augmentSC", "augmentOp", 
// and "augmentDG" must be defined.
// For convective BCs, (a type of mixed BC), the method "extract"
// which loads field data into field boundary storage, must also be
// defined (in addition to the "augmentXx" methods).
//
// See also boundary.h, edge.h, mesh.C.
//
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                        const=0;
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                          const
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                          const
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                          const
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                          const
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                          const
    { };
  virtual void describe  (char*)                                         const;
};


class Mixed : public Condition
// ===========================================================================
// Boundary condition class for mixed (a.k.a Robin) type BCs of form
//     dc/dn + K(c - C) = 0.
// Mixed BCs affect problem Helmholtz matrices, but only on the diagonal, 
// and element-boundary, terms. Syntax in session file is
//     <M> c = K, C </M>  or 
//     <M> c = K; C </M> 
// where 'c' is a field name and K and C can be evaluated as constants
// (perhaps using defined TOKENS). White space following the separators
// above (',', ';') preceding C is optional.
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                          const
    { };
  virtual void describe  (char*)                                         const;
private:
  real_t _K_;		// -- This is "K" above.
  real_t _C_;		// -- This is "C" above.
};


class Convective : public Condition
// ===========================================================================
// Mixed-type boundary condition class for mixed type BCs of form
//     dc/dt + V dc/dn   = 0, which is rearranged as
//     dc/dn + 1/V dc/dt = 0 and then has dc/dt approximated by backwards
// Euler: dc/dn^{n+1} + 1/(V D_T) (c^{n+1} - c^{n}) = 0.  This is a type of
// mixed BC where c^{n} is taken from past timestep data and V is a supplied
// convection speed.
//
// Syntax in session file is
//     <C> c = V </C>  or 
// where 'c' is a field name and V can be evaluated as a constant
// (perhaps using defined TOKENS).
//
// Note that evaluate does nothing; its action is replaced by extract.
// ===========================================================================
{
public:
  Convective             (const char*);
  virtual void evaluate  (const int_t, const int_t, const int_t,
			  const Element*, const int_t, const int_t,
			  const real_t*, const real_t*, real_t*)         const
    { };
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
  virtual void extract   (const int_t, const real_t*, const int_t,
			  const int_t, real_t*)                          const;
  virtual void describe  (char*)                                         const;
private:
  real_t _K_;		// -- This is "1/(V * D_T)" above.
};

#endif
