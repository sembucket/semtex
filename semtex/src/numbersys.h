#ifndef NUMBERSYS_H
#define NUMBERSYS_H

class NumberSys
// ===========================================================================
// This class is a holder for Field Element-boundary numbering
// schemes.  Different Fields can potentially have unique
// NumberSyss, or they may share them with other Fields, according
// to distribution of essential BCs.
//
// The numbering system is associated with the Geometry class, but
// augments it by allowing for connectivity and boundary conditions.
//
// Bmask and emask provide a hierarchy of descriptions of element-
// boundary-node essential boundary condition masks:
// 1. bmask is a vector of ints, 4*(np-1)*nel (i.e. nbndry) in length, which
//    describe the CCW traverse of element-boundary nodes.  Values
//    set to 1 indicate that the node has an imposed (essential) BC, values
//    set to 0 indicate that the node either has a natural or no BC.
// 2. emask is the next level of the hierarchy, a vector of ints nel long.
//    Values are set to 1 if the corresponding element has any bmask
//    values set to 1, otherwise set to 0.
//
// Nglobal specifies the number of global nodes, while nsolve ( <=
// nglobal) specifies the number of global nodes which have bmask
// values of zero, i.e. where the nodal values will be solved for
// instead of set explicitly.
//
// Nglobal and nsolve also effectively provide a top level in the
// bmask/emask hierarchy, since if their difference is non-zero then
// at least one bmask value (and emask value) will be 1.
// 
// Btog gives global node numbers to Element-boundary nodes, same
// length as bmask (nbndry).  The optimization level describes the
// scheme used to build btog.
//
// Inv_mass is the inverse of the values of the global mass matrix,
// but on Element boundaries only.  Used for Field smoothing
// operations after derivative operators (if required).  Length =
// nglobal.
//
// A NumberSys is uniquely identified by the entries in btog, or
// equivalently the entries of bmask, together with optimization
// level.
// ===========================================================================
{
friend class BCmgr;
public:
 ~NumberSys () { }; 

  integer        nGlobal () const { return _nglobal; }
  integer        nSolve  () const { return _nsolve;  }
  integer        nBand   () const { return _nbandw;  }

  const char*    fields  () const { return _fields;            }
  const integer* bmask   () const { return _bmask;             }
  const integer* emask   () const { return _emask;             }
  integer        fmask   () const { return _nglobal - _nsolve; }
  const integer* btog    () const { return _btog;              }
  const real*    imass   () const { return _imass;             }

private:
  integer  _optlev ;		// Optimization level used for btog.
  integer  _nglobal;		// Length of inv_mass.
  integer  _nsolve ;		// Number of non-masked global nodes.
  integer  _nbandw ;		// Bandwidth of btog (includes diagonal).

  char*    _fields;		// String with character labels for Fields.
  integer* _bmask ;		// 1 for essential-BC nodes, 0 otherwise.
  integer* _emask ;		// 1 if associated Element has any esstl set.
  integer* _btog  ;		// Gives numbers to all element-boundary knots.
  real*    _imass ;		// Inverse of global mass matrix;
};

#endif
