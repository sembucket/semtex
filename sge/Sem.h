#ifndef SEM_H
#define SEM_H
///////////////////////////////////////////////////////////////////////////////
// Sem.h: main header file for semtex spectral element solvers.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
//
// Conventions: 
// 1. Arrays are 0-offset.
// 2. Internal ident numbers id/ID start at 0.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <stdarg.h>		/* System C headers.  */
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include <iostream.h>		/* System C++ headers. */
#include <fstream.h>
#include <strstream.h>
#include <iomanip.h>

#include <femdef.h>		/* Semtex headers.     */
#include <List.h>
#include <Stack.h>
#include <Array.h>
#include <Blas.h>
#include <Lapack.h>
#include <Utility.h>
#include <Veclib.h>
#include <Femlib.h>

#include <Feml.h>
#include <Geometry.h>
#include <Mesh.h>

#define ROOTONLY if (Geometry::procID() == 0)
#define VERBOSE  ROOTONLY if (verbose)

class Element;
class Boundary;
class BoundarySys;
class AuxField;
class Field;
class Domain;
class BCmgr;
class PBCmgr;
class Statistics;
class MixPatch;


typedef struct patch {
  const Boundary* b1;
  const Boundary* b2;
  real*           x;
  real*           y;
  real*           val1;
  real*           val2;
} Patch;


class MixPatch
// ============================================================================
// Special class to implement patched BCs over 2 domains.
// ============================================================================
{
public:
  MixPatch::MixPatch (const Domain*, const Domain*);
  integer      np;
  integer      nz;
  List<Patch*> patches;
};


class Element
// ===========================================================================
// 2D quadrilateral element class with equal order in each direction.
//                                                
//   s ^                                    y ^     +4
//     |                                      |    / \
//     4      3                               |  1+   +3
//     +------+                               |    \  |     
//     |      |      <== Logical  space       |     \ |
//     |      |          Physical space ==>   |      \|
//     |      |                               |       +2  
//     +------+-> r                           +-----------> x
//     1      2
//
// Master element coordinates: -1 < r,s < 1; edge traverses CCW.
// ===========================================================================
{
private:
  const integer id  ;		// Element identifier.
  const integer np  ;		// Number of points on an edge.
  const integer npnp;		// Total number = np * np.
  const integer next;		// Number of points on periphery.
  const integer nint;		// Number of internal points.

  integer* emap;		// Indices of edges in nodal matrices.
  integer* pmap;		// Inversion of emap (pmap[emap[i]] = i).

  real*    xmesh;		// Physical space mesh.
  real*    ymesh;		// 2D row-major store.

  real*    drdx ;		// Partial derivatives (x, y) --> (r, s),
  real*    dsdx ;		//   evaluated at quadrature points.
  real*    drdy ;		//   (2D row-major storage.)
  real*    dsdy ;		//
  real*    G1   ;		// Geometric factor 1 at quadrature points.
  real*    G2   ;		// Geometric factor 2 at quadrature points.
  real*    G3   ;		// Geometric factor 3 at quadrature points.
  real*    G4   ;		// Geometric factor 4 at quadrature points.
  real*    delta;		// Local length scale.

  void map();
  void terminal (const integer side ,
		 integer&      start,
		 integer&      skip ) const
    // -- BLAS-conforming edge offsets & skips for element-edge traverses.
    {
      switch (side) {
      case 0: start = 0;           skip  = 1;   break;
      case 1: start = np - 1;      skip  = np;  break;
      case 2: start = np*(np - 1); skip  = -1;  break;
      case 3: start = 0;           skip  = -np; break;
      }
    }

  void helmRow    (const real**, const real**, const real, const real,
		   const integer, const integer, real*, real*)         const;
  void printMatSC (const real*, const real*, const real*)               const;

public:
  Element (const integer, const Mesh*, const real*, const integer);
 ~Element ();

  integer ID () const { return id; }

  // -- Elemental Helmholtz matrix construction, manipulation, operator.

  void HelmholtzSC   (const real,const real,real*,real*,
		      real*,real*,real*,integer*)                        const;
  void HelmholtzDg   (const real,const real,real*,real*)                 const;
  void Helmholtz     (const real,const real,real*,real*,real*)           const;
  void HelmholtzKern (const real, const real, real*, real*, real*,real*) const;

  // -- Element-boundary operators.

  void bndryDsSum  (const integer*, const real*, real*)                 const;
  void bndryMask   (const integer*, real*, const real*, const integer*) const;
  void bndryInsert (const integer*, const real*, real*)                 const;

  // -- Element-side operators.

  void sideGeom  (const integer, real*, real*, real*, real*, real*) const;
  void sideEval  (const integer, real*, const  char*)               const;
  void sideGrad  (const integer, const real*, real*, real*)         const;
  void sideGet   (const integer, const real*, real*)                const;
  void sideGetR  (const integer, real*)                             const;
  void sideDivR  (const integer, const real*, real*)                const;
  void sideDivR2 (const integer, const real*, real*)                const;

  // -- Element-operator functions.

  void grad  (real*, real*, const real**, const real**, real*) const;
  void gradX (const real* xr, const real* xs, real* dx)        const;
  void gradY (const real* yr, const real* ys, real* dy)        const;

  void divR (real*) const;
  void mulR (real*) const;

  void evaluate (const char*, real*) const;

  real integral (const char*)        const;
  real integral (const real*, real*) const;
  real area     ()                   const;
  void weight   (real*)              const;

  void lengthScale (real*)                                const;
  real CFL         (const real, const real*, const real*) const;
       
  real norm_inf (const real*) const;
  real norm_L2  (const real*) const;
  real norm_H1  (const real*) const;
       
  void e2g      (const real*, const integer*, real*, real*)       const;
  void e2gSum   (const real*, const integer*, real*, real*)       const;
  void g2e      (real*, const integer*, const real*, const real*) const;
  void g2eSC    (const real*, const integer*, real*,
		 real*, const real*, const real*, real*)          const;
  void e2gSumSC (real*, const integer*, real*, const real*)       const;

  // -- Probe functions.

  integer locate (const real, const real, real&, real&, const integer=0) const;
  real    probe  (const real, const real, const real*)                   const;

  // -- Debugging/informational routines.

  void printMesh  ()            const;
  void printBndry (const real*) const;
};


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
friend ifstream& operator >> (ifstream&, AuxField&);
friend ofstream& operator << (ofstream&, AuxField&);
friend class     Field;
friend class     PBCmgr;

public:
  AuxField (vector<Element*>&, const char = 0);

  char name     () const { return field_name; }
  void describe (char*) const;

  AuxField& operator  = (const real);
  AuxField& operator += (const real);
  AuxField& operator -= (const real);
  AuxField& operator *= (const real);
  AuxField& operator /= (const real);

  AuxField& operator  = (const AuxField&);
  AuxField& operator += (const AuxField&);
  AuxField& operator -= (const AuxField&);

  AuxField& operator  = (const char*);
  AuxField& axpy        (const real, const AuxField&);

  AuxField& innerProduct (const vector<AuxField*>&, const vector<AuxField*>&);
  AuxField& times        (const AuxField&, const AuxField&);
  AuxField& timesPlus    (const AuxField&, const AuxField&);
  AuxField& timesMinus   (const AuxField&, const AuxField&);

  AuxField& transform   (const integer);
  AuxField& transform32 (const integer, real*);
  AuxField& DLT2D       (const integer, const real* = 0);
  AuxField& DFfilt      (const real*);
  AuxField& addToPlane  (const integer, const real);
  AuxField& getPlane    (const integer, real*);
  AuxField& setPlane    (const integer, const real*);

  AuxField& gradient (const integer);
  AuxField& mulR     ();
  AuxField& divR     ();

  void gradient (const integer, const integer, real*, const integer) const;
  void mulR     (const integer, real*)                               const;
  void divR     (const integer, real*)                               const;

  void errors      (const Mesh*, const char*);
  void lengthScale (real*)                     const;
  real norm_inf    ()                          const;
  real mode_L2     (const integer)             const;
  real integral    ()                          const;
  real CFL         (const integer)             const;

  real probe (const Element*, const real, const real, const integer) const;
  real probe (const Element*, const real, const real, const real)    const;

  AuxField& reverse     ();
  AuxField& zeroNyquist ();

  static void swapData  (AuxField*, AuxField*);
  static void couple    (AuxField*, AuxField*, const integer);

protected:
  char              field_name ; // Identification tag.  '\0' by default.
  
  vector<Element*>& Elmt       ; // Quadrilateral elements.
  real*             data       ; // 2/3D data area, element x element x plane.
  real**            plane      ; // Pointer into data for each 2D frame.

private:
  AuxField& operator /= (const AuxField&) { return *this; }
  AuxField& operator *= (const AuxField&) { return *this; }
};


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
  virtual void evaluate  (const integer, const integer, const integer,
			  const Element*, const integer, const integer,
			  const real*, const real*, real*)           const = 0;
  virtual void set       (const integer, const integer*,
			  const real*, real*)                        const = 0;
  virtual void sum       (const integer, const integer*,
		          const real*, const real*, real*, real*)    const = 0;
  virtual void augmentSC (const integer, const integer, const integer,
			  const integer*, const real*, real*, real*) const = 0;
  virtual void augmentOp (const integer, const integer*,
			  const real*, const real*, real*)           const = 0;
  virtual void augmentDg (const integer, const integer*, 
			  const real*, real*)                        const = 0;
  virtual void describe  (char* tgt)                                 const = 0;

  const char* group      () const { return grp; }

  virtual ~Condition()   { }

protected:
  char* grp;
};


class Essential : public Condition
// ===========================================================================
// Essential BC applicator.  This one is for plain (constant value)
// Dirichlet/essential BCs.
// ===========================================================================
{
public:
  Essential              (const char*, const char*);
  virtual void evaluate  (const integer, const integer, const integer,
			  const Element*, const integer, const integer,
			  const real*, const real*, real*)              const;
  virtual void set       (const integer, const integer*,
			  const real*, real*)                           const;
  virtual void sum       (const integer, const integer*,
			  const real*, const real*, real*, real*)       const;
  virtual void augmentSC (const integer, const integer, const integer,
			  const integer*, const real*, real*, real*)    const;
  virtual void augmentOp (const integer, const integer*,
			  const real*, const real*, real*)              const;
  virtual void augmentDg (const integer, const integer*, 
			  const real*, real*)                           const;
  virtual void describe  (char*)                                        const;
private:
  real value;
};


class EssentialFunction : public Condition
// ===========================================================================
// Essential BC applicator for specified function Dirichlet/essential BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  EssentialFunction      (const char*, const char*);
  virtual void evaluate  (const integer, const integer, const integer,
			  const Element*, const integer, const integer,
			  const real*, const real*, real*)              const;
  virtual void set       (const integer, const integer*,
			  const real*, real*)                           const;
  virtual void sum       (const integer, const integer*,
			  const real*, const real*, real*, real*)       const;
  virtual void augmentSC (const integer, const integer, const intege,
			  const integer*, const real*, real*, real*)    const;
  virtual void augmentOp (const integer, const integer*,
			  const real*, const real*, real*)              const;
  virtual void augmentDg (const integer, const integer*, 
			  const real*, real*)                           const;
  virtual void describe  (char*)                                        const;
private:
  char* function;
};


class Natural : public Condition
// ===========================================================================
// Natural BC applicator.  This one is for plain (constant value)
// Neumann/natural BCs.
// ===========================================================================
{
public:
  Natural                (const char*, const char*);
  virtual void evaluate  (const integer, const integer, const integer,
			  const Element*, const integer, const integer,
			  const real*, const real*, real*)              const;
  virtual void set       (const integer, const integer*,
			  const real*, real*)                           const;
  virtual void sum       (const integer, const integer*,
			  const real*, const real*, real*, real*)       const;
  virtual void augmentSC (const integer, const integer, const intege,
			  const integer*, const real*, real*, real*)    const;
  virtual void augmentOp (const integer, const integer*,
			  const real*, const real*, real*)              const;
  virtual void augmentDg (const integer, const integer*, 
			  const real*, real*)                           const;
  virtual void describe  (char*)                                        const;
private:
  real value;
};


class NaturalFunction : public Condition
// ===========================================================================
// Natural BC applicator for specified function Neumann/natural BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  NaturalFunction        (const char*, const char*);
  virtual void evaluate  (const integer, const integer, const integer,
			  const Element*, const integer, const integer,
			  const real*, const real*, real*)              const;
  virtual void set       (const integer, const integer*,
			  const real*, real*)                           const;
  virtual void sum       (const integer, const integer*,
			  const real*, const real*, real*, real*)       const;
  virtual void augmentSC (const integer, const integer, const intege,
			  const integer*, const real*, real*, real*)    const;
  virtual void augmentOp (const integer, const integer*,
			  const real*, const real*, real*)              const;
  virtual void augmentDg (const integer, const integer*, 
			  const real*, real*)                           const;
  virtual void describe  (char*)                                        const;
private:
  char* function;
};


class NaturalHOPBC : public Condition
// ===========================================================================
// High-order pressure BC.  This is a natural BC, which has its value
// set at each timestep by an extrapolative process.  The value of dP/dn
// is obtained from the N-S equations.
// ===========================================================================
{
public:
  NaturalHOPBC           (const char*);
  virtual void evaluate  (const integer, const integer, const integer,
			  const Element*, const integer, const integer,
			  const real*, const real*, real*)              const;
  virtual void set       (const integer, const integer*,
			  const real*, real*)                           const;
  virtual void sum       (const integer, const integer*,
			  const real*, const real*, real*, real*)       const;
  virtual void augmentSC (const integer, const integer, const intege,
			  const integer*, const real*, real*, real*)    const;
  virtual void augmentDg (const integer, const integer*, 
			  const real*, real*)                           const;
  virtual void augmentOp (const integer, const integer*,
			  const real*, const real*, real*)              const;
  virtual void describe  (char*)                                        const;
private:

};


class Mixed : public Condition
// ===========================================================================
// Boundary condition class for mixed type BCs of form dc/dn + K(c - C) = 0.
// Mixed BCs affect problem Helmholtz matrices, but only on the diagonal, 
// and element-boundary, terms.
// ===========================================================================
{
public:
  Mixed                  (const char*, const char*);
  virtual void evaluate  (const integer, const integer, const integer,
			  const Element*, const integer, const integer,
			  const real*, const real*, real*)              const;
  virtual void set       (const integer, const integer*,
			  const real*, real*)                           const;
  virtual void sum       (const integer, const integer*,
			  const real*, const real*, real*, real*)       const;
  virtual void augmentSC (const integer, const integer, const intege,
			  const integer*, const real*, real*, real*)    const;
  virtual void augmentOp (const integer, const integer*,
			  const real*, const real*, real*)              const;
  virtual void augmentDg (const integer, const integer*, 
			  const real*, real*)                           const;
  virtual void describe  (char*)                                        const;
private:
  real K;		// -- This is "K" above.
  real C;		// -- This is "C" above.
};


class Boundary
// ===========================================================================
// Physical field element-wise boundary class.
// ===========================================================================
{
public:
  Boundary (const integer, const char*, const Condition*,
	    const Element*, const integer); 

  integer bOff  () const { return elmt -> ID() * Geometry::nExtElmt(); }
  integer dOff  () const { return doffset; }
  integer dSkip () const { return dskip; }

  integer ID       () const { return id; }
  void    print    () const;
  const char* group () const;
  void    geometry (real*, real*, real* = 0, real* = 0, real* = 0) const;


  void    evaluate  (const integer, const integer,  real*)        const;
  void    get       (const real*, real*)                          const;
  void    set       (const real*,   const integer*, real*)        const;
  void    sum       (const real*,   const integer*, real*, real*) const;
  void    augmentSC (const integer, const integer,
		     const integer*, real*, real*)                const;
  void    augmentOp (const integer*, const real*, real*)          const;
  void    augmentDg (const integer*, real*)                       const;
  void    curlCurl  (const integer,
		     const real*, const real*, const real*,
		     const real*, const real*, const real*,
		     real*, real*, real*, real*)                  const;

  void    addForGroup     (const char*, const real,  real*) const;
  void    setForGroup     (const char*, const real,  real*) const;
  real    flux            (const char*, const real*, real*) const;
  Vector  normalTraction  (const char*, const real*, real*) const;
  Vector  tangentTraction (const char*, const real*,
			   const real*, real*, real*)       const;

private:
  integer          id     ;	// Ident number.
  integer          np     ;	// Length of sides.
  const char*      bgroup ;	// Group string.
  const Condition* bcondn ;	// Boundary condition.

  const Element*   elmt   ;	// Corresponding element.
  integer          side   ;	// Corresponding side.

  integer          doffset;	// Offset in Field data plane (matches BLAS).
  integer          dskip  ;	// Skip   in Field data plane (matches BLAS).

  real*            x      ;	// Spatial location at nodes.
  real*            y      ;	// 
  real*            nx     ;	// Unit outward normal components at nodes.
  real*            ny     ;	// 
  real*            area   ;	// Weighted multiplier for parametric mapping.
};


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
  integer nGlobal () const { return ns_nglobal; }
  integer nSolve  () const { return ns_nsolve;  }
  integer nBand   () const { return ns_nbandw;  }

  const char*    fields () const { return ns_fields;              }
  const integer* bmask  () const { return ns_bmask;               }
  const integer* emask  () const { return ns_emask;               }
  const integer  fmask  () const { return ns_nglobal - ns_nsolve; }
  const integer* btog   () const { return ns_btog;                }
  const real*    imass  () const { return ns_imass;               }

private:
  integer  ns_optlev ;		// Optimization level used for btog.
  integer  ns_nglobal;		// Length of inv_mass.
  integer  ns_nsolve ;		// Number of non-masked global nodes.
  integer  ns_nbandw ;		// Bandwidth of btog (includes diagonal).

  char*    ns_fields ;		// String with character labels for Fields.
  integer* ns_bmask  ;		// 1 for essential-BC nodes, 0 otherwise.
  integer* ns_emask  ;		// 1 if associated Element has any esstl set.
  integer* ns_btog   ;		// Gives numbers to all element-boundary knots.
  real*    ns_imass  ;		// Inverse of global mass matrix;
};


typedef struct bctriple { char group; integer elmt; integer side; } BCtriple;


class BCmgr
// ===========================================================================
// This is a factory / retrieval service for classes derived from
// Condition, and maintains GROUP descriptors.  In addition, it reads
// and returns NumberSys objects from session.num.
// ===========================================================================
{
public:
  BCmgr (FEML*, vector<Element*>&);

  const char*      field        () const { return fields; }
  const char*      groupInfo    (const char);
  Condition*       getCondition (const char, const char, const integer = 0);
  NumberSys*       getNumberSys (const char, const integer = 0);
  List<BCtriple*>& getBCedges   () { return elmtbc; }
  integer          nBCedges     () const { return elmtbc.length(); }

  class CondRecd {
  public: 
    char        grp  ;
    char        fld  ;
    Condition*  bcn  ;
    char*       value;
  };
    
private:
  char*              fields  ;	// String containing field names.
  vector<char>       group   ;	// Single-character group tags.
  vector<char*>      descript;	// Group name strings.
  List<CondRecd*>    cond    ;	// Conditions in storage.
  List<BCtriple*>    elmtbc  ;	// Group tags for each element-side BC.
  vector<NumberSys*> numsys  ;	// Numbering schemes in storage.

  void buildnum  (const char*, vector<Element*>&);
  void buildsurf (FEML*, vector<Element*>&);
};


class PBCmgr
// ===========================================================================
// This class maintains internal storage for evaluation of "high-order"
// pressure boundary conditions & provides means for their evaluation.
// ===========================================================================
{
public:
  static void build      (const Field*);
  static void maintain   (const integer, const Field*, const AuxField***,
			  const AuxField***, const integer);
   static void evaluate   (const integer, const integer, const integer,
			  const integer, const real*, const real*, real*);
  static void accelerate (const Vector&, const Field*);

private:
  static real**** Pnx ;		// x component of dP / dn at domain  boundary.
  static real**** Pny ;		// y component of dP / dn at domain  boundary.
  static real**** Unx ;		// x component of normal velocity at boundary.
  static real**** Uny ;		// y component of normal velocity at boundary.
};


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

  const char               field () const { return field_name; }
  const integer            nSurf () const { return nbound; }
  const integer            mixBC () const { return mixed; }
  const vector<Boundary*>& BCs   (const integer = 0) const;
  const NumberSys*         Nsys  (const integer = 0) const;
  const real*              Imass (const integer = 0) const;

private:
  char               field_name;
  integer            nbound;	// Number of element edges with BCs.
  integer            mixed;	// Flags presence of mixed BC type.
  vector<Boundary*>* boundary;	// Boundary*'s           for modes 0, 1, 2.
  NumberSys**        number;	// NumberSys*'s          for modes 0, 1, 2.

  void buildbcs (const BCmgr*, const vector<Element*>&);
};


class MatrixSys
// ===========================================================================
// System of global and local Helmholtz matrices and partitions.
// Matrix factorisations are Cholesky, use LAPACK storage schemes.
// ===========================================================================
{
friend class Field;
friend ostream& operator << (ostream&, MatrixSys&);
friend istream& operator >> (istream&, MatrixSys&);
public:
  MatrixSys     (const real, const real, const integer,
		 const vector<Element*>&, const BoundarySys*);
  integer match (const real, const real,  const NumberSys*) const;

private:
  real HelmholtzConstant;	// Same for all modes.
  real FourierConstant  ;	// Varies with mode number.
 
  const vector<Boundary*>& BC;	// Internal copy of Boundary conditions.
  const NumberSys*         NS;	// Internal copy of NumberSys.

  integer  nel     ;		// Number of elemental matrices.
  integer  nband   ;		// Bandwidth of global matrix (incl. diagonal).
  integer  singular;		// If system is potentially singular.
  integer  nsolve  ;		// System-specific number of global unknowns.
  integer  npack   ;		// Number of reals for global matrix.

  real*    H     ;		// (Factored) packed global Helmholtz matrix.
  real**   hbi   ;		// Element external-internal coupling matrices.
  real**   hii   ;		// (Factored) internal-internal matrices.
  integer* bipack;		// Size of hbi for each element.
  integer* iipack;		// Size of hii for each element.
  
};

ostream& operator << (ostream&, MatrixSys&);
istream& operator >> (istream&, MatrixSys&);


class ModalMatrixSys
// ===========================================================================
// A way of organising MatrixSys*'s by Fourier mode.
// ===========================================================================
{
public:
  ModalMatrixSys (const real, const real, const integer, const integer,
		  const vector<Element*>&, const BoundarySys*);

  const MatrixSys* operator [] (const integer i) const { return Msys[i]; }

private:
  char*              fields;	// Character field tags for this system.
  vector<MatrixSys*> Msys  ;	// One MatrixSys for each Fourier mode.
};
  

class Field : public AuxField
// ===========================================================================
// Field adds boundary conditions and global numbering to AuxField.
// With the boundary conditions comes knowledge of which sections of
// the boundary hold essential (known) nodes, and which contain
// natural (unknown) nodes.
//
// A Field holds global node numbers and solve masks for element
// boundaries: where mesh value is given by an essential BC, solve
// mask is 0, for all other nodes have value 1.
//
// Helmholtz solution routines are also available.
// ===========================================================================
{
friend class PBCmgr;
public:
  Field  (BoundarySys*, vector<Element*>&, const char);
 ~Field  () { }

  Field& operator = (const AuxField& z) {AuxField::operator=(z); return *this;}
  Field& operator = (const real&     z) {AuxField::operator=(z); return *this;}
  Field& operator = (const char*     z) {AuxField::operator=(z); return *this;}

  Field& solve  (AuxField*, const ModalMatrixSys*);
  Field& solve  (AuxField*, const real);

  Field& smooth (AuxField* = 0);
  void   smooth (const int, real*) const;

  void evaluateBoundaries    (const integer);
  void evaluateM0Boundaries  (const integer);
  void addToM0Boundaries     (const real, const char*);
  
  void setPatch (MixPatch*) const;
  void getPatch (MixPatch*);

  static real   flux            (const Field*);
  static Vector normalTraction  (const Field*);
  static Vector tangentTraction (const Field*, const Field*, const Field* = 0);
  static void   normTractionV   (real*, real*, const Field*);
  static void   tangTractionV   (real*, real*, real*, const Field*,
				 const Field*, const Field* = 0);

  static void coupleBCs    (Field*, Field*, const integer);
  static real modeConstant (const char, const integer, const real);

  static void printBoundaries (const Field*);
  static void printConnect    (const Field*);

private:
  integer      nbound;		// Number of boundary edges.
  integer      n_line;		// Length of one boundary line.
  real*        sheet ;		// Wrap-around storage for data boundary.
  real**       line  ;		// Single plane's worth of sheet.
  BoundarySys* bsys  ;		// Boundary system information.

  void bTransform        (const integer);

  void getEssential      (const real*, real*,
			  const vector<Boundary*>&, const NumberSys*) const;
  void setEssential      (const real*, real*, const NumberSys*);

  void local2global      (const real*, real*, const NumberSys*)       const;
  void local2globalSum   (const real*, real*, const NumberSys*)       const;
  void global2local      (const real*, real*, const NumberSys*)       const;

  void constrain         (real*, const real, const real,
			  const real*, const NumberSys*)              const;
  void buildRHS          (real*, const real*, real*, real*,
			  const real**, const integer, const integer,
			  const vector<Boundary*>&, const NumberSys*) const;
  void HelmholtzOperator (const real*, real*, const real,
			  const real, real*, const integer)           const;
  void jacobi            (const real, const real, real*,
			  const NumberSys*)                           const;
};


class Domain
// ===========================================================================
// Physical domain storage class.
//
// Since users are expected to program with entities kept at the
// Domain level, all internal storage is exposed to view.
//
// Domain fields are written/read untransformed (in physical space).
// ===========================================================================
{
friend ifstream& operator >> (ifstream&, Domain&);
friend ofstream& operator << (ofstream&, Domain&);
public:
  Domain (FEML*, vector<Element*>&, BCmgr*);

  char*                 name;	// Session name.
  integer               step;	// Runtime step number.
  real                  time;	// Simulation time.
  vector<Element*>&     elmt;	// Shared for equal-order interpolations.
  vector<Field*>        u   ;	// Solution fields: velocities, pressure.
  vector<BoundarySys*>  b   ;	// Field boundary systems.

  integer  nField () const { return u.getSize(); }
  void report     ();
  void restart    ();
  void dump       ();
  void transform  (const integer);
  void setNumber  (const char, const NumberSys**) const;

private:
  char* field;		// Lower-case single character field names.
};


class Integration
// ===========================================================================
// Return coefficients for time-integration schemes.
// ===========================================================================
{
public:
  static const integer OrderMax;

  static void AdamsBashforth (const integer, real*);
  static void AdamsMoulton   (const integer, real*);
  static void StifflyStable  (const integer, real*);
  static void Extrapolation  (const integer, real*);
};


class FluidParticle
// ===========================================================================
// Class used to locate and integrate positions of massless particles.
// ===========================================================================
{
public:
  FluidParticle (Domain*, const integer, Point&);
  void           integrate (const integer); 
  integer        ID        () const { return id;     }
  const Element* inMesh    () const { return E;      }
  const Point&   location  () const { return P;      } 
  static integer IDMax     ()       { return ID_MAX; }

private:
  integer        id      ;	// Numeric tag.
  const Element* E       ;	// Pointer to the element particle is in.
  real           r       ;	// Corresponding "r" location within element.
  real           s       ;	// likewise for "s".
  Point          P       ;	// Physical space location.
  real*          u       ;	// Multilevel "x" velocity storage.
  real*          v       ;	// Multilevel "y" velocity storage.
  real*          w       ;	// Multilevel "z" velocity storage.

  static Domain* D       ;	// Velocity fields and class functions.
  static integer NDIM    ;	// Number of space dimensions.
  static integer NEL     ;	// Number of elements in mesh.
  static integer NZ      ;	// Number of z planes.
  static integer TORD    ;	// Order of N--S timestepping.
  static integer ID_MAX  ;	// Highest issued id.
  static real*   P_coeff ;	// Integration (predictor) coefficients.
  static real*   C_coeff ;	// Integration (corrector) coefficients.
  static real    DT      ;	// Time step.
  static real    Lz      ;	// Periodic length.
};


class HistoryPoint
// ===========================================================================
// Class used to provide x,y,z history information.
// ===========================================================================
{
public:
  HistoryPoint (const integer ID, const Element* e, const real R, 
		const real S, const real Z):
    id (ID), E (e), r (R), s (S), z (Z) { }

  integer               ID () const { return id; } 
  void                  extract (vector<AuxField*>&, real*) const;
  static const Element* locate (const real, const real,
				vector<Element*>&, real&, real&);

private:
  const integer  id ;		// Numeric identifier.
  const Element* E  ;		// Pointer to element.
  const real     r  ;		// Canonical-space r-location.
  const real     s  ;		// Canonical-space s-location.
  const real     z  ;		// Location in homogeneous direction.
};


class Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow
// solver.  This is designed to be overridden at implementation level
// if needed.
// ===========================================================================
{
public:
  Analyser  (Domain*, FEML*);
  ~Analyser () { }

  void analyse ();

protected:
  Domain*               src      ; // Source information.
  ofstream              par_strm ; // File for particle tracking.
  ofstream              his_strm ; // File for history points.
  ofstream              mdl_strm ; // File for modal energies.
  vector<HistoryPoint*> history  ; // Locations, etc. of history points.
};


class Statistics
// ===========================================================================
// Routines for statistical analysis of AuxFields.
// ===========================================================================
{
friend ifstream& operator >> (ifstream&, Statistics&);
friend ofstream& operator << (ofstream&, Statistics&);
public:
  Statistics (Domain*, vector<AuxField*>&);

  void update (AuxField***);
  void dump   ();

protected:
  const char*       name ;
  Domain*           base ;
  vector<AuxField*> src  ;
  vector<AuxField*> avg  ;
  integer           navg ;
};


template<class T> inline void roll (T* u, const integer n)
// ===========================================================================
// Stack roll template.  u is an array of type T, with at least n
// elements.  Roll up by one element.
// ===========================================================================
{
  if (n < 2) return;

  T tmp(u[n - 1]);

  for (register integer q(n - 1); q; q--)
    u[q] = u[q - 1];
  u[0] = tmp;
}

// -- Routines from misc.C:

ostream& printVector (ostream&, const char*, const integer, ... );
char*    upperCase   (char *);
void     writeField  (ofstream&, const char*, const int, const real,
		      vector<AuxField*>&);
#endif

