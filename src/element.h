#ifndef ELEMENT_H
#define ELEMENT_H

class Element
// ===========================================================================
// Virtual 2D quadrilateral element class, equal order in each direction.
//                                                
//         s                                y ^     +4
//         ^                                  |    / \
//     4   |                                  |  1+   +3
//     +------+                               |    \  |     
//     |   |  |      <== Logical  space       |     \ |
//     |   +--|-> r      Physical space ==>   |      \|
//     |      |                               |       +2  
//     +------+                               +-----------> x
//     1      2
//
// Master element coordinates: -1 < r,s < 1; edge traverses CCW.
//
// Elements in this classe will use Gauss--Lobatto--Legendre
// meshes, integration, and nodal basis functions.
//
// All 2D storage is row-major.
// ===========================================================================
{
public:
  Element (const integer, const integer, const Mesh*);
  ~Element();
  
  integer ID () const { return _id; }

  // -- Elemental Helmholtz matrix constructor, operator.

  void HelmholtzSC   (const real, const real, real*, real*,
		      real*, real*, real*, integer*)               const;
  void HelmholtzDiag (const real, const real, real*, real*)        const;
  void HelmholtzKern (const real, const real,
		      real*, real*, real*, real*)                  const;
  void HelmholtzOp   (const real, const real, real*, real*, real*) const;

  // -- Local/global projectors.

  void local2global      (const real*, const integer*, real*, real*)     const;
  void local2globalSum   (const real*, const integer*, real*, real*)     const;
  void local2globalSumSC (real*, const integer*,real*,const real*,real*) const;
  void global2local      (real*, const integer*,const real*,const real*) const;
  void global2localSC    (const real*, const integer*, real*,
			  real*, const real*, const real*, real*)        const;

  // -- Project from one interpolation order to another.

  void project (const integer, const real*, const integer, real*, real*) const;
  
  // -- Element-boundary operators.

  void bndryDsSum  (const integer*, const real*, real*)                 const;
  void bndryMask   (const integer*, real*, const real*, const integer*) const;
  void bndryInsert (const integer*, const real*, real*)                 const;

  // -- Element-side operators.

  void sideGeom  (const integer, real*, real*, real*, real*, real*) const;
  void sideEval  (const integer, real*, const  char*)               const;
  void sideGrad  (const integer, const real*, real*, real*, real*)  const;
  void sideGet   (const integer, const real*, real*)                const;
  void sideGetY  (const integer, real*)                             const;
  void sideMulY  (const integer, const real*, real*)                const;
  void sideDivY  (const integer, const real*, real*)                const;
  void sideDivY2 (const integer, const real*, real*)                const;

  // -- Element-operator functions.

  void grad  (real*, real*, real*)             const;
  void gradX (const real*, const real*, real*) const;
  void gradY (const real*, const real*, real*) const;

  void divY (real*) const;
  void mulY (real*) const;
  void mulX (real*) const;
  
  void evaluate (const char*, real*) const;

  real integral (const char*)        const;
  real integral (const real*, real*) const;
  real area     ()                   const;
  void weight   (real*)              const;

  void lengthScale (real*)                                       const;
  real CFL         (const real, const real*, const real*, real*) const;

  real norm_inf (const real*) const;
  real norm_L2  (const real*) const;
  real norm_H1  (const real*) const;

  // -- Probe functions.

  bool locate (const real, const real, real&, real&, const bool = false) const;
  real probe  (const real, const real, const real*, real*)               const;

  // -- Debugging/informational routines.

  void printMesh  ()            const;
  void printBndry (const real*) const;
  void printMatSC (const real*, const real*, const real*)       const;
  void Helmholtz  (const real, const real, real*, real*, real*) const;

protected:

  const integer _id  ;		// Element identifier.
  const integer _np  ;		// Number of points on an edge.
  const integer _npnp;		// Total number = np * np.
  const integer _next;		// Number of points on periphery.
  const integer _nint;		// Number of internal points.
  const bool    _cyl ;          // Cylindrical coordinate problem.

  const real* _zr   ;		// Master element mesh points on [-1, 1], r.
  const real* _zs   ;		// Master element mesh points on [-1, 1], s.
  const real* _wr   ;		// Master element quadrature weights, r.
  const real* _ws   ;		// Master element quadrature weights, s.
  const real* _DVr  ;		// Master element derivative operator, r.
  const real* _DTr  ;		// Transpose.
  const real* _DVs  ;		// Master element derivative operator, s.
  const real* _DTs  ;		// Transpose.

  integer*    _emap ;		// Indices of edges in nodal matrices.
  integer*    _pmap ;		// Inversion of emap (pmap[emap[i]] = i).

  real*       _xmesh;		// Physical space mesh.
  real*       _ymesh;		// 2D row-major store.

  real*       _drdx ;		// Partial derivatives (r, s) --> (x, y),
  real*       _dsdx ;		//   evaluated at quadrature points.
  real*       _drdy ;		//   (2D row-major storage.)
  real*       _dsdy ;		//

  real*       _delta;		// Local length scale.

  real*       _Q1   ;		// Geometric factor 1 at quadrature points.
  real*       _Q2   ;		// Geometric factor 2 at quadrature points.
  real*       _Q3   ;		// Geometric factor 3 at quadrature points.
  real*       _Q4   ;		// Geometric factor 4 at quadrature points.
  real*       _Q8   ;           // Like _Q4 but without possible factor of y.

  // -- Geometric and quadrature-specific internal functions.

  void mapping      ();
  void HelmholtzRow (const real, const real, const integer,
		     const integer, real*, real*) const;

  // -- BLAS-conforming edge offsets & skips for element-edge traverses.

  void terminal (const integer side, integer& start, integer& skip) const
  { switch (side) {
  case 0: start = 0;             skip  = 1;    break;
  case 1: start = _np - 1;       skip  = _np;  break;
  case 2: start = _np*(_np - 1); skip  = -1;   break;
  case 3: start = 0;             skip  = -_np; break;
  } }
};

#endif
