#ifndef MATRIX_H
#define MATRIX_H



class MatrixSys
// ===========================================================================
// System of global and local Helmholtz matrices and partitions.
// Matrix factorisations are Cholesky, use LAPACK storage schemes.
// ===========================================================================
{
friend class Field;
//friend ostream& operator << (ostream&, MatrixSys&);
//friend istream& operator >> (istream&, MatrixSys&);
public:
  MatrixSys     (const real, const real, const integer,
		 const vector<Element*>&, const BoundarySys*,
		 const SolverKind);
 ~MatrixSys     ();
  bool match    (const real, const real, const NumberSys*,
		 const SolverKind) const;

private:
  real  _HelmholtzConstant;	// Same for all modes.
  real  _FourierConstant  ;	// Varies with mode number.
 
  const vector<Boundary*>& _BC;	// Internal copy of Boundary conditions.
  const NumberSys*         _NS;	// Internal copy of NumberSys.

  integer    _nel     ;		// Number of elemental matrices.
  integer    _nglobal ;		// Number of unique element-boundary nodes.
  integer    _singular;		// If system is potentially singular.
  integer    _nsolve  ;		// System-specific number of global unknowns.
  SolverKind _method  ;		// Flag specifies direct or iterative solver.

  // -- For _method == DIRECT:

  integer    _nband ;		// Bandwidth of global matrix (incl. diagonal).
  integer    _npack ;		// Number of reals for global matrix.
  real*      _H     ;		// (Factored) packed global Helmholtz matrix.
  real**     _hbi   ;		// Element external-internal coupling matrices.
  real**     _hii   ;		// (Factored) internal-internal matrices.
  integer*   _bipack;		// Size of hbi for each element.
  integer*   _iipack;		// Size of hii for each element.

  // -- For _method == JACPCG:

  integer    _npts;		// Total number of unique meshpoints.
  real*      _PC  ;		// Diagonal preconditioner matrix.
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
		  const vector<Element*>&, const BoundarySys*, 
		  const SolverKind);
 ~ModalMatrixSys ();

  const MatrixSys* operator [] (const integer i) const { return _Msys[i]; }

private:
  char*              _fields;	// Character field tags for this system.
  vector<MatrixSys*> _Msys  ;	// One MatrixSys for each Fourier mode.
};

#endif
