#ifndef FIELD_H
#define FIELD_H
  
  

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
  Field  (BoundarySys*, real*, const int, vector<Element*>&, const char);
 ~Field  () { }

  Field& operator = (const AuxField& z) {AuxField::operator=(z); return *this;}
  Field& operator = (const real&     z) {AuxField::operator=(z); return *this;}
  Field& operator = (const char*     z) {AuxField::operator=(z); return *this;}


  Field& solve  (AuxField*, const ModalMatrixSys*);
  Field& solve  (AuxField*, const MatrixSys*);

  Field& smooth (AuxField* = 0);
  void   smooth (const int, real*) const;

  void evaluateBoundaries    (const int);
  void evaluateM0Boundaries  (const int);
  void addToM0Boundaries     (const real, const char*);

  static real   flux            (const Field*);
  static Vector normalTraction  (const Field*);
  static Vector moment          (const Field*, const Field*, const Field*);
  static Vector tangentTraction (const Field*, const Field*, const Field* = 0);
  static void   normTractionV   (real*, real*, const Field*);
  static void   tangTractionV   (real*, real*, real*, const Field*,
				 const Field*, const Field* = 0);

  static void coupleBCs    (Field*, Field*, const int);
  static real modeConstant (const char, const int, const real);

  static void printBoundaries (const Field*);
  static void printConnect    (const Field*);

private:
  int          _nbound;		// Number of boundary edges.
  int          _nline ;		// Length of one boundary line.
  real*        _sheet ;		// Wrap-around storage for data boundary.
  real**       _line  ;		// Single plane's worth of sheet.
  BoundarySys* _bsys  ;		// Boundary system information.

  void getEssential      (const real*, real*,
			  const vector<Boundary*>&, const NumberSys*) const;
  void setEssential      (const real*, real*, const NumberSys*);
  void local2global      (const real*, real*, const NumberSys*)       const;
  void global2local      (const real*, real*, const NumberSys*)       const;

  void constrain         (real*, const real, const real,
			  const real*, const NumberSys*, real*)       const;
  void buildRHS          (real*, const real*, real*, real*, const real**,
			  const int, const int, const vector<Boundary*>&,
			  const NumberSys*, real*)                    const;
  void HelmholtzOperator (const real*, real*, const real,
			  const real, const int, real*)               const;
};

#endif
