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
  Field  (BoundarySys*, real_t*, const int_t, vector<Element*>&, const char);
 ~Field  () { }

  Field& operator = (const AuxField& z) {AuxField::operator=(z); return *this;}
  Field& operator = (const real_t&   z) {AuxField::operator=(z); return *this;}
  Field& operator = (const char*     z) {AuxField::operator=(z); return *this;}

  Field& solve  (Domain*, AuxField*, const MatrixSys*, const vector<AuxField*>& , bool);

  Field& smooth (AuxField* = 0);

  void smooth (const int_t, real_t*) const;
  void evaluateBoundaries    (const int_t);
  void evaluateControl    (const int_t , real_t*, ofstream&);
  void evaluateM0Boundaries  (const int_t);
  void addToM0Boundaries     (const real_t, const char*);
   
  static real_t flux            (const Field*);
  static Vector normalTraction  (const Field*);
  static Vector moment          (const Field*, const Field*, const Field*);
  static Vector tangentTraction (const Field*, const Field*, const Field* = 0);
  static void   normTractionV   (real_t*, real_t*, const Field*);
  static void   tangTractionV   (real_t*, real_t*, real_t*, const Field*,
				 const Field*, const Field* = 0);

  static void   coupleBCs    (Field*, Field*, const int_t);
  static real_t modeConstant (const char, const int_t, const real_t);

  static void printBoundaries (const Field*);
  static void printConnect    (const Field*);
  real_t normc(real_t*);
  real_t normc_mixed(real_t*, real_t*);
  int_t  size_controlbc();
  void add_adjoint (real_t*, const int_t, const int_t);
  void add_adjoint_pressure(real_t*, real_t*, const int_t, const int_t);
  void controlmesh (real_t*,real_t*);
	  void control_normal_direction (real_t* ,real_t*);
private:
  int_t        _nbound;		// Number of boundary edges.
  int_t        _nline ;		// Length of one boundary line.
  real_t*      _sheet ;		// Wrap-around storage for data boundary.
  real_t**     _line  ;		// Single plane's worth of sheet.
  BoundarySys* _bsys  ;		// Boundary system information.

  void getEssential      (const real_t*, real_t*,
			  const vector<Boundary*>&, const NumberSys*)   const;
  void setEssential      (const real_t*, real_t*, const NumberSys*);
  void local2global      (const real_t*, real_t*, const NumberSys*)     const;
  void global2local      (const real_t*, real_t*, const NumberSys*)     const;

  void constrain         (real_t*, const real_t, const real_t,
			  const real_t*, const NumberSys*, real_t*)     const;
  void buildRHS          (real_t*, const real_t*, real_t*, real_t*,
			  const real_t**, const int_t, const int_t,
			  const vector<Boundary*>&, const NumberSys*,
			  real_t*)                                      const;
  void HelmholtzOperator (const real_t*, real_t*, const real_t,
			  const real_t, const int_t, real_t*, const vector<AuxField*>& , bool)           const;
};

#endif
