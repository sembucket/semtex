#ifndef FIELD_H
#define FIELD_H

class Field : public AuxField
//  ==========================================================================
/// Field adds boundary conditions and global numbering to AuxField.
/// With the boundary conditions comes knowledge of which sections of
/// the boundary hold essential (known) nodes, and which contain
/// natural (unknown) nodes.
///
/// A Field holds global node numbers and solve masks for element
/// boundaries: where domain boundary value is given by an essential
/// BC, the solve mask value is 0, for all other nodes the value is 1.
///
/// Helmholtz problem solution
/// --------------------------
/// Field solve routines provide solution to the discrete form of the
/// Helmholtz equation
/// \f[
///                 \nabla^2 u - \lambda^2 u = f
/// \f]
///
/// on domain \f$\Omega\f$, subject to essential BCs \f$u=g\f$ on
/// \f$\Gamma_g\f$ and natural BCs \f$\partial u/\partial n=h\f$ on
/// \f$\Gamma_h\f$, where the boundary \f$\Gamma\f$ of \f$\Omega\f$ is
/// the union of (non-overlapping) \f$\Gamma_g\f$ and \f$\Gamma_h\f$
/// and \f$\boldmath{n}\f$ is the unit outward normal vector on
/// \f$\Gamma\f$.  \f$\lambda^2\f$ is called the Helmholtz constant
/// below.
///
/// The Galerkin form, using integration by parts with weighting functions w
/// which are zero on \f$\Gamma_g\f$, is
/// \f[
///            (\nabla u, \nabla w) + \lambda^2 (u, w) = - (f, w) + <h, w>
/// \f]
/// where
/// + \f$(a, b)=\int a.b~d\Omega\f$ is an integration over the domain and
/// + \f$<a, b>=\int a.b~d\Gamma\f$ is an integration over the domain boundary.
///
/// The discrete (finite element) equivalent is to solve
/// \f[
///                   K.u + \lambda^2 M.u = - M.f + <h, w>
/// \f]
/// or
/// \f[
///                         H.u = - M.f + <h, w>
/// \f]
/// where \f$K, M\f$ and \f$H\f$ are respectively (assembled)
/// "stiffness", "mass" and Helmholtz matrices.
///
//  Some complications arise from dealing with essential boundary
//  conditions, since typically the elemental matrices K^e, M^e which
//  are assembled to form K and M do not account for the boundary
//  requirements on w.  There are a number of ways of dealing with this
//  issue: one approach is to partition H as it is formed (here F =
//  -M.f + <h, w>):
// 
//    +--------+-------------+ /  \     /  \
//    |        |             | |  |     |  |
//    |   Hp   |     Hc      | |u |     |F |   u : nodal values for solution.
//    |        |(constraint) | | s|     | s|    s
//    |        |             | |  |     |  |       (n_solve values)
//    +--------+-------------+ +--+     +--+
//    |        | H_ess: this | |  |  =  |  |
//    |        | partition   | |  |     |  |
//    |    T   | relates to  | |u |     |F |   u : are given essential BCs.
//    |  Hc    | essential   | | g|     | g|    g
//    |        | BC nodes    | |  |     |  |       (n_global - n_solve values)
//    |        | and is not  | |  |     |  |
//    |        | assembled.  | |  |     |  |
//    +--------+-------------+ \  /     \  /
//
//  Partition out the sections of the matrix corresponding to the known
//  nodal values (essential BCs), and solve instead the constrained
//  problem
//
//    +--------+               /  \     /  \     +-------------+ /  \
//    |        |               |  |     |  |     |             | |  |
//    |   Hp   |               |u |     |F |     |     Hc      | |  |
//    |        |               | s|  =  | s|  -  |             | |u |
//    |        |               |  |     |  |     |             | | g|
//    +--------+               \  /     \  /     +-------------+ |  |.
//                                                               |  |
//  Here n_global is the number of nodes that receive global node
//  numbers, typically those on the mesh edges.  N_solve is the number
//  of these nodes that have values that must be solved for,
//  i.e. n_global minus the number of global nodes situated on
//  essential-type boundaries.
// 
//  The action of Hc on u_g can be formed by a loop over the elements
//  which touch the essential boundaries and does not require storage
//  of the partition Hc.  M.f can also be formed on an
//  element-by-element basis and is cheap for nodal spectral elements
//  since M is diagonal.  Both tasks are performed by the
//  Field::constrain routine below.
//
/// Field names
/// -----------
/// The (one character) names of field variables are significant, and have
/// the following reserved meanings:
/// 
/// + u:  First velocity component.           (Cylindrical: axial     velocity.)
/// + v:  Second velocity component.          (Cylindrical: radial    velocity.)
/// + w:  Third velocity component.           (Cylindrical: azimuthal velocity.)
/// + p:  Pressure divided by density.
/// + c:  Scalar for transport or elliptic problems.
///
//  ==========================================================================
{
friend class BCmgr;
public:
  Field  (BoundarySys*, real_t*, const int_t, vector<Element*>&, const char);
 ~Field  () { }
 
  Field& operator = (const AuxField& z) {AuxField::operator=(z); return *this;}
  Field& operator = (const real_t&   z) {AuxField::operator=(z); return *this;}
  Field& operator = (const char*     z) {AuxField::operator=(z); return *this;}

  Field& solve  (AuxField*, const ModalMatrixSys*);

  Field& smooth (AuxField* = NULL);
  void   smooth (const int_t, real_t*) const;

  void evaluateBoundaries    (const Field*, const int_t, const bool = true);
  void evaluateM0Boundaries  (const Field*, const int_t);
  void addToM0Boundaries     (const real_t, const char*);
  void bTransform            (const int_t);

  void overwriteForGroup      (const char*, const AuxField*, AuxField*);

  static real_t scalarFlux    (const Field*);
  static Vector normTraction  (const Field*);
  static Vector tangTraction  (const Field*, const Field*, const Field* = 0);
  static void   normTractionV (real_t*, real_t*, const Field*);
  static void   tangTractionV (real_t*, real_t*, real_t*, const Field*,
			       const Field*, const Field* = 0);
  static void   traction      (real_t*, real_t*, real_t*,
			       const int_t, const int_t, const Field*,
			       const Field*, const Field*, const Field* = 0);
  static void   coupleBCs     (Field*, Field*, const int_t);
  static real_t modeConstant  (const char, const int_t, const real_t);

private:
  int_t        _nbound;		//!<  Number of boundary edges.
  int_t        _nline ;		//!<  Length of one boundary line.
  real_t*      _sheet ;		//!<  Wrap-around storage for data boundary.
  real_t**     _line  ;		//!<  Single plane's worth of sheet.
  BoundarySys* _bsys  ;		//!<  Boundary system information.

  void getEssential      (const real_t*, real_t*,
			  const vector<Boundary*>&, const NumberSys*) const;
  void setEssential      (const real_t*, real_t*, const NumberSys*);

  void local2global      (const real_t*, real_t*, const NumberSys*)   const;
  void global2local      (const real_t*, real_t*, const NumberSys*)   const;
  void local2globalSum   (const real_t*, real_t*, const NumberSys*)   const;

  void constrain         (real_t*, const real_t, const real_t,
			  const real_t*, const NumberSys*, real_t*)   const;
  void buildRHS          (real_t*,const real_t*, real_t*, real_t*, 
			  const real_t**, const int_t, const int_t,
			  const vector<Boundary*>&,
			  const NumberSys*, real_t*)                  const;
  void HelmholtzOperator (const real_t*, real_t*, const real_t,
			  const real_t, const int_t, real_t*)         const;
};

#endif
