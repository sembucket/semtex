#ifndef DOMAIN_H
#define DOMAIN_H


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

  char*                name;	// Session name.
  integer              step;	// Runtime step number.
  real                 time;	// Simulation time.
  vector<Element*>&    elmt;	// Shared for equal-order interpolations.
  vector<real*>        udat;	// Data storage area for solution fields.
  vector<Field*>       u   ;	// Solution fields: velocities, pressure.
  vector<BoundarySys*> b   ;	// Field boundary systems.

  integer nField    () const { return u.size(); }
  void    report    (ostream& stream = cout);
  void    restart   ();
  void    dump      ();
  void    setNumber (const char, const NumberSys**) const;

  // -- Required for base fields and stability analysis.

  vector<AuxField*> U       ; // -- Base velocity fields - no BCs.
  vector<real*>     Udat    ; // -- Data storage area for base auxfields.
  vector<real*>     baseFlow; // -- Fourier transformed base velocities.
  real              period  ; // -- total time for one period of base flow.

  void loadBase  ();
  void updateBase();

private:
  char* field;	      // Lower-case single character perturbation  field names.
  char* baseField;    // Upper-case single character base velocity field names.
};

#endif
