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
//
// NB: ua is added for nnewt.
// ===========================================================================
{
friend ifstream& operator >> (ifstream&, Domain&);
friend ofstream& operator << (ofstream&, Domain&);
public:
  Domain (FEML*, vector<Element*>&, BCmgr*);

  char*                name;	// Session name.
  int_t                step;	// Runtime step number.
  real_t               time;	// Simulation time.
  vector<Element*>&    elmt;	// Shared for equal-order interpolations.
  vector<real_t*>      udat;	// Data storage area for solution fields.
  vector<Field*>       u   ;	// Solution fields: velocities, pressure.
  vector<AuxField*>    ua  ;    // The above, plus any special AuxFields.
  vector<BoundarySys*> b   ;	// Field boundary systems.

  int_t nField    () const { return u.size(); }
  void  report    ();
  void  restart   ();
  void  dump      ();
  void  transform (const int_t);
  void  setNumber (const char, const NumberSys**) const;

private:
  char* field;		// Lower-case single character field names.
};

#endif
