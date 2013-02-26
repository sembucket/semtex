#ifndef FIELDFORCE_H
#define FIELDFORCE_H

const char secForce[] = "FORCE";
const char forcename = 'u';		// forcing fields are uvw

class VirtualForce
// ---------------------------------------------------------------------------
// virtual base class for various types of forcing
// ---------------------------------------------------------------------------
{

public:
  void         allocStorage (Domain *);
  AuxField*    allocAuxField (Domain *, char);
  void         readSteadyFromFile(char *, vector<AuxField*>);

  virtual void physical (AuxField*, const int, vector<AuxField*>) {};
  virtual void fourier  (AuxField*, const int, vector<AuxField*>) {};

  vector<AuxField*>     _a;		// -- storage for pre-processed part

protected:
  Domain*		_D;
  bool			_enabled;
};

class FieldForce
{
public:
  FieldForce            (Domain *, FEML *);
  void addPhysical      (AuxField*, AuxField*, const int, vector<AuxField*>);
//   void updateU		(real_t*, const int);
  //void updateC		(const real_t*, const int);
  void addFourier       (AuxField*, AuxField*, const int, vector<AuxField*>);
  void dump             ();
  void writeAux		(vector<AuxField *>);
protected:
  bool			_enabled;
  vector<VirtualForce*> _classes;   // -- vector of concrete forcing classes
  Domain*		_D;
  vector<AuxField*>	_u;         // -- storage for physical space velocity
};

class ConstForce : public VirtualForce
// ---------------------------------------------------------------------------
// A force constant in both space in time, applied in Fourier space
// ---------------------------------------------------------------------------
{
public:
  ConstForce            (Domain *, FEML *);
  void fourier  (AuxField*, const int, vector<AuxField*>);
protected:
  real_t                _v[3];	// Force components
};


class SteadyForce : public VirtualForce
{
public:
  SteadyForce           (Domain *, FEML *);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
};

class WhiteNoiseForce : virtual public VirtualForce
{
public:
  WhiteNoiseForce       (Domain *, FEML *);
  void fourier		(AuxField*, const int, vector<AuxField*>);
protected:
  real_t                _eps[3];
  int_t                 _mode;
  int_t                 _apply_step; // apply force every _apply_step'th step
};

class ModulatedForce : virtual public VirtualForce
{
public:
  ModulatedForce        (Domain *, FEML *);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  char                  _alpha[3][StrMax]; // -- temporally varying part
};

class SpatioTemporalForce : virtual public VirtualForce
{
public:
  SpatioTemporalForce   (Domain *, FEML *);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  char                  _alpha[3][StrMax]; // -- spatially and temporally varying part
};


class SpongeForce : virtual public VirtualForce
{
public:
  SpongeForce           (Domain *, FEML *);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  vector<AuxField*>     _Uref;
  AuxField*             _mask;
  char                  _mask_func[StrMax]; // -- mask function, f(x,y,z,t)
  int                   _update; // mask update frequency
};

class DragForce : virtual public VirtualForce
{
public:
  DragForce             (Domain *, FEML *);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  AuxField		*_mask;
  AuxField		*_umag;
};

class CoriolisForce : virtual public VirtualForce
{
public:
  CoriolisForce         (Domain *, FEML *);
  void physical         (AuxField*, const int, vector<AuxField*>);
  void OmegaTimesOmegaTimesX();
protected:
  char                  _omega[3][StrMax];    // -- angular velocity = f(t) ..
  char                  _DomegaDt[3][StrMax]; // -- and its time derivative
  vector<real_t>        _o;		      // -- evaluated at current time step
  vector<real_t>        _minus_o;	      // -- - omega
  vector<real_t>        _minus_2o;	      // -- - 2 * omega
  vector<real_t>        _DoDt;		      // -- evaluated at current time step
  int_t                 _unsteady;            // -- 1 if omega is unsteady
};


#endif
