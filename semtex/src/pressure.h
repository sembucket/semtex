#ifndef PRESSURE_H
#define PRESSURE_H


class PBCmgr
// ===========================================================================
// This class maintains internal storage for evaluation of "high-order"
// pressure boundary conditions & provides means for their evaluation.
// ===========================================================================
{
public:
  static void build      (const Field*);
  static void maintain   (const int_t, const Field*, const AuxField**,
			  const AuxField**, const bool = true);
  static void evaluate   (const int_t, const int_t, const int_t,
			  const int_t, const real_t*, const real_t*, real_t*);
  static void accelerate (const Vector&, const Field*);

private:
  static real_t**** _Pnx;	// x component of dP / dn at domain  boundary.
  static real_t**** _Pny;	// y component of dP / dn at domain  boundary.
  static real_t**** _Unx;	// x component of normal velocity at boundary.
  static real_t**** _Uny;	// y component of normal velocity at boundary.
};

#endif
