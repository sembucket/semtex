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
  static void maintain   (const integer, const Field*, const AuxField**,
			  const AuxField**, const bool = true);
  static void evaluate   (const integer, const integer, const integer,
			  const integer, const real*, const real*, real*);
  static void accelerate (const Vector&, const Field*);

private:
  static real**** _Pnx;		// x component of dP / dn at domain  boundary.
  static real**** _Pny;		// y component of dP / dn at domain  boundary.
  static real**** _Unx;		// x component of normal velocity at boundary.
  static real**** _Uny;		// y component of normal velocity at boundary.
};

#endif
