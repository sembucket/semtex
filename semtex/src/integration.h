#ifndef INTEGRATION_H
#define INTEGRATION_H


class Integration
// ===========================================================================
// Return coefficients for time-integration schemes.
// ===========================================================================
{
public:
  static const integer OrderMax;

  static void AdamsBashforth (const integer, real*);
  static void AdamsMoulton   (const integer, real*);
  static void StifflyStable  (const integer, real*);
  static void Extrapolation  (const integer, real*);
};

#endif
