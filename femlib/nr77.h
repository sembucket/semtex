#ifndef NR77_H
#define NR77_H
///////////////////////////////////////////////////////////////////////////////
// nr77.h: C++ header file for Numerical Recipes (2e) routines (FORTRAN-77).
//
// NB: Numerical Recipes routines are assumed compiled in double precision.
//
// Reference: Numerical Recipes in FORTRAN, 2nd Edition.  CUP 1992.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cfemdef.h>

#define F77NAME(x) x##_

extern "C" {
  void   F77NAME(spline) (const double*, const double*, const integer&,
			  const double&, const double&, double*);
  void   F77NAME(splint) (const double*, const double*, const double*, 
			  const integer&, const double&, double&);
  double F77NAME(rtsec)  (double(*)(const double&), const double&,
			  const double&, const double&);
  double F77NAME(golden) (const double&, const double&, const double&,
			  double(*)(const double&), const double&, double&);
  void   F77NAME(savgol) (double*, const integer&, const integer&,
			  const integer&, const integer&, const integer&);
  void   F77NAME(polint) (const double*, const double*, const integer&,
			  const double&, double&, double&);
  void   F77NAME(polcoe) (const double*, const double*, const integer&,
			  double*);
  void   F77NAME(rk4)    (const double*, const double*, const integer&,
			  const double&, const double&, double*,
			  void(*)(const double&, const double*, double*));
  void   F77NAME(svdcmp) (double*, const integer&, const integer&,
			  const integer&,
			  const integer&, double*, double*);
  void   F77NAME(svbksb) (const double*, const double*, const double*,
			  const integer&, const integer&,
			  const integer&, const integer&,
			  const double*, double*);
  double F77NAME(ratval) (const double&, const double*, const integer&,
			  const integer&);
  double F77NAME(brent)  (const double&, const double&, const double&,
			  double(*)(const double&), const double&, double&);
  void   F77NAME(mnbrak) (double&, double&, double&, double&, double&, double&,
			  double(*)(const double&));
  void   F77NAME(ludcmp) (double*, const integer&, const integer&,
			  integer*, double&);
  void   F77NAME(lubksb) (const double*, const integer&, const integer&,
			  const integer*, double*);
}


class Recipes {
public:
  static void spline   (const double* x, const double* y, const integer& n,
			const double& yp1, const double& yp2, double* y2) {
    F77NAME(spline) (x, y, n, yp1, yp2, y2);
  }
  static void splint   (const double* xa, const double* ya, const double* ya2,
			const integer& n, const double& x, double& y) {
    F77NAME(splint) (xa, ya, ya2, n, x, y);
  }
  static double rtsec  (double(*func)(const double&), const double& x1,
			const double& x2, const double& xacc) {
    return F77NAME(rtsec) (func, x1, x2, xacc);
  }
  static double golden (const double& ax, const double& bx, const double& cx,
			double(*f)(const double&), const double& tol,
			double& xmin) {
    return F77NAME(golden) (ax, bx, cx, f, tol, xmin);
  }
  static void savgol   (double* c, const integer& np, const integer& nl,
			const integer& nr, const integer& ld,
			const integer& m) {
    F77NAME(savgol) (c, np, nl, nr, ld, m);
  }
  static void polcoe   (const double* x, const double* y,
			const integer& n, double* cof) {
    F77NAME(polcoe) (x, y, n, cof);
  }
  static void polint   (const double* xa, const double* ya, const integer& n, 
			const double& x, double& y, double& dy) {
    F77NAME(polint) (xa, ya, n, x, y, dy);
  }
  static void rk4      (const double* y, const double* dydx, const integer& n,
			const double& x, const double& h, double* yout,
			void(*derivs)(const double&, const double*, double*)) {
    F77NAME(rk4) (y, dydx, n, x, h, yout, derivs);
  } 
  static void svdcmp   (double* a, const integer& m, const integer& n,
			const integer mp, const integer& np,
			double* w, double* v) {
    F77NAME(svdcmp) (a, m, n, mp, np, w, v); 
  }
  static void svbksb   (const double* u, const double* w, const double* v,
			const integer& m, const integer& n, const integer& mp,
                        const integer& np, const double* b, double* x) {
    F77NAME(svbksb) (u, w, v, m, n, mp, np, b, x);
  }
  static double ratval (const double& x, const double* cof, const integer& mm,
			  const integer& kk) {
    return F77NAME(ratval) (x, cof, mm, kk);
  }
  
  static double brent  (const double& ax, const double& bx, const double& cx,
			double(*f)(const double&),
			const double& tol, double& xmin) {
    return F77NAME(brent) (ax, bx, cx, f, tol, xmin);
  }
  static void mnbrak   (double& ax, double& bx, double& cx,
			double& fa, double& fb, double& fc,
			double(*func)(const double&)) {
    F77NAME(mnbrak) (ax, bx, cx, fa, fb, fc, func);
  }
  static void ludcmp   (double* a, const integer& n, const integer& np,
			integer* indx, double& d) {
    F77NAME(ludcmp) (a, n, np, indx, d);
  }
  static void lubksb (const double* a, const integer& n, const integer& np,
		      const integer* indx, double* b) {
    F77NAME(lubksb) (a, n, np, indx, b);
  }
};


#endif
