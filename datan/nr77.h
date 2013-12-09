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
  void   F77NAME(spline) (const double*, const double*, const int_t&,
			  const double&, const double&, double*);
  void   F77NAME(splint) (const double*, const double*, const double*, 
			  const int_t&, const double&, double&);
  double F77NAME(rtsec)  (double(*)(const double&), const double&,
			  const double&, const double&);
  double F77NAME(golden) (const double&, const double&, const double&,
			  double(*)(const double&), const double&, double&);
  void   F77NAME(savgol) (double*, const int_t&, const int_t&,
			  const int_t&, const int_t&, const int_t&);
  void   F77NAME(polint) (const double*, const double*, const int_t&,
			  const double&, double&, double&);
  void   F77NAME(polcoe) (const double*, const double*, const int_t&,
			  double*);
  void   F77NAME(rk4)    (const double*, const double*, const int_t&,
			  const double&, const double&, double*,
			  void(*)(const double&, const double*, double*));
  void   F77NAME(svdcmp) (double*, const int_t&, const int_t&,
			  const int_t&,
			  const int_t&, double*, double*);
  void   F77NAME(svbksb) (const double*, const double*, const double*,
			  const int_t&, const int_t&,
			  const int_t&, const int_t&,
			  const double*, double*);
  double F77NAME(ratval) (const double&, const double*, const int_t&,
			  const int_t&);
  double F77NAME(brent)  (const double&, const double&, const double&,
			  double(*)(const double&), const double&, double&);
  void   F77NAME(mnbrak) (double&, double&, double&, double&, double&, double&,
			  double(*)(const double&));
  void   F77NAME(ludcmp) (double*, const int_t&, const int_t&,
			  int_t*, double&);
  void   F77NAME(lubksb) (const double*, const int_t&, const int_t&,
			  const int_t*, double*);
  void   F77NAME(correl) (const double*, const double*,
			  const int_t&, double*);
}


class Recipes {
public:
  static void spline   (const double* x, const double* y, const int_t& n,
			const double& yp1, const double& yp2, double* y2) {
    F77NAME(spline) (x, y, n, yp1, yp2, y2);
  }
  static void splint   (const double* xa, const double* ya, const double* ya2,
			const int_t& n, const double& x, double& y) {
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
  static void savgol   (double* c, const int_t& np, const int_t& nl,
			const int_t& nr, const int_t& ld,
			const int_t& m) {
    F77NAME(savgol) (c, np, nl, nr, ld, m);
  }
  static void polcoe   (const double* x, const double* y,
			const int_t& n, double* cof) {
    F77NAME(polcoe) (x, y, n, cof);
  }
  static void polint   (const double* xa, const double* ya, const int_t& n, 
			const double& x, double& y, double& dy) {
    F77NAME(polint) (xa, ya, n, x, y, dy);
  }
  static void rk4      (const double* y, const double* dydx, const int_t& n,
			const double& x, const double& h, double* yout,
			void(*derivs)(const double&, const double*, double*)) {
    F77NAME(rk4) (y, dydx, n, x, h, yout, derivs);
  } 
  static void svdcmp   (double* a, const int_t& m, const int_t& n,
			const int_t mp, const int_t& np,
			double* w, double* v) {
    F77NAME(svdcmp) (a, m, n, mp, np, w, v); 
  }
  static void svbksb   (const double* u, const double* w, const double* v,
			const int_t& m, const int_t& n, const int_t& mp,
                        const int_t& np, const double* b, double* x) {
    F77NAME(svbksb) (u, w, v, m, n, mp, np, b, x);
  }
  static double ratval (const double& x, const double* cof, const int_t& mm,
			  const int_t& kk) {
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
  static void ludcmp   (double* a, const int_t& n, const int_t& np,
			int_t* indx, double& d) {
    F77NAME(ludcmp) (a, n, np, indx, d);
  }
  static void lubksb (const double* a, const int_t& n, const int_t& np,
		      const int_t* indx, double* b) {
    F77NAME(lubksb) (a, n, np, indx, b);
  }
  static void correl (const double* data1, const double* data2,
		      const int_t& n, double* ans) {
    F77NAME(correl) (data1, data2, n, ans);
  }
};


#endif
