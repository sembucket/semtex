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

#define F77NAME(x) x##_

extern "C" {
  void   F77NAME(spline) (double*, double*, const int&,
			  const double&, const double&, double*);
  void   F77NAME(splint) (double*, double*, double*, const int&,
			  const double&, double&);
  double F77NAME(rtsec)  (double(*)(const double&), const double&,
			  const double&, const double&);
  double F77NAME(golden) (const double&, const double&, const double&,
			  double(*)(const double&), const double&, double&);
  void   F77NAME(savgol) (double*, const int&, const int&,
			  const int&, const int&, const int&);
  void   F77NAME(polint) (double*, double*, const int&, const double&,
			  double&, double&);
}


class Recipes {
public:
  static void spline   (double* x, double* y, const int& n,
			const double& yp1, const double& yp2, double* y2) {
    F77NAME(spline) (x, y, n, yp1, yp2, y2);
  }
  static void splint   (double* xa, double* ya, double* ya2, const int& n,
			const double& x, double& y) {
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
  static void savgol   (double* c, const int& np, const int& nl,
			const int& nr, const int& ld, const int& m) {
    F77NAME(savgol) (c, np, nl, nr, ld, m);
  }
  static void polint   (double* xa, double* ya, const int& n, 
			const double& x, double& y, double& dy) {
    F77NAME(polint) (xa, ya, n, x, y, dy);
  }
};


#endif
