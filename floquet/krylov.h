#ifndef KRYLOV_H
#define KRYLOV_H

#include <Sem.h>
#define F77name(x) x ## _

extern "C" {
  void F77name(dgeev)		// -- Lapack eigensystem routine.
    (const char*    N    ,      // computes all eigenvalues of matrix ??
     const char*    V    ,      // and returns as ??
     const integer& dim1 ,
     double*        H    ,
     const integer& dim2 ,
     double*        wr   ,
     double*        wi   ,
     double* f1   ,
     const integer& f2   ,
     double*        Hvec ,
     const integer& dim3 ,
     double*          rwork,
     const integer& lwork,
     integer&       ier  );
}


class Krylov
//
//
//
{
public:
  Krylov(FEML*, int, int, int);         // constructor
  ~Krylov();        // destructor

  void report();
  int restart();
  int dump();
  int dim();
  void normalise(int);  // normalise selected column
  void setDomain(Domain*, int);
  void getDomain(Domain*, int);
  real norm(int);
  void roll();
  void scale(real);

private:
  char*          _session;    // name of current session.
  int            _kdim  ;
  int            _nvec  ;
  int            _length;    // length of Krylov column
  vector<real*>  _column;    // columns of Krylov matrix
  real*          _K     ;    // pointer to K matrix
  int            _Klen  ;    // total length of K matrix.

};

#endif
