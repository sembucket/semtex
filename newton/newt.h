#ifndef NEWT_H
#define NEWT_H
//////////////////////////////////////////////////////////////////////////////
// newt.h: header file for Newton's method steady state NS solver.
//
// http://www.netlib.org/templates/index.html
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "Sem.h"

typedef void   (*Advection) (Domain*, AuxField**, AuxField**);
void integrate (Domain*, Advection);
void nonlinear (Domain*, AuxField**, AuxField**);
void linear    (Domain*, AuxField**, AuxField**);

extern "C" {

  void F77NAME(bicgstab)	// -- Templates Bi-Conj-Grad-Stab solver.

    (const integer& N    ,
     const real*    B    ,
     real*          X    ,
     real*          WORK ,
     const integer& LDW  ,
     integer&       ITER ,
     real&          RESID,
     void (*MATVEC) (const real&, const real*, const real&, real*),
     void (*PSOLVE) (real*, const real*),
     integer&       INFO );

  void F77NAME (bcgsw)		// -- NSPCG Bi-Conj-Grad-Squared solver.

    (void (*SUBA)   (const real*, const integer*, const real*, const integer*,
	  	     const integer&, const real*, real*),
     void (*SUBQL)  (const real*, const integer*, const real*, const integer*,
		     const integer&, const real*, real*),
     void (*SUBQR)  (const real*, const integer*, const real*, const integer*,
		     const integer&, const real*, real*),
     const real*    COEFF ,
     const integer* JCOEFF,
     const real*    WFAC  ,
     const integer* JWFAC ,
     const integer& N     ,
     real*          U     ,
     const real*    UBAR  ,
     const real*    RHS   ,
     real*          WKSP  ,
     integer&       NW    ,
     integer*       IPARM ,
     real*          RPARM ,
     integer&       IER   );

  void F77NAME (dfault)		// -- NSPCG setup.
    (integer*, real*);

}

#endif
