#ifndef NEWT_H
#define NEWT_H
//////////////////////////////////////////////////////////////////////////////
// newt.h: header file for Newton's method steady state NS solver.
//
// http://www.netlib.org/templates/index.html
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "sem.h"

typedef void   (*Advection) (Domain*, AuxField**, AuxField**);
void integrate (Domain*, Advection);
void nonlinear (Domain*, AuxField**, AuxField**);
void linear    (Domain*, AuxField**, AuxField**);

extern "C" {

  void F77NAME(bicgstab)	// -- Templates Bi-Conj-Grad-Stab solver.

    (const int_t&   N    ,
     const real_t*  B    ,
     real_t*        X    ,
     real_t*        WORK ,
     const int_t&   LDW  ,
     int_t&         ITER ,
     real_t&        RESID,
     void (*MATVEC) (const real_t&, const real_t*, const real_t&, real_t*),
     void (*PSOLVE) (real_t*, const real_t*),
     int_t&         INFO );

  void F77NAME (bcgsw)		// -- NSPCG Bi-Conj-Grad-Squared solver.

    (void (*SUBA)   (const real_t*, const int_t*, const real_t*, const int_t*,
	  	     const int_t&, const real_t*, real_t*),
     void (*SUBQL)  (const real_t*, const int_t*, const real_t*, const int_t*,
		     const int_t&, const real_t*, real_t*),
     void (*SUBQR)  (const real_t*, const int_t*, const real_t*, const int_t*,
		     const int_t&, const real_t*, real_t*),
     const real_t* COEFF ,
     const int_t*  JCOEFF,
     const real_t* WFAC  ,
     const int_t*  JWFAC ,
     const int_t&  N     ,
     real_t*       U     ,
     const real_t* UBAR  ,
     const real_t* RHS   ,
     real_t*       WKSP  ,
     int_t&        NW    ,
     int_t*        IPARM ,
     real_t*       RPARM ,
     int_t&        IER   );

  void F77NAME (dfault)		// -- NSPCG setup.
    (int_t*, real_t*);

}

#endif
