//////////////////////////////////////////////////////////////////////////////
// helmholtz.C:  routines to solve elliptic problems in one variable.
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <Fem.h>


void Helmholtz (Domain* D, char* forcing, const real& lambda2)
// ---------------------------------------------------------------------------
// Solve Helmholtz's equation
//                                  2
//               div grad u - LAMBDA  u = f(x,y)
//
// subject to BCs.
// ---------------------------------------------------------------------------
{
  Field* Force = new Field (*D -> u[0]);

  if   (forcing) *Force = forcing;
  else           *Force = 0.0;

  D -> u[0] -> evaluateBoundaries (0);
  
  if (option ("ITERATIVE")) {
    *D -> u[0] = 0.0;
     D -> u[0] -> solve (Force, lambda2);
    
  } else {
    D -> u[0] -> assemble (lambda2);
    D -> u[0] -> solve    (Force);
  }

  D -> step ()++;
}
