/*****************************************************************************
 * elliptic.C:  routines to solve elliptic problems in one variable.
 *****************************************************************************/

static char RCSid[] = "$Id$";

#include <Fem.h>


void  Helmholtz (Domain* D, Mesh* M, char* forcing)
// ---------------------------------------------------------------------------
// Solve Helmholtz's equation
//                                  2
//               div grad u - LAMBDA  u = f(x,y)
//
// subject to BCs.  Parameter "LAMBDA2" needs to be set before entry.
// ---------------------------------------------------------------------------
{
  Field*  Force = new Field (*D -> u[0], *M);

  if   (forcing) *Force = forcing;
  else           *Force = 0.0;

  D -> u[0] -> evaluateBoundaries (0);  
  D -> u[0] -> buildSys (dparam ("LAMBDA2"));
  D -> u[0] -> solveSys (Force);

  D -> step ()++;
  D -> dump ();
}
