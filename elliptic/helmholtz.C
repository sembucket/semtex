//////////////////////////////////////////////////////////////////////////////
// helmholtz.C:  routines to solve elliptic problems in one variable.
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>


void Helmholtz (Domain*     D      ,
		const char* forcing)
// ---------------------------------------------------------------------------
// Solve Helmholtz's equation
//                                  2
//               div grad u - LAMBDA  u = f(x,y)
//
// subject to BCs.
// ---------------------------------------------------------------------------
{
  const real lambda2   = Femlib::value ("LAMBDA2"  );
  const real beta      = Femlib::value ("BETA"     );
  const int  iterative = Femlib::value ("ITERATIVE");
  const int  nmodes    = Geometry::nMode();
  const int  nZ        = Geometry::nZ();
  AuxField*  Force     = new AuxField (D -> Esys, nZ);

  if   (forcing) *Force = forcing;
  else           *Force = 0.0;
  
  if (iterative) {
   *D -> u[0] = 0.0;  
    D -> u[0] -> solve (Force, lambda2);

  } else {
    vector<Element*>&  E = D -> Esys;
    NumberSystem*      N = D -> Nsys[0];
    ModalMatrixSystem* M = new ModalMatrixSystem (lambda2, beta, nmodes, E, N);

    D -> u[0] -> solve (Force, M);
  }

  D -> step++;
}
