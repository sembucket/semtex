#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

void build_preconditioner(int nSlice, int nDofsSlice, int nDofsPlane, int localShift, vector<Element*> elmt, Mat P);
