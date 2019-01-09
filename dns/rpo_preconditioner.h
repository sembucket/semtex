#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

void build_preconditioner(int nSlice, int nDofsSlice, int nDofsPlane, int localSize, int localShift, vector<Element*> elmt, Mat P);
void build_preconditioner_ffs(int nSlice, int nDofsSlice, int nDofsPlane, int localSize, int** lShift, int* els, vector<Element*> elmt, Mat P);
