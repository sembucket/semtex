#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <fftw3.h>

void build_preconditioner_ffs(Context* context, Mat P);
void build_preconditioner_I(Context* context, Mat P);
void rpo_set_fieldsplits(Context* context);
