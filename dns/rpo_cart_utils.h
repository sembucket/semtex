#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

struct Context {
    int              nDofsSlice;
    int              nDofsPlane;
    Mesh*            mesh;
    vector<Element*> elmt;
    Domain*          domain;
    BCmgr*           bman;
    FEML*            file;
    DNSAnalyser*     analyst;
    FieldForce*      ff;
    vector<Field*>   ui;
    vector<Field*>   fi;
    vector<Field*>   u0;
    real_t           phi_i;
    real_t           tau_i;
    real_t           f_phi;
    real_t           f_tau;
    // parallel vector scattering data
    int              localSize;
    int*             lShift;
    IS               isl;  // local index set
    IS               isg;  // global index set
    VecScatter       global_to_semtex; // scatter from global data to semtex fields
    SNES             snes;
    Vec              x_delta; // for determining \delta x
    double           xmin;
    double           xmax;
    int              iteration;
    bool             travelling_wave;
    bool             build_dx;
    Domain*          write;   // additional fields for dumping at run time
};

void elements_to_logical(real_t* data_els, real_t* data_log);
void logical_to_elements(real_t* data_log, real_t* data_els);
void UnpackX(Context* context, vector<Field*> fields, real_t* phi, real_t* tau, Vec x);
void RepackX(Context* context, vector<Field*> fields, real_t  phi, real_t  tau, Vec x);
void assign_scatter_semtex(Context* context);
void phase_shift_z(Context* context, double phi,   double sign, vector<Field*> fields);
