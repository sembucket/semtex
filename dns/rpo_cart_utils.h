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
    vector<AuxField*>   ui;
    vector<AuxField*>   fi;
    vector<AuxField*>   u0;
    real_t           phi_i;
    real_t           tau_i;
    real_t           f_phi;
    real_t           f_tau;
    // parallel vector scattering data
    int              localSize;
    int*             lShift;
    IS               isl;              // local index set
    IS               isg;              // global index set
    VecScatter       global_to_semtex; // scatter from global data to semtex fields
    Vec              x_delta;          // for determining \delta x
    int              iteration;
    bool             travelling_wave;
    bool             build_dx;
    real_t**         coord_weights;
    int              n_mesh[3];
    int              n_mesh_max;
    int              n_mesh_sum;
    int*             n_mesh_sum_proc;
    int**            addToVector;      // semtex field data array to unique vector index
    double           u_scale[3];
    double           c_scale;
    AuxField*        uBar;
    double           beta;
};

void elements_to_logical(real_t* data_els, real_t* data_log);
void logical_to_elements(real_t* data_log, real_t* data_els);
void UnpackX(Context* context, vector<Field*> fields, real_t* phi, real_t* tau, Vec x);
void RepackX(Context* context, vector<Field*> fields, real_t  phi, real_t  tau, Vec x);
void assign_scatter_semtex(Context* context);
void phase_shift_z(Context* context, double phi, double sign, vector<Field*> fields);

void build_addToVector(Context* contex, vector<Field*> fields);
void build_coordWeights(Context* context);
void elements_to_vector(Context* context, int field_i, real_t* data_els, real_t* data_vec, bool fwd);
int LocalIndex(Context* context, int field_i, int plane_i, int point_i);
double GetScale(Context* context, int field_i, int mode_i, int mesh_i);
void _UnpackX(Context* context, vector<AuxField*> fields, real_t* phi, real_t* tau, Vec x, bool use_scale);
void _RepackX(Context* context, vector<AuxField*> fields, real_t phi, real_t tau, Vec x, bool use_scale);
void _phase_shift_z(Context* context, double phi, double sign, vector<Field*> fields);
void base_profile(Context* context, AuxField* ux, real_t scale, AuxField* uBar);
void velocity_scales(Context* context);
void _assign_scatter_semtex(Context* context);
