#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

struct Context {
    int              nSlice;
    int              nField;
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
    real_t*          theta_i;
    real_t*          phi_i;
    real_t*          tau_i;
    real_t           f_theta;
    real_t           f_phi;
    real_t           f_tau;
    // regular grid points in elements
    int_t*           el;
    real_t*          r;
    real_t*          s;
    // parallel vector scattering data
    int              localSize;
    int**            lShift;
    IS               isl;  // local index set
    IS               isg;  // global index set
    VecScatter       global_to_semtex; // scatter from global data to semtex fields
    bool             build_PC;
    // for the fieldsplit preconditioning
    IS*              is_s;
    IS*              is_u;
    IS*              is_p;
    SNES             snes;
    double           shift_scale; // for scaling the phase shifts
    double           dx_norm_prev;
    Vec              x_prev;
    Vec              x_delta;
    Vec              dx_test;
    int              nModesX;
    double           xmax;
    bool             x_fourier;
    bool             travelling_wave;
    int              nElsX;
    int              nElsY;
    int              prev_newton_it;
    int              iteration;
    Domain*          write_i; // additional fields for file writing (ui)
    bool             build_dx;
};

void data_transpose(real_t* data, int nx, int ny);
void elements_to_logical(int nex, int ney, real_t* data_els, real_t* data_log);
void logical_to_elements(int nex, int ney, real_t* data_log, real_t* data_els);
void SEM_to_Fourier(int plane_k, Context* context, Field* us, real_t* data_f);
void Fourier_to_SEM(int plane_k, Context* context, Field* us, real_t* data_f);
void UnpackX(Context* context, vector<Field*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x, double shift_scale);
void RepackX(Context* context, vector<Field*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x, double shift_scale);
void UnpackConstraints(Context* context, real_t* theta, real_t* phi, real_t* tau, Vec x, double shift_scale);
void RepackConstraints(Context* context, real_t* theta, real_t* phi, real_t* tau, Vec x, double shift_scale);
void assign_scatter_semtex(Context* context);
void phase_shift_x(Context* context, double theta, double sign, vector<Field*> fields);
void phase_shift_z(Context* context, double phi,   double sign, vector<Field*> fields);
