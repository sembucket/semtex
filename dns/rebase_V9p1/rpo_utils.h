#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

#include <fftw3.h>

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
    vector<AuxField*>   ui;
    vector<AuxField*>   fi;
    vector<AuxField*>   u0;
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
    VecScatter       global_to_semtex; // scatter from global data to semtex fields
    bool             build_PC;
    Vec              x_delta;
    int              nModesX;
    double           xmax;
    bool             x_fourier;
    int              travelling_wave;
    int              nElsX;
    int              nElsY;
    int              iteration;
    Domain*          write_i; // additional fields for file writing (ui)
    bool             build_dx;
    char*            session;
    fftw_complex*    data_s;
    fftw_complex*    data_f;
    fftw_plan        trans_fwd;
    fftw_plan        trans_bck;
    double*          rad_weights;
    double*          rad_coords;
    double           u_scale[3];
};

void data_transpose(real_t* data, int nx, int ny);
void elements_to_logical(int nex, int ney, real_t* data_els, real_t* data_log);
void logical_to_elements(int nex, int ney, real_t* data_log, real_t* data_els, int plane_i, int field_i);
void SEM_to_Fourier(int plane_k, Context* context, AuxField* us, real_t* data_r, real_t* data_i);
void Fourier_to_SEM(int plane_k, Context* context, AuxField* us, real_t* data_r, real_t* data_i, int field_i);
void UnpackX(Context* context, vector<AuxField*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x);
void RepackX(Context* context, vector<AuxField*> fields, real_t* theta, real_t* phi, real_t* tau, Vec x);
void assign_scatter_semtex(Context* context);
void phase_shift_x(Context* context, double theta, double sign, vector<Field*> fields);
void phase_shift_z(Context* context, double phi,   double sign, vector<Field*> fields);
void velocity_scales(Context* context);
int LocalIndex(Context* context, int field_i, int plane_i, int point_x, int point_y);
