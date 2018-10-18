typedef ModalMatrixSys Msys;

static Msys** preSolve (const Domain* D);
static void Solve (Domain*     D, const int_t i, AuxField*   F, Msys*       M);
static void waveProp (Domain*           D , const AuxField*** Us, const AuxField*** Uf);
static void setPForce (const AuxField** Us, AuxField**       Uf);
static void project (const Domain* D , AuxField**    Us, AuxField**    Uf);
void integrate (void (*advection) (Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
                Domain*      D , BCmgr*       B , FieldForce*  FF, int nStep);
