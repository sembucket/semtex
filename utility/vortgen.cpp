#include <sem.h>
#include <tensorcalcs.h>

#define FLAG_MAX 2
#define FLDS_MAX 2

//
// Globals
//
static char prog[] = "vortgen";
enum
{
  VORTICITY,
  VORTGEN
};

//
// Entrypoint
//
int main(int argc, char **argv)
{
  //
  // Parse args
  //
  char *session, *dump, *func;
  bool add[FLAG_MAX]{0}, need[FLAG_MAX]{0}; // Flag initialised as zero
  getargs(argc, argv, session, dump, func, add);

  ifstream file;
  file.open(dump);
  if (!file)
    message(prog, "can't open input field file", ERROR);

  //
  // Domain
  //
  FEML *F = new FEML(session);
  Mesh *M = new Mesh(F);
  int_t np = Femlib::ivalue("N_P");
  int_t nz = Femlib::ivalue("N_Z");
  int_t nel = M->nEl();
  Geometry::CoordSys system = (Femlib::ivalue("CYLINDRICAL")) ? Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set(np, nz, nel, system);
  int_t allocSize = Geometry::nTotal();
  int_t NDIM = Geometry::nDim();
  if (NDIM != 2)
    message(prog, "non 2D geometry");

  vector<Element *> element;
  element.resize(nel);
  for (int_t i = 0; i < nel; i++)
    element[i] = new Element(i, np, M);

  BCmgr *B = new BCmgr(F, element);
  Domain *D = new Domain(F, element, B);
  int_t NCOM;
  if (strstr(D->field, "uvw"))
    NCOM = 3;
  else if (strstr(D->field, "uv"))
    NCOM = 2;
  else
    message(prog, "lacking velocity components: is session valid?", ERROR);

  //
  // Initialise fields and gradients
  //
  vector<AuxField *> velocity;
  velocity.resize(NCOM);
  for (int_t i = 0; i < NCOM; i++)
    velocity[i] = D->u[i];

  AuxField *pressure;
  pressure = D->u[NCOM];

  vector<vector<AuxField *>> Vij;
  vector<vector<real_t *>> VijData;
  Vij.resize(3);
  VijData.resize(3);
  for (int_t i = 0; i < 3; i++)
  {
    Vij[i].resize(3);
    VijData[i].resize(3);
    for (int_t j = 0; j < 3; j++)
    {
      VijData[i][j] = new real_t[allocSize];
      Vij[i][j] = new AuxField(VijData[i][j], nz, element);
      *Vij[i][j] = 0.0;
    }
  }

  vector<AuxField *> addField(FLDS_MAX);
  int_t iAdd = 0;

  if (need[VORTICITY])
  {
    vector<AuxField *> vorticity;
    vector<real_t *> VorData;
    vorticity.resize(1);
    VorData.resize(1);
    VorData[0] = new real_t[allocSize];
    vorticity[0] = new AuxField(VorData[0], nz, element, 't');
    addField[iAdd++] = vorticity[0];
  }

  if (need[VORTGEN])
  {
    vector<AuxField *> vortgen;
    vector<real_t *> VortGenData;
    vortgen.resize(5);
    VortGenData.resize(5);
    for (int_t i = 0; i < 5; i++)
    {
      VortGenData[i] = new real_t[allocSize];
      vortgen[i] = new AuxField(VortGenData[i], nz, element, 'k' + i);
      addField[iAdd++] = vortgen[i];
    }
  }

  //
  // Cycle through field dump
  //
  while (getDump(D, file))
  {
    for (int_t i = 0; i < NDIM; i++)
      for (int_t j = 0; j < NCOM; j++)
      {
        *Vij[i][j] = *velocity[j];
        if (i == 2)
          Vij[i][j]->transform(FORWARD);
        Vij[i][j] -> gradient(i);
        if (i ==2)
          Vij[i][j] -> transform(INVERSE);
      }
  }

  file.close();
  return EXIT_SUCCESS;
}

//
// Parse command line arguments
//
static void getargs(int argc, char **argv, char *&session, char *&dump, char *&func, bool *flag)
{
  //
  // CLI
  //
  char usage[] =
      "vortgen.cpp: process semtex/NEKTON-type field files, computer and\n"
      "adding 2D vorticity and its derivatives\n"
      "\n"
      "Usage: %s [options] -s session dump.fld\n"
      "options:\n"
      "  -h        ... print this message \n"
      "  -v        ... add vorticity w=curl(u)\n"
      "  -G        ... add pressure and vorticity gradients\n"
      "\n"
      "Default behaviour if no option flags given: Calculate all\n";
  char buf[StrMax];
  while (--argc && **++argv == '-')
  {
    switch (*++argv[0])
    {
    case 'h':
      sprintf(buf, usage, prog);
      cout << buf;
      exit(EXIT_SUCCESS);
    case 's':
      if (*++argv[0])
        session = *argv;
      else
      {
        --argc;
        session = *++argv;
      }
    case 'v':
      flag[VORTICITY] = true;
      break;
    case 'G':
      flag[VORTICITY] = true;
      break;
    default:
      sprintf(buf, usage, prog);
      cout << buf;
      exit(EXIT_FAILURE);
      break;
    }
  }

  //
  // Check args
  //
  int_t sum;
  for (int_t i = 0; i < FLAG_MAX; i++)
  {
    sum += (flag[i]) ? 1 : 0;
  }
  if (!sum)
    flag[VORTGEN] = true;

  if (!session)
    message(prog, "no session file", ERROR);

  if (argc != 1)
    message(prog, "no field file", ERROR);
  else
    dump = *argv;
}

//
// Cycle through field dump
//
static bool getDump(Domain *D, ifstream &dump)
{
  dump >> *D;
  return dump.good();
}
