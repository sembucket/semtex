#include <sem.h>
#include <tensorcalcs.h>

#define FLAG_MAX 2

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
// Vorticity generation
//

void vortgen(void)
{
}