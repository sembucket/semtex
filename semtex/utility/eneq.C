//////////////////////////////////////////////////////////////////////////////
// eneq.C: from a field file containing correlations, compute terms of 
// the fluctuating flow energy equation, Tennekes & Lumley (3.2.1).
//
// Copyright (c) 2004 <--> $Date$, Hugh Blackburn
//
// NB: the input field file is assumed to contain only velocity and
// appropriate correlation data. The naming conventions employed for
// example in addfield.C are broken (names are here used for different
// variables).
//
// Usage:
// -----
// eneq [options] session session.avg
//   options:
//   -h        ... print this message
//
// Output to standard output.
//
// Input field names:
// ------------------
//
// u -- x velocity component (cylindrical: axial)
// v -- y velocity component (cylindrical: radial)
// w -- z velocity component (cylindrical: azimuthal)
// p -- pressure/density
//
// Names for components of the symmetric "Reynolds stress" correlations:
//
//                      / uu uv uw \     /  A  B  D \
//                      | .  vv vw |  =  |  .  C  E |
//                      \ .  .  ww /     \  .  .  F /
//
// Names for correlations specific to the energy equation:
// 
// a) Scalar: 
//    i) q = 0.5 [u^2 + v^2 (+ w^2)]
//   ii) d = sijsij
//
// b) Vector:
//    i) p u_i
//                      / pu \   / m \
//                      | pv | = | n |
//                      \ pw /   \ o /
//   ii) q u_i
//                      / qu \   / r \
//                      | qv | = | s |
//                      \ qw /   \ t /
//
//   iii) Sij u_j       / SxxU + SxyV + SxzW \   / a \
//                      | SyxU + SyyV + SxzW | = | b |
//                      \ SzxU + SzyV + SzzW /   \ c / 
// 
// c) Tensor: symmetric rate-of-strain tensor sij:
//
//                      / xx xy xz \     /  G  H  J \
//                      | .  yy yz |  =  |  .  I  K |
//                      \ .  .  zz /     \  .  .  L /
//
// Names for (output) terms in the fluctutating energy equation:
// -------------------------------------------------------------

// In the case that we consider non-Newtonian flow, there are ten
// terms to compute.  For Newtonian flow we only need six.  We will
// use arabic numeric characters to name these terms.
//
// 1: Uj dq/dxj
// 2: d(q.uj)/dxj
// 3: d(1/rho p.uj)/dxj
// 4: -2 kinvis d(sij.ui)/dxj
// 7: 2 kinvis sij.sij
// 0: ui.uj.Sij
//
// Note: the "missing" names 5,6,8,9 are for non-Newtonian flow.  The
// above terms should sum to zero for stationary turbulence, so we
// also write out
//
// S: sum of above terms.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>

static char prog[] = "eneq";
static void getargs (int,char**,const char*&,const char*&);
static void getMesh (const char*,vector<Element*>&);
static void makeBuf (vector<AuxField*>&,vector<AuxField*>&,vector<Element*>&);
static bool getDump (ifstream&,map<char, AuxField*>&,vector<Element*>&);
static bool doSwap  (const char*);
static void covary  (map<char, AuxField*>&,vector<AuxField*>&,
		     vector<AuxField*>&);

static const char* fieldNames(map<char, AuxField*>&);

int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver -- adapted from probe.C.
// ---------------------------------------------------------------------------
{
  const char           *session, *dump;
  ifstream             avgfile;
  vector<Element*>     elmt;
  vector<AuxField*>    outbuf, work;
  map<char, AuxField*> input;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, dump);

  avgfile.open (dump, ios::in);
  if (!avgfile) message (prog, "no field file", ERROR);
  
  getMesh (session, elmt);
  makeBuf (outbuf, work, elmt);
  
  while (getDump (avgfile, input, elmt)) {
    covary     (input, outbuf, work);
    writeField (cout, session, 0, 0.0, outbuf);
  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int          argc   ,
		     char**       argv   ,
		     const char*& session,
		     const char*& dump   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] =
    "Usage: %s [options] session dump.avg\n"
    "options:\n"
    "  -h ... print this message \n";
              
  char buf[StrMax];
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog); cerr << buf; exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog); cerr << buf; exit (EXIT_FAILURE);
      break;
    }

  if (argc != 2) {
    sprintf (buf, usage, prog); cerr << buf; exit (EXIT_FAILURE);
  } else {
    --argc; session = *++argv;
    --argc; dump    = *++argv;
  }
}


static void getMesh (const char*       session,
		     vector<Element*>& elmt   )
// ---------------------------------------------------------------------------
// Set up 2D mesh information. Note that parser tokens and Geometry
// are set here, too.
// ---------------------------------------------------------------------------
{
  FEML* F = new FEML (session);
  Mesh* M = new Mesh (F);
  
  const int_t nel = M -> nEl();  
  const int_t np  = Femlib::ivalue ("N_P");
  const int_t nz  = Femlib::ivalue ("N_Z");

  Geometry::CoordSys space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);
  elmt.resize   (nel);

  for (int_t k = 0; k < nel; k++) elmt[k] = new Element (k, np, M);
}


static void makeBuf (vector<AuxField*>& outbuf,
		     vector<AuxField*>& work  ,
		     vector<Element*>&  elmt  )
// ---------------------------------------------------------------------------
// Note that we only set up the output and work buffers here. The
// input buffers get created in getDump.
// ---------------------------------------------------------------------------
{
  const int_t nz   = Geometry::nZ();
  const int_t ntot = Geometry::nTotal();

  outbuf.resize (7);
  work.resize   (2);

  outbuf[0] = new AuxField (new real_t[ntot], nz, elmt, '0');
  outbuf[0] = new AuxField (new real_t[ntot], nz, elmt, '1');
  outbuf[0] = new AuxField (new real_t[ntot], nz, elmt, '2');
  outbuf[0] = new AuxField (new real_t[ntot], nz, elmt, '3');
  outbuf[0] = new AuxField (new real_t[ntot], nz, elmt, '4');
  outbuf[0] = new AuxField (new real_t[ntot], nz, elmt, '7');
  outbuf[0] = new AuxField (new real_t[ntot], nz, elmt, 'S');

  work[0]   = new AuxField (new real_t[ntot], nz, elmt, '\0');
  work[1]   = new AuxField (new real_t[ntot], nz, elmt, '\0');
}


static const char* fieldNames (map<char, AuxField*>& u)
// ---------------------------------------------------------------------------
// Return string containing single-character names of fields.
// ---------------------------------------------------------------------------
{
  int_t i;
  static char buf[StrMax];
  map<char, AuxField*>::iterator k;

  for (i = 0, k = u.begin(); k != u.end(); k++, i++) buf[i] = k -> first;
  buf[i] = '\0';

  return buf;
}



static bool getDump (ifstream&             file,
		     map<char, AuxField*>& u   ,
		     vector<Element*>&     elmt)
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  const int_t ntot = Geometry::nTotal();
  char        buf[StrMax], fields[StrMax];
  int_t       i, nf, np, nz, nel;
  bool        swab;

  if (file.getline(buf, StrMax).eof()) return false;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> np >> nz >> nz >> nel;
  file.getline (buf, StrMax);
  
  if (np != Geometry::nP() || nz != Geometry::nZ() || nel != Geometry::nElmt())
    message (prog, "size of dump mismatch with session file", ERROR);

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  nf = strlen  (fields);
  file.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);

  // -- Create AuxFields on first pass.

  if (u.size() == 0) {
    for (i = 0; i < nf; i++)
      u[fields[i]] = new AuxField (new real_t[ntot], nz, elmt, fields[i]);
  } else if (strcmp (fieldNames (u), fields) != 0)
    message (prog, "fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (map<char, AuxField*>::iterator k = u.begin(); k != u.end(); k++) {
    file >> *(k -> second);
    if (swab) k -> second -> reverse();
  }

  return file.good();
}


static bool doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping is required to make sense of binary input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}


static void covary  (map<char, AuxField*>& in  ,
		     vector<AuxField*>&    out ,
		     vector<AuxField*>&    work)
// ---------------------------------------------------------------------------
// This does the actual work of building energy equation terms.
// ---------------------------------------------------------------------------
{
  const char   list2d[] = "uvpABCqdmnrsabGHI";
  const char   list3d[] = "uvwpABCDEFqdmnorstabcGHIJKL";
  const char*  names    = fieldNames (in);
  const int_t  nvel     = (strchr (names, 'w')) ? 3 : 2;
  const real_t kinvis   = Femlib::value ("KINVIS");
  char         err[StrMax];

  if (nvel == 2 && !(strcmp (names, list2d) != 0)) {
    sprintf (buf,"list of names should match %s: have %s", list2d, names);
    message (prog, err, ERROR);
  } else if (!(strcmp (names, list3d) != 0)) {
    sprintf (buf,"list of names should match %s: have %s", list3d, names);
    message (prog, err, ERROR);
  }

  // -- Deal with all the correlations first, before any differentiation.

  // -- Turn ABC... into Reynolds stresses and q into a covariance.

  work[0] -> times (*in['u'], *in['u']); *in['A'] -= *work[0];
  *work[0] *= 0.5; *in['q'] -= *work[0]; *work[1]  = *work[0];
  work[0] -> times (*in['u'], *in['v']); *in['B'] -= *work[0];
  work[0] -> times (*in['v'], *in['v']); *in['C'] -= *work[0];
  *work[0] *= 0.5; *in['q'] -= *work[0]; *work[1] += *work[0];
  if (nvel == 3) {
    work[0] -> times (*in['u'], *in['w']); *in['D'] -= *work[0];
    work[0] -> times (*in['v'], *in['w']); *in['E'] -= *work[0];
    work[0] -> times (*in['w'], *in['w']); *in['F'] -= *work[0];
    *work[0] *= 0.5; *in['q'] -= *work[0]; *work[1] += *work[0];
  }

  // -- At this stage, work[1] holds the KE of the time-average flow.
  //    Use it to deal with q.ui: r, s (t).

  work[0] -> times (*work[1], *in['u']); *in['r'] -= *work[0];
  work[0] -> times (*work[1], *in['v']); *in['s'] -= *work[0];
  if (nvel == 3) {
    work[0] -> times (*work[1], *in['w']); *in['t'] -= *work[0];
  }

  // -- The pressure-velocity covariances: m, n (o).

  work[0] -> times (*in['p'], *in['u']); *in['m'] -= *work[0];
  work[0] -> times (*in['p'], *in['u']); *in['n'] -= *work[0];
  if (nvel == 3) {
    work[0] -> times (*in['p'], *in['u']); *in['o'] -= *work[0];
  }

  // -- Everything from now on involves derivatives in one way or another.

  // -- So let's get started by making term '1':

  *work[0]  = *in['q']; work[0]   -> gradient (0);
  *work[0] *= *in['u']; *out['1']  = *work[0];
  *work[0]  = *in['q']; work[0]   -> gradient (1);
  *work[0] *= *in['v']; *out['1'] += *work[0];
  if (nvel == 3) {
   (*work[0]  = *in['q']).transform(FORWARD).gradient(2).transform(INVERSE);
    *work[0] *= *in['w']; *out['1'] += *work[0];
  }

  // -- Term '2'. We can destroy 'r', 's', 't' along the way:
  
  in['r'] -> gradient (0); *out['2']  = *in['r'];
  in['s'] -> gradient (1); *out['2'] += *in['s'];
  if (nvel == 3) {
    (in['t']  -> transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['2'] += *in['t'];
  }

  // -- Term '3', destroying 'm', 'n', 'o' as we go:

  in['m'] -> gradient (0); *out['3']  = *in['m'];
  in['n'] -> gradient (1); *out['3'] += *in['n'];
  if (nvel == 3) {
    (in['o'] ->  transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['3'] *= *in['o'];
  }

  // -- Next we build the mean rate-of-strain tensor for further work.
  //    The '2d' components get stored in 'm', 'n', 'o' while 
  //    the '3d' ones, if needed,   go in 'r', 's', 't'.

  (*in['m']  = *in['u']) . gradient (0);
  (*in['o']  = *in['v']) . gradient (1);
  (*in['n']  = *in['u']) . gradient (1);
  (*work[0]  = *in['v']) . gradient (0);
  (*in['n'] += *work[0]) *= 0.5;

  if (nvel == 3) {
    (*in['r']  = *in['u']).transform(FORWARD).gradient(2).transform(INVERSE);
    (*work[0]  = *in['w']).gradient(0);
    (*in['r'] += *work[0]) *= 0.5;
    (*in['s']  = *in['v']).transform(FORWARD).gradient(2).transform(INVERSE);
    (*work[0]  = *in['w']).gradient(1);
    (*in['s'] += *work[0]) *= 0.5;
    (*in['t']  = *in['w']).transform(FORWARD).gradient(2).transform(INVERSE);
  }

  // -- Reduce the rate-of-strain correlations to covariances.

  *in['G'] -= *in['m'];
  *in['H'] -= *in['n'];
  *in['I'] -= *in['o'];

  if (nvel == 3) {
    *in['J'] -= *in['r'];
    *in['K'] -= *in['s'];
    *in['L'] -= *in['t'];
  }

  // -- Compute term '4':

  work[0] -> times (*in['m'], *in['u']);
  work[1] -> times (*in['n'], *in['v']);
  *work[0] += *work[1];
  if (nvel == 3) {
    work[1] -> times (*in['r'], *in['w']);
    *work[0] += *work[1];
  }
  (*in['a'] -= *work[0]) . gradient (0);
  *out['4']  = *in['a'];

  work[0] -> times (*in['n'], *in['u']);
  work[1] -> times (*in['o'], *in['v']);
  *work[0] += *work[1];
  if (nvel == 3) {
    work[1] -> times (*in['s'], *in['w']);
    *work[0] += *work[1];
  }
  (*in['b'] -= *work[0]) . gradient (1);
  *out['4'] += *in['b'];

  if (nvel == 3) {
    work[0] -> times (*in['r'], *in['u']);
    work[1] -> times (*in['s'], *in['v']);
    *work[0] += *work[1];
    work[1] -> times (*in['t'], *in['w']);
    *work[0] += *work[1];
    (*in['c'] -= *work[0]).transform(FORWARD).gradient(2).transform(INVERSE);
    *out['4'] += *in['c'];
  }
  
  *out['4'] *= -2.0 * kinvis;
  
  // -- Compute term '7':

  *in['H'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['J'] *= sqrt (2.0);
    *in['K'] *= sqrt (2.0);
  }

  *out['7'] = *in['d'];
  out['7'] -> timesMinus (*in['G'], *in['G']);
  out['7'] -> timesMinus (*in['H'], *in['H']);
  out['7'] -> timesMinus (*in['I'], *in['I']);
  if (nvel == 3) {
    out['7'] -> timesMinus (*in['J'], *in['J']);
    out['7'] -> timesMinus (*in['K'], *in['K']);
    out['7'] -> timesMinus (*in['L'], *in['L']);
  }
  *out['7'] *= 2.0 * kinvis;

  // -- Compute term '0':

  *in['B'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['D'] *= sqrt (2.0);
    *in['E'] *= sqrt (2.0);
  }
 
  out['0'] -> times     (*in['A'], *in['G']);
  out['0'] -> timesPlus (*in['B'], *in['H']);
  out['0'] -> timesPlus (*in['C'], *in['I']);
  if (nvel == 3) {
    out['0'] -> timesPlus (*in['D'], *in['J']);
    out['0'] -> timesPlus (*in['E'], *in['K']);
    out['0'] -> timesPlus (*in['F'], *in['L']);
  }
}

