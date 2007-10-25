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
// Also: the formulation has not been checked to work for cylindrical coords.
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
// -------------------------------------------------------
// a) Scalar: 
//    i) q = 0.5 [u^2 + v^2 (+ w^2)]
//   ii) d = SijSij
/// iii) e = l SijSij
//
// b) Vector: Naming:
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
//    iv) l u_i         / lu \   / f \
//                      | lv | = | g |
//                      \ lw /   \ h /
//
//     v) l Sij u_i     / la \   / i \
//                      | lb | = | j |
//                      \ lc /   \ k /
// 
// c) Tensor: Naming:
//     i) Sij           / xx xy xz \     /  G  H  J \
//                      | .  yy yz |  =  |  .  I  K |
//                      \ .  .  zz /     \  .  .  L /
//
//
//    ii) l Sij         / xx xy xz \     /  M  N  P \
//                      | .  yy yz |  =  |  .  O  Q |
//                      \ .  .  zz /     \  .  .  R /
//
// Names for (output) terms in the fluctutating energy equation:
// -------------------------------------------------------------
// In the case that we consider non-Newtonian flow, there are ten
// terms to compute.  We will use arabic numeric characters to name
// these terms.
//
// 1: Uj dq/dxj
// 2: -d(q.uj)/dxj
// 3: -d(1/rho p.uj)/dxj
// 4: 2 d(L sij.ui)/dxj
// 5: 2 d(l ui Sij)/dxj ** NEW
// 6: 2 d(l sij ui)/dxj ** NEW
// 7: -2 L sij.sij
// 8: -2 l sij.Sij      ** NEW
// 9: -2 l sij.sij      ** NEW
// 0: -ui.uj.Sij
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
static void makeBuf (map<char,AuxField*>&,vector<AuxField*>&,
		     vector<Element*>&);
static bool getDump (ifstream&,map<char, AuxField*>&,vector<Element*>&);
static bool doSwap  (const char*);
static void covary  (map<char,AuxField*>&,map<char,AuxField*>&,
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
  map<char, AuxField*> input, output;
  vector<AuxField*>    work,  outbuf;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, dump);

  avgfile.open (dump, ios::in);
  if (!avgfile) message (prog, "no field file", ERROR);
  
  getMesh (session, elmt);
  makeBuf (output, work, elmt);

  // -- Need to link outbuf to output so we can use writeField.
  //    Maybe in the longer term we should overload writeField.

  outbuf.resize (output.size()); int_t i = 0;
  for (map<char,AuxField*>::iterator k = output.begin();
       k != output.end(); k++, i++) outbuf[i] = k -> second;
  
  while (getDump (avgfile, input, elmt)) {
    covary     (input, output, work);
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
    --argc; session = *argv;
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


static void makeBuf (map<char, AuxField*>& output,
		     vector<AuxField*>&    work  ,
		     vector<Element*>&     elmt  )
// ---------------------------------------------------------------------------
// Note that we only set up the output and work buffers here. The
// input buffers get created in getDump.
// ---------------------------------------------------------------------------
{
  const int_t nz   = Geometry::nZ();
  const int_t ntot = Geometry::nTotal();

  work.resize (2);

  output['0'] = new AuxField (new real_t[ntot], nz, elmt, '0');
  output['1'] = new AuxField (new real_t[ntot], nz, elmt, '1');
  output['2'] = new AuxField (new real_t[ntot], nz, elmt, '2');
  output['3'] = new AuxField (new real_t[ntot], nz, elmt, '3');
  output['4'] = new AuxField (new real_t[ntot], nz, elmt, '4');
  output['5'] = new AuxField (new real_t[ntot], nz, elmt, '5');
  output['6'] = new AuxField (new real_t[ntot], nz, elmt, '6');
  output['7'] = new AuxField (new real_t[ntot], nz, elmt, '7');
  output['8'] = new AuxField (new real_t[ntot], nz, elmt, '8');
  output['9'] = new AuxField (new real_t[ntot], nz, elmt, '9');
  output['S'] = new AuxField (new real_t[ntot], nz, elmt, 'S');

  work[0]     = new AuxField (new real_t[ntot], nz, elmt, '\0');
  work[1]     = new AuxField (new real_t[ntot], nz, elmt, '\0');
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
		     map<char, AuxField*>& out ,
		     vector<AuxField*>&    work)
// ---------------------------------------------------------------------------
// This does the actual work of building energy equation terms.
// ---------------------------------------------------------------------------
{
  const char   list2d[] = "ABCGHIIMNOabdefgijmnpqrsuv";
  const char   list3d[] = "ABCDEFGHIJKLMNOPQRabcdefghijklmnopqrstuvw";
  const char*  names    = fieldNames (in);
  const int_t  nvel     = (strchr (names, 'w')) ? 3 : 2;
  const real_t kinvis   = Femlib::value ("KINVIS");
  char         err[StrMax];

  if (nvel == 2) {
    if (strcmp (names, list2d) != 0) {
      sprintf (err,"list of names should match %s: have %s", list2d, names);
      message (prog, err, ERROR);
    } 
  } else if (strcmp (names, list3d) != 0) {
    sprintf (err,"list of names should match %s: have %s", list3d, names);
    message (prog, err, ERROR);
  }

  // -- Deal with all the correlations first, before any differentiation.

  // -- Turn ABC... into Reynolds stresses.

  in['A'] -> timesMinus (*in['u'], *in['u']);
  in['B'] -> timesMinus (*in['u'], *in['v']);
  in['C'] -> timesMinus (*in['v'], *in['v']);
  if (nvel == 3) {
    in['D'] -> timesMinus (*in['u'], *in['w']);
    in['E'] -> timesMinus (*in['v'], *in['w']);
    in['F'] -> timesMinus (*in['w'], *in['w']);
  }

  // -- At this stage, in['q'] still holds the total kinetic energy.

  in['r'] -> timesMinus (*in['A'], *in['u']);
  in['r'] -> timesMinus (*in['B'], *in['v']);
  in['r'] -> timesMinus (*in['q'], *in['u']);
  in['s'] -> timesMinus (*in['B'], *in['u']);
  in['s'] -> timesMinus (*in['C'], *in['v']);
  in['s'] -> timesMinus (*in['q'], *in['v']);
  if (nvel == 3) {
    in['r'] -> timesMinus (*in['D'], *in['w']);
    in['s'] -> timesMinus (*in['E'], *in['w']);
    in['t'] -> timesMinus (*in['D'], *in['u']);
    in['t'] -> timesMinus (*in['E'], *in['v']);
    in['t'] -> timesMinus (*in['F'], *in['w']);
    in['t'] -> timesMinus (*in['q'], *in['w']);
  }

  // -- Rework in['q'] so it holds the fluctuating KE.

  *in['q'] = *in['A']; *in['q'] += *in['C'];
  if (nvel == 3) *in['q'] += *in['F'];
  *in['q'] *= 0.5;

  // -- The pressure-velocity covariances: m, n (o).

  in['m'] -> timesMinus (*in['p'], *in['u']);
  in['n'] -> timesMinus (*in['p'], *in['v']);
  if (nvel == 3) in['o'] -> timesMinus (*in['p'], *in['w']);

  // -- Everything from now on involves derivatives in one way or another.

  // -- So let's get started by making term '1' (advection):

  *work[0]  = *in['q']; work[0]   -> gradient (0);
  *work[0] *= *in['u']; *out['1']  = *work[0];
  *work[0]  = *in['q']; work[0]   -> gradient (1);
  *work[0] *= *in['v']; *out['1'] += *work[0];
  if (nvel == 3) {
   (*work[0]  = *in['q']).transform(FORWARD).gradient(2).transform(INVERSE);
    *work[0] *= *in['w']; *out['1'] += *work[0];
  }

  // -- Term '2' (turbulent transport). Destroy 'r', 's', 't' along the way:
  
  in['r'] -> gradient (0); *out['2']  = *in['r'];
  in['s'] -> gradient (1); *out['2'] += *in['s'];
  if (nvel == 3) {
    (in['t']  -> transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['2'] += *in['t'];
  }

  *out['2'] *= -1.0;

  // -- Term '3' (pressure work), destroying 'm', 'n', 'o' as we go:

  in['m'] -> gradient (0); *out['3']  = *in['m'];
  in['n'] -> gradient (1); *out['3'] += *in['n'];
  if (nvel == 3) {
    (in['o'] ->  transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['3'] += *in['o'];
  }

  *out['3'] *= -1.0;

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

  // -- Compute term '4' (viscous transport):

  work[0] -> times (*in['m'], *in['u']);
  work[1] -> times (*in['n'], *in['v']);
  *work[0] += *work[1];
  if (nvel == 3) {
    work[1] -> times (*in['r'], *in['w']);
    *work[0] += *work[1];
  }
  (*in['a'] -= *work[0]) *= *in['l'];
  (*work[0]  = *in['a']) . gradient (0);
  *out['4']  = *work[0];

  work[0] -> times (*in['n'], *in['u']);
  work[1] -> times (*in['o'], *in['v']);
  *work[0] += *work[1];
  if (nvel == 3) {
    work[1] -> times (*in['s'], *in['w']);
    *work[0] += *work[1];
  }
  (*in['b'] -= *work[0]) *= *in['l'];
  (*work[0]  = *in['b']) . gradient (1);
  *out['4'] += *work[0];

  if (nvel == 3) {
    work[0] -> times (*in['r'], *in['u']);
    work[1] -> times (*in['s'], *in['v']);
    *work[0] += *work[1];
    work[1] -> times (*in['t'], *in['w']);
    *work[0] += *work[1];
    (*in['c'] -= *work[0]) *= *in['l'];
    (*work[0]  = *in['c']).transform(FORWARD).gradient(2).transform(INVERSE);
    *out['4'] += *work[0];
  }
  
  *out['4'] *= 2.0;

  // -- Compute term '5':
  
  in['f'] -> timesMinus (*in['l'], *in['u']);
  in['g'] -> timesMinus (*in['l'], *in['v']);
  if (nvel == 3) in['h'] -> timesMinus (*in['l'], *in['w']);

  work[0] -> times     (*in['m'], *in['f']);
  work[0] -> timesPlus (*in['n'], *in['g']);
  if (nvel == 3) work[0] -> timesPlus (*in['r'], *in['h']);
  work[0] -> gradient (0);
  *out['5'] = *work[0];

  work[0] -> times     (*in['n'], *in['f']);
  work[0] -> timesPlus (*in['o'], *in['g']);
  if (nvel == 3) work[0] -> timesPlus (*in['s'], *in['h']);
  work[0] -> gradient (1);
  *out['5'] += *work[0];

  if (nvel == 3) {
    work[0] -> times     (*in['r'], *in['f']);
    work[0] -> timesPlus (*in['s'], *in['g']);
    work[0] -> timesPlus (*in['t'], *in['h']);
    work[0] -> transform(FORWARD).gradient(2).transform(INVERSE);
    *out['5'] += *work[0];
  }
  
  *out['5'] *= 2.0;

  // -- Now, "out of order", compute term '8':

  in['M'] -> timesMinus (*in['l'], *in['m']);
  in['N'] -> timesMinus (*in['l'], *in['n']);
  in['O'] -> timesMinus (*in['l'], *in['o']);
  if (nvel == 3) {
    in['P'] -> timesMinus (*in['l'], *in['r']);
    in['Q'] -> timesMinus (*in['l'], *in['s']);
    in['R'] -> timesMinus (*in['l'], *in['t']);
  }

  out['8'] -> times     (*in['M'], *in['m']);
  out['8'] -> timesPlus (*in['O'], *in['o']);
  if (nvel == 3) out['8'] -> timesPlus (*in['R'], *in['t']);
  work[0] -> times (*in['N'], *in['n']);
  if (nvel == 3) {
    work[0] -> timesPlus (*in['P'], *in['r']);
    work[0] -> timesPlus (*in['Q'], *in['s']);
  }
  *work[0] *= 2.0;
  (*out['8'] += *work[0]) *= -2.0;

  // -- Compute term '6':

  *in['i'] -= *in['a'];
  *in['j'] -= *in['b'];
  if (nvel == 3) *in['k'] -= *in ['c'];
  in['i'] -> timesMinus (*in['M'], *in['u']);
  in['i'] -> timesMinus (*in['N'], *in['v']);
  if (nvel == 3) in['i'] -> timesMinus (*in['P'], *in['w']);
  in['j'] -> timesMinus (*in['N'], *in['u']);
  in['j'] -> timesMinus (*in['O'], *in['v']);
  if (nvel == 3) in['j'] -> timesMinus (*in['Q'], *in['w']);
  if (nvel == 3) {
    in['k'] -> timesMinus (*in['P'], *in['u']);
    in['k'] -> timesMinus (*in['Q'], *in['v']);
    in['k'] -> timesMinus (*in['R'], *in['w']);
  }
  in['i'] -> timesMinus (*in['m'], *in['f']);
  in['i'] -> timesMinus (*in['n'], *in['g']);
  if (nvel == 3) in['i'] -> timesMinus (*in['r'], *in['h']);
  in['j'] -> timesMinus (*in['n'], *in['f']);
  in['j'] -> timesMinus (*in['o'], *in['g']);
  if (nvel == 3) in['j'] -> timesMinus (*in['s'], *in['h']);
  if (nvel == 3) {
    in['k'] -> timesMinus (*in['r'], *in['f']);
    in['k'] -> timesMinus (*in['s'], *in['g']);
    in['k'] -> timesMinus (*in['t'], *in['h']);
  }
  work[0] -> times     (*in['m'], *in['u']);
  work[0] -> timesPlus (*in['n'], *in['v']);
  if (nvel == 3) work[0] -> timesPlus (*in['r'], *in['w']);
  *work[0] *= *in['l'];
  *in['i'] -= *work[0];
  work[0] -> times     (*in['n'], *in['u']);
  work[0] -> timesPlus (*in['o'], *in['v']);
  if (nvel == 3) work[0] -> timesPlus (*in['s'], *in['w']);
  *work[0] *= *in['l'];
  *in['j'] -= *work[0];
  if (nvel == 3) {
      work[0] -> times     (*in['r'], *in['u']);
      work[0] -> timesPlus (*in['s'], *in['v']);
      work[0] -> timesPlus (*in['t'], *in['w']);
      *work[0] *= *in['l'];
      *in['k'] -= *work[0];
  }
  in['i'] -> gradient (0);
  in['j'] -> gradient (1);
  if (nvel == 3) in['k'] -> transform(FORWARD).gradient(2).transform(INVERSE);
  *out['6']  = *in['i'];
  *out['6'] += *in['j'];
  if (nvel == 3) *out['6'] += *in['k'];
  
  *out['6'] *= 2.0;
  
  // -- Compute term '7' (dissipation):

  *in['n'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['r'] *= sqrt (2.0);
    *in['s'] *= sqrt (2.0);
  }

  in['d'] -> timesMinus (*in['m'], *in['m']);
  in['d'] -> timesMinus (*in['n'], *in['n']);
  in['d'] -> timesMinus (*in['o'], *in['o']);
  if (nvel == 3) {
    in['d'] -> timesMinus (*in['r'], *in['r']);
    in['d'] -> timesMinus (*in['s'], *in['s']);
    in['d'] -> timesMinus (*in['t'], *in['t']);
  }
  (*out['7'] = *in['d']) *= *in['l'];
  *out['7'] *= -2.0;

  // -- Compute term 9:
  //    Note we have already premultiplied 'n', 'r', 's' by sqrt(2).

  work[0] -> times (*in['N'], *in['n']);
  if (nvel == 3) {
    work[0] -> timesPlus (*in['P'], *in['r']);
    work[0] -> timesPlus (*in['Q'], *in['s']);
  }
  *work[0] *= sqrt (2.0);
  work[0] -> timesPlus (*in['M'], *in['m']);
  work[0] -> timesPlus (*in['O'], *in['o']);
  if (nvel == 3) work[0] -> timesPlus (*in['R'], *in['t']);
  *work[0] *= 2.0;
  *in['e'] -= *work[0];
  
  in['e'] -> timesMinus (*in['d'], *in['l']);
  
  work[0] -> times (*in['n'], *in['n']);
  if (nvel == 3) {
    work[0] -> timesPlus (*in['r'], *in['r']);
    work[0] -> timesPlus (*in['s'], *in['s']);
  }
  work[0] -> timesPlus (*in['m'], *in['m']);
  work[0] -> timesPlus (*in['o'], *in['o']);
  if (nvel == 3) work[0] -> timesPlus (*in['t'], *in['t']);
  *work[0] *= *in['l'];

  *in['e'] -= *work[0];

  (*out['9'] = *in['e']) *= -2.0;
  
  // -- Compute term '0' (production):

  *in['n'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['r'] *= sqrt (2.0);
    *in['s'] *= sqrt (2.0);
  }
 
  out['0'] -> times     (*in['A'], *in['m']);
  out['0'] -> timesPlus (*in['B'], *in['n']);
  out['0'] -> timesPlus (*in['C'], *in['o']);
  if (nvel == 3) {
    out['0'] -> timesPlus (*in['D'], *in['r']);
    out['0'] -> timesPlus (*in['E'], *in['s']);
    out['0'] -> timesPlus (*in['F'], *in['t']);
  }

  *out['0'] *= -1.0;

  // -- Finally, compute the sum, 'S':

 (*out['S']  = *out['1']) *= -1.0;
  *out['S'] += *out['2'];
  *out['S'] += *out['3'];
  *out['S'] += *out['4'];
  *out['S'] += *out['5'];
  *out['S'] += *out['6'];
  *out['S'] += *out['7'];
  *out['S'] += *out['8'];
  *out['S'] += *out['9'];
  *out['S'] += *out['0'];
}

