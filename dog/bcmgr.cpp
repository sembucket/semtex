///////////////////////////////////////////////////////////////////////////////
// bcmgr.cpp: class functions for managing boundary conditions.
//
// Note for dog/lns: the new "robust" outflow boundary type
// <O>/outflow (see Ref 3 below) now included in semtex/dns is
// presently (as of April 2015) disallowed since it needs nonlinear
// products (such as u^2) for proper evaluation. For 3D modes, there
// would be no contribution at the wavenumber under
// consideration. Rather than only allow this BC type only for 2D
// analysis, it seems easier to remove it completely.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// SYNOPSIS
// --------
// BCmgr manufactures and stores instances of classes derived from the
// Condition virtual base class, which are subsequently used by the
// Boundary class to apply boundary conditions.
//
// INPUT
// -----
// The following sections are optionally present in FEML file: they
// are needed if the <SURFACES> section used any boundary conditions,
// but it is possible they could be absent if all boundaries are
// periodic.  The sections below are given as examples.
//
// <GROUPS NUMBER=5>

// #	tag	name	descriptor
// #    Note that name "axis" is reserved and must be associated with 
// #    the BC character tag 'A'.
// #    The name "velocity" also has a special meaning but has no required
// #    character tag.
// #    The one-character names should be unique.
//      1       a       axis
// 	2	v	velocity
// 	3	w	wall
//      4       h       radiation
//      5       e       exit 
// </GROUPS>
// 
// <BCS NUMBER=5>
// #	tag	group	number, followed by BCs.
//
// #    BCs for the "axis" group are internally assigned according to
//      field name and Fourier mode, see section below and ref [2].
//      1       a       4
//                      <A>     u                       </A>
//                      <A>     v                       </A>
//                      <A>     w                       </A>
//                      <A>     p                       </A>
// #    For all other groups, the user can select arbitrary combinations.
// 	2	v	4
// 			<D>	u = 1.0-4.0*(y-0.5)^2.0 </D>
// 			<D>	v = 0.0                 </D>
// 			<D>	w = 0.0                 </D>
// 			<H>	p                       </H>
// 	3	w	4
// 			<D>	u = 0.0                 </D>
// 			<D>	v = 0.0                 </D>
// 			<D>	w = 0.0                 </D>
// 			<H>	p                       </H>
//      4       h       4
//                      <M>     u = 1.0,2.0             </M>
//                      <M>     v = 1.0,0.5             </M>
//                      <M>     w = 0.5,2.0             </M>
//                      <H>     p                       </H>
//      5       e       4
//                      <N>     u = 0.                  </N>
//                      <N>     v = 0.                  </N>
//                      <N>     w = 0.                  </N>
//                      <D>     p = 0.                  </D>
// </BCS>
//
// <SURFACES NUMBER=6>
// # Surfaces associate a particular element face with a group of BCs (<B>)
// # OR nominate the face as periodic with another (<P>).
// #       tag     elmt    face    type
//         1       1       1       <P>     3       3       </P>
//         2       2       1       <P>     4       3       </P>
//         3       2       2       <B>     a       </B>
//         4       4       2       <B>     e       </B>
//         5       3       4       <B>     v       </B>
//         6       1       4       <B>     v       </B>
// </SURFACES>
// 
// GROUPS
// ------
// Groups provide the capability of associating a string description
// with a character tag, so that user routines can gain additional
// information about the tag.  Descriptions cannot contain whitespace.
//
// Typical (recognized) strings are "velocity", "axis", and "wall";
// the latter denotes that all velocity components are zero: this
// means e.g. that one Condition can gain information about behaviour
// of another Condition which shares the character tag (e.g. that it
// is also a zero-valued essential boundary).
//
// Each group which is used for specification of boundary conditions
// should have an associated descriptor set in the GROUPS section.
// This allows user routines to access boundary value storage areas by
// using the "addForGroup" and "zeroForGroup" routines.
//
// BCS
// ---
// BCs associate specific BC types and values with a particular group.
// Available types are:
//   <D> Dirichlet/essential,
//   <N> Neumann/natural,
//   <M> Mixed
//   <H> Natural pressure BC (no value specified, since it gets computed).
//   <A> Axis BCs for cylindrical coords.  Also, must belong to "axis" group.
//
// The character tags for variables as shown match those used
// internally as Field names, so that the order in which the BCs are
// supplied for each group is arbitrary, however the number must match
// or exceed the number of Field variables in the problem (this means
// that 3D BC specifications can also be used for 2D or 1D problems,
// but not vice-versa).
//
// For boundary value types where a function is supplied after the
// "=", it is necessary that no spaces appear in the function string.
//
// VARIABLE NAMES
// --------------
// The (one character) names of field variables are significant, and have
// the following reserved meanings:
// 
// u:  First velocity component.            (Cylindrical: axial     velocity.)
// v:  Second velocity component.           (Cylindrical: radial    velocity.)
// w:  Third velocity component.            (Cylindrical: azimuthal velocity.)
// p:  Pressure divided by density.
// c:  Scalar for transport or elliptic problems.
//
// CYLINDRICAL COORDINATE SYSTEM AND BCS, see ref [2]
// --------------------------------------------------
//
// Type <A> BCs are used for cylindrical coordinate systems when the
// edge of an element touches the axis (N.B. y (r) negative is
// illegal).  The BCs that are supplied for boundary type <A> depend
// on (a) the field in question (u, v, [w], p, [c]) and (b) the index
// of the Fourier mode.
//
// Fields u (axial velocity), p, & c are have the same treatment on
// the axis --- they have "scalar" type BCs there, due to the need for
// the fields to be single-valued at the axis.
//
// Fields v & w (radial and azimuthal velocities) have their BCs set
// after the fields have been coupled to give v~ and w~ (which
// *decouples* the Helmholtz problems for diffusion).  V & w must have
// a set phase relationship in the first Fourier mode resulting in a
// zero essential BC in the first Fourier mode for v~ and a zero
// natural BC in the first Fourier mode for w~.  In order to be able
// to set BCs for the coupled fields, v & w must have matching BC
// kinds on all remaining boundaries.  This is enforced in the
// enumeration phase.
//
// Summary:
//
//         +------------+------------------------------------------+
//         |  (Coupled) |       Axis BC for Fourier Mode Index     |
//         |   Variable |         0             1           2...   |
//         +------------+------------------------------------------+
//         |   u, p, c  |   du/dr = 0       u   = 0      u  = 0    |
//         |      v~    |     v~  = 0       v~  = 0      v~ = 0    |
//         |      w~    |     w~  = 0    dw~/dr = 0      w~ = 0    |
//         +------------+------------------------------------------+
//
// In order to deal with this modal dependence of BCs, three levels of
// Boundary pointers are maintained for cylindrical Fields,
// corresponding to Fourier modes 0, 1, 2 (and higher), even if there
// are no axial BCs to be applied.  These BC pointers are carried by
// all processors.
//
// The names of the numbering schemes that will be used for cylindrical
// 3D problems with axial BCs are (note case sensitivity):
//         +------------+------------------------------------------+
//         | (Coupled~) | Numbering scheme for Fourier mode index  |
//         |  Variable  |         0             1           2...   |
//         +------------+------------------------------------------+
//         |      u     |         u             U           U      |
//         |      v~    |         v             v           v      |
//         |      w~    |         w             W           w      |
//         |      p     |         p             P           P      |
//         |      c     |         c             C           C      |
//         +------------+------------------------------------------+
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag, JCP 97:414-443 (1991)
// [2] Blackburn & Sherwin, JCP 197:759-778 (2004)
// [3] Dong, Karniadakis & Chryssostomides JCP 261:83-105, (2014)
//
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>


BCmgr::BCmgr (FEML*             file,
	      vector<Element*>& elmt) :
// ---------------------------------------------------------------------------
// This constructor deals with <GROUPS>, <BCS> and <SURFACES> sections
// of FEML file, and loads internal tables for later use by 
// BCmgr::getCondition.
//
// In addition, it reads in prebuilt numbering schemes from file
// session.num for later retrieval.
// ---------------------------------------------------------------------------
  _axis (false)
{
  const char  routine[] = "BCmgr::BCmgr";
  char        buf[StrMax], err[StrMax], tag[StrMax], gat[StrMax];
  char        groupc, fieldc, testc, tagc;
  char*       trailer;
  const char* grpdesc;
  int_t       verbose = Femlib::ivalue ("VERBOSE");
  int_t       i, j, N, id, nbcs;
  Condition*  C;
  CondRecd*   R;

  // -- Load FIELDS.

  if (file -> seek ("FIELDS")) {
    file -> stream().ignore (StrMax, '\n');
    while (file -> stream().peek() == '#') // -- Skip comments.
      file -> stream().ignore (StrMax, '\n');
    i = 0;
    do
      file -> stream() >> buf[i++];
    while (buf[i - 1] != '<' && i < StrMax);
    if (buf[--i] == '<') {
      buf[i] = '\0';
      file -> stream() >> tag;
      if (!(strstr (tag,    "/FIELDS")))
	   message (routine, "FIELDS section not closed", ERROR);
    } else message (routine, "FIELDS section not closed", ERROR);
  } else   message (routine, "FIELDS section not found",  ERROR);

  strcpy ((_fields = new char [strlen (buf) + 1]), buf);

  VERBOSE cout << routine << ": Installing numbering systems ... ";

  buildnum (file -> root(), elmt);

  VERBOSE cout << "done" << endl;

  if (!file -> seek ("GROUPS")) {
    if (verbose)
      message (routine, "no GROUPS, assuming no boundary conditions", WARNING);
    return;
  }

  // -- Load GROUPS.

  N = file -> attribute ("GROUPS", "NUMBER");

  VERBOSE cout << "  Searching for " << N << " GROUPS ... ";

  _group   .resize (N);
  _descript.resize (N);
  
  for (i = 0; i < N; i++) {
    while (file -> stream().peek() == '#') // -- Skip comments.
      file -> stream().ignore (StrMax, '\n');
    file -> stream() >> id >> groupc >> buf;
    _group[i] = groupc;
    strcpy ((_descript[i] = new char [strlen (buf) + 1]), buf);
    if (!strcmp (_descript[i], "axis"))    _axis    = true;
  }
  
  VERBOSE cout << "done" << endl;

  // -- Load BCS.

  N = file -> attribute ("BCS", "NUMBER");

  VERBOSE cout << "  Searching for " << N << " BCS ... ";

  for (i = 0; i < N; i++) {

    while (file -> stream().peek() == '#') // -- Skip comments.
      file -> stream().ignore (StrMax, '\n');

    file -> stream() >> id >> groupc >> nbcs;

    grpdesc = groupInfo (groupc);	// -- Ensure matching group exists.

    for (j = 0; j < nbcs; j++) {

      // -- Open tag.

      file -> stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag: %s", tag);
	message (routine, err, ERROR);
      }

      // -- Decide if this is a value or function BC.

      file -> stream() >> fieldc;
      file -> stream() >> testc;
      if (testc == '=') { file -> stream() >> buf; strtod (buf, &trailer); }
      if (testc == '<') file -> stream().putback (testc);

      // -- Create appropriate derived Condition structure.
      //    We make at least one per BC group+field variable for 
      //    later retrieval.

      // -- The two primary kinds of BC are Dirichlet/Essential and
      //    Neumann/Natural.  Mixed is a linear combination of these.
      //    Those are all BCs defined by a declared value or function
      //    string that can be evaluated.

      // -- The Axis type amounts to either homogenous Dirichlet or
      //    Neumann depending on variable name and Fourier mode index,
      //    as explained in [2].

      // -- The remaining BC types are computed internally but again
      //    they come out to be either Dirichlet or Neumann types.

      // -- High-order pressure (H) is an internally-computed Neumann BC,
      //    explained in [1].

      // -- Outflow boundaries (O) comprise a set of computed Neumann
      //    (for velocity) or computed Dirichlet (for pressure) as
      //    outlined in [3], supplemented by homogeneous Neumann (for
      //    scalar). Since they appear together, the code allocates
      //    all of them and decides which one to apply for a given
      //    field variable.

      switch (tagc) {

      case 'D': case 'E':	// -- Dirichlet/Essential BC.
	if (testc != '=') {
	  sprintf (err, "expected an '=' in setting field '%c' BC", fieldc);
	  message (routine, err, ERROR);
	}
	if   (*trailer != 0) C = new EssentialFunction (buf);
	else                 C = new Essential         (buf);
	break;

      case 'N':			// -- Neumann/Natural BC.
	if (testc != '=') {
	  sprintf (err, "expected an '=' in setting field '%c' BC", fieldc);
	  message (routine, err, ERROR);
	}
	if   (*trailer != 0) C = new NaturalFunction (buf);
	else                 C = new Natural         (buf);
	break;

      case 'M':			// -- Mixed BC.
	if (testc != '=') {
	  sprintf (err, "expected an '=' in setting field '%c' BC", fieldc);
	  message (routine, err, ERROR);
	}
	if (!(strchr (buf, ';') || strchr (buf, ','))) {
	  sprintf (buf,"can't find multiplier and reference value in: %s",buf);
	  message (routine, buf, ERROR);
	}
	C = new Mixed (buf);
	break;

      case 'A':			// -- Axis BC.

	if (Geometry::system() != Geometry::Cylindrical)
	  message (routine, "axis BCs disallowed in Cartesian coords", ERROR);

	// -- Create two kinds to be retrieved later.
	//    Ensure that the group name is "axis" to aid retrieval.

	if (!strstr (groupInfo (groupc), "axis"))
	  message (routine, "type 'A' BC must belong to group \"axis\"",ERROR);

	strcpy (buf, "0.0");

	C = new Natural (buf);
	
	_cond.insert (_cond.end(), R = new CondRecd);
	R -> grp    = groupc;
	R -> fld    = fieldc;
	R -> bcn    = C;
	strcpy ((R -> value = new char [strlen (buf) + 1]), buf);
	
	C = new Essential (buf);
	break;

      case 'H':		// -- "High Order" computed natural pressure BC.
	if (fieldc != 'p') {
	  sprintf (err, "expected name 'p' with HOPBC, read '%c'", fieldc);
	  message (routine, err, ERROR);
	}
	C = new NaturalCBCp (this);
	break;

      default:
	sprintf (err, "unrecognized BC identifier: %c", tagc);
	message (routine, err, ERROR);
	break;
      }

      // -- Close tag.

      file -> stream() >> gat;

      if (strlen (gat) != 4
	  || gat[0] != '<'
	  || gat[1] != '/' 
	  || gat[2] != tagc
	  || gat[3] != '>') {
	sprintf (err, "close tag %s didn't match open tag %s", gat, tag);
	message (routine, err, ERROR);
      }

      // -- Install new Condition record in internal list.

      _cond.insert (_cond.end(), R = new CondRecd);
      R -> grp = groupc;
      R -> fld = fieldc;
      R -> bcn = C;
      strcpy ((R -> value = new char [strlen (buf) + 1]), buf);
    }
  }

  VERBOSE cout << "done" << endl;

  VERBOSE cout << "  Building internal list of BC edges ... ";

  buildsurf (file, elmt);

  VERBOSE cout << "done" << endl;
}


Condition* BCmgr::getCondition (const char  group,
				const char  field,
				const int_t mode )
// ---------------------------------------------------------------------------
// Search internal list of Conditions for appropriate one, such that
// the stored internal variables match input group & field names.
// Return pointer to the matching Condition.
//
// Have to get into more extended search to locate axis BCs.  Input
// variable 'mode' is only used for cylindrical coordinates: it is the
// Fourier mode number.  See file header above and ref [2] for more
// information.
// ---------------------------------------------------------------------------
{ 
  const char routine[] = "BCmgr::getCondition";
  char       buf[StrMax], err[StrMax], currgrp, currfld;
  CondRecd*  C;
  vector<CondRecd*>::iterator c;
  
  for (c = _cond.begin(); c != _cond.end(); c++) {
    C       = *c;
    currgrp = C -> grp;
    currfld = C -> fld;

    if (currgrp == group && currfld == field) {

      C -> bcn -> describe (buf);

      if (strstr (buf, "axis")) {

	switch (field) {
	case 'u': case 'p': case 'c':
	  if ((mode >  0 && strstr (buf, "essential")) ||
	      (mode == 0 && strstr (buf, "natural"    ))) return C -> bcn;
	  break;
	case 'v':
	  if               (strstr (buf, "essential"))    return C -> bcn;
	  break;
	case 'w':
	  if ((mode != 1 && strstr (buf, "essential")) ||
	      (mode == 1 && strstr (buf, "natural"    ))) return C -> bcn;
	  break;
	default:
	  sprintf (err, "unrecognised field '%c' on axis", field);
	  message (routine, err, ERROR);
	  break;
	} 

      } else 			// -- Default.
	return C -> bcn;
    }
  }

  sprintf (err, "can't find record for group '%c', field '%c'", group, field);
  message (routine, err, ERROR);
  return 0;
}


const char* BCmgr::groupInfo (const char name) const
// ---------------------------------------------------------------------------
// Given a group name, return pointer to string descriptor.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "BCmgr::groupInfo";
  const int_t N = _group.size();
  char        err[StrMax];
  int_t       i;

  for (i = 0; i < N; i++) if (name == _group[i]) return _descript[i];

  sprintf (err, "unknown group: %c", name);
  message (routine, err, WARNING);
  return 0;
}


NumberSys* BCmgr::getNumberSys (const char  name,
				const int_t mode)
// ---------------------------------------------------------------------------
// Return the NumberSys corresponding to the input field name and
// Fourier mode.  Mode numbers begin at zero, and are not expressed
// modulo number of modes per process.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "BCmgr::getNumberSys";
  const int_t nsys  = _numsys.size();
  const int_t cmode = clamp (mode,static_cast<int_t>(0),static_cast<int_t>(2));
  char        err[StrMax], selectname = name;
  NumberSys*  N = 0;

  // -- Switch the name for selected cylindrical coordinate modes.

  if (Geometry::cylindrical() && Geometry::nDim() == 3 && _axis)
    if      (name == 'u' && (cmode == 1 || cmode == 2)) selectname = 'U';
    else if (name == 'w' && (cmode == 1))               selectname = 'W';
    else if (name == 'p' && (cmode == 1 || cmode == 2)) selectname = 'P';
    else if (name == 'c' && (cmode == 1 || cmode == 2)) selectname = 'C';

  // -- Now search.
  
  for (int_t i = 0; i < nsys; i++)
    if (strchr (_numsys[i] -> fields(), selectname)) { N = _numsys[i]; break; }

  if (!N) {
    sprintf (err, "can't find scheme for field %c, mode %1d", name, mode);
    message (routine, err, ERROR);
  }
  
  return N;
}


void BCmgr::buildnum (const char*       session,
		      vector<Element*>& elmt   )
// ---------------------------------------------------------------------------
// Private member function.
//
// Retrieve numbering schemes (btog and bmsk values) from file
// "session.num": this is created by running the "enumerate" utility
// on root processor.
//
// The names of fields and their numbering schemes are significant.
// The convention employed is that the fields have lower-case
// single-character names.  Numbering schemes have the same names,
// *except* in the case of cylindrical coordinate systems where the
// domain includes the symmetry axis.  See top of this file for
// mode-related significance for upper-cased names of numbering
// schemes.
//
// After the numbering schemes have been set up, create the
// corresponding inverse mass matrices.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "BCmgr::buildnum";
  const int_t    np   = Geometry::nP();
  const int_t    NP   = Geometry::nPlane();
  const int_t    nel  = Geometry::nElmt();
  const int_t    npnp = Geometry::nTotElmt();
  const int_t    next = Geometry::nExtElmt();
  const int_t    ntot = Geometry::nBnode();
  char           buf[StrMax], err[StrMax], file[StrMax];
  ifstream       num;
  vector<real_t> work (npnp);
  real_t         *mass, *unity = &work[0];
  register int_t i, j, nset, nglobal;
  register int_t *gid, *q;

  // -- Read in NumberSystems.

  strcat   (strcpy (file, session), ".num");
  num.open (file);

  if (!num) {
    ROOTONLY {
      sprintf (buf, "enumerate -O1 %s > %s", session, file);
      if (system (buf)) {
        sprintf (err, "couldn't open session file %s, or %s", session, file);
        message (routine, err, ERROR);
      }
    }

    Femlib::synchronize();      // -- Ensure creation has completed.
    num.clear ();
    num.open  (file);

    if (!num) {
      sprintf (err, "couldn't find or create number file %s", file);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf >> buf >> buf;
  j = strlen (buf);
  for (i = 0; i < j; i++)
    if (!strchr (_fields, tolower (buf[i]))) {
      sprintf (err, "Fields nominated in %s.num (\"%s\") don't match \"%s\"",
	       session, buf, _fields);
      message (routine, err, ERROR);
    }

  num.getline(buf, StrMax).getline(buf, StrMax);

  num >> buf >> nset >> buf >> buf >> buf;

  _numsys.resize (nset);
  for (i = 0; i < nset; i++) {
    num >> buf;
    _numsys[i] = new NumberSys;
    _numsys[i] -> _fields = new char [strlen (buf) + 1];
    strcpy (_numsys[i] -> _fields, buf);
  }

  for (i = 0; i < nset; i++)
    for (j = 0; j < nset; j++) {
      if (i == j) continue;
      if (strpbrk (_numsys[i] -> _fields, _numsys[j] -> _fields)) {
	sprintf (err, "Field name duplication: %s <--> %s",
		 _numsys[i] -> _fields, _numsys[j] -> _fields);
	message (routine, err, ERROR);
      }
    }

  num >> buf >> buf;
  if (strcmp (buf, "NEL")) {
    sprintf (err, "expected \"NEL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != nel) {
      sprintf (err, "mismatch in number of elements: %1d vs. %1d", j, nel);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NP_MAX")) {
    sprintf (err, "expected \"NP_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != np) {
      sprintf (err, "mismatch in element order: %1d vs. %1d", j, np);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NEXT_MAX")) {
    sprintf (err, "expected \"NEXT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> j;

  num >> buf >> buf;
  if (strcmp (buf, "NINT_MAX")) {
    sprintf (err, "expected \"NINT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> j;

  num >> buf >> buf;
  if (strcmp (buf, "NTOTAL")) {
    sprintf (err, "expected \"NTOTAL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != NP) {
      sprintf (err, "mismatch in Field storage requirements: %1d vs %1d",j,NP);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NBOUNDARY")) {
    sprintf (err, "expected \"NBOUNDARY\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != ntot) {
      sprintf (err,"mismatch in number of boundary nodes: %1d vs. %1d",j,ntot);
      message (routine, err, ERROR);
    }
    _numsys[i] -> _btog  = new int_t [ntot];
    _numsys[i] -> _bmask = new int_t [ntot];
  }

  num >> buf >> buf;
  if (strcmp (buf, "NGLOBAL")) {
    sprintf (err, "expected \"NGLOBAL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _nglobal;

  num >> buf >> buf;
  if (strcmp (buf, "NSOLVE")) {
    sprintf (err, "expected \"NSOLVE\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _nsolve;

  num >> buf >> buf;
  if (strcmp (buf, "OPTIMIZATION")) {
    sprintf (err, "expected \"OPTIMIZATION\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _optlev;

  num >> buf >> buf;
  if (strcmp (buf, "BANDWIDTH")) {
    sprintf (err, "expected \"BANDWIDTH\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> _numsys[i] -> _nbandw;

  num.getline(buf, StrMax).getline(buf, StrMax).getline(buf, StrMax);

  // -- Consistency checks passed.  Now read in node numbers & BC masks.
  //    A bmask value of 1 implies an essential BC will be imposed at node.

  for (i = 0; i < ntot; i++) {
    num >> buf >> buf >> buf;
    for (j = 0; j < nset; j++) 
      num >> _numsys[j] -> _btog[i] >> _numsys[j] -> _bmask [i];
  }

  if (num.bad())
    message (routine, "failed reading to end of node-number file", ERROR);

  // -- Build emasks by inspecting bmasks for each element.
  //    An emask of 1 implies the element has at least one essential BC node.

  for (j = 0; j < nset; j++) {
    NumberSys* N = _numsys[j];
    N -> _emask = new int_t [nel];
    for (q = N -> _bmask, i = 0; i < nel; i++, q += next)
      N -> _emask[i] = Veclib::any (next, q, 1);
  }

  // -- Create inverse mass matrices, avoid division by zero on axis.

  for (j = 0; j < nset; j++) {
    
    // -- Create diagonal global mass matrix.

    nglobal = _numsys[j] -> _nglobal;
    mass    = _numsys[j] -> _imass = new real_t [nglobal];
    Veclib::zero (nglobal, mass, 1);
    for (gid = _numsys[j] -> _btog, i = 0; i < nel; i++, gid += next) {
      Veclib::fill (npnp, 1.0, unity, 1);
      elmt[i] -> bndryDsSum (gid, unity, mass);
    }

    // -- Invert.

    for (i = 0; i < nglobal; i++) mass[i] = 1.0 / mass[i];
  }
}


void BCmgr::buildsurf (FEML*             file,
		       vector<Element*>& Elmt)
// ---------------------------------------------------------------------------
// Private member function.
//
// Assemble list of element sides that have boundary conditions
// attached.  Element and side numbers are decremented by one, i.e are
// zero indexed in internal storage.
//
// As a part of internal checking, we want to ensure that the mesh for
// all "axis" group BCs has y=0.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "BCmgr::buildsurf";
  const int_t nsurf = file -> attribute ("SURFACES", "NUMBER");
  char        err[StrMax], tag[StrMax], group;
  int_t       i, t, elmt, side;
  BCtriple*   BCT;
 
  for (i = 0; i < nsurf; i++) {
    while ((file->stream().peek()) == '#') file->stream().ignore(StrMax, '\n');

    file -> stream() >> t >> elmt >> side >> tag;

    if (strcmp (tag, "<B>") == 0) {
      
      file -> stream() >> group;

      BCT = new BCtriple;

      BCT -> group = group;
      BCT -> elmt  = --elmt;
      BCT -> side  = --side;

      _elmtbc.insert (_elmtbc.end(), BCT);

      file -> stream() >> tag;
      if (strcmp (tag, "</B>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <B> with %s", t, tag);
	message (routine, err, ERROR);
      }
    } else
      file -> stream().ignore (StrMax, '\n');
  }

  if (Geometry::cylindrical()) {
    const int_t    np = Geometry::nP();
    vector<real_t> work (np);
    vector<BCtriple*>::iterator b;

    for (b = _elmtbc.begin(); b != _elmtbc.end(); b++) {
      BCT = *b;
      if (strstr (groupInfo (BCT -> group), "axis")) {
	Elmt[BCT -> elmt] -> sideGetY (BCT -> side, &work[0]);
	for (i = 0; i < np; i++)
	  if (::fabs (work[i]) > EPSDP) {
	    sprintf (err,
		     "elmt: %1d, side: %1d, offset: %1d, "
		     "y value (%g) too large on axis BC",
		     BCT -> elmt + 1,
		     BCT -> side + 1,
		     i, work[i]);
	    message (routine, err, ERROR);
	  }
      }
    }
  }
}


int_t BCmgr::nWall ()
// ---------------------------------------------------------------------------
// Count up the number of surfaces/element edges that have "wall" descriptor.
// ---------------------------------------------------------------------------
{
  vector<BCtriple*>::const_iterator b;
  int_t                             count = 0;
  BCtriple*                         BCT;

  for (b = _elmtbc.begin(); b != _elmtbc.end(); b++) {
    BCT = *b;
    if (strstr (groupInfo (BCT -> group), "wall")) count++;
  }

  return count;
}


///////////////////////////////////////////////////////////////////////////////
// Code following implements computed BC types.
///////////////////////////////////////////////////////////////////////////////


void BCmgr::buildComputedBCs (const Field* master)
// ---------------------------------------------------------------------------
// Build class-scope storage structures required for evaluation of various
// kinds of computed Boundary conditions for Navier--Stokes type problems.
//
// There is some wastage as memory is also allocated for locations that will
// not have computed BCs.
// ---------------------------------------------------------------------------
{
  int_t i, j;

  _nLine = master -> _nline;
  _nEdge = master -> _nbound;
  _nZ    = master -> _nz;
  _nP    = Geometry::nP();
  _nTime = Femlib::ivalue ("N_TIME");

  _work = new real_t [static_cast<size_t>
		      (5*sqr(_nP) + 7*_nP + Integration::OrderMax+1 + _nLine)];

  // -- The structure of each of the following arrays is
  //    _xx[time_level][z_plane] which evaluates to a real_t* that is
  //    a pointer to _nline real_t storage.  Thus _xx[0] is an
  //    equivalent pointer to Field->_line.

  _un    = new real_t** [static_cast<size_t>(_nTime)]; // -- u dot n.
  _hopbc = new real_t** [static_cast<size_t>(_nTime)]; // -- grad(p) dot n.

  for (i = 0; i < _nTime; i++) {
    _un   [i] = new real_t* [static_cast<size_t>(_nZ)];
    _hopbc[i] = new real_t* [static_cast<size_t>(_nZ)];
    
    _un   [i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    _hopbc[i][0] = new real_t [static_cast<size_t>(_nLine * _nZ)];
    
    Veclib::zero (_nLine * _nZ, _un   [i][0], 1);
    Veclib::zero (_nLine * _nZ, _hopbc[i][0], 1);

    for (j = 1; j < _nZ; j++) {
      _un   [i][j] = _un   [i][0] + j * _nLine;
      _hopbc[i][j] = _hopbc[i][0] + j * _nLine;
    }
  }
}


void BCmgr::maintainPhysical (const Field*             master,
			      const vector<AuxField*>& Us    ,
			      const int_t              nCom  )
// ---------------------------------------------------------------------------
// Just a stub for dog.  (Retained for compatibility with bcmgr.h.)
// ---------------------------------------------------------------------------
{ }


void BCmgr::maintainFourier (const int_t      step   ,
			     const Field*     master ,
			     const AuxField** Us     ,
			     const AuxField** Uf     ,
			     const bool       timedep)
// ---------------------------------------------------------------------------
// Update storage for evaluation of internally computed pressure boundary
// conditions.  Storage order for each edge represents a CCW traverse
// of element boundaries.
//
// If the velocity field varies in time on HOPB field boundaries
// (e.g. due to time-varying BCs) the local fluid acceleration will be
// estimated from input velocity fields by explicit extrapolation if
// timedep is true.  This correction cannot be carried out at the
// first timestep, since the required extrapolation cannot be done.
// If the acceleration is known, (for example, a known reference frame
// acceleration) it is probably better to leave timedep unset, and to
// use BCmgr::accelerate() to add in the accelerative term.  Note
// also that since grad P is dotted with n, the unit outward normal,
// at a later stage, timedep only needs to be set if there are
// wall-normal accelerative terms.  NB: The default value of timedep
// is true.
//
// Field* master gives a list of egdes with which to traverse storage
// areas (note this assumes equal-order interpolations).
//
// No smoothing is done to high-order spatial derivatives computed here.
// ---------------------------------------------------------------------------
{
  const real_t nu    = Femlib::value ("KINVIS");
  const real_t invDt = 1.0 / Femlib::value ("D_T");

  const AuxField* Ux = Us[0];
  const AuxField* Uy = Us[1];
  const AuxField* Uz = (Geometry::nPert() == 3) ? Us[2] : 0;

  const AuxField* Nx = Uf[0];
  const AuxField* Ny = Uf[1];

  const vector<Boundary*>& BC = master -> _bsys -> BCs (0);

  Boundary*       B;
  int_t           i, j, k, q, offset, skip, Je;

  // -- Roll grad P.n storage area up, load new level of nonlinear terms Uf.

  rollv (_hopbc, _nTime);

  for (i = 0; i < _nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();
    j      = i * _nP;

    for (k = 0; k < _nZ; k++) {
      Veclib::vmul  (_nP, Nx -> _plane[k] + offset, skip, B -> nx(), 1, 
                    _hopbc[0][k] + j, 1);
      Veclib::vvtvp (_nP, Ny -> _plane[k] + offset, skip, B -> ny(), 1, 
		     _hopbc[0][k] + j, 1, _hopbc[0][k] + j, 1);

      // -- For cylindrical coordinates, N_ are radius-premultiplied. Cancel.

      if (Geometry::cylindrical()) B -> divY (_hopbc[0][k] + j);
    }
  }

  // -- Add in -nu * curl curl u. There are 3 cases to deal with:
  //    perturbation is real_t, half-complex or full-complex.

  vector<real_t> work (5 * sqr(_nP) + 7 * _nP + Integration::OrderMax + 1);
  real_t         *UxRe, *UxIm, *UyRe, *UyIm, *UzRe, *UzIm, *tmp;
  real_t*        wrk   = &work[0];
  real_t*        xr    = wrk + 5*sqr(_nP) + 3*_nP;
  real_t*        xi    = xr + _nP;
  real_t*        yr    = xi + _nP;
  real_t*        yi    = yr + _nP;
  real_t*        alpha = yi + _nP;

  for (i = 0; i < _nEdge; i++) {
    B = BC[i];
    j = i * _nP;

    if (_nZ == 1) {

      UxRe = Ux -> _plane[0];
      UyRe = Uy -> _plane[0];

      if (Geometry::problem() == Geometry::O2_2D ||
	  Geometry::problem() == Geometry::SO2_2D ) { // -- Real perturbation.
	B->curlCurl(0,UxRe,0,UyRe,0,0,0,xr,0,yr,0,wrk);
      } else {			    // -- Half-complex perturbation.
	UzIm = Uz -> _plane[0];
	B->curlCurl(1,UxRe,0,UyRe,0,0,UzIm,xr,0,yr,0,wrk);
      }
      Veclib::svvttvp(_nP,-nu,xr,1,B->nx(),1,
		      _hopbc[0][0]+j,1,_hopbc[0][0]+j,1);
      Veclib::svvttvp(_nP,-nu,yr,1,B->ny(),1,
		      _hopbc[0][0]+j,1,_hopbc[0][0]+j,1); 

    } else {			    // -- Full complex perturbation.
      UxRe = Ux -> _plane[0];
      UxIm = Ux -> _plane[1];
      UyRe = Uy -> _plane[0];
      UyIm = Uy -> _plane[1];
      UzRe = Uz -> _plane[0];
      UzIm = Uz -> _plane[1];

      B->curlCurl(1,UxRe,UxIm,UyRe,UyIm,UzRe,UzIm,xr,xi,yr,yi,wrk);

      Veclib::svvttvp 
	(_nP, -nu, xr,1, B->nx(),1, _hopbc[0][0]+j,1, _hopbc[0][0]+j,1);
      Veclib::svvttvp 
	(_nP, -nu, xi,1, B->nx(),1, _hopbc[0][1]+j,1, _hopbc[0][1]+j,1);
      Veclib::svvttvp 
	(_nP, -nu, yr,1, B->ny(),1, _hopbc[0][0]+j,1, _hopbc[0][0]+j,1);
     Veclib::svvttvp 
	(_nP, -nu, yi,1, B->ny(),1, _hopbc[0][1]+j,1, _hopbc[0][1]+j,1);
    }
  }

  if (timedep) {

    // -- Estimate -d(u.n)/dt by backwards differentiation and add in.
    
    if (step > 1) {
      Je  = min (step - 1, _nTime);
      tmp = xr;
      Integration::StifflyStable (Je, alpha);
      
      for (i = 0; i < _nEdge; i++) {
	B      = BC[i];
	offset = B -> dOff ();
	skip   = B -> dSkip();
	j      = i * _nP;

	for (k = 0; k < _nZ; k++) {
	  Veclib::vmul  (_nP, Ux -> _plane[k] + offset, skip, B -> nx(), 1,
			 tmp, 1);
	  Veclib::vvtvp (_nP, Uy -> _plane[k] + offset, skip, B -> ny(), 1, 
			 tmp, 1, tmp, 1);
	  Blas::scal     (_nP, alpha[0], tmp, 1);
	  for (q = 0; q < Je; q++) 
	    Blas::axpy (_nP, alpha[q + 1], _un[q][k] + j, 1, tmp, 1);
	  Blas::axpy (_nP, -invDt, tmp, 1, _hopbc[0][k] + j, 1);

	}
      }
    }

    // -- Roll u.n storage area up, load new level.  It looks like we
    //    could avoid recomputing some things here by revising the
    //    loops above.

    rollv (_un, _nTime);
      
    for (i = 0; i < _nEdge; i++) {
      B      = BC[i];
      offset = B -> dOff ();
      skip   = B -> dSkip();
      j      = i * _nP;
    
      for (k = 0; k < _nZ; k++) {
	Veclib::vmul  (_nP, Ux -> _plane[k] + offset, skip, B -> nx(), 1,
		       _un[0][k] + j, 1);
	Veclib::vvtvp (_nP, Uy -> _plane[k] + offset, skip, B -> ny(), 1, 
		       _un[0][k] + j, 1, _un[0][k] + j, 1);
      }
    }
  }
}


void BCmgr::evaluateCNBCp (const int_t id   ,
			   const int_t plane,
			   const int_t step ,
			   real_t*     tgt  )
// ---------------------------------------------------------------------------
//  Refer KIO91.  Load pressure (p) BC value (tgt) with values
//  obtained from HOBC multi-level storage. Evaluation is confined to
//  a single element edge: parameter id tells us which this
//  corresponds to in our internal storage.
//
// The boundary condition for evaluation is
//
//   dP       /                                  du  \
//   -- = n . | N(u) - a + f + \nu*curlCurl(u) - --  |  =  n . grad P.
//   dn   ~   \ ~ ~    ~   ~                ~    dt  /     ~
//
// Grad P.n is estimated at the end of the current timestep using
// explicit extrapolation.  All the terms for this gradient were
// already stored, here we do the extrapolation only.
//
// This is basically the same as the old PBCmgr::evaluate method.
// ---------------------------------------------------------------------------
{
  if (step < 1) return;		// -- No evaluation during Field creation.

  Veclib::zero (_nP, tgt, 1);

  const int_t Je     = min (step, _nTime);
  const int_t offset = id * _nP;
  real_t*     beta   = _work;
  int_t       q;

  Integration::Extrapolation (Je, beta);
  
  for (q = 0; q < Je; q++)
    Blas::axpy (_nP, beta[q], _hopbc[q][plane] + offset, 1, tgt, 1);
}


void BCmgr::accelerate (const Vector& a,
			const Field*  u)
// ---------------------------------------------------------------------------
// Add in frame acceleration term on boundaries that correspond to
// essential velocity BCs (a = - du/dt).  Note that the acceleration
// time level should correspond to the time level in the most recently
// updated pressure gradient storage.  Work only takes place on zeroth
// Fourier mode.
//
// Yes, this is a HACK!
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = u -> _bsys -> BCs (0);
  int_t                    i;

  for (i = 0; i < u->_nbound; i++)
    BC[i] -> dotInForGroup ("velocity", a, _hopbc[0][0] + i*_nP);
}


void BCmgr::evaluateCEBCp (const Field* master,
			   const int_t  id    ,
			   const int_t  plane ,
			   const int_t  step  ,
			   real_t*      tgt   )
// ---------------------------------------------------------------------------
// Just a stub for dog.  (Prototype is needed by semtex/src/condition.cpp.)
// ---------------------------------------------------------------------------
{ }


void BCmgr::evaluateCNBCu (const Field* P    ,
			   const int_t  id   ,
			   const int_t  k    ,
			   const int_t  step ,
			   const char   cmpt ,
			   real_t*      tgt  )
// ---------------------------------------------------------------------------
// Just a stub for dog.  (Prototype is needed by semtex/src/condition.cpp.)
// ---------------------------------------------------------------------------
{ }
