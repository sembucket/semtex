#ifndef FEML_H
#define FEML_H
///////////////////////////////////////////////////////////////////////////////
// feml: Finite Element markup Language (FEML) header file.
//
// Copyright (c) 1997 <--> $Date$, Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

using namespace std;

#include <cfemdef.h>

#define KEYWORD_MAX 32


class FEML
// ===========================================================================
// Routines provide facilities to position input stream at location of
// a given keyword.
//
// This is the list of currently-defined (reserved) section keywords:
//   TOKENS
//   FIELDS
//   GROUPS
//   BCS
//   NODES
//   ELEMENTS
//   SURFACES
//   CURVES
//   USER
//   HISTORY
//   CUT
// Keywords are stored upper case, input is case-insensitive.
// The FEML class does not require that any of the above sections are actually
// used in an input file, it just treats them as reserved section tag-names.
// ===========================================================================
{

public:
  FEML  (const char*);
  ~FEML () { feml_file.close(); }

  bool        seek      (const char*);
  int_t       attribute (const char*, const char*);
  istream&    stream    ()       { return feml_file; }
  const char* root      () const { return feml_root; }  
  
private:
  char*     feml_root;		  // Name of FEML file, suffix removed.
  ifstream  feml_file;		  // Input stream.
  streampos keyPosn[KEYWORD_MAX]; // Locations corresponding to keywords.
  char*     keyWord[KEYWORD_MAX]; // Keywords used.

  bool      tokens ();
};

#endif
