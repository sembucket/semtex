#ifndef BCMGR_H
#define BCMGR_H

typedef struct bctriple { char group; int_t elmt; int_t side; } BCtriple;

class BCmgr
// ===========================================================================
// This is a factory / retrieval service for classes derived from
// Condition, and maintains GROUP descriptors.  In addition, it reads
// and returns NumberSys objects from session.num.
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
// ===========================================================================
{
public:
  BCmgr (FEML*, vector<Element*>&);

  const char*        field        () const { return _fields; }
  const char*        groupInfo    (const char) const;
  Condition*         getCondition (const char, const char, const int_t = 0);
  NumberSys*         getNumberSys (const char, const int_t = 0);
  vector<BCtriple*>& getBCedges   () { return _elmtbc; }
  int_t              nBCedges     () const { return _elmtbc.size(); }
  int_t              nWall        (); // Should be const: OSX compiler bug?

  class CondRecd {
  public: 
    char       grp  ;
    char       fld  ;
    Condition* bcn  ;
    char*      value;
  };
    
private:
  char*              _fields  ; // String containing field names.
  vector<char>       _group   ; // Single-character group tags.
  vector<char*>      _descript; // Group name strings.
  vector<CondRecd*>  _cond    ; // Conditions in storage.
  vector<BCtriple*>  _elmtbc  ; // Group tags for each element-side BC.
  vector<NumberSys*> _numsys  ; // Numbering schemes in storage.
  bool               _axis    ; // Session file declared and axis BC group.

  void buildnum  (const char*, vector<Element*>&);
  void buildsurf (FEML*, vector<Element*>&);
};

#endif
