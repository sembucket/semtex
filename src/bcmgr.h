#ifndef BCMGR_H
#define BCMGR_H

typedef struct bctriple { char group; int_t elmt; int_t side; } BCtriple;

class BCmgr
// ===========================================================================
// This is a factory / retrieval service for classes derived from
// Condition, and maintains GROUP descriptors.  In addition, it reads
// and returns NumberSys objects from session.num.
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
