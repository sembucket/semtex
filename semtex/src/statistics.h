#ifndef STATISTICS_H
#define STATISTICS_H

class Statistics
// ===========================================================================
// Routines for statistical analysis of AuxFields.
// ===========================================================================
{
friend ifstream& operator >> (ifstream&, Statistics&);
friend ofstream& operator << (ofstream&, Statistics&);
public:
  Statistics (Domain*, vector<AuxField*>&);

  void update (AuxField**);
  void dump   ();

protected:
  const char*       name;
  Domain*           base;
  vector<AuxField*> src ;
  vector<AuxField*> avg ;
  int_t             navg;
};

#endif
