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
  Statistics (Domain*);
  
  void initialise  ();

  void update      (AuxField**);
  void dump        ();

  void phaseUpdate (const int_t, AuxField**);

protected:
  const char*       _name;
  Domain*           _base;
  vector<AuxField*> _src ;	// -- Pointers to the base storage areas.
  vector<AuxField*> _avg ;	// -- Storage area for running averages.
  int_t             _navg;	// -- Number of averages so far.
  int_t             _iavg;	// -- Same as value of token "AVERAGE".
};

#endif
