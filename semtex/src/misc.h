#ifndef MISC_H
#define MISC_H


class Header
// ===========================================================================
// Nekton/Prism/Semtex-compatible header struct + I/O routines.
// ===========================================================================
{
public:
  Header();
 ~Header() { delete [] sess; delete [] sesd; delete [] flds; delete [] frmt; }

  bool swab() const;

  char*   sess;
  char*   sesd;
  integer nr  ;
  integer ns  ;
  integer nz  ;
  integer nel ;
  integer step;
  real    time;
  real    dt  ;
  real    visc;
  real    beta;
  char*   flds;
  char*   frmt;
};
istream& operator >> (istream&, Header&);
ostream& operator << (ostream&, Header&);

// -- Routines from misc.C:

ostream& printVector (ostream&, const char*, const integer, ... );
char*    upperCase   (char *);
void     writeField  (ofstream&, const char*, const integer, const real,
		      vector<AuxField*>&);

#endif
