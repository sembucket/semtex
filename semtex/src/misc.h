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

  char*  sess;
  char*  sesd;
  int_t  nr  ;
  int_t  ns  ;
  int_t  nz  ;
  int_t  nel ;
  int_t  step;
  real_t time;
  real_t dt  ;
  real_t visc;
  real_t beta;
  char*  flds;
  char*  frmt;
};
istream& operator >> (istream&, Header&);
ostream& operator << (ostream&, Header&);

// -- Routines from misc.C:

ostream& printVector (ostream&, const char*, const int_t, ... );
char*    upperCase   (char *);
void     writeField  (ostream&, const char*, const int_t, const real_t,
		      vector<AuxField*>&);

#endif
