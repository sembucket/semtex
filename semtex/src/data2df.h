#ifndef DATA2DF_H
#define DATA2DF_H

class Data2DF
// ============================================================================
// Canonical field class, each np X np element is defined on [-1,1] X [-1, 1].
// Data are arranged element-ordered in 2D planes to create a 3D scalar field.
// ============================================================================
{
friend istream& operator >> (istream&, Data2DF&);
friend ostream& operator << (ostream&, Data2DF&);

public:
  Data2DF  (const int_t nP, const int_t nZ, const int_t nEl,
	     const char Name='\0');
  ~Data2DF () { delete [] _data; delete [] _plane; }

  char getName () { return _name; }
  Data2DF& reverse    ();

  Data2DF& operator = (const Data2DF&);
  Data2DF& operator = (const real_t);

  Data2DF& DFT1D      (const int_t);
  Data2DF& DPT2D      (const int_t, const char);

  
  Data2DF& filter1D   (const real_t, const int_t);
  Data2DF& filter2D   (const real_t, const int_t);

  Data2DF& conjugate  (const bool);
  Data2DF& symmetrize (const bool);

private:
  const char  _name;
  const int_t _np, _nz, _nel, _np2;
  int_t       _nplane, _ntot;
  real_t*     _data;
  real_t**    _plane;
};

#endif
