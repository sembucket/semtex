#ifndef ARRAY_H
#define ARRAY_H
//////////////////////////////////////////////////////////////////////////////
// Array.h: templated 1D (vector) and 2D (matrix) array classes.
// 
// No subscript checking, and direct access to the underlying pointers
// can be obtained if desired using operator ().
//
// Matrix storage is row-major contiguous.
//////////////////////////////////////////////////////////////////////////////

// $Id$


template<class T>
class vector
// ===========================================================================
// 1D array class.
// ===========================================================================
{
public:
  // -- Creation/destruction.

  vector (const int n)          { num_elts = n; data = (n) ? new T[n] : 0; }
  vector ()                     { num_elts = 0; data = 0;                  } 
  vector (const vector<T>& src) { num_elts = src.num_elts;
				  data = new T[num_elts];
				  copy (src); } 
  ~vector ()                    { delete [] data; }

  // -- Assignment.

  vector<T>& operator = (const vector<T>& src) {
    if (data != src.data ) { setSize (src.num_elts); copy (src); }
    return *this;
  }

  vector<T>& operator = (const T& src) {
    register T* p = data + num_elts; while (p > data) *--p = src;
    return *this;
  }

  // -- Subscripting.

        T* operator () ()       { return data; }
  const T* operator () () const { return data; } 

        T& operator [] (const int i)       { return data[i]; }
  const T& operator [] (const int i) const { return data[i]; }
        T& operator () (const int i)       { return data[i]; }
  const T& operator () (const int i) const { return data[i]; } 

  // -- Size operators/information.

  int  getSize () const      { return num_elts; }
  void setSize (const int n) {
    if (n != num_elts) 
      {
	delete [] data;
	num_elts = n;
	data = (n) ? new T[n] : 0; 
      }
  }

private:
  int num_elts;
  T*  data;
  void copy (const vector<T>& src) {
    register T* p =     data + num_elts;
    register T* q = src.data + num_elts;
    while (p > data) *--p = *--q;
  }
};


template<class T> class matrix
// ==========================================================================
// 2D contiguous-storage row-major array.
// ==========================================================================
{
public:
  // -- Creation/destruction.

  matrix (const int n_rows, const int n_cols) {
    nr   = n_rows;
    nc   = n_cols;
    row  = new T* [nr];
    data = new T  [nr * nc];
    for (int i = 0; i < nr; i++) row[i] = data + i * nc;
  }
  matrix () { 
    nr   = nc = 0;
    row  = 0;
    data = 0;
  } 
  matrix (const matrix<T>& src) { 
    nr   = src.nr;
    nc   = src.nc;
    row  = new T* [nr];
    data = new T  [nr * nc];
    for (int i = 0; i < nr; i++) row[i] = data + i * nc;
    copy (src); 
  } 
  ~matrix () {
    delete [] data;
    delete [] row;
  }

  // -- Assignment.

  matrix<T>& operator = (const matrix<T>& src) {
    if (data != src.data ) {
      setSize (src.nr, src.nc);
      copy (src); }
    return *this;
  }

  matrix<T>& operator = (const T& src) {
    register T* p = data + nr * nc; while (p > data) *--p = src;
    return *this;
  }

  // -- Subscripting.

        T** operator () ()                          { return             row; }
  const T** operator () () const                    { return (const T**) row; }
  
        T* operator () (const int i)                    { return row[i]; }
  const T* operator () (const int i) const              { return row[i]; }
  
        T& operator () (const int i, const int j)       { return data[j+i*nc];}
  const T& operator () (const int i, const int j) const { return data[j+i*nc];}

  // -- Size operators/information.

  int  nRows   () const      { return nr;      }
  int  nCols   () const      { return nc;      }
  int  getSize () const      { return nr * nc; }
  void setSize (const int n_rows, const int n_cols) { 
    if (n_rows != nr && n_cols != nc) 
      {
	delete [] data;
	delete [] row;
	nr   = n_rows;
	nc   = n_cols;
	row  = (nr*nc) ? new T* [nr]    : 0;
	data = (nr*nc) ? new T  [nr*nc] : 0;
	for (int i = 0; i < nr; i++) row[i] = data + i * nc;
      }
  }

private:
  int  nr;
  int  nc;
  T**  row;
  T*   data;
  void copy (const matrix<T>& src) {
    register T* p =     data + nr * nc;
    register T* q = src.data + nr * nc;
    while (p > data) *--p = *--q;
  }
};

#endif
