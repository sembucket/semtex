//
// krylov.C
//

#include <Sem.h>
#include <krylov.h>
#include <time.h>

Krylov::Krylov(FEML*    F     ,
	       int      kdim  ,
	       int      nvec  ,
	       int      length)
// ---------------------------------------------------------------------------
// construct and allocate memory for 0..kdim vectors of length "_length"
//
//
{
  real*  alloc;
  int    i;


  strcpy ((_session = new char [strlen(F->root())+1]), F -> root());
 

  // set _kdim and _ndim & _length -> can do this better ?
  _kdim = kdim;
  _nvec = nvec;
  _length = length;

  _Klen = _length*(_kdim+1);

  // set number of columns
  _column.setSize(_kdim+1);

  // allocate column memory and assign to column's
  alloc = new real[(size_t) _Klen ];
  assert(alloc!=0);     // ensure memory allocated
  for( i = 0; i <= _kdim; i++)  _column[i] = alloc + i * _length;
  _K = alloc;           // redundant, but clearer.

  // catchall randomise or zero
  Veclib::vrandom (_Klen, _K, 1);
  Veclib::sadd    (_Klen, -0.5, _K, 1, _K, 1);
  //Veclib::zero ( _Klen, _K, 1);

}


Krylov::~Krylov()
{
}


int Krylov::restart()
// --------------------------------------------------------------------------- 
// if a restart file "session".kry can be found use it for input 
//
// returns result of search for restart file.
// ---------------------------------------------------------------------------
{
  char          restartfile[StrMax], s[StrMax];
  const char     routine[] = "Krylov::restart";  
  int          kdim, nvec, length, Klen;
  ifstream strm (strcat (strcpy (restartfile, _session), ".kry"));

  if (strm) {
    cout << "-- Initial Krylov Matrix: read from file " << restartfile << endl;

    // session name
    strm.getline(s, StrMax);

    // creation date
    strm.getline(s, StrMax);

    // K dimensions
    strm.getline(s, StrMax);
    istrstream(s, strlen(s)) >> kdim; 
    if( kdim != _kdim) message (routine, "K dimension mismatch", ERROR);

    // N vec dimension
    strm.getline(s, StrMax);
    istrstream(s, strlen(s)) >> nvec; 
    if( nvec != _nvec) {
      cout << " changing N vec from " << _nvec << " to " << nvec << endl;
      _nvec = nvec;
    }

    // Vector length
    strm.getline(s, StrMax);
    istrstream(s, strlen(s)) >> length; 
    if( length != _length) 
      message (routine, "(No. fields * Planesize) mismatch", ERROR);

    // Matrix length
    strm.getline(s, StrMax);
    istrstream(s, strlen(s)) >> Klen; 
    if( Klen != _Klen) 
      message (routine, "Krylov matrix length mismatch", ERROR);

    // format
    strm.getline(s, StrMax);
    if (!strstr (s, "binary"))
      message (routine, "input field file not in binary format", ERROR);


    // K matrix data
    strm.read((char *) _K, (_length*(_kdim+1)*sizeof(real)));

    strm.close();

    return (EXIT_SUCCESS);
  }

  // if we get here, then file doesn't exist.
  return (EXIT_FAILURE);

}

int Krylov::dump()
// ---------------------------------------------------------------------------
// Write a krylov dump file.
// ---------------------------------------------------------------------------
{
  char      outputfile[StrMax];
  const char *kry_format[] = {
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25d "    "K dimension\n",
    "%-25d "    "N dimension\n",
    "%-25d "    "Column Length\n",
    "%-25d "    "Matrix Length\n",    
    "%-25s "    "Format\n"
  };

  ofstream strm(strcat( strcpy(outputfile, _session), ".kry"));

  if (strm){

    cout << "writing Krylov Matrix file" << endl;

    char      s1[StrMax], s2[StrMax];
    time_t    tp (time (0));

    sprintf( s1, kry_format[0], _session);
    strm << s1;

    strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
    sprintf  (s1, kry_format[1], s2);
    strm << s1;

    sprintf( s1, kry_format[2], _kdim);
    strm << s1;

    sprintf( s1, kry_format[3], _nvec);
    strm << s1;

    sprintf( s1, kry_format[4], _length);
    strm << s1;

    sprintf( s1, kry_format[5], _Klen);
    strm << s1;

    sprintf (s2, "binary ");
    Veclib::describeFormat (s2 + strlen (s2));
    sprintf (s1, kry_format[6], s2);
    strm << s1;

    // binary write of memory areas.
    strm.write((char *) _K, _Klen*sizeof(real));
    
    strm.close();

    return (EXIT_SUCCESS);
  }

  message("Krylov::dump", "Unable to open output file", ERROR);

  return (EXIT_FAILURE);
}

int Krylov::dim(){
  return _kdim;
}


real Krylov::norm(int col)
{
  return Blas::nrm2(_length, _column[col], 1);
}

void Krylov::normalise(int col)
{
  real norm;

  norm = Blas::nrm2 (_length, _column[col], 1);
  Blas::scal (_length, 1.0/norm, _column[col], 1);

}

void Krylov::scale(real norm)
{
  Blas::scal (_Klen, 1.0/norm, _K, 1);
}


void Krylov::roll()
{
  register integer i;

  for (i = 1; i <=_kdim; i++) 
    Veclib::copy (_length, _column[i], 1, _column[i-1], 1);
}


void Krylov::setDomain(Domain* D,
		       int     col)
{
  int j;

  for (j = 0; j < D->nField(); j++)
	D -> u[j] -> setPlane(0, _column[col]+j*Geometry::planeSize() );

  // convert to fourier space -> not required
  // D -> transform(FORWARD);

}

void Krylov::getDomain(Domain* D,
		       int     col)
{
  int j;

  // convert to physical space -> not required
  // D -> transform(INVERSE);  

  for (j = 0; j < D->nField(); j++)
    D -> u[j] -> getPlane(0, _column[col]+j*Geometry::planeSize() );
}


void Krylov::report()
{
}




