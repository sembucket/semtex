/*****************************************************************************
 *                  HEADER FILE FOR 1-D FFT UTILITIES
 *****************************************************************************/

#define SQR(a)     ((a) * (a))
#define MIN(a, b)  ((a) < (b) ?     (a) : (b))
#define MAX(a, b)  ((a) > (b) ?     (a) : (b))
#define STR_MAX 256

typedef double   real;
typedef struct  {real Re, Im;} complex;
enum    err_lev {WARNING, ERROR, REMARK};

void      message     (const char*, const char*, int);
complex*  cvector     (long, long);
int*      ivector     (long, long);
real*     rvector     (long, long);
void      freeCvector (complex*, long);
void      freeIvector (int*, long);
void      freeRvector (real*, long);

int   ispow2    (int);
int   roundpow2 (int);
void  preFFT    (complex*, const int, const int);
void  cFFT      (complex*, const int, const complex*, const int, const int);
void  rcFFT     (complex*, const int, const complex*, const int, const int);
void  pcFFT     (complex*, const int, const int);
