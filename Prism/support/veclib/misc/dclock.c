#if !defined(i860) && !defined(dclock)

/*
 * Double-precision timing routine
 */

#include <time.h>

double dclock(void)
{
  static double tps = 1.0 / CLOCKS_PER_SEC;
  return (double) clock() * tps;
}

#endif
