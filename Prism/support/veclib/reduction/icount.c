/*
 * Number of logical true values
 */

int icount(int n, int *x, int incx)
{
  int count = 0;
  while (n--) {
    if (*x) count++;
    x += incx;
  }
  return count;
}
