/*---------------------------------------------------------------------
Sorts an array a[1..n] into ascending numerical order by Shell's method
(diminishing increment sort). n is input; a is replaced on output by its
sorted rearrangement.
----------------------------------------------------------------------*/
#include "stdio.h"
#include "stdlib.h"

void shell(unsigned long n, double a[]) {
  unsigned long i, j, inc;
  double v;
  inc = 1;                       // Determine the starting increment.
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do {                         // Loop over the partial sorts.
    inc /= 3;
    for (i = inc + 1; i<=n; i++) {   // Outer loop of straight insertion.
      v=a[i];
      j=i;
      while (a[j-inc] > v) {   // Inner loop of straight insertion.
	a[j] = a[j-inc];
	j -= inc;
	if (j <= inc) break;
      }
      a[j]=v;
    }
  } while (inc > 1);
}
