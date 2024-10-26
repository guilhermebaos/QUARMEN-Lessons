// - the function "random" returns pseudo-random integers uniformly distributed 
//   in [0,RAND_MAX]. In the GNU C library RAND_MAX=2147483647.
// - the function "srandom" initializes the pseudo-random number sequence
//   using the integer returned by "time"
// - transform to a) integers uniformly distributed in [1,6]
//                b) floating numbers uniformly distributed in (0,1)


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main () {
  int mini = 6;
  int maxi = 12;

  int i, j, k, n;
  double x;

  srand(time(NULL));                                // initialize random sequence

  n = 20000;                                        // number of samples
  for (i = 0; i < n; i++) {
    j = rand();                                     // integers in [0, RAND_MAX]

    k = j % (maxi - mini) + mini;                   // integers in [mini, maxi]
    x = ((double) j + 1.)/((double) RAND_MAX + 2.); // floating in (0,1)

    printf("%d\n", k);                                // choose output (j, k or x)
  }
}
