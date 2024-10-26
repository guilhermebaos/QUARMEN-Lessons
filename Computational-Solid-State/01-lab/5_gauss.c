// sample a gaussian with mean 0 and variance 1 using the central limit theorem
// through a sum of random numbers uniformly distributed in (0,1).
// CLT says that the sum of N random variables (with arbitrary
// distribution but finite variance) has a gaussian distribution
//
// this uniform distribution has mean 0.5 and variance 1/12;
// the sum of k such numbers has mean 0.5*k and variance k/12.
// we need to shift and rescale the sum to get a gaussian with
// zero mean and unit variance

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,j,k,n;
  double x;

  srandom((unsigned)time(NULL));

  n=2000;                        // number of gaussian samples
  k=5;                           // number of terms in the sum
  for (i=0; i<n ; i++) {
    x=0.0;
    for (j=0; j<k; j++) {
      x+=rnd();                  // this is the sum
    }
    x=(x-0.5*k)*pow(12.0/k,0.5); // shift mean and rescale standard deviation
    printf("%10.5f\n",x);        // exp(-0.5*pow(x,2))/pow(2*pi,0.5) (CLT)
  }

}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
