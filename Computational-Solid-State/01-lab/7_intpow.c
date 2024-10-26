// this is an example of importance sampling: modify p(x) to get lower variance
//
// - we again calculate int_0^1 f(x) dx with f(x)=x^k and k=-0.8
// - this time, we choose p(x)=0.3*x^(-0.7) (see 3_pow.c)
// - as a result, g(x)=f(x)/p(x)=x^(-0.1)/0.3 and the integral of
//   p(x)g(x)^2 = x^(-0.9)/0.3 is finite.
// - we get results statistically consistent with the exact value 5
//   and the expected behavior (~1/sqrt(n)) for the standard deviation

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,n;
  double x,r,k,kinv;

  srandom((unsigned)time(NULL));

  n=10000;          // number of samples
  k=-0.7;           // the exponent of x in p(x)
  kinv=1.0/(1.0+k);
  for (i=0; i<n ; i++) {
    r=rnd();        // pick a reandom number in (0,1)
    x=pow(r,kinv);  // sample the distribution pi(x)=0.3*x^(-0.7)
    printf("%f\n",pow(x,-0.1)/0.3); // g(x)=x^(-0.1)/0.3
  }
}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
