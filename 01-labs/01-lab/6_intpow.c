// this is an integral, not just a sampling.
//
// - we calculate int_0^1 f(x)dx = int_0^1 p(x)g(x)dx
//   where p(x) is a probability density and g(x)=f(x)/p(x) is the
//   remaining factor in the integrand.
// - the integral is estimated by a) sampling n values {x_i} from p(x), and
//                                b) averaging g(x_i) over the sampled points
// - this code outputs {g(x_i)}, the average is done with statfor.x afterwards.
//
// - here f(x)=x^k, p(x)=1 and g(x)=f(x), for the cases k=-0.4 and k=-0.8
//
// - both integrals are finite, but for k=-0.8 the variance diverges
//   because the integral of p(x)*g(x)^2=x^(-1.6) from 0 to 1 is infinity.
// - running the code a few times, you can note that for k=-0.8:
//     a) the estimate of the integral differs from the exact value 5 
//        by more than one standard deviation more often than expected;
//     b) the standard deviation itself is somewhat erratic,
//        and it does not scale as 1/sqrt(n).

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,n;
  double x,k1,k2;

  srandom((unsigned)time(NULL));

  n= 5000;                // number of samples
  k1=-0.4;                // exponent in g(x) for k=-0.4
  k2=-0.8;                // exponent in g(x) for k=-0.8
  for (i=0; i<n ; i++) {
    x=rnd();              // sample uniform distribution in (0,1)
    printf("%f %f\n",pow(x,k1),pow(x,k2)); // g(x)=x^k, k=-0.4,-0.8
  }
}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
