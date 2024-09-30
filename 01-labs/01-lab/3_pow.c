// sample p(x)=(k+1)*x^k in (0,1):
// - as in 3_exp.c, pick a random number r in (0,1) and select s^(-1)(r)
// - here, s(x)=\int_0^x p(y)dy=x^(k+1) and s^(-1)(r)=r^[1/(1+k)]
//
// (this example shows that p(x) can diverge)

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,n;
  double x,r,k,kinv;

  srandom((unsigned)time(NULL));

  n=20000;                 // number of samples
  k=-0.7;                  // exponent in p(x)
  kinv=1.0/(1.0+k);
  for (i=0; i<n ; i++) {
    r=rnd();               // pick a random number in (0,1)
    x=pow(r,kinv);         // random numbers distributed as 0.3*x^(-0.7)
    printf("%10.5f\n",x);
  }

}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
