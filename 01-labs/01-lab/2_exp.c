// generalization of the tower sampling to a continuous distribution:
// - given a probability density p(x), the equivalent of the sum s[i] 
//   is s(x)=\int_0^x p(y)dy
// - pick a random number r in (0,1); the equivalent of sampling the i-th
//   event is sampling x=s^(-1)(r)
// - for p(x)=exp(-x) in (0,infinite), s(x)=1-exp(-x) and s^(-1)(r)=-log(1-r)
  
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,n;
  double x,r;

  srandom((unsigned)time(NULL));

  n=2000;                  // number of samples
  for (i=0; i<n ; i++) {
    r=rnd();               // pick a random number r in (0,1)
    x = -log(1.0-r);       // random numbers distributed as exp(-x)
    printf("%10.5f\n",x);
  }

}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
