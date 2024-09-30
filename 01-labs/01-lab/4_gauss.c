// - sample a gaussian with zero mean and variance one
//   using (one version of) the Box-Muller algorithm
// - shift and scale to get a gaussian with mean 1 and variance 4

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,n;
  double x,y,mean,variance;

  srandom((unsigned)time(NULL));

  mean=1.0;
  variance=4.0;
  n=2000;
  for (i=0; i<n ; i++) {
    x=cos(3.14159265358979*rnd())*pow(-2.0*log(rnd()),0.5); // Box-Muller
    y=mean+sqrt(variance)*x;                                // shift and scale
    printf("%10.5f %10.5f\n",x,y);
  }

}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
