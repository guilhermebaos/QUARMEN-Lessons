// use tower sampling to sample 6 events with randomly assigned probability:
// - s[i] sum of the assigned probabilities for events {1,...,i}
// - f[i] fraction of samples corresponding to the i-th event
// - on output, the first column is the assigned probability of the i-th event,
//   and the second column is its estimate f[i] 
//
// no postprocessing required, you may just want to check that more 
// sampling, i.e. larger n, gives more precise estimates of the assigned 
// probabilities.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double rnd();

int main () {
  int i,j,n;
  double x,s[6],f[6];

  srandom((unsigned)time(NULL));

  n=2000;                                  // number of samples
  s[0]=0.0;
  for (i=1; i<7; i++) {s[i]=s[i-1]+rnd();} // construct s[i], i=1,...,6
  for (i=1; i<7; i++) {s[i]=s[i]/s[6];}    // normalize so that s[6]=1

  for (i=1; i<7; i++) {f[i]=0.0;}          // initialize f[i]
  for (j=0; j<n; j++) {                    // loop on samples
    x = rnd();                             // pick a random number x in (0,1)
    i = 1;
    while (s[i]<x) {i++;}                  // pick largest i such that s[i] < x
    f[i]++;                                // increment event counter for that i
  }

  for (i=1; i<7; i++) {printf("%f %f\n",s[i]-s[i-1],f[i]/n);}

}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
