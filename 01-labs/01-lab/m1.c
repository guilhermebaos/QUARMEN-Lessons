// back to sampling
//
// - sample p(x) with the metropolis algorithm:
//    0 start from a chosen point x
//    1 propose a move to xp=x+r, with r uniformly sampled in [-delta/2,delta/2]
//      (delta controls the size of the move)
//    2 accept the move with probability p(xp)/p(x)
//    3 if the move is accepted update x<--xp
//    4 go back to 1 
//   after a number of "equilibation steps" (to be identified), the random 
//   walk samples p(x).
//
// - comparison with direct sampling: 
//   drawbacks: a) identify and exclude from sampling the equilibration steps
//              b) samples are in general correlated, which adversely affects
//                 the efficiency.
//   advantage: can sample any probability distribution
//
// - here p(x)=exp(-2*beta x^2)/sqrt(pi*0.5/beta)
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,n;
  double x,xp,acc,p,delta,beta,lnp,lnpp;

  n=20000;    // number of monte carlo steps
  beta=0.3;   // parameter related to the variance (sigma2=0.25/beta)
  delta=5.0;  // size of the move
  x=10.0;     // initial position
  lnp=-2*beta*pow(x,2); // value of (the logarithm of) p(x)

  srandom((unsigned)time(NULL));

  for (i=0; i<n ; i++) {
    xp=x+delta*(rnd()-0.5);   // trial move
    lnpp=-2*beta*pow(xp,2.0); // new value of log(p)
    p=exp(lnpp-lnp);          // probability of accepting the move
    if(p>rnd()){              // acceptance test
      acc=1;                  // move is accepted
      x=xp;                   // update x
      lnp=lnpp;               // update log(p)
    }
    else{
      acc=0;                  // move is rejected; do not update x 
    }
    printf("%f %f\n",x,acc);  // print x (after equilibration sample p(x))
  }
}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
