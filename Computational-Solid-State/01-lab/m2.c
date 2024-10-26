// use the metropolis algorithm to do integrals
//
// - harmonic oscillator:
//
//           h = -1/2 nabla^2 + 1/2 x^2, psi=exp(-beta x^2)
//
// - the variational energy on the state |psi> is
//   <psi|h|psi>/<psi|h|psi> = int psi^2(x)E_L(x)dx / int psi^2(x)dx
//                           = int p(x)E_L(x)dx
//   where p(x)=psi^2(x) / int psi^2(x)dx is a probability density, and
//   the "local energy" is
//         E_L(x)=-1/2[nabla^2 psi(x)]/psi(x) + V(x)
//               =-1/2[(2*beta*x)^2-2*beta] + 1/2 x^2 (kinetic + potential)
//
// - sample {x_i} from p(x) (see m1.c) and calculate E_L(x_i);
//
// - the variational energy is estimated as the average of {E_L(x_i)};
//   this average is done afterwards with statfor.x, excluding equilibration 
//   steps (by hand) and taking into account the correlation of the random 
//   walk (automatically)
//
// - wave functions, variational energy, local energy etc. will be discussed
//   in more detail. Here we want to see that integrals can be computed by
//   the metropolis algorithm, with due care for equilibration steps and 
//   correlation of the random walk.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double rnd();

int main () {
  int i,n;
  double x,xp,acc,p,delta,beta,lnpsi,lnpsip,ekin,epot,elocal;

  n=20000;    // number of mnte carlo steps
  beta=0.3;   // variational parameter
  delta=5.0;  // size of the move
  x=10.0;     // initial position

  srandom((unsigned)time(NULL));

  lnpsi=-beta*pow(x,2);               // value of (the logarithm of) psi(x)
  ekin=-0.5*(pow(2*beta*x,2)-2*beta); // kinetic energy term
  epot=+0.5*pow(x,2);                 // potential energy term
  elocal=ekin+epot;                   // local energy E_L(x)

  for (i=0; i<n ; i++) {
    xp=x+delta*(rnd()-0.5);           // trial move
    lnpsip=-beta*pow(xp,2.0);         // new value of log(psi)
    p=exp(2*(lnpsip-lnpsi));          // acceptance prob. [psi(xp)/psi(x)]^2
    if(p>rnd()){                      // acceptance test
      acc=1;                          // move accepted
      x=xp;                           // update x
      lnpsi=lnpsip;                   // update log(psi)
      ekin=-0.5*(pow(2*beta*x,2)-2*beta); // calculate kinetic energy term
      epot=+0.5*pow(x,2);             // calculate potential energy term
      elocal=ekin+epot;               // calculate local energy
    }
    else{
      acc=0;
    }
    printf("%f %f %f %f %f\n",x,ekin,epot,elocal,acc); // print local energy etc
  }
}

double rnd() {
  double r=((double)random()+1.)/((double)RAND_MAX+2.);
  return r;
}
