// hydrogen atom with the metropolis method
//
//       	h = -1/2 \nabla^2 - 1/r, psi=exp(-\beta r)
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main () {
  int i,n,seed;
  double r,rp,x,y,z,xp,yp,zp,p,delta,beta,lnpsi,lnpsip,ekin,epot,etot;

  srand48((unsigned)time(NULL));        // initialization of rng

  n=10000;                              // number of monte carlo steps
  beta=1.3;                             // variational parameter
  delta=0.4;                            // size of the move
  x=10.0;                               // initial position x
  y=0.1;                                // initial position y
  z=0.01;                               // initial position z

  r=sqrt(x*x+y*y+z*z);                  // distance e-n at (x,y,z)
  lnpsi=-beta*r;                        // log psi
  ekin=-0.5*(beta*beta-2*beta/r);       // kinetic part of local energy at r
  epot=-1.0/r;                          // potential energy at r
  etot=ekin+epot;                       // local energy at r

  for (i=0; i<n ; i++) {
    xp=x+delta*(drand48()-0.5);         //
    yp=y+delta*(drand48()-0.5);         // propose a move to (xp,yp,zp)
    zp=z+delta*(drand48()-0.5);         //
    rp=sqrt(xp*xp+yp*yp+zp*zp);         // distance e-n at (xp,yp,zp)
    lnpsip=-beta*rp;                    // log psi
    p=fmin(1.0,exp(2*(lnpsip-lnpsi)));  // probability to accept the move
    if(p>drand48()){
      x=xp;                             //  *****************************
      y=yp;                             //  * update electron position, *
      z=zp;                             //  * distance e-n,             *
      r=rp;                             //  * and log psi               *
      lnpsi=lnpsip;                     //  *****************************
      ekin=-0.5*(pow(beta,2)-2*beta/r); //  +++++++++++++++++++++++++++++
      epot=-1.0/r;                      //  + compute other properties  +
      etot=ekin+epot;                   //  +++++++++++++++++++++++++++++
    }
    printf("%f %f %f %f %f %f\n",x,ekin,epot,etot,p,r);
  }
}
