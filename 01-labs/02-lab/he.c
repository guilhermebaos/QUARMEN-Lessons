// helium atom with the metropolis method
// (atomic units: energies in hartree, distances in bohr)
//
//    h = -1/2 (\nabla_1^2+\nabla_2^2 - 2/|r_1| - 2/|r_2| + 1/|r_1-r_2|
//
//    psi_0 = exp(-\beta |r_1|) * exp(-\beta |r_2|)
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main () {
  int i,n,seed;
  double p,r,delta,beta,lnpsi,lnpsip,ekin,epot,etot;
  double x1,y1,z1,r1,xp1,yp1,zp1,rp1;
  double x2,y2,z2,r2,xp2,yp2,zp2,rp2;

  srand48((unsigned)time(NULL));        // initialization of rng

  n=10000;                              // number of monte carlo steps
  beta=2.1;                             // variational parameter
  delta=1.0;                            // size of the move
  x1=10.0;                              // initial position x1
  y1=0.1;                               // initial position y1
  z1=0.01;                              // initial position z1
  x2=0.1;                               // initial position x2
  y2=0.2;                               // initial position y2
  z2=0.3;                               // initial position z2

  r1=sqrt(x1*x1+y1*y1+z1*z1);           // distance e1-n
  r2=sqrt(x2*x2+y2*y2+z2*z2);           // distance e2-n
  lnpsi=-beta*(r1+r2);                  // log psi at R=(x1,y1,z1,x2,y2,z2)
  ekin=-0.5*(2*pow(beta,2)-2*beta*(1.0/r1+1.0/r2)); // kinetic part of elocal at R
  r=sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2)); // distance e-e
  epot=-2.0/r1-2.0/r2+1.0/r;            // potential energy at R
  etot=ekin+epot;                       // local energy at R

  for (i=0; i<n ; i++) {
    xp1=x1+delta*(drand48()-0.5);       //
    yp1=y1+delta*(drand48()-0.5);       //
    zp1=z1+delta*(drand48()-0.5);       // propose a move to (xp1,yp1,zp1,
    xp2=x2+delta*(drand48()-0.5);       //                    xp2,yp2,zp2)
    yp2=y2+delta*(drand48()-0.5);       //
    zp2=z2+delta*(drand48()-0.5);       //
    rp1=sqrt(xp1*xp1+yp1*yp1+zp1*zp1);  // distance e-n at (xp1,yp1,zp1)
    rp2=sqrt(xp2*xp2+yp2*yp2+zp2*zp2);  // distance e-n at (xp2,yp2,zp2)
    lnpsip=-beta*(rp1+rp2);             // log psi at Rp
    p=fmin(1.0,exp(2*(lnpsip-lnpsi)));  // probability to accept the move
    if(p>drand48()){
      x1=xp1;                             //  *****************************
      y1=yp1;                             //  * update electron positions,*
      z1=zp1;                             //  * distances e-n,            *
      r1=rp1;                             //  * and log psi               *
      x2=xp2;                             //  *                           *
      y2=yp2;                             //  *                           *
      z2=zp2;                             //  *                           *
      r2=rp2;                             //  *                           *
      lnpsi=lnpsip;                       //  *****************************
      ekin=-0.5*(2*pow(beta,2)-2*beta*(1.0/r1+1.0/r2)); //  +++++++++++++++
      r=sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2)); // +                +
      epot=-2.0/r1-2.0/r2+1.0/r;          //  + compute other properties  +
      etot=ekin+epot;                     //  +++++++++++++++++++++++++++++
    }
    printf("%f %f %f %f %f %f %f %f\n",x1,x2,ekin,epot,etot,p,r1,r2);
  }
}
