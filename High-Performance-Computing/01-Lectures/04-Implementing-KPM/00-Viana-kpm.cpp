#define _Alignof alignof
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cmath>
#include <chrono>
#include <fstream>
#include <vector>
#include <sstream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include <random>
#include <numeric>
#include <unistd.h>

struct kpm {
  const unsigned L;
  const double W;
  const double EnergyScale; 
  const double t;
  unsigned index;
  Eigen::Array<double, -1, -1> v;
  Eigen::Array<double, -1, -1> phi0;
  Eigen::Array<double, -1, 1> U;
  std::mt19937 rnd;
  std::uniform_real_distribution<> dist;

  kpm (unsigned ll, double jj, double w) : L(ll),  W(w), EnergyScale( (2 * jj + W/2)*1.0001 ), t(jj/EnergyScale), v(L, 2), phi0(1, L), U(L) {
    
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    rnd.seed(seed);	
  };
  
  void initialize() {
    index = 0;
    for(unsigned i = 0; i < L; i++)
      {
	v(i,index) = (dist(rnd) - 0.5) * 2 * sqrt(3);
	phi0(0,i) = v(i,index);
	U(i) = W * (dist(rnd) - 0.5)/EnergyScale;
      }  
  };
  
  void hamiltonian()
  {
    index = 1;
    for(unsigned i = 0; i < L; i++)
      v(i, 1) = -t * (v( (i + 1)%L , 0)  + v( (i + L - 1)%L , 0)) + U(i) *v(i, 0) ;
  };

  void kpm_iteration()
  {
    index++;
    unsigned i0 = (index)%2;
    unsigned i1 = (index - 1)%2;
    for(unsigned i = 0; i < L; i++)
      v(i, i0 ) = -2 * t * (v( (i + 1)%L , i1  )  + v( (i + L - 1)%L , i1 ))  +2*U(i)*v(i,i1)  - v(i, i0 );
  };
};

int main ()
{
  unsigned NAverages = 1;
  unsigned NMoments = 1024;
 
  struct kpm v(16384*64, 1., 1.);
  Eigen::Array<double, -1, 1> mu (NMoments);
  
  mu.setZero();
  
  for(unsigned av = 0; av < NAverages; av++)
    {
      v.initialize();
      v.hamiltonian();
      mu.segment(0, 2) += ( (v.phi0.matrix() * v.v.matrix()).array().transpose() - mu.segment(0, 2))/double( av + 1);
      
      for(unsigned m = 2; m < NMoments ; m += 2)
	{
	  for(unsigned j = 0; j < 2; j++)
	    v.kpm_iteration();
	  
	  mu.segment(m, 2) += ( (v.phi0.matrix() * v.v.matrix()).array().transpose() - mu.segment(m, 2))/(double(av) + 1); 
	}
    }

  std::cout << mu/double(v.L) << std::endl;
  
}
