#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cmath>
#include <complex>

#include <random>
#include <Eigen/Core>

#include <omp.h>


class kpm {
public:

    // Parameters
    const unsigned L;
    const double t;

    // Computing mu
    unsigned index;
    Eigen::Array<double, -1, -1> v;
    Eigen::Array<double, -1, -1> phi0;

    // Disorder (W is the strength of the disorder)
    const double W;
    const double EnergyScale;
    Eigen::Array<double, -1, 1> U;

    // Random number generation
    std::mt19937 rnd;
    std::uniform_real_distribution<> dist;


    // Initialize RNG with a seed (8 numbers for maximum entropy according to the prof)
    kpm(unsigned ll, double tt, double w) : L(ll), W(w), EnergyScale((2 * tt + W/2) * 1.0001), t(tt / EnergyScale), v(L, 2), phi0(1, L), U(L) {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        rnd.seed(seed);
    }

    // The starting vector is just a random vector because the Chebyshev polynomial of order 0 is 1
    void initialize() {
        index = 0;
        for(unsigned i = 0; i < L; i++) {
            v(i, index) = (dist(rnd) - 0.5) * 2 * sqrt(3);
            phi0(0, index) = v(i, index);
            U(i) = W * (dist(rnd) - 0.5) / EnergyScale;
        }
    }

    // Multiply by the the Hamiltonian for a 1D tight-binding model
    // We write the multiplication by hand because we know the structure
    // If we didn't, we could try to use sparse matrices, but in this case is just adding complexity
    void hamiltonian() {
        index = 1;
        for(unsigned i = 0; i < L; i++) {
            v(i, 1) = -t * (v((i + 1) % L, 0) + v((i - 1 + L) % L, 0)) + U(i) * v(i, 0);
        }
    }


    // Do an iteration of the recursive kpm formula
    void kpm_iteration() {
        index++;
        unsigned i0 = index % 2;
        unsigned i1 = index % 2;

        for(unsigned i = 0; i < L; i++) {
            v(i, i0) = -2 * t * (v((i + 1) % L, i1) + v((i - 1 + L) % L, i1)) + 2 * U(i) * v(i, i1) - v(i, i0);
        }
    }
};

// There is something missing below, there might be something missing above!
int main () {
    // Parameters
    unsigned NAverages = 2;
    unsigned NMoments = 16;

    // We choose t < 0.5 because the spectrum is 2t
    // We want the spectrum to be between -1 and 1 because that is the domain of the Chebyshev polynomials 
    kpm v(16, 1., 1.);

    // Compute the moments 
    Eigen::Array<double, -1, 1> mu(NMoments, 1);
    mu.setZero();

    // Compute moments
    for (unsigned av = 0; av < NAverages; av++) {
        v.initialize();
        v.hamiltonian();

        // Do the two inner products (its memory efficient because we only need to access phi once)
        mu.segment(0, 2) = ((v.phi0.matrix() * v.v.matrix()).array().transpose() - mu.segment(0, 2)) / double(av + 1);

        for(unsigned m = 2; m < NMoments; m += 2) {
            for(unsigned j = 0; j < 2; j++) 
                v.kpm_iteration();            
            
            mu.segment(m, 2) += ((v.phi0.matrix() * v.v.matrix()).array().transpose() - mu.segment(m, 2)) / double(av + 1.);
        }
    }
}


// This program is very fast, hence we can go to very large system sizes
// In the limit of a very large system the majority of the time might be spent moving memory to the CPU cache