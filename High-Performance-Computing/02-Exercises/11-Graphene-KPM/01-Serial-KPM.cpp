#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cmath>
#include <complex>

#include <random>
#include <Eigen/Core>

#include <omp.h>


int get_index(int n0, int n1, int alpha, unsigned N) {
    // Row-Major Ordering
    return alpha + 2 * (n1 + N * n0);

    // Column-Major Ordering
    // return n0 + N * (n1 + N * alpha);
}

std::array<int, 3> get_coords(int index, unsigned N) {
    // Row-Major Ordering
    int alpha = index % 2;
    int n1 = ((index - alpha) / 2) % N;
    int n0 = (index - alpha - 2 * n1) / (2 * N);
    return {n0, n1, alpha};
}

// Remember that class members are initialized in order of their declaration in the class
// Not in the order stated by the constructor
class kpm {
public:

    // Disorder and Energy Scale (W is the strength of the disorder)
    const double EnergyScale;
    const double W;
    Eigen::Array<double, -1, 1> U;

    // Parameters
    const unsigned L;
    const double t;

    // Computing mu
    unsigned index;
    Eigen::Array<double, -1, -1> v;
    Eigen::Array<double, -1, -1> phi0;

    // Random number generation
    std::mt19937 rnd;
    std::uniform_real_distribution<> dist;


    // Initialize RNG with a seed (8 numbers for maximum entropy according to the prof)
    kpm(unsigned ll, double tt, double ww) :
        EnergyScale((4 * tt + ww/2) * 1.0001),
        W(ww / EnergyScale),
        U(2*ll*ll),

        L(ll),
        t(tt / EnergyScale),

        v(2*ll*ll, 2),
        phi0(1, 2*ll*ll)
    {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        rnd.seed(seed);
    }

    // The starting vector is just a random vector because the Chebyshev polynomial of order 0 is 1
    void initialize() {
        index = 0;
        for(unsigned i = 0; i < 2*L*L; i++) {
            v(i, index) = (dist(rnd) - 0.5) * 2 * sqrt(3);
            v(i, index+1) = 0;
            phi0(0, i) = v(i, index);
            U(i) = W * (dist(rnd) - 0.5);
        }
    }

    // Multiply by the the Hamiltonian for a 1D tight-binding model
    // We write the multiplication by hand because we know the structure
    // If we didn't, we could try to use sparse matrices, but it would make it less efficient
    // Factor is 1 after initialization and 2 when doing KPM iteration
    void hamiltonian_kpm(int factor) {
        if (factor == 1) index = 1;
        else             index++;
    
        unsigned i0 = index % 2;
        unsigned i1 = (index - 1) % 2;
    
        for (unsigned n0 = 0; n0 < L; n0++) {
          for (unsigned n1 = 0; n1 < L; n1++) {
            for (unsigned alpha = 0; alpha < 2; alpha++) {
              unsigned beta = (alpha + 1) % 2;
              int me  = get_index(n0, n1, alpha, L);
    
              // Disorder and Chebyshev factor
              v(me, i0) = factor * U(me) * v(me, i1) - v(me, i0);
    
              // Same‐cell hop:  A <--> B in (n0,n1)
              int same_cell = get_index(n0, n1, beta, L);
              v(me, i0) += -factor * t * v(same_cell, i1);
    
              // 3) Now split by which sublattice alpha is:
              if (alpha == 0) {
                // For A‐sites, the two “inter‐cell” neighbors are:
                //   B at (n0, n1−1)   and  B at (n0−1, n1)
                int down  = get_index(n0, (n1 - 1 + L) % L, beta, L);
                int left  = get_index((n0 - 1 + L) % L, n1, beta, L);
                v(me, i0) += -factor * t * v(down, i1);
                v(me, i0) += -factor * t * v(left, i1);
              }
              else {
                // For B‐sites, the two “inter‐cell” neighbors are:
                //   A at (n0, n1+1)   and  A at (n0+1, n1)
                int up    = get_index(n0, (n1 + 1) % L, beta, L);
                int right = get_index((n0 + 1) % L, n1, beta, L);
                v(me, i0) += -factor * t * v(up, i1);
                v(me, i0) += -factor * t * v(right, i1);
              }
            }
          }
        }
    }
};


int main(int argc, char* argv[]) {
    // Parse inputs
    int ll = std::stoi(argv[1]);
    double tt = std::stod(argv[2]);
    double ww = std::stod(argv[3]);
    unsigned NAverages = std::stoi(argv[4]);
    unsigned NMoments = std::stoi(argv[5]);
    unsigned save = std::stoi(argv[6]);

    // Make sure the number of moments is even
    NMoments += NMoments % 2;

    // Initialize the class
    kpm kpmobj(ll, tt, ww);

    // Compute the moments 
    Eigen::Array<double, -1, 1> mu(NMoments, 1);
    mu.setZero();

    // Compute moments
    for (unsigned av = 0; av < NAverages; av++) {
        kpmobj.initialize();
        kpmobj.hamiltonian_kpm(1);

        // Do the two inner products and compute iterative average (its memory efficient because we only need to access phi once)
        mu.segment(0, 2) += ((kpmobj.phi0.matrix() * kpmobj.v.matrix()).array().transpose() - mu.segment(0, 2)) / double(av + 1);

        for(unsigned m = 2; m < NMoments; m += 2) {
            for(unsigned j = 0; j < 2; j++) 
                kpmobj.hamiltonian_kpm(2);            
            
            mu.segment(m, 2) += ((kpmobj.phi0.matrix() * kpmobj.v.matrix()).array().transpose() - mu.segment(m, 2)) / double(av + 1.);
        }
    }

    // Normalize
    mu = mu/double(kpmobj.L * kpmobj.L);
    
    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(ll) + "-" + std::to_string(tt) + "-" + std::to_string(ww) + "-" + std::to_string(NAverages) + "-" + std::to_string(NMoments)  + ".txt");
        for (unsigned i = 0; i < NMoments; i++) {
            outfile << (mu[i]) << " ";
            outfile << std::endl;
        }
        outfile.close(); 
    }
}


// This program is very fast, hence we can go to very large system sizes
// In the limit of a very large system the majority of the time might be spent moving memory to the CPU cache