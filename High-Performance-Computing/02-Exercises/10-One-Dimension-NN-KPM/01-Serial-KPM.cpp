#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cmath>
#include <complex>

#include <random>
#include <Eigen/Core>

#include <omp.h>


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
    Eigen::Array<double, -1, -1> KetXi;
    Eigen::Array<double, -1, -1> BraXi;

    // Random number generation
    std::mt19937 rnd;
    std::uniform_real_distribution<> dist;


    // Initialize RNG with a seed (8 numbers for maximum entropy according to the prof)
    kpm(unsigned ll, double tt, double ww) :
        EnergyScale((2 * tt + ww/2) * 1.01),
        W(ww / EnergyScale),
        U(ll),

        L(ll),
        t(tt / EnergyScale),

        KetXi(ll, 2),
        BraXi(1, ll)
    {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        rnd.seed(seed);
    }

    // The starting vector is just a random vector because the Chebyshev polynomial of order 0 is 1
    void initialize() {
        index = 0;
        for(unsigned i = 0; i < L; i++) {
            KetXi(i, index) = (dist(rnd) - 0.5) * 2 * sqrt(3);
            BraXi(0, i) = KetXi(i, index);
            U(i) = W * (dist(rnd) - 0.5);
        }
    }

    // Multiply by the the Hamiltonian for a 1D tight-binding model
    // We write the multiplication by hand because we know the structure
    // If we didn't, we could try to use sparse matrices, but it would make it less efficient
    void hamiltonian_kpm(int factor) {

        // factor == 1 is initialization, factor == 2 is a recursive KPM step
        if (factor == 1) index = 1;
        else             index++;
    
        unsigned i0 = index % 2;
        unsigned i1 = (index - 1) % 2;
    
        for(unsigned i = 0; i < L; i++) {
            KetXi(i, i0) = -factor * t * (KetXi((i + 1) % L, i1) + KetXi((i - 1 + L) % L, i1)) + factor * U(i) * KetXi(i, i1) - KetXi(i, i0);
        }
    }
};


// ./a.out 1000 1 0 2 1000 2000 2 1 0
int main(int argc, char* argv[]) {
    // Parse inputs
    int ll = std::stoi(argv[1]);
    double tt = std::stof(argv[2]);
    double ww = std::stof(argv[3]);
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
        mu.segment(0, 2) += ((kpmobj.BraXi.matrix() * kpmobj.KetXi.matrix()).array().transpose() - mu.segment(0, 2)) / double(av + 1);

        for(unsigned m = 2; m < NMoments; m += 2) {
            for(unsigned j = 0; j < 2; j++) 
                kpmobj.hamiltonian_kpm(2);            
            
            mu.segment(m, 2) += ((kpmobj.BraXi.matrix() * kpmobj.KetXi.matrix()).array().transpose() - mu.segment(m, 2)) / double(av + 1.);
        }
    }

    // Normalize
    mu = mu/double(kpmobj.L);
    
    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(ll) + "-" + std::to_string(tt) + "-" + std::to_string(ww) + "-" + std::to_string(NAverages) + "-" + std::to_string(NMoments)  + ".txt");
        for (unsigned i = 0; i < NMoments; i++) {
            outfile << mu[i] << " ";
            outfile << std::endl;
        }
        outfile.close(); 
    }
}