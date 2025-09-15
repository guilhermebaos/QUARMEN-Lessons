#include <cstdlib>
#include <chrono>

#include <iostream>
#include <fstream>

#include <cmath>
#include <complex>

#include <random>

#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using Eigen::SelfAdjointEigenSolver;


// Energy scale is t
#define t 1




// Nascent Dirac delta
double nascent_delta(double x, double eps, int func) {
    switch(func) {
        case 1:
            return eps / (M_PI * (x * x + eps * eps));

        case 2:
            return exp(-x * x / (eps * eps)) / (eps * sqrt(M_PI));
        
        default:
            return -1;
    }
}


// ./a.out 512 8 0 0.05 1000 1 0 0
int main(int argc, char* argv[]) {
    // Parse inputs
    unsigned L = std::stoi(argv[1]);
    unsigned Nsc = std::stoi(argv[2]);
    double w = std::stod(argv[3]);
    double gam = std::stod(argv[4]);
    unsigned points = std::stoi(argv[5]);
    unsigned tnum = std::stoi(argv[6]);
    unsigned save = std::stoi(argv[7]);
    unsigned timing = std::stoi(argv[8]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Start counting time
    auto start = std::chrono::high_resolution_clock::now();
    
    // Define imaginary constant
    std::complex<double> I(0, 1);

    // Set the seed for RNG
    std::seed_seq seed{
        static_cast<uint32_t>(34070037),
        static_cast<uint32_t>(3141791178),
        static_cast<uint32_t>(97172939),
        static_cast<uint32_t>(92441404),
        static_cast<uint32_t>(4106087501),
        static_cast<uint32_t>(2981147808),
        static_cast<uint32_t>(264707443),
        static_cast<uint32_t>(2165443747)
    };
    std::mt19937 rnd;
    rnd.seed(seed);

    // Generate a single random seed
    // std::random_device r;
    // for (unsigned i = 0; i < 8; i++)
    //     std::cout << std::to_string(r()) << ", ";
    
    // Random number generation (dist is by default from -1 to +1 which is what we want)
    std::uniform_real_distribution<> dist;

    // Build the Hamiltonian
    unsigned sizeH = L;

    MatrixXcd realH(sizeH, sizeH);
    realH.fill(0.);

    for (unsigned n = 0; n < L-1; n++) {
        // Hopping to next atom in the lattice
        realH(n, n+1) = -t;
        realH(n+1, n) = -t;
        
        // Disorder term
        realH(n, n) = w * (dist(rnd) - 0.5);
    }
    realH(L-1, L-1) = w * (dist(rnd) - 0.5);


    // Build the energy vector
    VectorXd ee(points);
    double energyScale = (2 * t + w/2) * 1.1;
    for (unsigned id = 0; id < points; id++) {
        ee[id] = -energyScale + (2. * id / double(points-1)) * energyScale;
    }

    
    // Create the DOS vector
    VectorXd dos(points);
    dos.fill(0.);

    // K-Dependent part
    for (unsigned ell = 0; ell < Nsc; ell++) {
        double k = 2 * M_PI * ell / Nsc;

        // Twisting terms
        realH(L-1, 0) = std::complex<double> (-t, 0) * std::exp(-I * k);
        realH(0, L-1) = std::complex<double> (-t, 0) * std::exp(I * k);
        
        // Diagonalize the Hamiltonian
        SelfAdjointEigenSolver<MatrixXcd> solver(realH);
        VectorXd vals = solver.eigenvalues();
        // MatrixXd vecs = solver.eigenvectors();

        // Sum the Dirac deltas centered at each eigenenergy
        for (unsigned ie = 0; ie < sizeH; ie++) {
            for (unsigned id = 0; id < points; id++) {
                dos[id] += nascent_delta(ee[id] - vals[ie], gam, 1);
            }
        }
    }
    
    dos /= double(L * Nsc);


    // Output time
    if (timing) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Running time for tnum = " << tnum << ": " << elapsed.count() << " seconds\n";
    }
    

    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(L) + "-" + std::to_string(Nsc) + "-" + std::to_string(w) + "-" + std::to_string(gam) + ".txt");
        
        for (unsigned i = 0; i < points; i++) {
            outfile << ee[i] << " ";
        }
        outfile << std::endl;

        for (unsigned i = 0; i < points; i++) {
            outfile << dos[i] << " ";
        }
        outfile << std::endl;
        outfile.close(); 
    }
}

