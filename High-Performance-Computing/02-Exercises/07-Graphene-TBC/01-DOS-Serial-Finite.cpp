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


// Number of atoms in the basis
#define basis 2

// Energy scale is t
#define t 1

// Number of averages over disorder
#define R 1


// We will choose below n0 to be x-axis and n1 to be y-axis and row-major ordering
// This means that we have
// ^ n1
// |   
// |   3  7  11 15
// |   2  6  10 14
// |   1  5  9  13
// |   0  4  8  12
// O-------------------> n0



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


unsigned get_index(unsigned n0, unsigned n1, unsigned alpha, unsigned N) {
    // Row-Major Ordering
    return alpha + basis * (n1 + N * n0);

    // Column-Major Ordering
    // return n0 + N * (n1 + N * alpha);
}

// std::array<int, 3> get_coords(int index, unsigned N) {
//     // Row-Major Ordering
//     int alpha = index % basis;
//     int n1 = ((index - alpha) / basis) % N;
//     int n0 = (index - alpha - basis * n1) / (basis * N);
//     return {n0, n1, alpha};
// }


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
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};

    std::mt19937 rnd;
    rnd.seed(seed);
    
    // Random number generation (dist is by default from -1 to +1 which is what we want)
    std::uniform_real_distribution<> dist;
    
    // Build the energy vector
    VectorXd ee(points);
    double energyScale = (3 * t + w/2) * 1.2;
    for (unsigned id = 0; id < points; id++) {
        ee[id] = -energyScale + (2. * id / double(points-1)) * energyScale;
    }

    // Build the Hamiltonian
    unsigned sizeH = basis * L * L;

    MatrixXcd realH(sizeH, sizeH);
    realH.fill(0.);

    for (unsigned n0 = 0; n0 < L; n0++) {
        for (unsigned n1 = 0; n1 < L; n1++) {
            // Local energy
            unsigned here = get_index(n0, n1, 1, L);
            unsigned sameCell = get_index(n0, n1, 0, L);
            
            // Same cell hopping
            realH(here, sameCell) = -t;
            realH(sameCell, here) = -t;
            
            // Make sure we are not at the boundary
            if (n1 != L-1) {
                unsigned up = get_index(n0, n1 + 1, 0, L);
                realH(here, up) = -t;
                realH(up, here) = -t;
            }

            // Make sure we are not at the boundary
            if (n0 != L-1) {
                unsigned right = get_index(n0 + 1, n1, 0, L);
                realH(here, right) = -t;
                realH(right, here) = -t;
            }
        }
    }

    // Averaging over disorder
    VectorXd dos(points);
    dos.fill(0.);
    for (unsigned r = 0; r < R; r++) {
        for (unsigned n = 0; n < sizeH; n++) {
            realH(n, n) = w * (dist(rnd) - 0.5);
        }

        // K-Dependent part
        for (unsigned ellx = 0; ellx < Nsc; ellx++) {
            // Get k values
            double kx = 2 * M_PI * ellx / Nsc;
            for (unsigned elly = 0; elly < Nsc; elly++) {
                double ky = 2 * M_PI * elly / Nsc;

                // Twisting terms
                for (unsigned n0 = 0; n0 < L; n0++) {
                    unsigned here = get_index(n0, L-1, 1, L);

                    unsigned up = get_index(n0, 0, 0, L);
                    realH(here, up) = std::complex<double> (-t, 0) * std::exp(-I * ky);
                    realH(up, here) = std::complex<double> (-t, 0) * std::exp(I * ky);
                }
                for (unsigned n1 = 0; n1 < L; n1++) {
                    unsigned here = get_index(L-1, n1, 1, L);

                    unsigned right = get_index(0, n1, 0, L);
                    realH(here, right) = std::complex<double> (-t, 0) * std::exp(-I * kx);
                    realH(right, here) = std::complex<double> (-t, 0) * std::exp(I * kx);
                }
                
                // Diagonalize the Hamiltonian
                SelfAdjointEigenSolver<MatrixXcd> solver(realH);
                VectorXd vals = solver.eigenvalues();
                // MatrixXd vecs = solver.eigenvectors();

                // Sum the Dirac deltas centered at each eigenenergy
                // TODO: Iterative average
                for (unsigned ie = 0; ie < sizeH; ie++) {
                    for (unsigned id = 0; id < points; id++) {
                        dos[id] += nascent_delta(ee[id] - vals[ie], gam, 1);
                    }
                }
            }
        }
    }
    
    dos /= double(L * L * Nsc * Nsc * R);


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

