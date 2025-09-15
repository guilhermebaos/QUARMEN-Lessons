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
#define t 1.

// Helper definition
#define basis 2

// Number of averages over disorder
#define R 32

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


// Order the atoms of the same cell next to each other
unsigned get_index(unsigned n, unsigned alpha) {
    return alpha + basis * n;
}


// ./a.out 64 100 0 0.05 0.5 1000 1 1 0
int main(int argc, char* argv[]) {
    // Parse inputs
    unsigned L = std::stoi(argv[1]);
    unsigned K = std::stoi(argv[2]);
    double w = std::stod(argv[3]);
    double gam = std::stod(argv[4]);
    double Delta = std::stod(argv[5]);
    unsigned points = std::stoi(argv[6]);
    unsigned tnum = std::stoi(argv[7]);
    unsigned save = std::stoi(argv[8]);
    unsigned timing = std::stoi(argv[9]);

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
    double energyScale = (4 * t + w/2) * 1.1;
    for (unsigned id = 0; id < points; id++) {
        ee[id] = -energyScale + (2. * id / double(points-1)) * energyScale;
    }

    // Build the Hamiltonian
    unsigned sizeH = basis * L;

    MatrixXcd realH(sizeH, sizeH);
    realH.fill(0.);

    // Hoppings inside the same super cell
    for (unsigned n = 0; n < L; n++) {
        unsigned hereA = get_index(n, 0);
        unsigned hereB = hereA + 1;

        // Same cell hopping
        realH(hereA, hereB) = -t/2.;
        realH(hereB, hereA) = -t/2.;

        // Make sure we are not at the boundary
        if (n != L-1) {
            unsigned rightA = get_index(n+1, 0);
            unsigned rightB = rightA + 1;

            // Hoppings t0
            realH(hereA, rightA) = -t;
            realH(rightA, hereA) = -t;
            
            realH(hereB, rightB) = -t;
            realH(rightB, hereB) = -t;

            // Hopping t1
            realH(hereB, rightA) = -t/2.;
            realH(rightA, hereB) = -t/2.;

            // Hopping t2
            realH(hereA, rightB) = -t/4.;
            realH(rightB, hereA) = -t/4.;
        }
    }

    // Averaging over disorder
    VectorXd dos(points);
    dos.fill(0.);
    
    #pragma omp parallel
    {
        // Thread-specific vector
        VectorXd dosThread(points);
        dosThread.fill(0.);

        #pragma omp for schedule(dynamic, 1)
            for (unsigned r = 0; r < R; r++) {
                for (unsigned n = 0; n < L; n++) {

                    // Same atom hopping
                    unsigned hereA = get_index(n, 0);
                    unsigned hereB = hereA + 1;
                    
                    realH(hereA, hereA) = - Delta / 2.;
                    realH(hereB, hereB) = + Delta / 2.;
                    
                    // Add Anderson disorder
                    realH(hereA, hereA) += w * (dist(rnd) - 0.5);
                    realH(hereB, hereB) += w * (dist(rnd) - 0.5);
                }
                
                // Twisting terms
                unsigned hereA = get_index(L-1, 0);
                unsigned hereB = hereA + 1;

                unsigned rightA = get_index(0, 0);
                unsigned rightB = rightA + 1;

                // K-Dependent part (sample K points k)
                for (unsigned ell = 0; ell < K; ell++) {
                    // Sample from -pi to +pi
                    double k = M_PI * 2 * (dist(rnd) - 0.5);

                    // Hoppings t0
                    realH(hereA, rightA) = std::complex<double> (-t, 0) * std::exp(-I * k);
                    realH(rightA, hereA) = std::complex<double> (-t, 0) * std::exp(I * k);
                    
                    realH(hereB, rightB) = std::complex<double> (-t, 0) * std::exp(-I * k);
                    realH(rightB, hereB) = std::complex<double> (-t, 0) * std::exp(I * k);

                    // Hopping t1
                    realH(hereB, rightA) = std::complex<double> (-t/2., 0) * std::exp(-I * k);
                    realH(rightA, hereB) = std::complex<double> (-t/2., 0) * std::exp(I * k);

                    // Hopping t2
                    realH(hereA, rightB) = std::complex<double> (-t/4., 0) * std::exp(-I * k);
                    realH(rightB, hereA) = std::complex<double> (-t/4., 0) * std::exp(I * k);
                    
                    // Diagonalize the Hamiltonian
                    SelfAdjointEigenSolver<MatrixXcd> solver(realH);
                    VectorXd vals = solver.eigenvalues();
                    // MatrixXd vecs = solver.eigenvectors();

                    // Sum the Dirac deltas centered at each eigenenergy
                    for (unsigned ie = 0; ie < sizeH; ie++) {
                        for (unsigned id = 0; id < points; id++) {
                            dosThread[id] += nascent_delta(ee[id] - vals[ie], gam, 1);
                        }
                    }
                }
            }

        #pragma omp critical
        {
            dos += dosThread;
        }
    }
    
    
    dos /= double(L * K * R);


    // Output time
    if (timing) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Running time for tnum = " << tnum << ": " << elapsed.count() << " seconds\n";
    }
    

    if (save) {
        // Construct the filename
        std::string name = "out2";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(L) + "-" + std::to_string(K) + "-" + std::to_string(w) + "-" + std::to_string(gam) + "-" + std::to_string(R) + ".txt");
        
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

