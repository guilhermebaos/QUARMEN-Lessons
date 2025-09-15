#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cmath>
#include <complex>

#include <omp.h>
#include <Eigen/Core>


// Eigen definitions
using Eigen::MatrixXcd;
using Eigen::VectorXcd;


int main(int argc, char* argv[]) {
    // Parse inputs
    int tnum = std::stoi(argv[1]);
    int timesteps = std::stoi(argv[2]);
    double dt = std::stof(argv[3]);
    int N = std::stoi(argv[4]);
    int n0 = std::stoi(argv[5]);
    int sig = std::stoi(argv[6]);
    double k0 = std::stof(argv[7]);
    int save = std::stoi(argv[8]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Define imaginary constant
    std::complex<double> i(0, 1);

    // Create an array with the wave function and its discrete Fourier transform
    MatrixXcd psi(timesteps, N);
    MatrixXcd psiFourier(timesteps, N);

    psi.fill(0.);
    psiFourier.fill(0.);
    
    // Set up the initial state
    for (int n = 0; n < N; n++) {
        psi(0, n) = std::exp(-1.0 * (n - n0) * (n - n0) / (2.0 * sig * sig) - i * k0 * (1.0 * n));
    }

    // TEST our transform
    for (int t = 0; t < timesteps-1; t++) {

        // Change into Bloch eigenbasis
        for (int ell = 0; ell < N; ell++) {
            for (int n = 0; n < N; n++) {
                psiFourier(t, ell) += std::exp(-i * (2. * M_PI * ell * n / N)) * psi(0, n);
            }
    
            psiFourier(t, ell) /= sqrt(N);
        }


        // Change into real space basis
        #pragma omp parallel firstprivate(psiFourier)
        {   
            // Temporary function
            MatrixXcd psiTemp(timesteps, N);
            psiTemp.fill(0.);

            #pragma omp for collapse(1) schedule(dynamic, 10)
            for (int n = 0; n < N; n++) {
                for (int ell = 0; ell < N; ell++) {
                    psiTemp(t+1, n) += std::exp(i * (2. * M_PI * ell * n / N)) * psiFourier(t, ell);
                }

                psiTemp(t+1, n) /= sqrt(N);
            }

            #pragma omp critical
            {   
                // Add to the final variable our intermediate result
                psi += psiTemp;
            }
        }

    }

    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile1(name + "-real-" + std::to_string(N) + "-" + std::to_string(n0) + "-" + std::to_string(sig) + "-" + std::to_string(k0) + ".txt");
        for (int i = 0; i < timesteps; i++) {
            for (int j = 0; j < N; j++) {
                outfile1 << psi(i, j) << " ";
            }
            outfile1 << std::endl;
        }

        outfile1.close(); 
        
        // Write the output to a file
        std::ofstream outfile2(name + "-fourier-" + std::to_string(N) + "-" + std::to_string(n0) + "-" + std::to_string(sig) + "-" + std::to_string(k0) + ".txt");
        for (int i = 0; i < timesteps; i++) {
            for (int j = 0; j < N; j++) {
                outfile2 << psiFourier(i, j) << " ";
            }
            outfile2 << std::endl;
        }

        outfile2.close();
    }
}

