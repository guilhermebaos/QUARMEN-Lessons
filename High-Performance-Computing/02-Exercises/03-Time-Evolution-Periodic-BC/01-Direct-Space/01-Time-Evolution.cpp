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
    double Fimag = std::stof(argv[3]);
    int N = std::stoi(argv[4]);
    int n0 = std::stoi(argv[5]);
    int sig = std::stoi(argv[6]);
    double k0 = std::stof(argv[7]);
    int save = std::stoi(argv[8]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Define F number
    std::complex<double> F(0, Fimag);

    // Define imaginary constant
    std::complex<double> i(0, 1);

    // Create an array with the wave function
    MatrixXcd psi(timesteps, N);

    // Set up the initial state
    for (int n = 0; n < N; n++) {
        psi(0, n) = std::exp(-1.0 * (n - n0) * (n - n0) / (2.0 * sig * sig) - i * k0 * (1.0 * n));
        psi(1, n) = std::exp(-1.0 * (n - n0) * (n - n0) / (2.0 * sig * sig) - i * k0 * (1.0 * n));
    }

    // Do the time evolution
    for (int q = 1; q < timesteps-1; q++) {
        for (int n = 0; n < N; n++) {
            psi(q+1, n) = F * (psi(q, (n+1) % N) + psi(q, (n-1 + N) % N)) + psi(q-1, n);
        }
    }

    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(tnum) + "-" + std::to_string(N) + "-" + std::to_string(n0) + "-" + std::to_string(sig) + "-" + std::to_string(k0) + ".txt");
        for (int i = 0; i < timesteps; i++) {
            for (int j = 0; j < N; j++) {
                outfile << psi(i, j) << " ";
            }
            outfile << std::endl;
        }

        outfile.close();   
    }
}

