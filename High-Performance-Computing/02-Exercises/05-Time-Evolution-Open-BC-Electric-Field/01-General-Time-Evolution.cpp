#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cmath>
#include <complex>
#include <chrono>

#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

// Energy scale is t
#define t 1


// Eigen definitions
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;


// ./a.out 1 200 0.1 256 128 4 4 0 1
int main(int argc, char* argv[]) {
    // Parse inputs
    int tnum = std::stoi(argv[1]);
    int timesteps = std::stoi(argv[2]);
    double dt = std::stod(argv[3]);
    int N = std::stoi(argv[4]);
    int n0 = std::stoi(argv[5]);
    int sig = std::stoi(argv[6]);
    double k0 = std::stod(argv[7]);
    double E = std::stod(argv[8]);
    int save = std::stoi(argv[9]);
    int time = std::stoi(argv[10]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Start counting time
    auto start = std::chrono::high_resolution_clock::now();
    
    // Define imaginary constant
    std::complex<double> i(0, 1);

    // Create an array with the wave function and its reciprocal-space transform
    Eigen::MatrixXcd psi(timesteps, N);
    psi.fill(0.);

    Eigen::MatrixXcd psiBloch(timesteps, N);
    psiBloch.fill(0.);
    
    // Set up the initial state
    for (int n = 0; n < N; n++) {
        psi(0, n) = std::exp(-1.0 * (n - n0) * (n - n0) / (2.0 * sig * sig) - i * k0 * (1.0 * n));
    }

    // Build the Hamiltonian
    MatrixXd realH(N, N);
    realH.fill(0.);

    for (int n = 0; n < N-1; n++) {
        // Electric field term
        realH(n, n) = - E * n;

        // Hopping to next atom in the lattice
        realH(n, n+1) = -t;
        realH(n+1, n) = -t;
    }

    // Last electric field term
    realH(N-1, N-1) = -E * (N-1);

    
    // Diagonalize the Hamiltonian
    SelfAdjointEigenSolver<MatrixXd> solver(realH);
    VectorXd vals = solver.eigenvalues();
    MatrixXd vecs = solver.eigenvectors();


    // Change into the Hamiltonian's eigenbasis
    // We compute the inner product between psi and the eigenvector alpha and that is the component alpha of psiBloch
    // If the Hamiltonian was complex we should use the complex conjugate of vecs, because it is real its eigenvectors are also real
    for (int alpha = 0; alpha < N; alpha++) {
        for (int n = 0; n < N; n++) {
            psiBloch(0, alpha) += psi(0, n) * vecs(n, alpha);
        }
    }

    // Time evolution
    for (int alpha = 0; alpha < N; alpha++) {
        double ek = vals[alpha];

        for (int timeIndex = 1; timeIndex < timesteps; timeIndex++) {
            psiBloch(timeIndex, alpha) = std::exp(-i * (ek * 1.) * (timeIndex * dt)) * psiBloch(0, alpha);
        }
    }

    // Parallelization only works if timesteps >> N
    // This is because diagonalization is O(N^3) while this part is O(N^2 * timesteps)
    // Moreover, this takes a lot of memory because we copy psiBloch into each thread, so beware
    #pragma omp parallel shared(psi) firstprivate(psiBloch)
    {
        
        // Temporary function
        Eigen::MatrixXcd psiTemp(timesteps, N);
        psiTemp.fill(0.);

        // Change back into position basis
        // For each n we go through all components alpha of psiBloch and take the n value of the corresponding eigenvector
        #pragma omp for collapse(2) schedule(dynamic, 10)
        for (int timeIndex = 1; timeIndex < timesteps; timeIndex++) {
            for (int n = 0; n < N; n++) {
                for (int alpha = 0; alpha < N; alpha++) {
                    psiTemp(timeIndex, n) += psiBloch(timeIndex, alpha) * vecs(n, alpha);
                }
            }
        }

        #pragma omp critical
        {   
            // Add to the final variable our intermediate result
            psi += psiTemp;
        }
    }
    
    // Output time
    if (time) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Running time for tnum = " << tnum << ": " << elapsed.count() << " seconds\n";
    }
    

    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile1(name + "-real-" + std::to_string(N) + "-" + std::to_string(n0) + "-" + std::to_string(sig) + "-" + std::to_string(k0) + "-" + std::to_string(E) + ".txt");
        for (int i = 0; i < timesteps; i++) {
            for (int j = 0; j < N; j++) {
                outfile1 << psi(i, j) << " ";
            }
            outfile1 << std::endl;
        }

        outfile1.close(); 
        
        // Write the output to a file
        std::ofstream outfile2(name + "-fourier-" + std::to_string(N) + "-" + std::to_string(n0) + "-" + std::to_string(sig) + "-" + std::to_string(k0) + "-" + std::to_string(E) + ".txt");
        for (int i = 0; i < timesteps; i++) {
            for (int j = 0; j < N; j++) {
                outfile2 << psiBloch(i, j) << " ";
            }
            outfile2 << std::endl;
        }

        outfile2.close();
    }
}

