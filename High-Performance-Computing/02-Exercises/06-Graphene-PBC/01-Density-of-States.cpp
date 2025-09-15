#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cmath>
#include <complex>

#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>


// Eigen definitions
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;


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



int get_index(int n0, int n1, int alpha, int N) {
    // Row-Major Ordering
    return alpha + 2 * (n1 + N * n0);

    // Column-Major Ordering
    // return n0 + N * (n1 + N * alpha);
}

std::array<int, 3> get_coords(int index, int N) {
    // Row-Major Ordering
    int alpha = index % 2;
    int n1 = ((index - alpha) % N) / 2;
    int n0 = (index - alpha - 2 * n1) / (2 * N);
    return {n0, n1, alpha};
}


int main(int argc, char* argv[]) {
    // Parse inputs
    int tnum = std::stoi(argv[1]);
    int N = std::stoi(argv[2]);
    double gam = std::stof(argv[3]);
    double eA = std::stof(argv[4]);
    double eB = std::stof(argv[5]);
    double t = std::stof(argv[6]);
    double dosMax = std::stof(argv[7]);
    int dosPos = std::stoi(argv[8]);
    int save = std::stoi(argv[9]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Build the energy vector
    VectorXd ee(dosPos);
    for (int id = 0; id < dosPos; id++) {
        ee[id] = -dosMax + (2. * id / dosPos) * dosMax;
    }


    // Size of the matrix
    int sizeH = get_index(N - 1, N - 1, 1, N) + 1;

    // Build the Hamiltonian
    MatrixXd realH(sizeH, sizeH);
    realH.fill(0.);

    for (int n0 = 0; n0 < N; n0++) {
        for (int n1 = 0; n1 < N; n1++) {
            // Local energy
            int temp0 = get_index(n0, n1, 0, N);
            int temp1 = get_index(n0, n1, 1, N);

            realH(temp0, temp0) = eA;
            realH(temp1, temp1) = eB;
            
            // Same cell hopping
            realH(temp0, temp1) = -t;
            realH(temp1, temp0) = -t;

            // Change in n1
            int temp2 = get_index(n0, (n1 - 1 + N) % N, 1, N);
            realH(temp0, temp2) = -t;
            realH(temp2, temp0) = -t;

            // Change in n0
            int temp3 = get_index((n0 - 1 + N) % N, n1, 1, N);
            realH(temp0, temp3) = -t;
            realH(temp3, temp0) = -t;
        }
    }

    // Diagonalize the Hamiltonian
    SelfAdjointEigenSolver<MatrixXd> solver(realH);
    VectorXd vals = solver.eigenvalues();
    MatrixXd vecs = solver.eigenvectors();

    // Sum the Dirac deltas centered at each eigenenergy
    VectorXd dos(dosPos);
    dos.fill(0.);
    for (int ie = 0; ie < sizeH; ie++) {
        for (int id = 0; id < dosPos; id++) {
            dos[id] += nascent_delta(ee[id] - vals[ie], gam, 1);
        }
    }
    dos /= N*N;

    // Save the results
    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile(name + "-"  + std::to_string(N) + "-" + std::to_string(gam) + ".txt");
        for (int id = 0; id < dosPos; id++) {
            outfile << ee[id] << " ";
        }
        outfile << std::endl;
        for (int id = 0; id < dosPos; id++) {
            outfile << dos[id] << " ";
        }
        outfile.close();
    }
}
