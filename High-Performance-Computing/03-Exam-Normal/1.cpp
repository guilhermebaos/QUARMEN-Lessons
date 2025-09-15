#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <Eigen/Core>

// Parameters
#define points 2048


// Energy boundaries for t = 1
#define emin -4.
#define emax 4.


using Eigen::MatrixXd;
using Eigen::VectorXd;


double nascent_delta(double x, double eps, int func) {
    switch(func) {
        case 1:
            return eps / (M_PI * (x * x + eps * eps));
            break;

        case 2:
            return exp(-x * x / (eps * eps)) / (eps * sqrt(M_PI));
            break;
        
        default:
            return -1;
    }
}


// ./a.out 1 1 100000 0.05 1 0
int main(int argc, char* argv[]) {
    // Parse inputs
    int tnum = std::stoi(argv[1]);
    double Delta = std::stod(argv[2]);
    int N = std::stoi(argv[3]);
    float gam = std::stof(argv[4]);
    int save = std::stoi(argv[5]);
    int timing = std::stoi(argv[6]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Start counting time
    auto start = std::chrono::high_resolution_clock::now();

    // Imaginary constant
    std::complex<double> I(0, 1);

    // Step in energy
    double estep = (emax - emin) / points;

    // Energy values used
    VectorXd ee(points + 1);
    ee.fill(0.);
    for (int i = 0; i <= points; i++) {
        ee[i] = emin + estep * i;
    }


    // Create an array with the final result
    VectorXd dos(points + 1);

    #pragma omp parallel
    {
        // Private intermediate result
        VectorXd dosPrivate(points + 1);
        dosPrivate.fill(0.);

        // Initial and final indexes

        // Calculate the dos at many points
        #pragma omp for schedule(dynamic, 2)
            
            // Calculate DOS
            for (int m = 0; m < N; m++) {
                // Calculate energy for this k
                double k = 2 * M_PI * m / N;
                double cosk = cos(k);
                double ecommon = -2 * cosk;
                double tab2 = std::norm(0.5 + 0.5 * std::exp(-I * std::complex<double> (k, 0)) + 0.25 * std::exp(I * std::complex<double> (k, 0)));
                double ediff = sqrt(pow(Delta, 2) / 4 + tab2);
                
                // Calculate the effect of this energy on the DOS
                for (int i = 0; i <= points; i++) {
                    dosPrivate[i] += nascent_delta(ee[i] - (ecommon + ediff), gam, 1) / N;
                    dosPrivate[i] += nascent_delta(ee[i] - (ecommon - ediff), gam, 1) / N;
                }
            }
        
        #pragma omp critical
        {
            // Here we sum the results of each thread one at a time
            dos += dosPrivate;
        }
    }

    // Output time
    if (timing) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Running time for tnum = " << tnum << ": " << elapsed.count() << " seconds\n";
    }

    if (save) {
        // Construct the filename
        std::string name = "out1";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(tnum) + "-" + std::to_string(Delta) + "-" + std::to_string(N) + "-" + std::to_string(gam) + ".txt");
        for (int i = 0; i <= points; i++) {
            outfile << ee[i] << " ";
            outfile << dos[i] << " ";
            outfile << std::endl;
        }

        outfile.close();   
    }
}

