#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <Eigen/Core>

// Parameters
#define points 2048


// Energy boundaries for t = 1
#define emin -2.
#define emax 2.


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


int main(int argc, char* argv[]) {
    // Parse inputs
    int tnum = std::stoi(argv[1]);
    int pbc = std::stoi(argv[2]);
    int N = std::stoi(argv[3]);
    float gam = std::stof(argv[4]);
    int save = std::stoi(argv[5]);
    float emaxfrac = std::stof(argv[6]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Step in energy
    double erange = emax - emin;
    double erangelim = erange * emaxfrac;
    double estep = erange / (points * 1.0 / emaxfrac);

    // Energy values used
    VectorXd ee(points + 1);
    ee.fill(0.);
    for (int i = 0; i <= points; i++) {
        ee[i] = -erangelim / 2 + estep * i;
    }


    // Create an array with the final result
    MatrixXd dos(N, points + 1);
    MatrixXd dosExact(N, points +1);


    #pragma omp parallel
    {
        // Private intermediate result
        MatrixXd dosPrivate(N, points + 1);
        dosPrivate.fill(0.);

        // Private intermediate result
        MatrixXd dosExactPrivate(N, points + 1);
        dosExactPrivate.fill(0.);

        // Initial and final indexes
        int mMin = pbc ? 0 : 1;
        int mMax = pbc ? N-1 : N;

        // Calculate the dos at many points
        #pragma omp for collapse(2) schedule(dynamic, 5)
            
            // Calculate DOS
            for (int m = mMin; m <= mMax; m++) {
                // Calculate energy for this k
                double k = pbc ? 2 * M_PI * m / N : m * M_PI / (N + 1);
                double ek = -2 * cos(k);
                
                // Calculate the effect of this energy on the DOS
                for (int i = 0; i <= points; i++) {
                    
                    // Calculate LDOS for each n
                    for (int n = 0; n < N; n++) {
                        dosPrivate(n, i) += pbc ? nascent_delta(ee[i] - ek, gam, 1) / (N * N) : (2 * std::pow(sin(k * (n + 1)), 2)) * nascent_delta(ee[i] - ek, gam, 1) / (N + 1);
                    }
                }
            }
        
        #pragma omp for schedule(dynamic, 5)
            
            // Calculate exact DOS
            for (int i = 0; i <= points; i++) {
                double ek = ee[i];

                // Calculate LDOS for each n
                for (int n = 0; n < N; n++) { 
                    if (ek * ek != 4) {
                        dosExactPrivate(n, i) += pbc ? 1 / (2 * N * M_PI * sqrt(1 - ek * ek / 4)) : (2 * std::pow(sin(acos(- ek / 2) * (n + 1)), 2)) / (2 * (N + 1) * M_PI * sqrt(1 - ek * ek / 4));
                    }
                }                
            }
        
        #pragma omp critical
        {
            // Here we sum the results of each thread one at a time
            dos += dosPrivate;
            dosExact += dosExactPrivate;
        }
    }

    if (save) {
        // Construct the filename
        std::string name = "lout";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(tnum) + "-" + std::to_string(pbc) + "-" + std::to_string(N) + "-" + std::to_string(gam) + ".txt");
        for (int i = 0; i <= points; i++) {
            outfile << ee[i] << " ";

            for (int n = 0; n < N; n++) {
                outfile << dos(n, i) << " ";
            }

            for (int n = 0; n < N; n++) {
                outfile << dosExact(n, i) << " ";
            }

            outfile << std::endl;
        }

        outfile.close();   
    }
}

