#include <iostream>
#include <fstream>
#include <cstdlib>
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

        case 2:
            return exp(-x * x / (eps * eps)) / (eps * sqrt(M_PI));
        
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

    // Square lattice of size N x N
    int N2 = N * N;

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

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
    VectorXd dosExact(points + 1);


    #pragma omp parallel
    {
        // Private intermediate result
        VectorXd dosPrivate(points + 1);
        dosPrivate.fill(0.);

        // Private intermediate result
        VectorXd dosExactPrivate(points + 1);
        dosExactPrivate.fill(0.);

        // Initial and final indexes
        int mMin = pbc ? 0 : 1;
        int mMax = pbc ? N-1 : N;

        // Calculate the dos at many points
        #pragma omp for collapse(3) schedule(dynamic, 5) 
            
            // Calculate DOS
            for (int mx = mMin; mx <= mMax; mx++) {
                
                for (int my = mMin; my <= mMax; my++) {
                    // Calculate energy for this k
                    double ek = pbc ? -2 * (cos(2 * M_PI * mx / N) + cos(2 * M_PI * my / N)) : -2 * (cos(mx * M_PI / (N + 1)) + cos(my * M_PI / (N + 1)));
                    
                    // Calculate the effect of this energy on the DOS
                    for (int i = 0; i <= points; i++) {
                        dosPrivate[i] += nascent_delta(ee[i] - ek, gam, 1) / N2;
                    }
                }
            }
        
        #pragma omp for schedule(dynamic, 5)
            
            // Calculate exact DOS
            for (int i = 0; i <= points; i++) {
                double ek = ee[i];
                
                // Calculate the effect 
                if (ek * ek != 4) {
                    // dosExactPrivate[i] += 1 / (2 * M_PI * sqrt(1 - ek * ek / 4));
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
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(tnum) + "-" + std::to_string(pbc) + "-" + std::to_string(N) + "-" + std::to_string(gam) + ".txt");
        for (int i = 0; i <= points; i++) {
            outfile << ee[i] << " ";
            outfile << dos[i] << " ";
            outfile << dosExact[i] << " ";
            outfile << std::endl;
        }

        outfile.close();   
    }
}

