#include <iostream>
#include <fstream>
#include <cstdlib>

#include <utility>  // for std::pair
#include <chrono>

#include <cmath>
#include <complex>

#include <random>
#include <Eigen/Core>

#include <omp.h>


// Remember that class members are initialized in order of their declaration in the class
// Not in the order stated by the constructor
class kpm {
public:

    // Disorder and Energy Scale (W is the strength of the disorder)
    const double EnergyScale;
    const double W;
    Eigen::Array<double, -1, 1> U;

    // Parameters
    const unsigned L;
    const double t;

    // Ghost points 
    Eigen::Array<double, 1, 2> GhostLeft;
    Eigen::Array<double, 1, 2> GhostRight;

    // Computing mu
    unsigned index;
    Eigen::Array<double, -1, -1> KetXi;
    Eigen::Array<double, -1, -1> BraXi;

    // Random number generation
    std::mt19937 rnd;
    std::uniform_real_distribution<> dist;


    // Initialize RNG with a seed (8 numbers for maximum entropy according to the prof)
    kpm(unsigned ll, double tt, double ww) :
        EnergyScale((2 * tt + ww/2) * 1.01),
        W(ww / EnergyScale),
        U(ll),

        L(ll),
        t(tt / EnergyScale),

        GhostLeft(Eigen::Array<double, 1, 2>::Zero()),
        GhostRight(Eigen::Array<double, 1, 2>::Zero()),

        KetXi(ll, 2),
        BraXi(1, ll)
    {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        rnd.seed(seed);
    }

    // The starting vector is just a random vector because the Chebyshev polynomial of order 0 is 1
    void initialize() {
        index = 0;
        for(unsigned i = 0; i < L; i++) {
            KetXi(i, index) = (dist(rnd) - 0.5) * 2 * sqrt(3);
            BraXi(0, i) = KetXi(i, index);
            U(i) = W * (dist(rnd) - 0.5);
        }
    }

    // Multiply by the the Hamiltonian for a 1D tight-binding model
    // We write the multiplication by hand because we know the structure
    // If we didn't, we could try to use sparse matrices, but it would make it less efficient
    void hamiltonian_kpm(int factor) {

        // factor == 1 is initialization, factor == 2 is a recursive KPM step
        if (factor == 1) index = 1;
        else             index++;
    
        unsigned i0 = index % 2;
        unsigned i1 = (index - 1) % 2;
    
        // We need to store two vectors, i0 is the n index and i1 is the n-1 index
        KetXi(0, i0) = -factor * t * (KetXi(1, i1) + GhostLeft(i1)) + factor * U(0) * KetXi(0, i1) - KetXi(0, i0);
        KetXi(L-1, i0) = -factor * t * (GhostRight(i1) + KetXi(L - 2, i1)) + factor * U(L - 1) * KetXi(L - 1, i1) - KetXi(L-1, i0);
        
        // Notice that compared to the serial code we don't use % to get the correct indexes
        // This is because the indexes don't wrap around, they stop at the ghosts
        for(unsigned i = 1; i < L-1; i++) {
            KetXi(i, i0) = -factor * t * (KetXi(i+1, i1) + KetXi(i-1, i1)) + factor * U(i) * KetXi(i, i1) - KetXi(i, i0);
        }
    }


    // Update ghost points with data from other threads
    void pull_ghosts(double ValueLeft, double ValueRight) {
        GhostLeft(index % 2) = ValueLeft;
        GhostRight(index % 2) = ValueRight;
    }


    // Update ghost points with data from other threads
    double left() {
        return KetXi(0, index % 2);
    }
    double right() {
        return KetXi(L-1, index % 2);
    }
};


// ./a.out 1000 1 0 2 1000 2000 2 1 0
int main(int argc, char* argv[]) {
    // Parse inputs
    unsigned ll = std::stoi(argv[1]);
    double tt = std::stof(argv[2]);
    double ww = std::stof(argv[3]);
    unsigned NAverages = std::stoi(argv[4]);
    unsigned NMoments = std::stoi(argv[5]);
    unsigned tnum = std::stoi(argv[6]);
    unsigned save = std::stoi(argv[7]);
    unsigned time = std::stoi(argv[8]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Start counting time
    auto start = std::chrono::high_resolution_clock::now();

    // Make sure the number of moments is even
    NMoments += NMoments % 2;

    // Compute the moments 
    Eigen::Array<double, -1, 1> mu(NMoments, 1);
    mu.setZero();

    // Shared array of ghost points
    Eigen::Array<double, -1, 1> ghosts(tnum * 2, 1);
    ghosts.setZero();

    // Each thread will have their own instance of the kpm class
    #pragma omp parallel firstprivate(ll, tt, ww, NAverages, NMoments) shared(ghosts)
    {
        // Get thread number
        unsigned posThread = omp_get_thread_num();

        // Divide the system size by the number of threads
        ll /= tnum;

        // Initialize the class
        kpm kpmobj(ll, tt, ww);

        // Thread-specific moments
        Eigen::Array<double, -1, 1> muThread(NMoments, 1);
        muThread.setZero();

        // Compute moments
        for (unsigned av = 0; av < NAverages; av++) {
            kpmobj.initialize();

            // Push ghost points (we store left1, right1, left2, right2, ...)
            ghosts[2 * posThread] = kpmobj.left();
            ghosts[2 * posThread + 1] = kpmobj.right();
            #pragma omp barrier

            // Get ghost points (right of thread+1 and left of thread-1)
            kpmobj.pull_ghosts(ghosts[2 * ((posThread - 1 + tnum) % tnum) + 1], ghosts[2 * ((posThread + 1) % tnum)]);
            #pragma omp barrier
            
            // Do first iteration, push, sync and pull
            kpmobj.hamiltonian_kpm(1);
            ghosts[2 * posThread] = kpmobj.left();
            ghosts[2 * posThread + 1] = kpmobj.right();
            #pragma omp barrier

            kpmobj.pull_ghosts(ghosts[2 * ((posThread - 1 + tnum) % tnum) + 1], ghosts[2 * ((posThread + 1) % tnum)]);
            #pragma omp barrier
            

            // Do the two inner products and compute average
            // Notice that we only access BraXi once
            muThread.segment(0, 2) += ((kpmobj.BraXi.matrix() * kpmobj.KetXi.matrix()).array().transpose() - muThread.segment(0, 2)) / double(av + 1);

            for(unsigned m = 2; m < NMoments; m += 2) {
                for(unsigned j = 0; j < 2; j++) {
                    kpmobj.hamiltonian_kpm(2);   
                    
                    ghosts[2 * posThread] = kpmobj.left();
                    ghosts[2 * posThread + 1] = kpmobj.right();
                    #pragma omp barrier

                    kpmobj.pull_ghosts(ghosts[2 * ((posThread - 1 + tnum) % tnum) + 1], ghosts[2 * ((posThread + 1) % tnum)]);
                    #pragma omp barrier
                }
                
                muThread.segment(m, 2) += ((kpmobj.BraXi.matrix() * kpmobj.KetXi.matrix()).array().transpose() - muThread.segment(m, 2)) / double(av + 1);
            }
        }

        #pragma omp critical
        {
            mu += muThread;
        }
    }

    // Normalize
    mu = mu / double(ll);

    // Output time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if (time) {
        std::cout << "Running time for tnum = " << tnum << ": " << elapsed.count() << " seconds\n";
    }
    
    if (save) {
        // Construct the filename
        std::string name = "out";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(ll) + "-" + std::to_string(tt) + "-" + std::to_string(ww) + "-" + std::to_string(NAverages) + "-" + std::to_string(NMoments)  + ".txt");
        for (unsigned i = 0; i < NMoments; i++) {
            outfile << mu[i] << " ";
            outfile << std::endl;
        }
        outfile.close(); 
    }
}


// This program is very fast, hence we can go to very large system sizes
// In the limit of a very large system the majority of the time might be spent moving memory to the CPU cache