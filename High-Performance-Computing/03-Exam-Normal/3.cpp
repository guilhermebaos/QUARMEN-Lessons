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


// Helper definitions
#define basis 2
#define Delta 0.125


// We will choose below n0 to be x-axis and n1 to be y-axis and row-major ordering
// This means that we have, inside each thread
// ^ n1 (L1)
// |   
// |   2  5  8  11
// |   1  4  7  10
// |   0  3  6  9
// O-------------------> n0 (L0)


unsigned get_index(unsigned n0, unsigned n1, unsigned alpha, unsigned L1) {
    // Row-Major Ordering
    return alpha + basis * (n1 + L1 * n0);
}


// Remember that class members are initialized in order of their declaration in the class
// Not in the order stated by the constructor
class kpm {
public:

    // Disorder and Energy Scale (W is the strength of the disorder)
    const double EnergyScale;
    const double W;
    Eigen::Array<double, -1, 1> U;

    // Parameters
    const unsigned L0;
    const unsigned L1;
    const unsigned L;
    const double t;

    // Ghost points 
    Eigen::Array<double, -1, -1> GhostLeft;
    Eigen::Array<double, -1, -1> GhostRight;

    // Computing mu
    unsigned index;
    Eigen::Array<double, -1, -1> KetXi;
    Eigen::Array<double, -1, -1> BraXi;

    // Random number generation
    std::mt19937 rnd;
    std::uniform_real_distribution<> dist;


    // Initialize RNG with a seed (8 numbers for maximum entropy according to the prof)
    kpm(unsigned ll0, unsigned ll1, double tt, double ww) :
        EnergyScale((4 * tt + ww/2) * 1.01),
        W(ww / EnergyScale),
        U(ll0 * ll1 * basis),

        L0(ll0),
        L1(ll1),
        L(basis * ll0 * ll1),
        t(tt / EnergyScale),

        GhostLeft(basis * ll1, 2),
        GhostRight(basis * ll1, 2),

        KetXi(L, 2),
        BraXi(1, L)
    {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        rnd.seed(seed);
    }

    // The starting vector is just a random vector because the Chebyshev polynomial of order 0 is 1
    void initialize() {
        index = 0;
        for(unsigned i = 0; i < L; i++) {
            KetXi(i, 0) = (dist(rnd) - 0.5) * 2 * sqrt(3);
            KetXi(i, 1) = 0;
            BraXi(0, i) = KetXi(i, index);
            U(i) = W * (dist(rnd) - 0.5);
        }
    }

    // Multiply by the the Hamiltonian for the system model
    // We write the multiplication by hand because we know the structure
    // If we didn't, we could try to use sparse matrices, but it would make it less efficient
    void hamiltonian_kpm(int factor) {

        // factor == 1 is initialization, factor == 2 is a recursive KPM step
        if (factor == 1) index = 1;
        else             index++;
    
        unsigned i0 = index % 2;
        unsigned i1 = (index - 1) % 2;
    
        for (unsigned n0 = 0; n0 < L0; n0++) {
            for (unsigned n1 = 0; n1 < L1; n1++) {
                for (unsigned alpha = 0; alpha < basis; alpha++) {
                    // Trick because basis = 2
                    unsigned beta = (alpha + 1) % 2;

                    // Get current index
                    int me = get_index(n0, n1, alpha, L1);
            
                    // Disorder and Chebyshev factor (notice that KetXi(me, i0) = 0 in the first iteration, as desired)
                    KetXi(me, i0) = factor * U(me) * KetXi(me, i1) - KetXi(me, i0);

                    // Self-energy
                    KetXi(me, i0) += factor * (alpha == 0 ? -1. : 1.) * (Delta / 2.) * KetXi(me, i1);

                    // Same‐cell t1 hop
                    int same_cell = get_index(n0, n1, beta, L1);
                    KetXi(me, i0) += -factor * (t / 2.) * KetXi(same_cell, i1);
            
                    // Same-orbital t0 hop to the left, same n1
                    if (n0 != 0) {
                        int left = get_index(n0 - 1, n1, alpha, L1);
                        KetXi(me, i0) += -factor * t * KetXi(left, i1);
                    } else {
                        // Get left ghosts at height n1
                        KetXi(me, i0) += -factor * t * GhostLeft(basis*n1 + alpha, i1);
                    }

                    // Same-orbital t0 hop to the right, same n1
                    if (n0 != L0-1) {
                        int right  = get_index(n0 + 1, n1, alpha, L1);
                        KetXi(me, i0) += -factor * t * KetXi(right, i1);
                    } else {
                        // Get right ghosts at height n1
                        KetXi(me, i0) += -factor * t * GhostRight(basis*n1 + alpha, i1);
                    }
            
                    // Split by which sublattice alpha is:
                    if (alpha == 0) {

                        // Hopping t1 with sites to the left, same n1
                        if (n0 != 0) {
                            int left  = get_index(n0 - 1, n1, beta, L1);
                            KetXi(me, i0) += -factor * (t/2.) * KetXi(left, i1);
                        } else {
                            // Get left ghosts at height n1
                            KetXi(me, i0) += -factor * (t/2.) * GhostLeft(basis*n1 + beta, i1);
                        }

                        // Hopping t2 with sites to the right, same n1
                        if (n0 != L0-1) {
                            int right = get_index(n0 + 1, n1, beta, L1);
                            KetXi(me, i0) += -factor * (t/4.) * KetXi(right, i1);
                        } else {
                            // Get right ghosts at height n1
                            KetXi(me, i0) += -factor * (t/4.) * GhostRight(basis*n1 + beta, i1);
                        }

                        if (n1 != 0) {
                            // For A‐sites, we have t1 hopping with B-sites below
                            int down = get_index(n0, n1 - 1, beta, L1);
                            KetXi(me, i0) += -factor * (t/2.) * KetXi(down, i1);

                            // Hopping t2 with sites to the left
                            if (n0 != 0) {
                                int left  = get_index(n0 - 1, n1 - 1, beta, L1);
                                KetXi(me, i0) += -factor * (t/4.) * KetXi(left, i1);
                            } else {
                                // Get left ghosts at height n1 - 1
                                KetXi(me, i0) += -factor * (t/4.) * GhostLeft(basis*(n1-1) + beta, i1);
                            }
                        }
                    }
                    else {

                        // Hopping t1 with sites to the right, same n1
                        if (n0 != L0-1) {
                            int right  = get_index(n0 + 1, n1, beta, L1);
                            KetXi(me, i0) += -factor * (t/2.) * KetXi(right, i1);
                        } else {
                            // Get left ghosts at height n1
                            KetXi(me, i0) += -factor * (t/2.) * GhostRight(basis*n1 + beta, i1);
                        }

                        // Hopping t2 with sites to the left, same n1
                        if (n0 != 0) {
                            int left = get_index(n0 - 1, n1, beta, L1);
                            KetXi(me, i0) += -factor * (t/4.) * KetXi(left, i1);
                        } else {
                            // Get left ghosts at height n1
                            KetXi(me, i0) += -factor * (t/4.) * GhostLeft(basis*n1 + beta, i1);
                        }
                        
                        if (n1 != L1-1) {
                            // For B‐sites, we have t1 hopping with A-sites above
                            int up = get_index(n0, n1 + 1, beta, L1);
                            KetXi(me, i0) += -factor * (t/2.) * KetXi(up, i1);

                            // Hopping t2 with sites to the right
                            if (n0 != L0-1) {
                                int right = get_index(n0 + 1, n1 + 1, beta, L1);
                                KetXi(me, i0) += -factor * (t/4.) * KetXi(right, i1);
                            } else {
                                // Get right ghosts at height n1 + 1
                                KetXi(me, i0) += -factor * (t/4.) * GhostRight(basis*(n1+1) + beta, i1);
                            }
                        }
                    }
                }
            }
        }
    }


    // Update ghost points with data from other threads
    void pull_ghosts(int posThread, int tnum, Eigen::Array<double, -1, 1>& GhostsPublic) {
        int border = basis * L1;

        GhostLeft.col(index % 2) = GhostsPublic.segment(basis * border * (posThread - 1 + tnum) % tnum + border, border);
        GhostRight.col(index % 2) = GhostsPublic.segment(basis * border * (posThread + 1) % tnum, border);
    }

    void push_ghosts(int posThread, Eigen::Array<double, -1, 1>& GhostsPublic) {
        int border = basis * L1;

        GhostsPublic.segment(2 * border * posThread, border) = KetXi.col(index % 2).segment(0, border);
        GhostsPublic.segment(2 * border * posThread + border, border) = KetXi.col(index % 2).segment(L - border, border);
    }
};


// ./a.out 128 1 0 2 1000 1 1 0
int main(int argc, char* argv[]) {
    // Parse inputs
    unsigned ll0 = std::stoi(argv[1]);
    unsigned ll1 = std::stoi(argv[2]);
    double ww = std::stof(argv[3]);
    unsigned NAverages = std::stoi(argv[4]);
    unsigned NMoments = std::stoi(argv[5]);
    unsigned tnum = std::stoi(argv[6]);
    unsigned save = std::stoi(argv[7]);
    unsigned time = std::stoi(argv[8]);

    // Parameter
    double tt = 1.0;

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);

    // Start counting time
    auto start = std::chrono::high_resolution_clock::now();

    // Make sure the number of moments is even
    NMoments += NMoments % 2;

    // Make sure system size is divisible by tnum
    ll0 -= ll0 % tnum;

    // Compute the moments 
    Eigen::Array<double, -1, 1> mu(NMoments, 1);
    mu.setZero();

    // Shared array of ghost points (left0, right0, left1, right1, ...)
    Eigen::Array<double, -1, 1> ghosts(tnum * 2 * basis * ll1, 1);
    ghosts.setZero();

    // Each thread will have their own instance of the kpm class
    #pragma omp parallel firstprivate(ll0, ll1, tt, ww, NAverages, NMoments) shared(ghosts)
    {
        // Get thread number
        unsigned posThread = omp_get_thread_num();

        // Divide the system size by the number of threads
        ll0 /= tnum;

        // Initialize the class
        kpm kpmobj(ll0, ll1, tt, ww);

        // Thread-specific moments
        Eigen::Array<double, -1, 1> muThread(NMoments, 1);
        muThread.setZero();

        // Compute moments
        for (unsigned av = 0; av < NAverages; av++) {
            // Initialize, push ghost points, sync and get ghost points from other threads
            kpmobj.initialize();
            kpmobj.push_ghosts(posThread, ghosts);
            #pragma omp barrier

            kpmobj.pull_ghosts(posThread, tnum, ghosts);
            #pragma omp barrier
            
            
            // Do first iteration
            kpmobj.hamiltonian_kpm(1);
            kpmobj.push_ghosts(posThread, ghosts);
            #pragma omp barrier

            kpmobj.pull_ghosts(posThread, tnum, ghosts);
            #pragma omp barrier

            // Do the two inner products and compute average
            // Notice that we only access BraXi once
            muThread.segment(0, 2) += ((kpmobj.BraXi.matrix() * kpmobj.KetXi.matrix()).array().transpose() - muThread.segment(0, 2)) / double(av + 1);

            for(unsigned m = 2; m < NMoments; m += 2) {
                for(unsigned j = 0; j < 2; j++) {
                    kpmobj.hamiltonian_kpm(2);
                    kpmobj.push_ghosts(posThread, ghosts);
                    #pragma omp barrier

                    kpmobj.pull_ghosts(posThread, tnum, ghosts);
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
    mu = mu / double(ll0 * ll1);

    // Output time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if (time) {
        std::cout << "Running time for tnum = " << tnum << ": " << elapsed.count() << " seconds\n";
    }
    
    if (save) {
        // Construct the filename
        std::string name = "out3";

        // Write the output to a file
        std::ofstream outfile(name + "-" + std::to_string(ll0) + "-" + std::to_string(ll1) + "-" + std::to_string(ww) + "-" + std::to_string(NAverages) + "-" + std::to_string(NMoments)  + ".txt");
        for (unsigned i = 0; i < NMoments; i++) {
            outfile << mu[i] << " ";
            outfile << std::endl;
        }
        outfile.close(); 
    }
}


// This program is very fast, hence we can go to very large system sizes
// In the limit of a very large system the majority of the time might be spent moving memory to the CPU cache