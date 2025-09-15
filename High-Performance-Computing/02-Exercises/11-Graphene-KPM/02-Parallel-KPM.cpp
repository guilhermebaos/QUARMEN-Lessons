#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cmath>
#include <complex>
#include <chrono>

#include <random>
#include <Eigen/Core>

#include <omp.h>


// Remeber that corner points have *two neighbors*
// This means that each thread has to store its entire perimeter as ghosts
// If a thread is a L by L square then it has 4L * basis neighbors

// Number of elements in the basis
#define basis 2


// We will choose below n0 to be x-axis and n1 to be y-axis and row-major ordering
// This means that we have, inside each thread
// ^ n1
// |   
// |   3  7  11 15
// |   2  6  10 14
// |   1  5  9  13
// |   0  4  8  12
// O-------------------> n0


// Moreover, our threads are ordered row-major as well, but in the usual order
// For example, if we have 9 threads they will be placed in the 2D material as:
// O--------------> tx   
// |   0  1  2
// |   3  4  5 
// |   6  7  8 
// |   
// v ty


int get_index(int n0, int n1, int alpha, unsigned N) {
    // Row-Major Ordering
    return alpha + basis * (n1 + N * n0);

    // Column-Major Ordering
    // return n0 + N * (n1 + N * alpha);
}

// std::array<int, 3> get_coords(int index, unsigned N) {
//     // Row-Major Ordering
//     int alpha = index % basis;
//     int n1 = ((index - alpha) / basis) % N;
//     int n0 = (index - alpha - basis * n1) / (basis * N);
//     return {n0, n1, alpha};
// }

// Variable `side` is top=0, right=1, bottom=2 and left=3
// Variable `posThread` starts counting at 0, not 1
unsigned ghost_index(unsigned L, unsigned posThread, unsigned side) {
    return (4 * posThread + side) * basis * L;
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
    const unsigned L;
    const double t;

    // Computing mu
    unsigned index;
    Eigen::Array<double, -1, -1> KetXi;
    Eigen::Array<double, -1, -1> BraXi;

    // Ghost points 
    Eigen::Array<double, -1, -1> GhostsPrivate;

    // Random number generation (dist is by default from -1 to +1 which is what we want)
    std::mt19937 rnd;
    std::uniform_real_distribution<> dist;


    // Initialize RNG with a seed (8 numbers for maximum entropy according to the prof)
    // Notice that ll is the linear size for the lattice, meaning we have a ll*ll lattice
    kpm(unsigned ll, double tt, double ww) :
        EnergyScale((4 * tt + ww/2) * 1.0001),
        W(ww / EnergyScale),
        U(2*ll*ll),

        L(ll),
        t(tt / EnergyScale),

        KetXi(basis*ll*ll, 2),
        BraXi(1, basis*ll*ll),

        GhostsPrivate(4 * ll * basis, 2)
    {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        rnd.seed(seed);
    }

    // The starting vector is just a random vector because the Chebyshev polynomial of order 0 is 1
    void initialize() {
        index = 0;
        for(unsigned i = 0; i < basis*L*L; i++) {
            KetXi(i, 0) = (dist(rnd) - 0.5) * 2 * sqrt(3);
            KetXi(i, 1) = 0;
            BraXi(0, i) = KetXi(i, index);
            U(i) = W * (dist(rnd) - 0.5);
        }
    }

    // Multiply by the the Hamiltonian for a 1D tight-binding model
    // We write the multiplication by hand because we know the structure
    // If we didn't, we could try to use sparse matrices, but it would make it less efficient
    // Factor is 1 after initialization and 2 when doing KPM iteration
    void hamiltonian_kpm(int factor) {
        if (factor == 1) {
            index = 1;
        } else {
            index++;
        }
    
        unsigned i0 = index % 2;
        unsigned i1 = (index - 1) % 2;
    
        for (unsigned n0 = 0; n0 < L; n0++) {
            for (unsigned n1 = 0; n1 < L; n1++) {
                for (unsigned alpha = 0; alpha < basis; alpha++) {
                    // Special trick because basis = 2
                    unsigned beta = (alpha + 1) % 2;

                    // Get current index
                    int me  = get_index(n0, n1, alpha, L);
            
                    // Disorder and Chebyshev factor (notice that KetXi(me, i0) = 0 in the first iteration, as desired)
                    KetXi(me, i0) = factor * U(me) * KetXi(me, i1) - KetXi(me, i0);
            
                    // Same‐cell hop:  A <--> B in (n0,n1)
                    int same_cell = get_index(n0, n1, beta, L);
                    KetXi(me, i0) += -factor * t * KetXi(same_cell, i1);
            
                    // Split by which sublattice alpha is:
                    unsigned border = basis * L;
                    if (alpha == 0) {
                        // For A‐sites, the two “inter‐cell” neighbors are:
                        //   B at (n0, n1−1)   and  B at (n0−1, n1)
                        if (n0 != 0) {
                            int left  = get_index(n0 - 1, n1, beta, L);
                            KetXi(me, i0) += -factor * t * KetXi(left, i1);
                        } else {
                            // Get ghosts at position 3=left
                            KetXi(me, i0) += -factor * t * GhostsPrivate(3*border + basis * n1 + beta, i1);
                        }

                        if (n1 != 0) {
                            int down  = get_index(n0, n1 - 1, beta, L);
                            KetXi(me, i0) += -factor * t * KetXi(down, i1);
                        } else {
                            // Get ghosts at position 2=down
                            KetXi(me, i0) += -factor * t * GhostsPrivate(2*border + basis * n0 + beta, i1);
                        }
                        
                    }
                    else {
                        // For B‐sites, the two “inter‐cell” neighbors are:
                        //   A at (n0, n1+1)   and  A at (n0+1, n1)
                        
                        if (n0 != L-1) {
                            int right = get_index(n0 + 1, n1, beta, L);
                            KetXi(me, i0) += -factor * t * KetXi(right, i1);
                        } else {
                            // Get ghosts at position 1=right
                            KetXi(me, i0) += -factor * t * GhostsPrivate(1*border + basis * n1 + beta, i1);
                        }

                        if (n1 != L-1) {
                            int up    = get_index(n0, n1 + 1, beta, L);
                            KetXi(me, i0) += -factor * t * KetXi(up, i1);
                        } else {
                            // Get ghosts at position 0=top
                            KetXi(me, i0) += -factor * t * GhostsPrivate(basis * n0 + beta, i1);
                        }
                    }
                }
            }
        }
    }



    // Update ghost points with data from other threads
    // We order the ghosts clockwise top -> right -> bottom -> left
    void pull_ghosts(const Eigen::Array<unsigned, -1, 1>& GhostStart, const Eigen::Array<double, -1, 1>& GhostsPublic) {
        // Where to store on our private vector
        unsigned i0 = index % 2;

        // Store in order: Top, right, bottom, left
        unsigned border = basis * L;
        for (unsigned i = 0; i < 4; i++) {
            GhostsPrivate.col(i0).segment(i * border, border) = GhostsPublic.segment(GhostStart[i], border);
        }
    }

    void push_ghosts(int posThread, Eigen::Array<double, -1, 1>& GhostsPublic) {
        // Where to store on our private vector
        unsigned i0 = index % 2;
        
        // Retrieve in order: Top, right, bottom, left
        unsigned border = basis * L;
        for (unsigned i = 0; i < 4; i++) {

            // For left and right we can just copy the segment in order because of the ordering we chose
            if (i == 1) {
                GhostsPublic.segment(ghost_index(L, posThread, i), border) = KetXi.col(i0).segment(basis*L*L - border, border);
            }
            else if (i == 3) {
                GhostsPublic.segment(ghost_index(L, posThread, i), border) = KetXi.col(i0).segment(0, border);
            }

            // For top and bottom we have to copy them in order
            else if (i == 0) {
                for (unsigned n0 = 0; n0 < L; n0++) {
                    for (unsigned alpha = 0; alpha < basis; alpha++) {
                        GhostsPublic(ghost_index(L, posThread, i) + basis * n0 + alpha) = KetXi(get_index(n0, L-1, alpha, L), i0);
                    }
                }
            }
            else if (i == 2) {
                for (unsigned n0 = 0; n0 < L; n0++) {
                    for (unsigned alpha = 0; alpha < basis; alpha++) {
                        GhostsPublic(ghost_index(L, posThread, i) + basis * n0 + alpha) = KetXi(get_index(n0, 0, alpha, L), i0);
                    }
                }
            }
        }
    }
};


// ./a.out 100 1 0 2 1024 1 1 0
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

    // Check that tnum is a perfect square
    int tnumSqrt = static_cast<int>(sqrt(tnum));
    if (tnumSqrt * tnumSqrt != tnum) {
        return -1;
    }

    // Divide the system size by the square root of the number of threads
    unsigned llThread = ll / tnumSqrt;

    // Start counting time
    auto start = std::chrono::high_resolution_clock::now();

    // Make sure the number of moments is even
    NMoments += NMoments % 2;

    // Compute the moments 
    Eigen::Array<double, -1, 1> mu(NMoments, 1);
    mu.setZero();

    // Shared array of ghost points
    Eigen::Array<double, -1, 1> ghosts(tnum * (basis * 4 * llThread), 1);
    ghosts.setZero();

    // Each thread will have their own instance of the kpm class
    #pragma omp parallel firstprivate(ll, tt, ww, NAverages, NMoments) shared(ghosts)
    {
        // Get thread number
        unsigned posThread = omp_get_thread_num();

        // Aray with the index of where to get the ghosts from
        Eigen::Array<unsigned, -1, 1> ghostIndexes(4, 1);

        // Get this threads' position on the 2D map of threads (View above the ordering)
        unsigned tx = posThread % tnumSqrt;
        unsigned ty = posThread / tnumSqrt;

        // Get the neighboring threads
        unsigned upThread    = ((ty - 1 + tnumSqrt) % tnumSqrt) * tnumSqrt + tx;
        unsigned downThread  = ((ty + 1) % tnumSqrt) * tnumSqrt + tx;
        unsigned leftThread  = ty * tnumSqrt + ((tx - 1 + tnumSqrt) % tnumSqrt);
        unsigned rightThread = ty * tnumSqrt + ((tx + 1)     % tnumSqrt);

        // Get index of the neighbor threads in the public ghosts array
        // Top: Get bottom from the thread on top of this one
        ghostIndexes[0] = ghost_index(llThread, upThread, 2);
        // Right: Get left from the thread on the right of this one
        ghostIndexes[1] = ghost_index(llThread, rightThread, 3);
        // Bottom: Get top from the thread below this one
        ghostIndexes[2] = ghost_index(llThread, downThread, 0);
        // Left: Get right from the thread on the left of this one
        ghostIndexes[3] = ghost_index(llThread, leftThread, 1);

        // Initialize the class
        kpm kpmobj(llThread, tt, ww);

        // Thread-specific moments
        Eigen::Array<double, -1, 1> muThread(NMoments, 1);
        muThread.setZero();

        // Compute moments
        for (unsigned av = 0; av < NAverages; av++) {
            // Initialize, push ghost points, sync and get ghost points from other threads
            kpmobj.initialize();
            kpmobj.push_ghosts(posThread, ghosts);
            #pragma omp barrier

            kpmobj.pull_ghosts(ghostIndexes, ghosts);
            #pragma omp barrier
            
            // Do first iteration
            kpmobj.hamiltonian_kpm(1);
            kpmobj.push_ghosts(posThread, ghosts);
            #pragma omp barrier

            kpmobj.pull_ghosts(ghostIndexes, ghosts);
            #pragma omp barrier
            

            // Do the two inner products and compute average
            // Notice that we only access BraXi once
            muThread.segment(0, 2) += ((kpmobj.BraXi.matrix() * kpmobj.KetXi.matrix()).array().transpose() - muThread.segment(0, 2)) / double(av + 1);

            for(unsigned m = 2; m < NMoments; m += 2) {
                for(unsigned j = 0; j < 2; j++) {
                    kpmobj.hamiltonian_kpm(2);
                    kpmobj.push_ghosts(posThread, ghosts);
                    #pragma omp barrier
                    
                    kpmobj.pull_ghosts(ghostIndexes, ghosts);
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
    mu = mu/double(ll * ll);

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