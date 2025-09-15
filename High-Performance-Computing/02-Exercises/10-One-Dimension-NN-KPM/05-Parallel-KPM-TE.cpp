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


// Initial state parameters (function of lattice size)
#define n0factor 0.75
#define sigfactor 0.04
#define k0 1.5


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
    Eigen::Array<std::complex<double>, 1, 2> GhostLeft;
    Eigen::Array<std::complex<double>, 1, 2> GhostRight;

    // Computing mu
    unsigned index;
    Eigen::Array<std::complex<double>, -1, -1> KetPsi;

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

        KetPsi(ll, 2)
    {
        std::random_device r;
        std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        rnd.seed(seed);
    }

    // The starting vector is just a random vector because the Chebyshev polynomial of order 0 is 1
    void initialize(Eigen::Array<std::complex<double>, -1, 1> psi) {
        index = 0;
        for(unsigned i = 0; i < L; i++) {
            KetPsi(i, 0) = psi(i);
            KetPsi(i, 1) = 0;
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
        KetPsi(0, i0) = -factor * t * (KetPsi(1, i1) + GhostLeft(i1)) + factor * U(0) * KetPsi(0, i1) - KetPsi(0, i0);
        KetPsi(L-1, i0) = -factor * t * (GhostRight(i1) + KetPsi(L - 2, i1)) + factor * U(L - 1) * KetPsi(L - 1, i1) - KetPsi(L-1, i0);
        
        // Notice that compared to the serial code we don't use % to get the correct indexes
        // This is because the indexes don't wrap around, they stop at the ghosts
        for(unsigned i = 1; i < L-1; i++) {
            KetPsi(i, i0) = -factor * t * (KetPsi(i+1, i1) + KetPsi(i-1, i1)) + factor * U(i) * KetPsi(i, i1) - KetPsi(i, i0);
        }
    }


    // Update ghost points with data from other threads
    void pull_ghosts(std::complex<double> ValueLeft, std::complex<double> ValueRight) {
        GhostLeft(index % 2) = ValueLeft;
        GhostRight(index % 2) = ValueRight;
    }


    // Update ghost points with data from other threads
    std::complex<double> left() {
        return KetPsi(0, index % 2);
    }
    std::complex<double> right() {
        return KetPsi(L-1, index % 2);
    }
};


// ./a.out 100000 1 0 0.1 200 1000 2 1 0
int main(int argc, char* argv[]) {
    // Parse inputs
    unsigned ll = std::stoi(argv[1]);
    double tt = std::stof(argv[2]);
    double ww = std::stof(argv[3]);
    double dt = std::stof(argv[4]);
    unsigned timesteps = std::stoi(argv[5]);
    unsigned NCheb = std::stoi(argv[6]);
    unsigned tnum = std::stoi(argv[7]);
    unsigned save = std::stoi(argv[8]);
    unsigned time = std::stoi(argv[9]);

    // When a parallel block appears, use tnum threads
    omp_set_num_threads(tnum);
    
    // Start counting time
    auto start = std::chrono::high_resolution_clock::now();

    // Make NCheb even
    NCheb += NCheb % 2;

    // Define imaginary constant
    std::complex<double> I(0, 1);

    // Precompute powers of I
    Eigen::Array<std::complex<double>, -1, 1> powI(4, 1);
    powI[0] = 1;
    powI[1] = I;
    powI[2] = -1;
    powI[3] = -I;

    // Ensure ll is divisible by tnum
    ll -= ll % tnum;

    // Set up final array with all psi
    Eigen::Array<std::complex<double>, -1, -1> psiTime(ll, timesteps);

    // Set up the initial state
    Eigen::Array<std::complex<double>, -1, 1> psi0(ll, 1);
    double n0 = ll * n0factor;
    double sig = ll * sigfactor;

    for (unsigned n = 0; n < ll; n++) {
        psi0(n, 0) = std::exp(-1.0 * (n - n0) * (n - n0) / (2.0 * sig * sig) - I * k0 * (1.0 * n));
    }

    // Pre-compute Bessel functions
    Eigen::Array<double, -1, -1> Bessel(NCheb, timesteps);
    Bessel.setZero();

    #pragma omp parallel shared(Bessel)
    {   
        Eigen::Array<double, -1, -1> BesselThread(NCheb, timesteps);

        #pragma omp for schedule(static, 5) nowait
            for (unsigned i = 0; i < NCheb; i++) {
                for (unsigned j = 0; j < timesteps; j++) {
                    BesselThread(i, j) = jn(i, j * dt);
                }
            }

        #pragma omp critical
        {
            Bessel += BesselThread;
        }
    }


    // Shared array of ghost points
    Eigen::Array<std::complex<double>, -1, 1> ghosts(2*tnum, 1);
    ghosts.setZero();

    // Each thread will have their own instance of the kpm class
    // TODO: Figure out if there is a way to not copy in "Bessel" because it is a huge matrix
    #pragma omp parallel firstprivate(ll, tt, ww, NCheb, psi0, Bessel, I) shared(ghosts)
    {
        // Get thread number
        unsigned posThread = omp_get_thread_num();


        // Divide the system size by the number of threads
        unsigned llThread = ll / tnum;

        // Create thread-specific vector to store psi
        Eigen::Array<std::complex<double>, -1, -1> psiTimeThread(llThread, timesteps);
        psiTimeThread.setZero();

        // Initialize the class
        kpm kpmobj(llThread, tt, ww);

        // Compute T_n(H)|psi>

        // Set up the initial state with the correct slice from the vector psi0
        kpmobj.initialize(psi0.segment(posThread * llThread, llThread));

        // Push ghost points (we store left0, right0, left1, right1, ...)
        ghosts[2 * posThread] = kpmobj.left();
        ghosts[2 * posThread + 1] = kpmobj.right();
        #pragma omp barrier

        // Get ghost points (right of thread-1 and left of thread+1)
        kpmobj.pull_ghosts(ghosts[2 * ((posThread - 1 + tnum) % tnum) + 1], ghosts[2 * ((posThread + 1) % tnum)]);
        #pragma omp barrier
        
        // Do first iteration, push, sync and pull
        kpmobj.hamiltonian_kpm(1);
        ghosts[2 * posThread] = kpmobj.left();
        ghosts[2 * posThread + 1] = kpmobj.right();
        #pragma omp barrier

        kpmobj.pull_ghosts(ghosts[2 * ((posThread - 1 + tnum) % tnum) + 1], ghosts[2 * ((posThread + 1) % tnum)]);
        #pragma omp barrier

        // Compute time evolution using first two Chebychev polynomials
        for (unsigned i = 0; i < llThread; i++) {
            for (unsigned j = 0; j < timesteps; j++) {
                psiTimeThread(i, j) += Bessel(0, j) * kpmobj.KetPsi(i, 0);
                psiTimeThread(i, j) += 2. * (-I) * Bessel(1, j) * kpmobj.KetPsi(i, 1);
            }
        }

        // Do KPM loop
        for(unsigned m = 2; m < NCheb; m += 2) {
            for(unsigned j = 0; j < 2; j++) {
                kpmobj.hamiltonian_kpm(2);   
                
                ghosts[2 * posThread] = kpmobj.left();
                ghosts[2 * posThread + 1] = kpmobj.right();
                #pragma omp barrier

                kpmobj.pull_ghosts(ghosts[2 * ((posThread - 1 + tnum) % tnum) + 1], ghosts[2 * ((posThread + 1) % tnum)]);
                #pragma omp barrier
            }

            // Compute time evolution using two new Chebychev polynomials
            for (unsigned i = 0; i < llThread; i++) {
                for (unsigned j = 0; j < timesteps; j++) {
                    psiTimeThread(i, j) += 2. * pow(-1, m) * powI[m % 4] * Bessel(m, j) * kpmobj.KetPsi(i, 0);
                    psiTimeThread(i, j) += 2. * pow(-1, m+1) * powI[(m+1) % 4] * Bessel(m+1, j) * kpmobj.KetPsi(i, 1);
                }
            }
        }

        // Write to the general psi
        #pragma omp critical
        {   
            psiTime.block(posThread * llThread, 0, llThread, timesteps) = psiTimeThread;
        }
    }

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
        std::ofstream outfile(name + "-" + std::to_string(ll) + "-" + std::to_string(tt) + "-" + std::to_string(ww) + "-" + std::to_string(NCheb) + ".txt");
        for (unsigned i = 0; i < ll; i++) {
            for (unsigned j = 0; j < timesteps; j++) {
                outfile << abs(psiTime(i, j)) << " ";
            }
            outfile << std::endl;
        }
        outfile.close(); 
    }
}


// This program is very fast, hence we can go to very large system sizes
// In the limit of a very large system the majority of the time might be spent moving memory to the CPU cache