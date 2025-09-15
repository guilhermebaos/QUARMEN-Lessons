#include <iostream>
#include <random>
#include <omp.h>

#define M 100000000


int main() {
    // When a parallel block appears, use four threads
    unsigned tnum = 4;
    omp_set_num_threads(tnum);

    // Create shared array
    unsigned count = 0;

    // We use a private variable for each thread to prevent false sharing
    #pragma omp parallel reduction(+:count)
    {   
        // Make sure it's 0
        count = 0;

        // Use the random device as a seed for the pseudo-random number generator (we use 8 numbers because the prof read somewhere that it maximizes the entropy of the generated numbers)
        std::random_device r;
        std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
        std::mt19937 e2(seed2);

        // Use the pseudo-random number generator to get a uniform distribution 
        std::uniform_real_distribution dist;

        // The chunk size should be equal to the size of the cache
        // This way whenever the thread accesses the memory it gets the highest amount of useful data 
        // This is very hardware-dependent so it's a last-minute worry!
        #pragma omp for schedule(static, 5)

            // Use a Monte-Carlo method to get PI (we generate random points inside a square of size 1 in the first quadrant)
            for(int i = 0; i < M; i++) {
            double x = dist(e2);
            double y = dist(e2);
            
            // Determine how many points are inside the quarter circle in the first quadrant
            count += (x * x + y * y < 1 ? 1 : 0);
        }

    }

    // PI is the fraction of points inside the quarter circle multiplied by 4
    std::cout << 4 * (count * 1.0 / M) << std::endl;

    return 0;
}