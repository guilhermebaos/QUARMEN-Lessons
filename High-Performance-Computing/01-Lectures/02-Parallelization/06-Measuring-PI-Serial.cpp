#include <iostream>
#include <chrono>
#include <random>
#include <omp.h>

#define M 10000000


int main() {
    // Use the random device as a seed for the pseudo-random number generator (we use 8 numbers because the prof read somewhere that it maximizes the entropy of the generated numbers)
    std::random_device r;
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 e2(seed2);

    // Use the pseudo-random number generator to get a uniform distribution 
    std::uniform_real_distribution dist;

    // Use a Monte-Carlo method to get PI (we generate random points inside a square of size 1 in the first quadrant)
    unsigned count = 0;
    for(int i = 0; i < M; i++) {
        double x = dist(e2);
        double y = dist(e2);
        
        // Determine how many points are inside the quarter circle in the first quadrant
        count += (x * x + y * y < 1 ? 1 : 0);
    }

    // PI is the fraction of points inside the quarter circle multiplied by 4
    std::cout << 4 * (count * 1.0 / M) << std::endl;

    return 0;
}