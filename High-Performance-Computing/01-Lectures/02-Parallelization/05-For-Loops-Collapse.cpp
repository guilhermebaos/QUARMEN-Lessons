#include <omp.h>
#include <cmath>
#include <chrono>
#include <thread>
#include <iostream>


int main() {
    // We create an array of four unsigned variables, all of them are 0
    unsigned count[4] = {0};

    // When a parallel block appears, use four threads
    omp_set_num_threads(4);

    #pragma omp parallel shared(count)
    {
        unsigned tid = omp_get_thread_num();

        #pragma omp for schedule(static, 20) collapse(2)
            for(unsigned i = 0; i < 11; i++) {
                for(unsigned i = 0; i < 11; i++) {
                    for(unsigned i = 0; i < 11; i++) {
                        std::this_thread::sleep_for(std::chrono::milliseconds(10 * (tid + 1)));
                        count[tid]++;
                    }
                }
            }
    }
    
    std::cout << std::endl;
    std::cout << "Total Workload: ";

    for(unsigned i = 0; i < 4; i++) {
        std::cout << count[i] << " ";
    }

    std::cout << std::endl;
}


// If we run with schedule(static, 20) then we can see the effect of collapsing loops

// With collapse(3) we collapse all loops, and we are able to distribute the workload evenly
// Total Workload: 340 340 331 320

// With collapse(2) we collapse two loops, meaning that the threads with higher number do less work due to the order of the assignment
// Total Workload: 440 440 231 220 