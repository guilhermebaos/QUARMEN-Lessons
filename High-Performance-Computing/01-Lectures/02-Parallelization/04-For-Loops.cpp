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

        #pragma omp for schedule(static, 2)
            for(unsigned i = 0; i < 11; i++) {
                std::cout << "tid: " << tid << "   index: " << i << std::endl;
                
                std::this_thread::sleep_for(std::chrono::milliseconds(10 * (tid + 1)));
                count[tid]++;
        }
    }
    
    std::cout << std::endl;
    std::cout << "Total Workload: ";

    for(unsigned i = 0; i < 4; i++) {
        std::cout << count[i] << " ";
    }

    std::cout << std::endl;
}


// If we use schedule(static, 20) then thread 0 gets all the indexes

// If we use schedule(static, 2) then the indices are divided in chunks of 2 and assigned to each thread:
// Thread 0 gets 0, 1
// Thread 1 gets 2, 3
// Thread 2 gets 4, 5
// Thread 3 gets 6, 7
// Thread 0 gets 8, 9
// Thread 1 gets 10



// If we use schedule(dynamic, 2) we have that:
// First we divide the indeces into chunks of a given size.
// Secondly we assign the first index to the first thread that arrives (that is, that finishes the previous job).





