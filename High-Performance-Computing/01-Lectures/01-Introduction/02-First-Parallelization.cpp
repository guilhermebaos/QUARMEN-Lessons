#include <iostream>
#include <omp.h>

int main () {
    std::cout << "A" << std::endl;

    // Execute the code in all threads at the same time
    #pragma omp parallel {
        std::cout << "B" << std::endl;
    }


    // Execute the code in one thread at a time (no parallelization)
    #pragma omp critical {
        int thread = omp_get_thread_num()
        std::cout << "B " << thread << std::endl;

        if (thread < 8) {
            std::cout << "THREAD < 8" << std::endl;
            std::cout << std::endl;
        }
    }

    return 0;
}