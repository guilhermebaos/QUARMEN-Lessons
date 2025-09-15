#include <omp.h>
#include <cmath>
#include <unistd.h>
#include <iostream>


int main() {
    double sum = 0;

    omp_set_num_threads(4);

    std::cout << "Global variable: " << sum << std::endl;


    #pragma omp parallel reduction(+:sum)
    {
        sum = omp_get_thread_num();

        #pragma omp critical 
        {
            std::cout << "Local variable: " << sum << std::endl;

            sleep(omp_get_thread_num() + 1);
        }
    }

    std::cout << "Global variable: " << sum << std::endl;
}


// In this program we have a global variable sum
// When we do reduction(+:sum) we define a local variable sum
// such that at the end of the parallel block the value of the local variable is added to the globaL one

// In this case at the end we get 0 + 1 + 2 + 3 = 6 (the sum over thread numbers)
