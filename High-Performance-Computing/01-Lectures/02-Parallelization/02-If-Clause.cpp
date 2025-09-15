#include <omp.h>
#include <stdio.h>
#include <iostream>

int main(int argc, char *argv[]) {
    unsigned M = atoi(argv[1]);
    double x = atof(argv[2]);

    omp_set_num_threads(11);

    #pragma omp parallel if(M < 10)
    {
        #pragma omp critical
        {
            printf("M = %d \n", M);
            printf("x = %f \n", x);
            printf("Number of threads: %d \n", omp_get_num_threads());
            printf("Current thread: %d \n\n", omp_get_thread_num());
        }
    }
}