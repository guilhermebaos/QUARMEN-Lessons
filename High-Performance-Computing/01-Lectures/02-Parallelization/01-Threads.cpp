#include <omp.h>
#include <stdio.h>

int main() {
    int a = 0, tid = 1;
    float x = 1.0;

    printf("1st Parallel Region:\n\n");

    #pragma omp parallel shared(x, a) private(tid) num_threads(2)
    {
        printf("X1: a = %d, tid = %d, x = %f \n\n", a, tid, x);
        tid = omp_get_thread_num();

        printf("X1: a = %d, tid = %d, x = %f \n\n", a, tid, x);
        a = tid;

        printf("X1: a = %d, tid = %d, x = %f \n\n", a, tid, x);
    }

    printf("X1: a = %d, tid = %d, x = %f \n\n", a, tid, x);

    printf("2nd Parallel Region:\n\n");

    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
    }

    return 0;
}