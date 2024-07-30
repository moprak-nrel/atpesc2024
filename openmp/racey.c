#include <omp.h>
#include <stdio.h>

int main() {
    printf("I think");
    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp task
            printf(" car");
            #pragma omp task
            printf(" race");
        }
    }
    printf("s are fun\n");
}
