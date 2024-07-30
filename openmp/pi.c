/*
   Using openmp to compute pi using integral of 4/(1+x*x) from 0 to 1.
   */

#include <omp.h>
#include <stdio.h>
static long num_steps = 10000000;
/* static long num_steps = 4000; */

double pi_serial(double step) {
    double sum = 0;
    for (int i = 0; i < num_steps; ++i) {
        double x = (i - 0.5) * step;
        sum += 4.0 / (1.0 + x * x);
    }
    return sum;
}

double pi_omp_critical(double step) {
    double sum = 0;
    #pragma omp parallel
    {
        int nt = omp_get_num_threads();
        int tid = omp_get_thread_num();
        double tsum = 0;
        for (int i = tid; i < num_steps; i += nt) {
            double x = (i - 0.5) * step;
            tsum += 4.0 / (1.0 + x * x);
        }
        #pragma omp critical
        sum += tsum;
    }
    return sum;
}

double pi_modern_omp(double step) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_steps; ++i) {
        double x = (i - 0.5) * step;
        sum += 4.0 / (1.0 + x * x);
    }
    return sum;
}

double pi_recursive(double step, int start, int end){
    double sum1, sum2;
    if (end - start <= 1000000) {
        double sum = 0;
        for( int i = start; i < end; ++i) {
            double x = (i - 0.5) * step;
            sum += 4.0 / (1.0 + x * x);
        }
        return sum;
    } else {
        #pragma omp task shared(sum1)
        sum1 = pi_recursive(step, start, start + (int)(end - start)/2);
        #pragma omp task shared(sum2)
        sum2 = pi_recursive(step, start + (int)(end-start)/2, end);
        #pragma omp taskwait
        return (sum1 + sum2);
    }
}

int main() {
    double x, sum = 0.0, step = 1.0 / (double)num_steps;
    double start_time, run_time_serial, run_time_critical, run_time_modern;
    double run_time;
    double pi;

    start_time = omp_get_wtime();
    sum = pi_serial(step);
    pi = sum * step;
    run_time = omp_get_wtime() - start_time;
    printf("\t serial: \t\t pi: %f, time: %f \n", pi, run_time);

    start_time = omp_get_wtime();
    sum = pi_omp_critical(step);
    pi = sum * step;
    run_time = omp_get_wtime() - start_time;
    printf("\t omp critical: \t\t pi: %f, time: %f \n", pi, run_time);

    start_time = omp_get_wtime();
    sum = pi_modern_omp(step);
    pi = sum * step;
    run_time = omp_get_wtime() - start_time;
    printf("\t omp reduce: \t\t pi: %f, time: %f \n", pi, run_time);

    start_time = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        {
            sum = pi_recursive(step, 0, num_steps);
        }
    }
    pi = sum*step;
    run_time = omp_get_wtime() - start_time;
    printf("\t omp task: \t\t pi: %f, time: %f \n", pi, run_time);
}
