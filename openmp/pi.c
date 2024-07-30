/*
   Using openmp to compute pi using integral of 4/(1+x*x) from 0 to 1.
*/

#include <omp.h>
#include <stdio.h>
static long num_steps = 10000000;

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

int main() {
  double x, sum = 0.0, step = 1.0 / (double)num_steps;
  double start_time, run_time_serial, run_time_critical, run_time_modern;
  start_time = omp_get_wtime();
  sum = pi_serial(step);
  double pi_s = sum * step;
  run_time_serial = omp_get_wtime() - start_time;

  start_time = omp_get_wtime();
  sum = pi_omp_critical(step);
  double pi_c = sum * step;
  run_time_critical = omp_get_wtime() - start_time;

  start_time = omp_get_wtime();
  sum = pi_modern_omp(step);
  double pi_m = sum * step;
  run_time_modern = omp_get_wtime() - start_time;
  printf("\tserial:\t\t pi: %f, time: %f \n\tomp critical:\t pi: %f, time: %f "
         "\n\tomp reduce:\t pi: %f, time: %f \n",
         pi_s, run_time_serial, pi_c, run_time_critical, pi_m, run_time_modern);
}
