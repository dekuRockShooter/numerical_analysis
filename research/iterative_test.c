#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "iterative.h"
#include "lu.h"

//void print_vec(double* a, int n) {
    //for (int k = 0; k < n; ++k) {
        //printf("%f ", a[k]);
    //}
//}
void print_vec(double* a, int n);
void mat_vec_mul(double** a, double* b, double* c, int n);

void write_csv(int method, convergence_stats* cstats, FILE* f) {
    // method, converged, iterations, infinity norm, tolerance.
    fprintf(f,
            "%d,%d,%d,%f,%f\n",
            method,
            cstats->converged,
            cstats->iterations,
            cstats->max_residual,
            cstats->tolerance);
}

void init_matrices(double** a, double* x, double* b, int n) {
    // The random matrix test.
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            a[j][k] = 1 + 10.0 * ((double) rand()) / RAND_MAX;
        }
        b[j] = 1 + 20.0 * ((double) rand()) / RAND_MAX;
        x[j] = 0.0;
    }
    // The diagonally dominant matrix test.
    //for (int j = 0; j < n; ++j) {
        //for (int k = 0; k < n; ++k) {
            //a[j][k] = 10.0 * ((double) rand()) / RAND_MAX;
        //}
        //a[j][j] = 10.0*n + 10.0*((double) rand()) / RAND_MAX;
        //b[j] = 100.0*((double) rand()) / RAND_MAX;
        //x[j] = 0.0;
    //}
}

void init_vec(double* x, int n) {
    for (int j = 0; j < n; ++j) {
        x[j] = 0.0;
    }
}

int main(int argc, char** argv) {
    int n = 1000;
    int max_iter = (argc == 2) ? (int) strtod(argv[1], NULL) : 10;
    double tol = 1e-10;
    FILE* fout = fopen("./iter_refine_sub.csv", "a");
    convergence_stats* cstats =
        (convergence_stats*) malloc(sizeof(convergence_stats));
    int j = 0;
    srand(time(NULL));
    double** a = (double**) calloc(n, sizeof(double*));
    double* x = (double*) calloc(n, sizeof(double));
    double* b = (double*) calloc(n, sizeof(double));
    for (int j = 0; j < n; ++j) {
        a[j] = (double*) calloc(n, sizeof(double));
    }
    //fprintf(fout, "method,converged,iterations,infinity_norm,tolerance");
    LU lu;
    int converged_count = 0;
    while (j < 1) {
        init_matrices(a, x, b, n);
        //double** a_copy = (double**) calloc(n, sizeof(double*));
        //for (int j = 0; j < n; ++j) {
            //a_copy[j] = (double*) calloc(n, sizeof(double));
            //for (int k = 0; k < n; ++k) {
                //a_copy[j][k] = a[j][k];
            //}
        //}

        //init_vec(x, n);
        //gauss_seidel(a, x, b, n, max_iter, tol, cstats);
        //write_csv(1, cstats, fout);
//
        //init_vec(x, n);
        //jacobian(a, &x, b, n, max_iter, tol, cstats);
        //write_csv(2, cstats, fout);

        //init_vec(x, n);
        //iter_refine(a, x, b, n, max_iter, tol, cstats);
        //write_csv(3, cstats, fout);
        //if (cstats->converged == 1) {
            //++converged_count;
        //}
        lu_factor(&lu, a, n);
        lu_solve(&lu, x, b);
        ++j;
    }
    printf("%d\n", converged_count);
    fclose(fout);
    return 0;
}
