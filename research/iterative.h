#ifndef ITERATIVE_H
#define ITERATIVE_H


typedef struct convergence_stats {
    int iterations;
    int converged;
    double max_residual;
    double tolerance;
} convergence_stats;


void gauss_seidel(
        double** a,
        double* x_approx,
        double* b,
        int n,
        int max_iter,
        double tol,
        convergence_stats* c_stats);

void jacobian(
        double** a,
        double** x_approx,
        double* b,
        int n,
        int max_iter,
        double tol,
        convergence_stats* c_stats);

void iter_refine(
        double** a,
        double* x_approx,
        double* b,
        int n,
        int max_iter,
        double tol,
        convergence_stats* c_stats);

#endif
