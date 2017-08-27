#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "iterative.h"
#include "lu.h"




//void print_vec(double* a, int n) {
    //for (int k = 0; k < n; ++k) {
        //printf("%f ", a[k]);
    //}
//}
void print_vec(double* a, int n);


void mat_mul(double** a, double** b, double** c, int n, int m) {
    int j = 0;
    int k = 0;
    int i = 0;
    double dot_prod = 0.0;
    while (j < n) {
        k = 0;
        while (k < m) {
            i = 0;
            dot_prod = 0.0;
            while (i < n) {
                dot_prod = dot_prod + a[j][i]*b[i][k];
                ++i;
            }
            c[j][k] = dot_prod;
            ++k;
        }
        ++j;
    }
}
void mat_vec_mul(double** a, double* b, double* c, int n) {
    int j = 0;
    int k = 0;
    int i = 0;
    int m = 1;
    double dot_prod = 0.0;
    while (j < n) {
        k = 0;
        while (k < m) {
            i = 0;
            dot_prod = 0.0;
            while (i < n) {
                dot_prod = dot_prod + a[j][i]*b[i];
                ++i;
            }
            c[j] = dot_prod;
            ++k;
        }
        ++j;
    }
}

void residual(double** a, double* x, double* b, double* r, int n) {
    int j = 0;
    int k = 0;
    int i = 0;
    int m = 1;
    double dot_prod = 0.0;
    while (j < n) {
        k = 0;
        while (k < m) {
            i = 0;
            dot_prod = 0.0;
            while (i < n) {
                dot_prod = dot_prod + a[j][i]*x[i];
                ++i;
            }
            r[j] = dot_prod - b[j];
            ++k;
        }
        ++j;
    }
}

double inf_norm(double* a, int n) {
    double max = fabs(a[0]);
    double elem = 0.0;
    int j = 1;
    while (j < n) {
        elem = fabs(a[j]);
        if (elem > max) {
            max = elem;
        }
        ++j;
    }
    return max;
}

void gauss_seidel(
        double** a,
        double* x_approx,
        double* b,
        int n,
        int max_iter,
        double tol,
        convergence_stats* c_stats) {
    int k = 0;
    int j = 0;
    int count = 0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double* r = (double*) calloc(n, sizeof(double));
    double inf_norm_val = 0.0;
    while (count < max_iter) {
        residual(a, x_approx, b, r, n);
        inf_norm_val = inf_norm(r, n);
        if (inf_norm_val < tol) {
            break;
        }
        k = 0;
        while (k < n) {
            sum1 = 0.0;
            j = 0;
            while (j < k) {
                sum1 = sum1 + a[k][j]*x_approx[j];
                ++j;
            }
            sum2 = 0.0;
            j = k + 1;
            while (j < n) {
                sum2 = sum2 + a[k][j]*x_approx[j];
                ++j;
            }
            x_approx[k] = (b[k] - sum1 - sum2)/a[k][k];
            ++k;
        }
        ++count;
    }
    if (c_stats != NULL) {
        c_stats->iterations = count;
        c_stats->converged = count < max_iter;
        c_stats->max_residual = inf_norm_val;
        c_stats->tolerance = tol;
    }
}

void jacobian(
        double** a,
        double** x_approx,
        double* b,
        int n,
        int max_iter,
        double tol,
        convergence_stats* c_stats) {
    int k = 0;
    int j = 0;
    int count = 0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double* x_new = (double*) calloc(n, sizeof(double));
    double* x_old = *x_approx;
    double r[n];
    double inf_norm_val = 0.0;
    while (count < max_iter) {
        residual(a, x_old, b, r, n);
        inf_norm_val = inf_norm(r, n);
        if (inf_norm_val < tol) {
            break;
        }
        k = 0;
        while (k < n) {
            sum1 = 0.0;
            j = 0;
            while (j < k) {
                sum1 = sum1 + a[k][j]*x_old[j];
                ++j;
            }
            sum2 = 0.0;
            j = k + 1;
            while (j < n) {
                sum2 = sum2 + a[k][j]*x_old[j];
                ++j;
            }
            x_new[k] = (b[k] - sum1 - sum2)/a[k][k];
            ++k;
        }
        x_old = x_new;
        ++count;
    }
    *x_approx = x_old;
    if (c_stats != NULL) {
        c_stats->iterations = count;
        c_stats->converged = count < max_iter;
        c_stats->max_residual = inf_norm_val;
        c_stats->tolerance = tol;
    }
}

void iter_refine(
        double** a,
        double* x_approx,
        double* b,
        int n,
        int max_iter,
        double tol,
        convergence_stats* c_stats) {
    double** a_lu = (double**) calloc(n, sizeof(double*));
    for (int j = 0; j < n; ++j) {
        a_lu[j] = (double*) calloc(n, sizeof(double));
        for (int k = 0; k < n; ++k) {
            a_lu[j][k] = a[j][k];
        }
    }
    LU lu;
    int count = 0;
    lu_factor(&lu, a_lu, n);
    double* dx = calloc(n, sizeof(double));
    double* r = calloc(n, sizeof(double));
    double inf_norm_val_prev = 0.0;
    double inf_norm_val_cur = 0.0;
    int j = 0;
    while (count < max_iter) {
        residual(a, x_approx, b, r, n);
        inf_norm_val_cur = inf_norm(r, n);
        if (inf_norm_val_cur < tol) {
            break;
        }
        lu_solve(&lu, dx, r);
        j = 0;
        while (j < n) {
            x_approx[j] = x_approx[j] - dx[j];
            ++j;
        }
        inf_norm_val_prev = inf_norm_val_cur;
        ++count;
    }
    if (c_stats != NULL) {
        c_stats->iterations = count;
        c_stats->converged = count < max_iter;
        c_stats->max_residual = inf_norm_val_cur;
        c_stats->tolerance = tol;
    }
}
