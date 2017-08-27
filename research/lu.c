#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lu.h"
#include "pivot.h"


void print_vec(double* a, int n) {
    for (int k = 0; k < n; ++k) {
        printf("%f ", a[k]);
    }
}

void print_veci(int* a, int n) {
    for (int k = 0; k < n; ++k) {
        printf("%d ", a[k]);
    }
}

void print_mat(double** a, int n) {
    for (int j = 0; j < n; ++j) {
        print_vec(a[j], n);
        printf("\n");
    }
}

void swap_elem(double* a, int j, int k) {
    double tmp = a[j];
    a[j] = a[k];
    a[k] = tmp;
}

void swap_elemi(int* a, int j, int k) {
    int tmp = a[j];
    a[j] = a[k];
    a[k] = tmp;
}

void swap_rows(double** a, int j, int k) {
    double* tmp = a[j];
    a[j] = a[k];
    a[k] = tmp;
}

void lu_init(LU* lu, double** a, int n) {
    lu->n = n;
    lu->lu = a;
}

void backward_sub(double** a, double* x, double* b, int n) {
    int k = n - 1;
    int j = 0;
    double prod = 0.0;
    while (k > -1) {
        j = k + 1;
        prod = 0.0;
        while (j < n) {
            prod = prod + a[k][j]*x[j];
            ++j;
        }
        x[k] = (b[k] - prod)/a[k][k];
        --k;
    }
}

void forward_sub(double** a, double* x, double* b, int n) {
    int k = 0;
    int j = 0;
    double prod = 0.0;
    while (k < n) {
        j = 0;
        prod = 0.0;
        while (j < k) {
            prod = prod + a[k][j]*x[j];
            ++j;
        }
        x[k] = (b[k] - prod)/a[k][k];
        ++k;
    }
}

void lu_factor(LU* lu, double** a, int n) {
    // Step 1.  Initialize scale factors and permutation vector.
    int row_idx = 0;
    int col_idx = 0;
    double max = 0.0;
    int* perm_vec = calloc(n, sizeof(int));
    double* scale_factors = calloc(n, sizeof(double));
    double cur_elem = 0.0;
    while (row_idx < n) {
        col_idx = 0;
        max = 0.0;
        while (col_idx < n) {
            cur_elem = fabs(a[row_idx][col_idx]);
            if (cur_elem > max) {
                max = cur_elem;
            }
            ++col_idx;
        }
        perm_vec[row_idx] = row_idx;
        scale_factors[row_idx] = max;
        ++row_idx;
    }
    row_idx = 0;
    int j = 0; // Row index.
    int k = 0; // Column index.
    int i = 0; // Index used for the dot product.
    int p = 0; // Pivot row index.
    double dot_prod = 0.0;
    while (j < n) {
        // Step 2.  Find and swap the pivot row.
        //scaled_partial_pivot(a, j; n, perm_vec, scale_factors);
        max = fabs(a[j][j]) / scale_factors[j];
        cur_elem = 0.0;
        p = j;
        i = j + 1;
        while (i < n) {
            cur_elem = fabs(a[i][j]) / scale_factors[i];
            if (cur_elem > max) {
                max = cur_elem;
                p = i;
            }
            ++i;
        }
        if (p != j) {
            swap_rows(a, p, j);
            swap_elemi(perm_vec, p, j);
            swap_elem(scale_factors, p, j);
        }
        // Step 3.  Calculate the j'th row of L.
        k = 0;
        while (k < j) {
            // Partial dot product.
            dot_prod = 0.0;
            i = 0;
            while (i < k) {
                dot_prod = dot_prod + a[j][i]*a[i][k];
                ++i;
            }
            a[j][k] = (a[j][k] - dot_prod) / a[k][k];
            ++k;
        }
        // Step 4.  Calculate the j'th row of U.
        while (k < n) {
            // Partial dot product.
            dot_prod = 0.0;
            i = 0;
            while (i < j) {
                dot_prod = dot_prod + a[j][i]*a[i][k];
                ++i;
            }
            a[j][k] = a[j][k] - dot_prod;
            ++k;
        }
        ++j;
    }
    // Step 5.  Find the transpose of the permutation vector.
    int* perm_vec_transpose = (int*) calloc(n, sizeof(int));
    j = 0;
    while (j < n) {
        perm_vec_transpose[perm_vec[j]] = j;
        ++j;
    }
    lu->lu = a;
    lu->n = n;
    lu->perm_vec_transpose = perm_vec_transpose;
    lu->perm_vec= perm_vec;
}

void lu_solve(LU* lu_struct, double* x, double* b) {
    int n = lu_struct->n;
    double b_tpose[n];
    double lu_diag[n];
    double** lu = lu_struct->lu;
    int* perm_tpose = lu_struct->perm_vec_transpose;
    int* perm = lu_struct->perm_vec;
    int j = 0;
    while (j < n) {
        //b_tpose[j] = b[perm_tpose[j]];
        b_tpose[j] = b[perm[j]];
        lu_diag[j] = lu[j][j];
        lu[j][j] = 1.0;
        ++j;
    }
    forward_sub(lu, x, b_tpose, n);
    j = 0;
    while (j < n) {
        lu[j][j] = lu_diag[j];
        ++j;
    }
    backward_sub(lu, x, x, n);
}
