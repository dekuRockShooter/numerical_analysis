#include "interpolation_1d.h"

double lagrange(double x, double* X, double* Y, int n) {
    double prod = 1.0;
    double sum = 0.0;
    int j = 0;
    int k = 0;
    while (j < n) {
        k = 0;
        prod = 1.0;
        // L_j[x]
        while (k < n) {
            if (k == j) {
                ++k;
            }
            prod = prod*((x - X[k])/(x[j] - x[k]));
            ++k;
        }
        sum = sum + Y[j]*prod;
        ++j;
    }
    return sum;
}

// 'grid' must be an n by n matrix.  The entries below the diagonal are not
// used.
void neville(double x, double* X, double* Y, int n, double** grid) {
    int k = 0;
    int j = 0;
    while (k < n) {
        grid[k][k] = Y[k];
        j = k - 1;
        while (j > -1) {
            grid[j][k] = 
                (grid[j][k-1]*(x - X[k]) + grid[j+1][k]*(X[j] - x))
                / (X[j] - X[k]);
            --j;
        }
        ++k;
    }
}

double cubic_spline(double x, double* X, double* Y, int n) {
    double a[n][n];
    int j = 1;
    int m = n - 1;
    // Set a[0][k] according to natural, clamped, or other.
    while (j < m) {
        a[j][j-1] = x[j];
        a[j][j] = x[j] + x[j+1];
        a[j][j+1] = x[j+1];
        ++j
    }
    // Set a[m][k] according to natural, clamped, or other.
    //lu_factor_tri_diag(lu, a, n);
    //lu_solve_tri_diag(lu, s, b, n);
    //s*x;
}

double linear_spline(double x, double* X, double* Y, int n) {
    int j = bsearch(X, x, 0, n - 1, n);
    return (Y[j]*(X[j+1] - x) + Y[j+1]*(x - X[j]))/(X[j+1] - X[j]);
}

// Return the index j whose element in a is the largest number 
// smaller than 'x'.  This means a[j] <= x <= a[j + 1].
static int bsearch(double* a, int x, int beg, int end, int n) {
    if (!(end < n)) {
        return -1;
    }
    int mid = 0;
    int less_than = 0;
    while (beg < end) {
        mid = (beg + end) / 2;
        if (a[mid] < x) {
            beg = mid + 1;
        }
        else if (a[mid] > x) {
            end = mid - 1;
        }
        else {
            break;
        }
    }
    return (a[beg] > x) ? beg - 1 : beg;
}
