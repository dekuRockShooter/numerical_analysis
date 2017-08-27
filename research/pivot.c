#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pivot.h"


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

void partial_pivot(double** a, int piv_row, int n, int* perm) {
    double max = abs(a[piv_row][piv_row]);
    double cur_elem = 0.0;
    int piv_idx = piv_row;
    // TODO: compare with while when inlining.
    for (int j = piv_row + 1; j < n; ++j) {
        cur_elem = abs(a[j][piv_row]); // TODO: compare with handmade abs.
        if (cur_elem > max) {
            max = cur_elem;
            piv_idx = j;
        }
    }
    if (piv_idx != piv_row) {
        swap_rows(a, piv_idx, piv_row);
        swap_elemi(perm, piv_idx, piv_row);
    }
}

void scaled_partial_pivot(
        double** a,
        int piv_row,
        int n,
        int* perm_vec,
        double* scale_factors) {
    double max = abs(a[piv_row][piv_row]) / scale_factors[piv_row];
    double cur_elem = 0.0;
    int piv_idx = piv_row;
    for (int j = piv_row + 1; j < n; ++j) {
        cur_elem = abs(a[j][piv_row]) / scale_factors[j];
        if (cur_elem > max) {
            max = cur_elem;
            piv_idx = j;
        }
    }
    if (piv_idx != piv_row) {
        swap_rows(a, piv_idx, piv_row);
        swap_elemi(perm_vec, piv_idx, piv_row);
        swap_elem(scale_factors, piv_idx, piv_row);
    }
}
