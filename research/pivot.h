#ifndef PIVOT_H
#define PIVOT_H


void partial_pivot(double** a, int piv_row, int n, int* perm);
void scaled_partial_pivot(double** a, int piv_row, int n, int* perm,
        double* scale_factors);

#endif
