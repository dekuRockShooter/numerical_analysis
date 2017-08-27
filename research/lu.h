#ifndef LU_H
#define LU_H

typedef struct LU {
    double** lu;
    int n;
    int* perm_vec_transpose;
    int* perm_vec;
} LU;

void lu_factor(LU* lu, double** a, int n);
void lu_solve(LU* lu, double* x, double* b);

#endif
