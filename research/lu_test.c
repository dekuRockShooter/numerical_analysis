#include <stdlib.h>
#include <stdio.h>
#include "lu.h"

//void print_vec(double* a, int n) {
    //for (int k = 0; k < n; ++k) {
        //printf("%f ", a[k]);
    //}
//}
//
//void print_mat(double** a, int n) {
    //for (int j = 0; j < n; ++j) {
        //print_vec(a[j], n);
        //printf("\n");
    //}
//}
void print_vec(double* a, int n);
void print_mat(double** a, int n);

int main() {
    LU* lu = (LU*) malloc(sizeof(LU));
    int n = 3;
    //double** a = (double**) calloc(n, sizeof(double*));
    //for (int j = 0; j < n; ++j) {
        //a[j] = (double*) calloc(n, sizeof(double));
        //for (int k = 0; k < n; ++k) {
            //a[j][k] = j + k;
        //}
    //}
    double** a = (double**) calloc(n, sizeof(double*));
    double* x = (double*) calloc(n, sizeof(double));
    double* b = (double*) calloc(n, sizeof(double));
    a[0] = (double*) calloc(n, sizeof(double));
    a[1] = (double*) calloc(n, sizeof(double));
    a[2] = (double*) calloc(n, sizeof(double));
    a[0][0] = 3.0; a[0][1] = 1.0; a[0][2] = 0.0;
    a[1][0] = -1.0; a[1][1] = 2.0; a[1][2] = 2.0;
    a[2][0] = 5.0; a[2][1] = 0.0; a[2][2] = -1.0;
    b[0] = 6.0; b[1] = -7.0; b[2] = 10.0;
    x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
    //print_mat(a, n);
    lu_factor(lu, a, n);
    print_mat(a, n);
    printf("\n");
    lu_solve(lu, x, b);
    print_vec(x, n);
    printf("\n");
    return 0;
}
