#ifndef INTERPOLATION_1D_H
#define INTERPOLATION_1D_H

double lagrange(double x, double* X, double* Y, int n);
void neville(double x, double* X, double* Y, int n, double** grid);
double cubic_spline(double x, double* knots, double* Y, int n);
double linear_spline(double x, double* knots, double* Y, int n);

static int bsearch(double* a, int beg, int end, int n);

#endif
