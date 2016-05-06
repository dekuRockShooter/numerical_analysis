"""
gauss_elim_naive_2.py

Includes code for functions that do basic vector and
matrix arithmetic.  Most of these functions support
the function ge_1(aug) which takes an n by n+1
augmented matrix and returns a row-reduced matrix and
an approximate solution vector of the corresponding linear
system.  It uses gaussian elimination with a naive pivot
strategy.  That is, at each stage in the row reduction it
chooses, as the pivot, the first nonzero entry that lies
on or below the diagonal in the current column.

revision 03/12/12  changed filename, commented out tests [ds]
revision 02/17/15  fixed test for zero entry in findPivotrow1 [ds]

matrix examples: [[28, 12, 20, 28], [4, 32, 28, 16]] , 2 by 4
                 [[28, 12, 20, 28]] , 1 by 4
                 [[[28], [12], [20], [28]] , 4 by 1


vector example [28, 12, 20, 28] 

"""
import math
import matplotlib.pyplot as plt
from diff_int import gaussian_quadrature2, gaussian_quadrature

### Supporting functions for Naive Gaussian Elimination function ge_1


def rows(mat):
    "return number of rows"
    return(len(mat))

def cols(mat):
    "return number of cols"
    return(len(mat[0]))
 
def zero(m,n):
    "Create zero matrix"
    new_mat = [[0 for col in range(n)] for row in range(m)]
    return new_mat
 
def transpose(mat):
    "return transpose of mat"
    new_mat = zero(cols(mat),rows(mat))
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            new_mat[col][row] = mat[row][col]
    return(new_mat)

def dot(A,B):
    "vector dot product"
    if len(A) != len(B):
        print("dot: list lengths do not match")
        return()
    dot=0
    for i in range(len(A)):
        dot = dot + A[i]*B[i]
    return(dot)

def getCol(mat, col):
    "return column col from matrix mat"
    return([r[col] for r in mat])

def getRow(mat, row):
    "return row row from matrix mat"
    return(mat[row])

def matMult(mat1,mat2):
    "multiply two matrices"
    if cols(mat1) != rows(mat2):
        print("matMult: mismatched matrices")
        return()
    prod = zero(rows(mat1),cols(mat2))
    for row in range(rows(mat1)):
        for col in range(cols(mat2)):
            prod[row][col] = dot(mat1[row],getCol(mat2,col))
    return(prod)

def vectorQ(V):
    "mild test to see if V is a vector"
    if type(V) != type([1]):
        return(False)
    if type(V[0]) == type([1]):
        return(False)
    return(True)

def scalarMult(a,mat):
    "multiply a scalar times a matrix"
    if vectorQ(mat):
        return([a*m for m in mat])
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            mat[row][col] = a*mat[row][col]
    return(mat)

def addVectors(A,B):
    "add two vectors"
    if len(A) != len(B):
        print("addVectors: different lengths")
        return()
    return([A[i]+B[i] for i in range(len(A))])

def swaprows(M,i,j):
    "swap rows i and j in matrix M"
    N=copyMatrix(M)
    T = N[i]
    N[i] = N[j]
    N[j] = T
    return N

def copyMatrix(M):
    return([[M[row][col] for col in range(cols(M))]for row in
            range(rows(M))])

def addrows(M, f, t, scale=1):
    "add scale times row f to row t"
    N=copyMatrix(M)
    T=addVectors(scalarMult(scale,N[f]),N[t])
    N[t]=T
    return(N)
    
def show(mat):
    "Print out matrix"
    for row in mat:
        print(row)

### vectors vs rowVectors and colVectors
### the latter are matrices

def vec2rowVec(vec):
    "[a,b,c] -> [[a,b,c]]"
    return([vec])

def vec2colVec(vec):
    return(transpose(vec2rowVec(vec)))

def colVec2vec(col_mat):
    rowVec = transpose(col_mat)
    return(rowVec[0])

def augment(mat,vec):
    "given nxn mat and n length vector return augmented matrix"
    amat = []
    for row in range(rows(mat)):
        amat.append(mat[row]+[vec[row]])
    return(amat)



### The naive gaussian elimination code begins here.

def findPivotrow1(mat,col):
#    Finds index of the first row with nonzero entry on or
#    below diagonal.  If there isn't one return(-1).

    epsilon = 10**(-17)
    for row in range(col, rows(mat)):
#        if mat[row][col] != 0:
        if abs(mat[row][col]) > epsilon:
            return(row)
    return(-1)


def partialPivot(mat, row, num_of_rows, scale_factors=None):
    """Find the index of the pivot row.

    This function uses either partial pivoting or scaled partial
    pivoting to find the pivot row.  To use scaled partial pivoting,
    the function must be passed a list of scale factors.

    Args:
        mat: An nxn matrix.
        row: the index of the row to start from.  Pivot rows will
             be searched from [row, n].
        num_of_rows: The number of rows in mat.
        scale_factors (list): A list of the scale factors to use in
                              scaled partial pivoting.  If this is not
                              given, then partial pivoting is used.
                              Defaults to None.

    Returns:
        The index of the pivot row, or -1 if no pivot row could be
        found.
    """
    piv_row = -1
    piv_val = -1
    for cur_row in range(row, num_of_rows):
        entry = mat[cur_row][row] # partial pivoting.
        if scale_factors is not None:
            entry = entry / scale_factors[cur_row] # scaled partial pivoting.
        if abs(entry) > piv_val:
            piv_val = abs(entry)
            piv_row = cur_row
    return piv_row

def identity(n):
    """Return an n by n identity matrix."""
    mat = []
    for row in range(n):
        mat.append([])
        for col in range(n):
            if row == col:
                mat[row].append(1)
            else:
                mat[row].append(0)
    return mat

def lu_factor(mat):
    """Fatorize a matrix into its lower and upper parts.

    No check is done to ensure that the given matrix is indeed LU
    factorable.  The caller must know in advance that the matrix
    can safely be factored.  Matrices that cause this function to
    raise an error are those that need to interchange rows during
    row reduction.

    Returns:
        A list of mat's upper and lower triangular matrices.
    """
    n = len(mat)
    L = identity(n)
    U = matMult(L, mat)
    for col_idx in range(n - 1):
        upper = identity(n)
        lower = identity(n)
        for row_idx in range(col_idx + 1, n):
            multiplier = U[row_idx][col_idx] / U[col_idx][col_idx]
            upper[row_idx][col_idx] = -multiplier
            lower[row_idx][col_idx] = multiplier
        U = matMult(upper, U)
        L = matMult(L, lower)
    return [L, U]

def swap(vec, j, k):
    """Swap vec[j] and vec[k]."""
    temp = vec[j]
    vec[j] = vec[k]
    vec[k] = temp

def rowReduce(M):
    "returns a row reduced version of M"

    N = copyMatrix(M)
    cs = cols(M)-2   # no need to consider last two cols
    rs = rows(M)
    scale_factors = []
    # Initialize the scales for each row
    for row in N:
        max = -1
        for col in row:
            col_abs = abs(col)
            if col_abs > max:
                max = col_abs
        scale_factors.append(max)
    # Reduce.
    for col in range(rs-1):
        #j = findPivotrow1(N,col)
        #j = partialPivot(N,col, len(N))
        j = partialPivot(N,col, len(N), scale_factors)
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swaprows(N,col,j)
                swap(scale_factors,col,j)
            scale = -1.0 / N[col][col]
            for row in range(col+1,rs):                
                N=addrows(N, col, row, scale * N[row][col])
    return(N)


def backSub(M, rhs_idx=0):
    """Solves a row reduced matrix.

    Args:
        M: an nx(n+m) augmented matrix.  The left matrix is an nxn
           matrix.  The right matrix is an nxm matrix.
        rhs_idx: The index of the first column of the right matrix.

    Returns:
        The solution vector.
    """

#   given a row reduced augmented matrix with nonzero 
#   diagonal entries, returns a solution vector
    

    if rhs_idx == 0:
        rhs_idx = cols(M)
    cs = rhs_idx-1 # cols not counting augmented col
    sol = [0 for i in range(cs)] # place for solution vector
    for i in range(1, rhs_idx):
        row = cs-i # work backwards
        sol[row] = ((M[row][cs] - sum([M[row][j]*sol[j] for
                    j in range(row+1,cs)])) / M[row][row]) 
    return(sol)

def append_col(mat, vec):
    """Append a column to a matrix."""
    for mat_row, vec_row in zip(mat, vec):
        mat_row.append(vec_row)

def inverse(mat):
    """Find the inverse of a matrix.

    Args:
        mat: An nxn matrix.

    Returns:
        The inverse of mat, if it exists.
        An empty list, otherwise.
    """
    n = len(mat)
    if n != len(mat[0]):
        return []
    identity_n = identity(n)
    inverse = []
    # Create the augmented matrix [A | I].
    for mat_row, identity_row in zip(mat, identity_n):
        mat_row.extend(identity_row)
    mat = rowReduce(mat)
    for k in range(n):
        inverse.append([])
    for col_idx in range(n, n + n):
        # Create the augmented matrix [A' | B[col_idx]].  The solution to this
        # is the next column vector of the inverse.
        for row_idx in range(n):
            mat[row_idx][n] = mat[row_idx][col_idx]
        sol = backSub(mat, n + 1)
        append_col(inverse, sol)
    show(inverse)
    return inverse


def diag_test(mat):

#   Returns True if no diagonal element is zero, False
#   otherwise.
    

    for row in range(rows(mat)):
        if mat[row][row] == 0:
            return(False)
    else:
        return(True)


def ge_1(aug):    

#   Given an augmented matrix it returns a list.  The [0]
#   element is the row-reduced augmented matrix, and 
#   ge_1(aug)[1] is a solution vector.  The solution
#   vector is empty if there is no unique solution.
    

    aug_n = rowReduce(aug)
    if diag_test(aug_n):
        sol = backSub(aug_n)
    else:
        print("\nge_1(): There is no unique solution")
        sol = []
    results = [aug_n, sol]
    return(results)


### The next two functions support checking a solution.

def getAandb(aug):
    "Returns the coef. matrix A and the vector b of Ax=b"
    m = rows(aug)
    n = cols(aug)
    A = zero(m,n-1)
    b = zero(m,1)
    for i in range(m):
        for j in range(n-1):
            A[i][j] = aug[i][j]
            
    for i in range(m):
        b[i] = aug[i][n-1]
    Aandb = [A,b]
    return(Aandb)

def checkSol_1(aug,x):
    "For aug=[A|b], returns Ax, b, and b-Ax as vectors"
    A  = getAandb(aug)[0]
    b  = getAandb(aug)[1]
    x_col_vec = vec2colVec(x)
    Ax = matMult(A,x_col_vec)
    r  = addVectors(b,scalarMult(-1.0,colVec2vec(Ax)))
    L  = [Ax,b,r]
    return(L)

def lsp_discrete(x, y, n):
    """Calculate the coefficients of a least squares polynomial.

    This function uses the normal power function basis for the
    polynomials.

    Args:
        x: the list of independent variables.
        y: the list of dependent variables such that x[k] = y[k].
        n: the degree of the least squares polynomial for the data.

    Returns:
        A list of coefficients for the least squares polynomial of
        degree n.  The k'th element is the coefficient of the k'th
        degree term.
    """
    m = len(x)
    n = n + 1
    sum = 0
    A = []
    # Initialize matrix.
    for j in range(n):
        A.append([0 for _ in range(n + 1)])
    for j in range(n):
        for k in range(n):
            sum = 0
            #g = lambda x: x**(j + k)
            #sum = gaussian_quadrature2(g, a, b)
            for i in range(m):
                sum = sum + (x[i]**(j + k))
            A[j][k] = sum
        sum = 0
        #g = lambda x: x**(j) * f(x)
        #sum = gaussian_quadrature2(g, a, b)
        for i in range(m):
            sum = sum + ((x[i]**j)*y[i])
        A[j][-1] = sum # This is the b vector of the equation.
    return ge_1(A)[1]

def lsp_continuous(f, a, b, n):
    """Calculate the coefficients of a least squares polynomial.

    This function uses the normal power function basis for the
    polynomials.

    Args:
        f: the function to approximate.
        a: the lower bound of the interval of approximation.
        b: the upper bound of the interval of approximatio.
        n: the degree of the least squares polynomial for the data.

    Returns:
        A list of coefficients for the least squares polynomial of
        degree n.  The k'th element is the coefficient of the k'th
        degree term.
    """
    n = n + 1
    sum = 0
    A = []
    # Initialize matrix.
    for j in range(n):
        A.append([0 for _ in range(n + 1)])
    for j in range(n):
        for k in range(n):
            sum = 0
            g = lambda x: x**(j + k)
            sum = gaussian_quadrature2(g, 5, a, b)
            A[j][k] = sum
        sum = 0
        g = lambda x: x**(j) * f(x)
        sum = gaussian_quadrature2(g, 5, a, b)
        A[j][-1] = sum # This is the b vector of the equation.
    #show(A)
    show(ge_1(A)[1])

def lsp_continuous_legendre(f, a, b, n):
    """Calculate the coefficients of a least squares polynomial.

    This function uses the Legendre polynomials as the basis for the
    polynomials.

    Args:
        f: the function to approximate.
        a: the lower bound of the interval of approximation.
        b: the upper bound of the interval of approximatio.
        n: the degree of the least squares polynomial for the data.

    Returns:
        A list of coefficients for the least squares polynomial of
        degree n.  The k'th element is the coefficient of the k'th
        degree term (which in this case is the k'th Legendre
        polynomial).  For example, c0*P0(x) + ... + cn*Pn(x). If a
        is not -1 and b is not 1, then the Legendre polynomials should
        be shifted to be of the form Pk((2x - a - b)/(b - a)) in order
        to see the correct results.

    """
    n = n + 1
    A = []
    #legendre = [lambda x: 1.0, lambda x: x, lambda x: 0.5*(3*(x**2) - 1),
                #lambda x: 0.5*(5*(x**3) - (3*x)),
                #lambda x: (1.0/8)*(35*(x**4) - (30*(x**2)) + 3),
               #]
    legendre = [lambda x: 1.0, lambda x: x, lambda x: 0.5*(3*(x**2) - 1),
                lambda x: 0.5*(5*(x**3) - (3*x)),
                lambda x: (1.0/8)*(35*(x**4) - (30*(x**2)) + 3),
               ]
    # Initialize matrix.
    for j in range(n):
        A.append([0 for _ in range(n + 1)])
    for j in range(n):
        for k in range(n):

            def g(x):
                x = ((2.0*x) - a - b)/(b - a)
                return legendre[j](x) * legendre[k](x)

            A[j][k] = gaussian_quadrature2(g, 5, a, b)
        #g = lambda x: f(x) * legendre[j](x)

        def g(x):
            x1 = ((2.0*x) - a - b)/(b - a)
            return f(x) * legendre[j](x1)

        # This is the b vector of the equation.
        A[j][-1] = gaussian_quadrature2(g, 5, a, b) 
    return ge_1(A)[1]

def plot2(x_vals, y_vals, x2_vals, y2_vals, x_min, x_max, y_min, y_max, title):
    plt.axvline(x=0,color='black')
    plt.axhline(y=0,color='black')
    plt.axis([x_min, x_max, y_min, y_max], 1/6)
    plt.scatter(x_vals, y_vals, color='red', label='data points')
    plt.plot(x2_vals, y2_vals, color='blue', label='LSE')
    plt.title(title)
    plt.legend(loc='top right')
    plt.show()

def plot(x_vals, y_vals, x_min, x_max, y_min, y_max, title):
    plt.axvline(x=0,color='black')
    plt.axhline(y=0,color='black')
    plt.axis([x_min, x_max, y_min, y_max], 1/6)
    plt.plot(x_vals, y_vals[0], label='f(x)')
    plt.plot(x_vals, y_vals[1], label='1/2')
    for k in range(2, len(y_vals)):
        plt.plot(x_vals, y_vals[k], label='cos(pi*x*' + str(k - 1) + ')')
        #plt.plot(x_vals, y_vals[k], label='degree ' + str(k - 1))
    plt.legend(loc='lower center')
    plt.title(title)
    plt.show()


#degree: 1
#329.013193034
#
#degree: 2
#0.00144291288593
#
#degree: 3
#0.000527341203057
#
#degree: 4
#4.24597021092e-05
#
def problem_1():
    x_vals = [0.0, 0.2, 0.5, 0.7, 1.1, 1.5, 1.9, 2.3, 2.8, 3.1]
    y_vals = [2.56, 13.18, 30.11, 42.05, 67.53, 95.14, 124.87, 156.73, 199.50,
         226.72]
    x_best_fit = [k * 0.01 for k in range(320)]
    for degree in range(1, 5):
        print
        print "degree:", degree
        coeffs = lsp_discrete(x_vals, y_vals, degree)
        
        def best_fit(x):
            sum = 0
            for k in range(len(coeffs)):
                sum = sum + (coeffs[k] * (x)**k)
            return sum

        y_best_fit = [best_fit(x) for x in x_best_fit]
        plot2(x_vals, y_vals, x_best_fit, y_best_fit, 0, 3.3, 0, 200.0,
              'Least Squares Polynomial of Degree ' + str(degree))
        print sum([(y_vals[k] - best_fit(x_vals[k]))**2 for k in range(len(x_vals))])

def problem_2():
    t = [0.0100, 0.9998, 2.1203, 3.0023, 3.9892, 5.0017]
    y = [1.9262, 1.0042, 0.4660, 0.2496, 0.0214, 0.0130]
    exp = lambda k, x2: math.e**(t[k] * x2)
    diff = lambda k, x1, x2: y[k] - (x1*exp(k, x2))

    def dE_dx1(x1, x2):
        sum = 0
        for k in range(len(t)):
            sum = sum + (diff(k, x1, x2) * -exp(k, x2))
        return 2.0 * sum

    def dE_dx2(x1, x2):
        sum = 0
        for k in range(len(t)):
            sum = sum + (diff(k, x1, x2) * x1 * t[k] * -exp(k, x2))
        return 2.0 * sum

    def dE_dx1dx1(x1, x2):
        sum = 0
        for k in range(len(t)):
            sum = sum + (exp(k, 2.0*x2))
        return 2.0 * sum

    def dE_dx1dx2(x1, x2):
        sum = 0
        for k in range(len(t)):
            term1 = x1 * t[k] * exp(k, 2.0*x2)
            term2 = (diff(k, x1, x2) * t[k] * exp(k, x2))
            sum = sum + (term1 - term2)
        return 2.0 * sum

    def dE_dx2dx2(x1, x2):
        sum = 0
        for k in range(len(t)):
            term1 = ((x1*t[k])**2.0) * exp(k, 2.0*x2)
            term2 = (diff(k, x1, x2) * x1 * (t[k]**2.0) * exp(k, x2))
            sum = sum + (term1 - term2)
        return 2.0 * sum

    x1 = 2.145046416783
    x2 = -0.720054160509
    for k in range(10):
        system = [
                     [dE_dx1dx1(x1, x2), dE_dx1dx2(x1, x2), -dE_dx1(x1, x2)],
                     [dE_dx1dx2(x1, x2), dE_dx2dx2(x1, x2), -dE_dx2(x1, x2)],
                 ]
        dx = ge_1(system)[1]
        x1 = x1 + dx[0]
        x2 = x2 + dx[1]
        print x1, x2
        print dE_dx1(x1, x2)
        print dE_dx2(x1, x2)
        print
    t_best_fit = [k*0.01 for k in range(540)]
    y_best_fit = [x1 * math.e**(x2 * x) for x in t_best_fit]
    #print str((len(t) - len(y)))
    plot2(t, y, t_best_fit, y_best_fit, 0, 5.1, 0, 2.0,
          'Least Squares Exponential')
    print sum(iter((y[x] - (x1 * math.e**(x2 * t[x])))**2
              for x in range(len(t))))
### Test examples.

#problem_1()
#problem_2()

# Example 1, Ch.8.1 in Burden & Faires 9E.
x = [j for j in range(1, 11)]
y = [1.3, 3.5, 4.2, 5.0, 7.0, 8.8, 10.1, 12.5, 13.0, 15.6]
#lsp(x, y, 1)

print
# Example 2, Ch.8.1 in Burden & Faires 9E.
x = [0.0, 0.25, 0.5, 0.75, 1.0]
y = [1.0000, 1.2840, 1.6487, 2.1170, 2.7183]
#lsp(x, y, 2)

print
# Example 1, Ch.8.2 in Burden & Faires 9E.
#lsp_continuous(lambda x: math.sin(math.pi * x), 0, 1, 2)

# To use the Legendre polynomials for their orthogonality, first do a linear
# transformation to change the interval of integration to [-1, 1].
#lsp_continuous_legendre(lambda x: 0.5*((0.5*(x + 1))**(0.5)), -1, 1, 3)
#print lsp_continuous_legendre(lambda x: x**(3), -1, 1, 4)

def problem_3():

    def f(x):
        if x >= 0.5 and x <= 1:
            return 0
        if x <= -0.5 and x >= -1:
            return 0
        if (x < -1) or (x > 1):
            return
        return 0.5

    a = -1.0
    b = 1.0
    legendre = [lambda x: 1.0, lambda x: x, lambda x: 0.5*(3*(x**2) - 1),
                lambda x: 0.5*(5*(x**3) - (3*x)),
                lambda x: (1.0/8)*(35*(x**4) - (30*(x**2)) + 3),
               ]
    x_best_fit = [k * 0.01 for k in range(-310, 310)]
    y_vals = [f(x) for x in x_best_fit]
    y_best_fits = [y_vals]
    coeffs = []
    for k in range(0, 5):
        f1 = lambda x: f(x) * legendre[k](x)
        f2 = lambda x: (legendre[k](x))**2
        numer = gaussian_quadrature2(f1, 5, -0.5, 0.5)
        denom = gaussian_quadrature2(f2, 5, -1.0, 1.0)
        coeffs.append(numer / denom)

        def best_fit(x):
            sum = 0
            for k in range(len(coeffs)):
                sum = sum + (coeffs[k] * legendre[k](x))
            return sum

        y_best_fit = [best_fit(x) for x in x_best_fit]
        y_best_fits.append(y_best_fit)
        #plot(x_best_fit, y_vals, x_best_fit, y_best_fit, -1, 1, -1, 1.0, 'title')
    plot(x_best_fit, y_best_fits, -1, 1, -1, 1,
             'Least Squares Polynomials Using Legendre Polynomials')
    plot(x_best_fit, y_best_fits, -3, 3, -5, 5,
             'Least Squares Polynomials Using Legendre Polynomials')

def problem_4():

    def f(x):
        if x >= 0.5 and x <= 1:
            return 0
        if x <= -0.5 and x >= -1:
            return 0
        if (x < -1) or (x > 1):
            return
        return 0.5

    a = -0.5
    b = 0.5
    x_best_fit = [k * 0.01 for k in range(-310, 310)]
    y_best_fits = []
    y_vals = [f(x) for x in x_best_fit]
    coeffs = []
    y_best_fits.append(y_vals)
    for k in range(0, 6):
        if k == 0:
            f1 = lambda x: 0.5 * f(x)
            f2 = lambda x: 0.25
        else:
            f1 = lambda x: f(x) * math.cos(k * math.pi * x)
            f2 = lambda x: math.cos(k * math.pi * x)**2
        numer = gaussian_quadrature2(f1, 5, -0.5, 0.5)
        denom = gaussian_quadrature2(f2, 5, -1.0, 1.0)
        coeffs.append(numer / denom)

        def best_fit(x):
            sum = 0
            sum = sum + (coeffs[0] * 0.5)
            for k in range(1, len(coeffs)):
                sum = sum + (coeffs[k] * math.cos(x*k*math.pi))
            return sum

        y_best_fit = [best_fit(x) for x in x_best_fit]
        y_best_fits.append(y_best_fit)
    print coeffs
    plot(x_best_fit, y_best_fits, -1, 1, -1.0, 1.4,
             'Least Squares Polynomials Using Cosines')
    plot(x_best_fit, y_best_fits, -3, 3, -1.0, 1.4,
             'Least Squares Polynomials Using Cosines')

#problem_3()
problem_4()
