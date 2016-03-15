from math import *
import matplotlib.pyplot as plt
import numpy as np

def plotterActual():
    """Return the lists x_k and y_k, calculated as 2sin(3pix)e^-2x."""
    A = [0., 1./24, 1./12, 1./6, 1./4, 1./3, 5./12, 1./2, 2./3, 5./6, 1.]
    B = []
    for i in A:
        B.append(yValue(i))
    return A,B


    
  
def valuecurve():
    """Return the lists .005k and y_k, calculated as 2sin(3pix)e^-2x."""
    A = []
    B = []
    for x in np.arange (0,1.005,0.005):
        A.append(x)
        B.append(yValue(x))
    return A,B

def diff(A,B):
    E = []
 
    for i in range(len(A)):
        E.append(abs(A[i] - B[i]))
 
    return E

def valueInter(A,B):
    """Return the lists .005k and y_k, approx to 2sin(3pix)e^-2x."""
    C = []
    D = []
    for x in np.arange (0,1.005,0.005):
        C.append(x)
        D.append(lagrange_interpolate(A,B,x))
    return C,D
           
def yValue(x):
    """Returns 2sin(3pix)e^-2x"""
    yVal = (2*sin(3*pi*x))/(e**(2*x))
    return yVal
def derivative(x):
    """Returns the derivative of 2sin(3pix)e^-2x"""
    
    deriv = -(4*sin(3*pi*x)-6*pi*cos(3*pi*x))/(e**(2*x))
    return deriv

def lagrange_interpolate(A,B, current):
    """Uses Lagrange interpolation with the nodes and values in A, B."""
    C = []
    
    sum = 0
    for i in range(len(A)):
        temp = 1
        k = i
        temp = _lagrange(k, current, A)
        C.append(B[i] * temp)
        
    for i in range(len(C)):
        sum = sum + C[i]
       
    return sum

def _lagrange(k, x, nodes):
    """Computes Lk(x).

    The k'th Lagrange polynomial is evaluated using the given nodes.

    Args:
        k: the index of the node for which the function returns 1.
        x: a real number.
        nodes: the nodes (given x values) for which the function
               values are known.

    Returns:
        Lk(x).
    """
    numerator = 1;
    denomenator = 1;
    result = 0
    for j in range(len(nodes)):
        if nodes[j] == nodes[k]:
            continue
        numerator = numerator * (x - nodes[j])
        denomenator = denomenator * (nodes[k] - nodes[j])
    return 1.0 * numerator / denomenator

def _lagrange_derivative(k, nodes):
    """Evaluate the derivative of the k'th Lagrange function at x_k.

    Args:
        k: the index of the node at which to evaluate the derivative.
        nodes: the known set of x values.
    """
    sum = 0
    for j in range(len(nodes)):
        if j != k:
            sum = sum + 1.0 * 1 / (nodes[k] - nodes[j])
    return sum

def _hermite_left(k, x, nodes):
    """Evaluates the k'th Hermite function for function values.

    Args:
        k: the index of the node for which the function returns 1.
        x: a real number.
        nodes: the nodes (given x values) for which the function
               values are known.

    Returns:
        The function value evaluated at x.
    """
    term_1 = _lagrange(k, x, nodes)
    term_1 = term_1 * term_1
    derivative = _lagrange_derivative(k, nodes)
    return term_1 * (1 - (2 * derivative * (x - nodes[k])))

def _hermite_right(k, x, nodes):
    """Evaluates the k'th Hermite function for derivative values.

    Args:
        k: the index of the node for which the function returns 1.
        x: a real number.
        nodes: the nodes (given x values) for which the function
               values are known.

    Returns:
        The function value evaluated at x.
    """
    term_1 = _lagrange(k, x, nodes)
    term_1 = term_1 * term_1
    return term_1 * (x - nodes[k])

def hermite_interp(x, nodes, values, deriv_vals):
    """Interpolate a set of data using Hermite interpolation.

    Args:
        x: a real number to evaluate.
        nodes: the known set of x values.
        values: the known function values at each node.
        deriv_vals: the known function derivative values at each node.

    Returns:
        f(x), where f is the unknown function that is being
        approximated, such that f(nodes[k]) = values[k] and
        f'(nodes[k]) = deriv_vals[k]
    """
    sum = 0
    for k in range(len(nodes)):
        left_term = _hermite_left(k, x, nodes) * values[k]
        right_term = _hermite_right(k, x, nodes) * deriv_vals[k]
        sum = sum + left_term + right_term
    return sum

def _pos_linear(x, last_node=False):
    """Implement the positive linear spline basis function.

    Args:
        x: a real number.
        last_node: boolean for changing the bounds.

    Returns:
        x if x is in [0, 1) and last_node is False.
        x if x is in [0, 1] and last_node is True.
        0, otherwise.
    """
    if last_node:
        if (x >= 0) and (x <= 1):
            return x
    elif (x >= 0) and (x < 1):
        return x
    return 0

def _neg_linear(x):
    """Implement the negative linear spline basis function.

    Args:
        x: a real number.

    Returns:
        1 - x if x is in [0, 1).
        0, otherwise.
    """
    if (x >= 0) and (x < 1):
        return 1 - x
    return 0

def linear_spline_interp(x, nodes, values):
    """Interpolate a set of data using linear spline interpolation.

    Args:
        x: a real number to evaluate.
        nodes: the known set of x values.
        values: the known function values at each node.

    Returns:
        f(x), where f is the unknown function that is being
        approximated, such that f(nodes[k]) = values[k].  The values
        in between nodes are approximated with a line that connects
        the boundary nodes.
    """
    sum = _neg_linear(1.0*(x - nodes[0]) / (nodes[1] - nodes[0])) * values[0]
    for k in range(1, len(nodes) - 1):
        pos = _pos_linear(1.0*(x - nodes[k - 1]) / (nodes[k] - nodes[k - 1]))
        neg = _neg_linear(1.0*(x - nodes[k]) / (nodes[k + 1] - nodes[k]))
        sum = sum + (1.0*(pos + neg) * values[k])
    last_term = _pos_linear(1.0*(x - nodes[k]) / (nodes[k + 1] - nodes[k]),
                            True)
    return sum + (last_term * values[k + 1])

def _pos_cubic(x, last_node=False):
    """Positive cubic basis function.

    Args:
        x: a real number.
        last_node: boolean for changing the bounds.

    Returns:
        1 - _neg_cubic(x) if x is in [0, 1) and last_node is False.
        1 if x is 1 and last_node is True.
        0, otherwise.
    """
    if last_node:
        if (x >= 0) and (x <= 1):
            return 1 - _neg_cubic(x)
    if (x >= 0) and (x < 1):
        return 1 - _neg_cubic(x)
    return 0

def _neg_cubic(x):
    """Negative cubic basis function.

    Args:
        x: a real number.

    Returns:
        ((x - 1)**2) * (2*x + 1) if x is in [0, 1)
        0, otherwise.
    """
    if (x >= 0) and (x < 1):
        return ((x - 1)**2) * (2*x + 1)
    return 0

def _pos_cubic_deriv(x):
    """Positive cubic derivative basis function.

    Args:
        x: a real number.

    Returns:
        x * ((x - 1)**2) if x is in [0, 1)
        0, otherwise.
    """
    if (x >= 0) and (x < 1):
        return x * ((x - 1)**2)
    return 0

def _neg_cubic_deriv(x, last_node=False):
    """Negative cubic derivative basis function.

    Args:
        x: a real number.
        last_node: boolean for changing the bounds.

    Returns:
        (x**2) * (x - 1) if x is in [0, 1) and last_node is False.
        1 if x is 1 and last_node is True.
        0, otherwise.
    """
    if last_node and (x == 1):
        return 1
    if (x >= 0) and (x < 1):
        return (x**2) * (x - 1)
    return 0

def _cubic_basis(x, nodes, values):
    """Approximate input using the cubic basis functions.

    This function very roughly approximates f(x), where f is an unknown
    function whose values are known only at the nodes.  To get a good
    approximation, this must be used in coordination with
    _cubic_deriv_basis.

    Args:
        x: the input to approximate the function value of.
        nodes: the known set of x values.
        values: the known function values at each node.

    Returns:
        A crude approximation to f(x)
    """
    sum = values[0]*_neg_cubic(1.0*(x - nodes[0]) / (nodes[1] - nodes[0]))
    for k in range(1, len(nodes) - 1):
        pos = _pos_cubic(1.0*(x - nodes[k - 1]) / (nodes[k] - nodes[k - 1]))
        neg = _neg_cubic(1.0*(x - nodes[k]) / (nodes[k + 1] - nodes[k]))
        sum = sum + (1.0*(pos + neg) * values[k])
    last_term = values[k + 1]*_pos_cubic(1.0*(x - nodes[k]) /
                                         (nodes[k + 1] - nodes[k]), True)
    return sum + last_term

def _cubic_deriv_basis(x, nodes, deriv_vals):
    """Correct the approximation of _cubic_basis.

    This functions corrects over and under approximations from
    _cubic_basis.  It should be called using the same x values, nodes,
    and values as those used when calling _cubic_basis so that the
    approximations will be good.

    Args:
        x: the input to approximate the function value of.
        nodes: the known set of x values.
        deriv_vals: the known function derivative values at each node.

    Returns:
        A correction to the approximation of f(x), as approximated by
        _cubic_basis.
    """
    denom = (nodes[1] - nodes[0])
    sum = deriv_vals[0] * denom * _pos_cubic_deriv(1.0*(x - nodes[0]) / denom)
    for k in range(1, len(nodes) - 1):
        denom = (nodes[k + 1] - nodes[k])
        pos = _pos_cubic_deriv(1.0*(x - nodes[k]) / denom) * denom
        denom = (nodes[k] - nodes[k - 1])
        neg = _neg_cubic_deriv(1.0*(x - nodes[k - 1]) / denom) * denom
        sum = sum + (1.0*(pos + neg) * deriv_vals[k])
    denom = (nodes[k + 1] - nodes[k])
    last_term = deriv_vals[k + 1]*_neg_cubic_deriv(1.0*(x - nodes[k]) / denom,
                                                   False) * denom
    return sum + last_term

def cubic_spline_interp(x, nodes, values, deriv_vals):
    """Interpolate a set of data using cubic spline interpolation.

    Args:
        x: a real number to evaluate.
        nodes: the known set of x values.
        values: the known function values at each node.
        deriv_vals: the known function derivative values at each node.

    Returns:
        f(x), where f is the unknown function that is being
        approximated, such that f(nodes[k]) = values[k].  The values
        in between nodes are approximated with a cubic curve that
        connects the boundary nodes.
    """
    return _cubic_basis(x, nodes, values) +\
           _cubic_deriv_basis(x, nodes, deriv_vals)

def show(vec):
    for row in vec:
        print row

def plot(xValActual,yValActual,xValIP,yValIP, diffY, title):
    plt.title(title + r' of $(2*sin(3*pi*x))/(e**(2*x))$')
    plt.plot(xValActual,yValActual,color='red', label = 'Actual values')
    plt.plot(xValIP,yValIP,color='blue', label = 'Interpolated Values')
    plt.plot(xValActual,diffY, color = 'black', label = 'Absolute difference')

    plt.axis([0,1,-1,1.5], 1/6)
    plt.axhline(y=0,color='black')
    plt.axvline(x=0,color='black')
    plt.legend(loc = 'top')

    plt.show()
    plt.savefig("example.pdf")

def plotDiff(xValActual,diffY, title):
    plt.title(r'$Plot Difference Between Actual Values and '+ title +
              ' Interpolated Values$')
    plt.plot(xValActual,diffY, color = 'red', label = 'Difference')
    plt.axis([-.25,1.1,.0, max(diffY)])
    plt.xticks(np.arange(min(xValActual)-.25, max(xValActual)+.15, .25))
    plt.axhline(y=0,color='black')
    plt.axvline(x=0,color='black')
    plt.legend(loc = 'lower right')
    plt.show()


nodes, values = plotterActual()
deriv_vals = [derivative(x) for x in nodes]
x_005, y_005 = valuecurve()
y_005_approx_lagrange = valueInter(nodes, values)[1] # Lagrange approx
y_005_approx_hermite = [hermite_interp(x, nodes, values, deriv_vals)
                        for x in x_005]
y_005_approx_linear = [linear_spline_interp(x, nodes, values) for x in x_005]
y_005_approx_cubic = [cubic_spline_interp(x, nodes, values, deriv_vals)
                      for x in x_005]
# Lagrange plots
diff_lagrange = diff(y_005, y_005_approx_lagrange)
plot(x_005, y_005, x_005, y_005_approx_lagrange, diff_lagrange,
     'Lagrange interpolation')
## Hermite plots
diff_hermite = diff(y_005, y_005_approx_hermite)
diff_hermite_scaled = [10**7*x for x in diff_hermite]
plot(x_005, y_005, x_005, y_005_approx_hermite, diff_hermite_scaled,
     'Hermite interpolation')
# Linear spline
diff_linear = diff(y_005, y_005_approx_linear)
plot(x_005, y_005, x_005, y_005_approx_linear, diff_linear,
     'Linear spline interpolation')
# Cubic spline
diff_cubic = diff(y_005, y_005_approx_cubic)
diff_cubic_scaled = [10**2*x for x in diff_cubic]
plot(x_005, y_005, x_005, y_005_approx_cubic, diff_cubic_scaled,
     'Cubic spline interpolation')
