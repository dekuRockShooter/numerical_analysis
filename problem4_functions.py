from math import exp, sin, pi, e
from operator import mul
import numpy as np
import matplotlib.pyplot as plt

list1 = [0., 1./24., 1./12., 1./6., 1./4., 1./3., 5./12., 1./2., 2./3., 5./6., 1.]
list2 = [0., 1./10., 2./10., 3./10., 4./10., 5./10., 6./10., 7./10., 8./10., 9./10., 1.]

def g_func(x):
    "function g(x)=(2e^(-2x))*sin(3(pi)x)"
    ans = 2 * exp(-2 * x) * sin(3 * pi * x)
    return ans

def get_g_results(nodelist):
    ans = []
    for element in nodelist:
        ans.append(g_func(element))
    return nodelist, ans

def frange(start, end, jump):
    "Creates a range that jump as its step. Use it for floats"
    ans = []
    while start <= end*1.01:
        ans.append(start)
        start += jump
    return ans

def pi_func(x, nodelist):
    x_minus_nodes = []

    for element in nodelist:
        x_minus_nodes.append(x - element)

    #multiply all elements in list
    ans = reduce(mul, x_minus_nodes)
    return ans

def get_pi_results(start, end, nodelist):
    #range from start to end with step of 0.005
    #x_values = frange(start, end, 0.005)
    x_values = np.arange(start, end, 0.005)
    ans = []
    for element in x_values:
        ans.append(pi_func(element, nodelist))
    return x_values, ans

def lagrange_6(x, nodelist):
    #all elements from nodelist except nodelist[6]
    temp = nodelist[:6] + nodelist[(7):]
    num = []
    denom = []
    #append to numerator and denominator lists
    for element in temp:
        num.append(x - element)
        denom.append(nodelist[6] - element)
    #multiply all elements in numerator and denominator
    #finally divide
    ans = reduce(mul, num) / reduce(mul, denom)
    return ans

def get_lagrange_results(start, end, nodelist):
    #range from start to end with step of 0.005
    x_values = np.arange(start, end, 0.005)
    ans = []
    for element in x_values:
        ans.append(lagrange_6(element, nodelist))
        #ans.append(_lagrange(6, element, nodelist))
    return x_values, ans

def plotLagrange(lagrange, xmin, xmax, ymin, ymax, xticks):
    plt.title(r'$L_6(x)$')
    plt.plot(lagrange[0],lagrange[1],color='blue', label = 'Lagrange_6 Values')
    
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xticks(xticks)
    plt.grid(True)
    plt.axhline(y=0,color='black')
    plt.axvline(x=0,color='black')
    plt.legend(loc = 'upper left')

    plt.show()

def plotPi(piValues, xmin, xmax, ymin, ymax):
    plt.title(r'$Pi(x)$')
    plt.plot(piValues[0],piValues[1],color='red', label = 'Pi Values')
    
    plt.axis([xmin, xmax, ymin, ymax])
    plt.grid(True)
    plt.axhline(y=0,color='black')
    plt.axvline(x=0,color='black')
    plt.legend(loc = 'upper left')

    plt.show()

#prints the return lists of the get_*_results() functions
#Usage: print_results(get_*_results(...))
def print_results(two_lists):
    for x,y in zip(two_lists[0], two_lists[1]):
        print(str(x),str(y))

#********************************LAGRANGE***********************************
#***************************************************************************
# FIRST NODE SET:
# L6 with first node set
approx_coords = get_lagrange_results(-2.0, 2.0, list1)

# [x0, xn]
plotLagrange(approx_coords, -0.1, 1.1, -2, 2, list1)
# [x0, xn] zoomed out
plotLagrange(approx_coords, -0.1, 1.1, -160, 30, list1)

# Show before x0
plotLagrange(approx_coords, -1, 1, -20, 200, list1)
# Show after xn
plotLagrange(approx_coords, 0.0, 2, -20, 200, list1)

# Everything zoomed out
plotLagrange(get_lagrange_results(-0.25, 1.1, list1), -0.5, 1.5, -210, 5000,
             list1)
#***************************************************************************
# SECOND NODE SET:
# L_6 with second node set
approx_coords = get_lagrange_results(-2.0, 2.0, list2)

# [x0, xn]
plotLagrange(approx_coords, -0.1, 1.1, -2, 2, list2)
# [x0, xn] zoomed out
plotLagrange(approx_coords, -0.1, 1.1, -8, 2, list2)

# Show before x0
plotLagrange(approx_coords, -1, 1, -20, 200, list2)
# Show after xn
plotLagrange(approx_coords, 0.0, 2, -20, 200, list2)

# Everything zoomed out
plotLagrange(get_lagrange_results(-0.25, 1.1, list2), -0.5, 1.5, -210, 5000,
             list2)
#***************************************************************************

#********************************** PI *************************************
#***************************************************************************
# FIRST NODE SET:
# pi with first node set
approx_coords = get_pi_results(-2., 2.01, list1)

# [x0, xn]
plotPi(approx_coords , -0.1, 1.1, -0.000010, 0.000010)
# [x0, xn] zoomed in
plotPi(approx_coords , -0.1, 1.1, -0.0000003, 0.0000003)
# [x0, xn] zoomed out
plotPi(approx_coords , -0.2, 1.1, -0.00011, 0.00011)

# Show before x0
plotPi(approx_coords , -0.4, 1.1, -0.00013, 0.00002)
# Show after xn
plotPi(approx_coords , -0.2, 1.1, -0.00011, 0.0003)

# Everything zoomed out
plotPi(approx_coords, -0.4, 1.2, -0.006, 0.006)
#***************************************************************************
# SECOND NODE SET:
# pi with second node set
approx_coords = get_pi_results(-2., 2.01, list2)

# [x0, xn]
plotPi(approx_coords , -0.1, 1.1, -0.000010, 0.000010)
# [x0, xn] zoomed in
plotPi(approx_coords , -0.1, 1.1, -0.0000003, 0.0000003)
# [x0, xn] zoomed out
plotPi(approx_coords , -0.2, 1.2, -0.00011, 0.00011)

# Show before x0
plotPi(approx_coords , -0.4, 1.1, -0.00013, 0.00002)
# Show after xn
plotPi(approx_coords , -0.1, 1.2, -0.00005, 0.00035)

# Everything zoomed out
plotPi(approx_coords, -0.4, 1.4, -0.010, 0.010)
