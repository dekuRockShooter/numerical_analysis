from math import exp, sin, pi, e
from operator import mul

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
    ans = []
    while start <= end*1.01:
        ans.append(start)
        start += jump
    return ans

def pi_func(x, nodelist):
    x_minus_nodes = []
    for element in nodelist:
        x_minus_nodes.append(x - element)
    ans = reduce(mul, x_minus_nodes)
    return ans

def get_pi_results(start, end, nodelist):
    x_values = frange(start, end, 0.1)
    ans = []
    for element in x_values:
        ans.append(pi_func(element, nodelist))
    return x_values, ans

def print_results(two_lists):
    for x,y in zip(two_lists[0], two_lists[1]):
        print(str(x),str(y))
