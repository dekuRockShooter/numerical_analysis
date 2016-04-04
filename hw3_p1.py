from math import *

def g_x(x):
    return (exp((-x**2)/2))

def gx_prime(x):
    return ((exp((-x**2)/2))*(-x))

def central_diff(x, h):
    return ((g_x(x+h) - g_x(x-h))/(2*h)) 

def problem_1():
    true = gx_prime(1.4)
    print("g'(1.4) = %1.20f" %(true))
    print("            h                          g(1.4)                        r.error         ")
    print("  ------------------------     -------------------------     ------------------------")
    for i in range(1,21):
        value = central_diff(1.4,(10.**(-i)))
        print("   %.20f       %2.20f       %.20f    " %((10.**(-i)), value, abs((true-value))/abs(true)))
    return

problem_1()
