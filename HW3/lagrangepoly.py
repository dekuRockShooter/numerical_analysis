import scipy
import matplotlib.pyplot as plt
import numpy as np
import math 
from math import pow



""" Computes the Larrange interpolation polynomial given a list of x and g(x) values using a poly1d array to hold the poly
   Args: 
      x: a list of x values
      g(x): a list of y values
   Returns:
      result:  a poly1d array that holds the Lagrange polynomial. A poly1d array is well suited for this task
               because printing the array shows the degree of x. We can also find the derivative in one line of code.
"""
def lagrangeInterpolation(x, gx):
   result = scipy.poly1d([0.0]) 
   for i in range(0,len(x)): 
      temp_numerator = scipy.poly1d([1.0]) 
      denumerator = 1.0 
      for j in range(0,len(x)):
          if i != j:
              temp_numerator *= scipy.poly1d([1.0,-x[j]]) #numerator for Li
              denumerator *= x[i]-x[j] #denumerator for Li
      result += (temp_numerator/denumerator) * gx[i] #create the linear combo
   return result;




# array of x values 
x = scipy.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8])
# array of g(x) values 
gx = scipy.array([1.0,.9058372,.81873075,.74081822,.67032005,.60653066,.54881164,.49658530,.44932896])

# gets the larrange interpolating polynomial result
lagrangePoly = lagrangeInterpolation(x, gx)
print "The Lagrange Polynomial: "
print lagrangePoly
print ""

# gets the derivative of the lagrange polynomial using polyder, so handy!
derivativePoly = np.polyder(lagrangePoly)
print "The first derivative of the Lagrange Polynomial: "
print derivativePoly
print ""

#evaluates the derivative of the lagrange polynomial at x value 1.4: p'(1.4)
evaluatedDerivative = derivativePoly(1.4)
print "The derivative evaluated at p'(1.4): "
print evaluatedDerivative


