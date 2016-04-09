import scipy
import matplotlib.pyplot as plt
import numpy as np
import math 
from math import pow



# Computes Larrange Interpolation Polynomial Given list of points x and g(x)  
def lagrangeInterpolation(x, gx):
   #setting result = 0
   result = scipy.poly1d([0.0]) 
   for i in range(0,len(x)): #number of polynomials L_k(x).
      temp_numerator = scipy.poly1d([1.0]) # resets temp_numerator such that a new numerator can be created for each i.
      denumerator = 1.0 #resets denumerator such that a new denumerator can be created for each i.
      for j in range(0,len(x)):
          if i != j:
              temp_numerator *= scipy.poly1d([1.0,-x[j]]) #finds numerator for L_i
              denumerator *= x[i]-x[j] #finds denumerator for L_i
      result += (temp_numerator/denumerator) * gx[i] #linear combination
   return result;

   
            
    



# Main Test

# array of x values 
x = scipy.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8])
# array of g(x) values 
gx = scipy.array([1.0,.9058372,.81873075,.74081822,.67032005,.60653066,.54881164,.49658530,.44932896])

# gets the larrange interpolating polynomial result
lagrangePoly = lagrangeInterpolation(x, gx)
print "The Lagrange Polynomial: "
print lagrangePoly
print ""

# gets the derivative of the lagrange polynomial
derivativePoly = np.polyder(lagrangePoly)
print "The first derivative of the Lagrange Polynomial: "
print derivativePoly
print ""

#evaluates the derivative of the lagrange polynomial at x value 1.4: p'(1.4)
evaluatedDerivative = derivativePoly(1.4)
print "The derivative evaluated at p'(1.4): "
print evaluatedDerivative


