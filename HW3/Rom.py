import numpy as np
from math import *

#romberg integration over a function, not a set of points

# initial ,final , function to pass, number of devisions
def romberg_integration(x_i, x_f, funct, max_divs):
	#max_divs should be a power of 2
	estimates = []
	temp = []
	iteration = 1

	ret_val = 0
	#trapezoids
	for i in range(int(log(max_divs,2)) + 1):
		num_divs = int(pow(2., i))
		inc = (x_f - x_i) / num_divs
		x_init = x_i
		temp_sum = 0
		for j in range(num_divs):
			x_next = x_init + inc
			temp_sum += (x_next - x_init) * (funct(x_init) + funct(x_next)) / 2
			x_init = x_next
		temp.append(temp_sum)
	#append estimates
	estimates.append(temp)

	#reducing to the approximation
	while len(temp) != 1:
		temp = []
		mult = pow(4, iteration)
		trape_ests = estimates[iteration - 1]
		for i in range(len(trape_ests) - 1):
			temp.append( (mult * trape_ests[i+1] - 
					trape_ests[i]) / (mult - 1 ))
		iteration += 1
		estimates.append(temp)
	for i in estimates:
		print i
	return temp[0]

#
#
#

def sin(x):
        return np.sin(np.pi*x)
def exp(x):
        return np.e**(-(x**2)/2)
def sinc(x):
        if(x != 0):
                return (np.sin(2*pi*x))/((2*pi*x))
        else:
                return 1

# lower, upper , functions , n where n is number of xi - x0 
romberg_integration(0,1,sinc, 32)
