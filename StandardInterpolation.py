from math import *
import matplotlib.pyplot as plt
import numpy as np

def plotterActual():
    A = [0., 1./24, 1./12, 1./6, 1./4, 1./3, 5./12, 1./2, 2./3, 5./6, 1.]
    B = []
    for i in A:
        B.append(yValue(i))
    return A,B


    
  
def valuecurve():
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
    C = []
    D = []
    for x in np.arange (0,1.005,0.005):
        C.append(x)
        D.append(interpolate(A,B,x))
    return C,D
           
def yValue(x):
    yVal = (2*sin(3*pi*x))/(e**(2*x))
    return yVal
def derivative(x):
    
    deriv = -(4*sin(3*pi*x)-6*pi*cos(3*pi*x))/(e**(2*x))
    return deriv

def interpolate(A,B, current):
    C = []
    
    sum = 0
    for i in range(len(A)):
        temp = 1
        k = i
        for j in range(len(A)):
            if(k!=j):
                temp = temp * (current - A[j])/(A[k]-A[j])
        
        C.append(B[i] * temp)
        
    for i in range(len(C)):
        sum = sum + C[i]
       
    return sum


def show(vec):
    for row in vec:
        print row
def plot(xValActual,yValActual,xValIP,yValIP):
    plt.title(r'$(2*sin(3*pi*x))/(e**(2*x))$')
    plt.plot(xValActual,yValActual,color='red', label = 'Actual values')
    plt.plot(xValIP,yValIP,color='blue', label = 'Interpolated Values')
    
    plt.axis([-1,1.5,-2,2], 1/6)
    plt.axhline(y=0,color='black')
    plt.axvline(x=0,color='black')
    plt.legend(loc = 'lower right')

    plt.show()
    plt.savefig("example.pdf")
def plotDiff(xValActual,diffY):
    plt.title(r'$Plot Difference Between Actual Values and Interpolated Values$')
    plt.plot(xValActual,diffY, color = 'red', label = 'Difference')
    plt.axis([-.25,1.1,-.025,.10])
    plt.xticks(np.arange(min(x)-.25, max(x)+.15, .25))
    plt.axhline(y=0,color='black')
    plt.axvline(x=0,color='black')
    plt.legend(loc = 'lower right')
    plt.show()
    
xValActual,yValActual = valuecurve()
#plot(xValActual,yValActual)
x, y = plotterActual()
C = interpolate(x,y, x[1])


xVal, yVal = valueInter(x,y)

plot(xValActual,yValActual,xVal, yVal)

diffY = diff(yValActual, yVal)
#show(diffX)
#show(diffY)
plotDiff(xValActual, diffY)

#fromip = interpolate(xVa, Yva)


