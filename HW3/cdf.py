fx = [1.00000000, 0.90483742, 0.81873075, 0.74081822, 0.67032005, 0.60653066, 0.54881164, 0.49658530, 0.44932896]
"""Centered Difference using Richardson's Extrapolation Method to find the derivative of a function at a point
    Args:
        FofX: a list of y-values
    Returns:
        f'(h): a list containing approximations of the derivative that increase in accuracy
"""
def cdfRichardson(FofX):
    fPrime =[]
    #D1(h)
    fPrime.append((FofX[8] - FofX[0])/(2*0.4))
    #D1(h/2)
    fPrime.append((FofX[6] - FofX[2])/(2*0.2))
    #D2(h)
    fPrime.append(fPrime[1] + (fPrime[1] - fPrime[0])/3)
    #D1(h/4)
    fPrime.append((FofX[5] - FofX[3])/(2*0.1))
    #D2(h/2)
    fPrime.append(fPrime[3] + (fPrime[3] - fPrime[1])/3)
    #D3(h/4)
    fPrime.append(fPrime[4] + (fPrime[4] - fPrime[2])/15)

    return fPrime
"""A method to show the contents of the list"""
def show(array):
    for i in array:
        print i
        


test = cdfRichardson(fx)
show(test)
