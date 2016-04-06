fx = [1.00000000, 0.90483742, 0.81873075, 0.74081822, 0.67032005, 0.60653066, 0.54881164, 0.49658530, 0.44932896]

def cdfRichardson(FofX):
    fPrime =[]
    fPrime.append((FofX[8] - FofX[0])/(2*0.4))
    fPrime.append((FofX[6] - FofX[2])/(2*0.2))
    fPrime.append(fPrime[1] + (fPrime[1] - fPrime[0])/3)
    fPrime.append((FofX[5] - FofX[3])/(2*0.1))
    fPrime.append(fPrime[3] + (fPrime[3] - fPrime[1])/3)
    fPrime.append(fPrime[4] + (fPrime[4] - fPrime[2])/15)

    return fPrime

def show(array):
    for i in array:
        print i
        


test = cdfRichardson(fx)
show(test)
