import math

def f(x, y):
    return y-2*x/y

def Y(x):
    return math.sqrt(1+2*x)

def RK4(l, r, h):
    n = int((r-l)/h+1)
    x = [l+i*h for i in range(n)]
    y = [1]*n
    realy = y
    for i in range(n-1):
        K1 = f(x[i], y[i])
        K2 = f(x[i]+h/2, y[i]+K1*h/2)
        K3 = f(x[i]+h/2, y[i]+K2*h/2)
        K4 = f(x[i]+h, y[i]+K3*h)
        y[i+1] = y[i]+(K1+2*K2+2*K3+K4)*h/6
        realy[i+1] = Y(x[i+1])
        print(y[i], realy[i])

    return y
