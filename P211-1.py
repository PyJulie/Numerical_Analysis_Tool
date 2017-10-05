from numpy import *
import matplotlib.pyplot as plt

def jacobi(A, b):
    n = len(b)
    x = [0]*n
    for step in range(30):
        for i in range(n):
            tmp = 0
            for j in range(n):
                if i != j:
                    tmp += A[i][j]*x[j]
            x[i] = (b[i]-tmp)/A[i][i]
    return x

def SOR(A, b, w):
    n = len(b)
    x = [0]*n
    for step in range(30):
        for i in range(n):
            deltaX = b[i]
            for j in range(n):
                deltaX -= A[i][j]*x[j]
            deltaX /= A[i][i]
            deltaX *= w
            x[i] = x[i] + deltaX
    return x

n = 6
Hilbet = [[1 / (i + 1 + j) for i in range(n)] for j in range(n)]
x = [1]*n
Xaxis = [i for i in range(n)]
b = (mat(Hilbet)*mat(x).T).T.tolist()[0]
print(x)
xj6 = jacobi(Hilbet, b)
plt.plot(Xaxis, xj6, 'b')
'''
xs1 = SOR(Hilbet, b, 1)
xs125 = SOR(Hilbet, b, 1.25)
xs15 = SOR(Hilbet, b, 1.5)
print(xs1, '\n', xs125, '\n', xs15)
plt.plot(Xaxis, x, 'r')
plt.plot(Xaxis, xs1, 'g')
plt.plot(Xaxis, xs125, 'blue')
plt.plot(Xaxis, xs15, 'yellow')
plt.legend(['x', 'SOR:w = 1', 'SOR:w = 1.25', 'SOR:w = 1.5'])
plt.title("SOR n = 6")
plt.show()
'''

n = 8
Hilbet = [[1 / (i + 1 + j) for i in range(n)] for j in range(n)]
x = [1]*n
Xaxis = [i for i in range(n)]
b = (mat(Hilbet)*mat(x).T).T.tolist()[0]
print(x)
xj8 = jacobi(Hilbet, b)
plt.plot(Xaxis, xj8, 'g')
'''
xs1 = SOR(Hilbet, b, 1)
xs125 = SOR(Hilbet, b, 1.25)
xs15 = SOR(Hilbet, b, 1.5)
print(xs1, '\n', xs125, '\n', xs15)
plt.plot(Xaxis, x, 'r')
plt.plot(Xaxis, xs1, 'g')
plt.plot(Xaxis, xs125, 'blue')
plt.plot(Xaxis, xs15, 'yellow')
plt.legend(['x', 'SOR:w = 1', 'SOR:w = 1.25', 'SOR:w = 1.5'])
plt.title("SOR n = 8")
plt.show()
'''

n = 10
Hilbet = [[1 / (i + 1 + j) for i in range(n)] for j in range(n)]
x = [1]*n
Xaxis = [i for i in range(n)]
b = (mat(Hilbet)*mat(x).T).T.tolist()[0]
print(x)
xj10 = jacobi(Hilbet, b)
plt.plot(Xaxis, xj10, 'y')
plt.plot(Xaxis, x, 'r')
plt.legend(['jac:n = 6', 'jac:n = 8', 'jac:n = 10', 'x'])
plt.show()
'''
xs1 = SOR(Hilbet, b, 1)
xs125 = SOR(Hilbet, b, 1.25)
xs15 = SOR(Hilbet, b, 1.5)
print(xs1, '\n', xs125, '\n', xs15)
plt.plot(Xaxis, x, 'r')
plt.plot(Xaxis, xs1, 'g')
plt.plot(Xaxis, xs125, 'blue')
plt.plot(Xaxis, xs15, 'yellow')
plt.legend(['x', 'SOR:w = 1', 'SOR:w = 1.25', 'SOR:w = 1.5'])
plt.title("SOR n = 10")
plt.show()
'''