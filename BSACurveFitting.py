'''
{1,x,x^2,...,x^n}作基
'''

import matplotlib.pyplot as plt





def Gauss(A, x):
    eps = 1e-9
    k = 0
    equ = len(A)
    var = len(A[0])
    for col in range(var):
        max_r = k
        for i in range(k + 1, equ):
            if abs(A[i][col]) > abs(A[max_r][col]):
                max_r = i
        if abs(A[max_r][col]) < eps:
            return 0
        if max_r != k:
            for j in range(col, var):
                A[k][j], A[max_r][j] = A[max_r][j], A[k][j]
            x[k], x[max_r] = x[max_r], x[k]
        x[k] /= A[k][col]
        for j in range(col + 1, var):
            A[k][j] /= A[k][col]
        A[k][col] = 1
        for i in range(equ):
            if i != k and A[i][k] != 0:
                x[i] -= x[k] * A[i][k]
                for j in range(col + 1, var):
                    A[i][j] -= A[k][j] * A[i][col]
                A[i][col] = 0
        k += 1
        if k >= equ:
            break


def IP(i, j,x):
    ans = 0
    for k in x:
        ans += k**(i + j)
    return ans


def dn(n,x,yi):
    ans = 0
    for i in range(len(x)):
        ans += yi[i] * x[i]**n
    return ans


def Sstar(s, X, n):
    ans = s[0]
    tmp = X
    for i in range(1, n):
        ans += s[i] * tmp
        tmp *= X
    return ans


