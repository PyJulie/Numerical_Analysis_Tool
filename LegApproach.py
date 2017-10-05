'''
勒让德正交多项式
'''

import math
import numpy as np
import matplotlib.pyplot as plt



def chengji(fx1, fx2, n1, n2):
    ans = []
    for n in range(n1 + n2 - 1):
        temp = 0
        for i in range(n1):
            if (n - i < 0):
                break
            elif (n - i >= n2):
                continue
            else:
                temp = temp + fx1[i] * fx2[n - i]
        ans.append(temp)
    return ans


def neiji(ans, a, b):
    ans.insert(0, 0)
    for n in range(len(ans)):
        if (n == 0):
            continue
        ans[n] /= float(n)

#    print ans
    computeB = [(b**i) * ans[i] for i in range(len(ans))]
    computeA = [(a**i) * ans[i] for i in range(len(ans))]
    ans = sum(computeB) - sum(computeA)
    return ans


def an(fn, n, a, b):  # n为fn的系数个数
    fnn = fn[:]
    fnn.insert(0, 0)
    fz = neiji(chengji(fnn, fn, n + 1, n), a, b)
    fm = neiji(chengji(fn, fn, n, n), a, b)
    ans = fz / fm
    return ans


def bn(fn, fnn, n, a, b):
    if (n == 1):
        return 0
    fz = neiji(chengji(fn, fn, n, n), a, b)
    fm = neiji(chengji(fnn, fnn, n - 1, n - 1), a, b)
    return fz / fm


def fn3(fn2, fn1, n3, a, b):  # n3为fn3的系数个数
    a = -1 * an(fn2, n3 - 1, a, b)
    b = bn(fn2, fn1, n3 - 1, a, b)

    L = chengji([a, 1], fn2, 2, n3 - 1)
    R = chengji([b], fn1, 1, n3 - 2)

    for i in range(len(R)):
        L[i] -= R[i]
    return L


def xishu(n):
    ans = float(math.factorial(2 * n)) / ((2**n) * ((math.factorial(n))**2))
    return ans




def f(x):
    return math.sin(2 * x**4 * math.pi)


def P(x, n, Legendra):
    ans = Legendra[n][0]
    tmp = x
    for i in range(1, n + 1):
        ans += Legendra[n][i] * tmp
        tmp *= x
    return ans


def IP(n,h,Legendra):  # Inner product
    x = [i for i in np.arange(-1, 1, h)]
    ans = 0
    for i in x:
        ans += h * (f(i) * P(i, n,Legendra) + 4 * f(i + h / 2) * P(i + h / 2, n,Legendra) +
                    f(i + h) * P(i + h, n,Legendra)) / 6
    return ans


def norm(n,h):
    x = [i for i in np.arange(-1, 1, h)]
    ans = 0
    for i in x:
        ans += h * (P(i, n,h) * P(i, n,h) + 4 * P(i + h, n,h) * P(i + h, n,h) + P(
            i + h, n,h) * P(i + h, n,h)) / 6
    return ans




def Sstar(x,aStar,h):
    ans = 0
    n = len(aStar)
    for i in range(n):
        ans += aStar[i] * P(x, i,h)
    return ans



