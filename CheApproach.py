'''
切比雪夫正交多项式
'''

import math
import numpy as np
import matplotlib.pyplot as plt




def qiemixuefu(Tn, Tn_1):  # compute T(n+1)
    Tn.insert(0, 0)
    for i in range(len(Tn)):
        Tn[i] *= 2
    for i in range(len(Tn_1)):
        Tn[i] -= Tn_1[i]
    return Tn





def f(x):
    return math.sin(2 * x**4 * math.pi)


def P(x, n,T):
    ans = 0
    tmp = 1
    for i in range(len(T[n])):
        ans += T[n][i] * tmp
        tmp *= x
    return ans


def fc(x, k):
    return f(math.cos(x)) * math.cos(k * x)


def Cstar(k,h):
    x = [i for i in np.arange(0, math.pi, h)]
    ans = 0
    for i in x:
        ans += h * (fc(i, k) + 4 * fc(i + h / 2, k) + fc(i + h, k)) / 6
    return 2 / math.pi * ans





def Sstar(x,aStar,T):
    ans = aStar[0] / 2
    n = len(aStar)
    for i in range(1, n):
        ans += aStar[i] * P(x, i,T)
    return ans



