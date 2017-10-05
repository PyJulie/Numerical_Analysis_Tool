from numpy import *
import numpy as np
import matplotlib.pyplot as plt

def IP(a, b):
    if len(a) != len(b):
        return -1
    else:
        n = len(a)
    ans = 0
    for i in range(n):
        ans += a[i]*b[i]
    return ans

def steepest(A, b):
    x = [0 for i in range(len(b))]
    step = 10
    Xaxis = [i for i in range(step+1)]
    x1 = [0]
    x2 = [0]
    y1 = [1]*(step+1)
    y2 = [2]*(step+1)
    while step > 0:
        r = (mat(b).T - mat(A)*mat(x).T).T.tolist()[0]
        alpha = IP(r, r)/IP((mat(A)*mat(r).T).T.tolist()[0], r)
        x = (mat(x)+alpha*mat(r)).tolist()[0]
        x1.append(x[0])
        x2.append(x[1])

        step-=1
    return x

