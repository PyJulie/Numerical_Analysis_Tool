from numpy import *
import numpy as np

def PowerMethod(A):
    n = len(A)
    v = [1]*n
    u = [1]*n
    exp = 1e-6
    Lambda = 10
    Uk = 1
    while(abs(Lambda - Uk)>exp):
        v = (mat(A)*mat(u).T).T.tolist()[0]
        Lambda = Uk
        Uk = max(v)
        u = [i/Uk for i in v]
    return Lambda


