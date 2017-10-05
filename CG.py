from numpy import *


def IP(a, b):
    if len(a) != len(b):
        return -1
    else:
        n = len(a)
    ans = 0
    for i in range(n):
        ans += a[i]*b[i]
    return ans

def CGmethod(A, b):
    exp = 1e-10
    n = len(b)
    x = [0 for i in range(n)]
    r = (mat(b).T - mat(A) * mat(x).T).T.tolist()[0]
    p = copy(r)
    step = 10
    x1 = [0]
    x2 = [0]
    y1 = [1]
    y2 = [2]
    while step > 0:
        step-=1
        tmp = IP((mat(A)*mat(p).T).T.tolist()[0], p)
        if tmp == 0:
            break
        alpha = IP(r, r)/tmp
        x = (mat(x)+alpha*mat(p)).tolist()[0]
        x1.append(x[0])
        x2.append(x[1])
        y1.append(1)
        y2.append(2)
        r_nx = (mat(r).T - alpha*mat(A)*mat(p).T).T.tolist()[0]
        beta = IP(r_nx, r_nx)/IP(r, r)
        r = r_nx
        for i in range(n):
            if r[i] > exp:
                break
        else:
            break
        p = (mat(r) + beta*mat(p)).tolist()[0]
    while step > 0:
        x1.append(x[0])
        x2.append(x[1])
        y1.append(1)
        y2.append(2)
        step-=1
    return x


