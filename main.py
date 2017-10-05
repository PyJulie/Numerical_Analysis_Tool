import math
from numpy import *
from math import sqrt
from ui import Ui_title
from PyQt5 import QtWidgets,QtGui
from CG import *
from rk import *
from power_method import *
from steepest_descent_method import *
import BSA
import CheApproach
import LegApproach
import LegCurveFitting
import BSACurveFitting
exp = 1e-5

def f(x):
    C = 2
    return x*x-C

def dichotomy(num):
    x = sqrt(num)
    y = num / 2.0
    low = 0.0
    up = num * 1.0
    count = 1
    while abs(y - x) > 0.00000001:
        print
        count, y
        count += 1
        if (y * y > num):
            up = y
            y = low + (y - low) / 2
        else:
            low = y
            y = up - (up - y) / 2
    return y

def sqrt_newton(num):
    x=sqrt(num)
    y=num/2.0
    count=1
    while abs(y-x)>0.00000001:
        count+=1
        y=((y*1.0)+(1.0*num)/y)/2.0000
    return y

class mywindow(QtWidgets.QMainWindow,Ui_title):
    def __init__(self):
        super(mywindow,self).__init__()
        self.setupUi(self)


    def erfen(self):
        answer = dichotomy(int(self.lineEdit.text()))
        answer_str = str(answer)
        self.lineEdit_2.setText(answer_str)

    def newton(self):
        answer = sqrt_newton(int(self.lineEdit.text()))
        answer_str = str(answer)
        self.lineEdit_2.setText(answer_str)

    def sr(self):
        A = [[int(self.lineEdit_3.text()),int(self.lineEdit_4.text())],[int(self.lineEdit_5.text()),int(self.lineEdit_6.text())]]
        b = [int(self.lineEdit_7.text()),int(self.lineEdit_8.text())]
        x = CGmethod(A,b)
        self.lineEdit_9.setText(str(x[0]))
        self.lineEdit_10.setText(str(x[1]))
    def sdm(self):
        A = [[int(self.lineEdit_3.text()), int(self.lineEdit_4.text())],
             [int(self.lineEdit_5.text()), int(self.lineEdit_6.text())]]
        b = [int(self.lineEdit_7.text()), int(self.lineEdit_8.text())]
        x = steepest(A, b)
        self.lineEdit_9.setText(str(x[0]))
        self.lineEdit_10.setText(str(x[1]))

    def rk(self):
        s = self.lineEdit_11.text().split()
        a = []
        for i in s:
            a.append(double(i))
        answer = RK4(a[0],a[1],a[2])
        j=a[0]
        for answersingle in answer:
            self.textEdit.append(str(j)+" "+str(answersingle)+'\n')
            j+=a[2]

    def pm(self):
        A=[[double(self.lineEdit_12.text()),double(self.lineEdit_13.text()),double(self.lineEdit_14.text())],
           [double(self.lineEdit_15.text()),double(self.lineEdit_16.text()),double(self.lineEdit_17.text())],
           [double(self.lineEdit_18.text()),double(self.lineEdit_19.text()),double(self.lineEdit_20.text())]]
        answer = PowerMethod(A)
        self.lineEdit_21.setText(str(answer))

    def bsa(self):
        n = int(self.lineEdit_22.text())
        n += 1
        d = [BSA.dn(i) for i in range(n)]
        Hilbet = [[1 / (i + 1 + j) for i in range(n)] for j in range(n)]
        BSA.Gauss(Hilbet, d)
        x = [i for i in np.arange(0, 1, 0.01)]
        fx = [BSA.f(i, 0) for i in x]
        y = [BSA.Sstar(d, i, n) for i in x]
        err = [y[i] - fx[i] for i in range(len(fx))]

        plt.plot(x, fx, 'r')
        plt.plot(x, y, 'g')
        plt.plot(x, err, 'grey')
        plt.legend(['f(x)', 'S(x)', 'error'])
        plt.show()

    def che(self):
        power = int(self.lineEdit_24.text())
        h = double(self.lineEdit_25.text())
        power += 1

        T = [[1], [0, 1]]

        for i in range(2, power, 1):
            temp = T[i - 1][:]
            T.append(CheApproach.qiemixuefu(temp, T[i - 2]))
        aStar = [CheApproach.Cstar(i,h) for i in range(power)]

        x = [i for i in np.arange(-1, 1, 0.01)]
        fx = [CheApproach.f(i) for i in x]
        y = [CheApproach.Sstar(i,aStar,T) for i in x]
        err = [y[i] - fx[i] for i in range(len(fx))]

        plt.plot(x, fx, 'r')
        plt.plot(x, y, 'g')
        plt.plot(x, err, 'grey')
        plt.legend(['f(x)', 'S(x)', 'error'])
        plt.show()

    def leg(self):
        power = int(self.lineEdit_24.text())
        h = double(self.lineEdit_25.text())
        power += 1
        f_1 = [0]
        f0 = [1]
        f1 = LegApproach.fn3(f0, f_1, 2, -1, 1)
        f2 = LegApproach.fn3(f1, f0, 3, -1, 1)
        f3 = LegApproach.fn3(f2, f1, 4, -1, 1)
        Legendra = [f0, f1]  # 结果保存在x里面 从P0开始
        for i in range(1, power, 1):  # 生成多项式的个数
            Legendra.append(LegApproach.fn3(Legendra[i], Legendra[i - 1], i + 2, -1, 1))

        for i in range(len(Legendra)):
            for j in range(len(Legendra[i])):
                Legendra[i][j] *= LegApproach.xishu(i)
        aStar = [LegApproach.IP(i,h,Legendra) / LegApproach.norm(i,h) for i in range(power)]

        x = [i for i in np.arange(-1, 1, 0.01)]
        fx = [LegApproach.f(i) for i in x]
        y = [LegApproach.Sstar(i,aStar,h) for i in x]
        err = [y[i] - fx[i] for i in range(len(fx))]

        plt.plot(x, fx, 'r')
        plt.plot(x, y, 'g')
        plt.plot(x, err, 'grey')
        plt.legend(['f(x)', 'S(x)', 'error'])
        plt.show()

    def bsa2(self):
        n = int(self.lineEdit_23.text())
        x = [1, 2, 3]  # 输入x值
        yi = [1, 0.5, 0.33]  # 输入y值
        n += 1
        d = [BSACurveFitting.dn(i,x,yi) for i in range(n)]
        Hilbet = [[BSACurveFitting.IP(i, j,x) for i in range(n)] for j in range(n)]
        BSACurveFitting.Gauss(Hilbet, d)
        y = [BSACurveFitting.Sstar(d, i, n) for i in x]

        plt.plot(x, y, 'g')
        plt.plot(x, yi, 'y')
        plt.legend(['S(x)', 'yi'])
        plt.show()

    def leg2(self):
        power = int(self.lineEdit_23.text())
        xi = [i for i in range(100)]  # 输入x值
        yi = [i * i for i in range(100)]  # 输入y值

        power += 1
        x = xi
        maxXi = max(x)
        minXi = min(x)

        for i in range(len(x)):
            x[i] = x[i] / (maxXi - minXi) * 2 - 1 - minXi

        power += 1
        x = xi
        maxXi = max(x)
        minXi = min(x)

        f_1 = [0]
        f0 = [1]
        f1 = LegCurveFitting.fn3(f0, f_1, 2, -1, 1)
        f2 = LegCurveFitting.fn3(f1, f0, 3, -1, 1)
        f3 = LegCurveFitting.fn3(f2, f1, 4, -1, 1)
        Legendra = [f0, f1]  # 结果保存在x里面 从P0开始
        for i in range(1, power, 1):  # 生成多项式的个数
            Legendra.append(LegCurveFitting.fn3(Legendra[i], Legendra[i - 1], i + 2, -1, 1))

        for i in range(len(Legendra)):
            for j in range(len(Legendra[i])):
                Legendra[i][j] *= LegCurveFitting.xishu(i)
        aStar = [LegCurveFitting.IP(i,x,yi,Legendra) / LegCurveFitting.norm(i,yi,Legendra,x) for i in range(power)]
        y = [LegCurveFitting.Sstar(i,aStar,Legendra) for i in x]

        plt.plot(xi, y, 'g')
        plt.plot(xi, yi, 'y')
        plt.legend(['S(x)', 'yi'])
        plt.show()

if __name__=="__main__":
    import sys
    app=QtWidgets.QApplication(sys.argv)
    myshow = mywindow()
    myshow.show()
    sys.exit(app.exec())