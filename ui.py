# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'D:\ui.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_title(object):
    def setupUi(self, title):
        title.setObjectName("title")
        title.resize(713, 525)
        self.centralwidget = QtWidgets.QWidget(title)
        self.centralwidget.setObjectName("centralwidget")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(0, 0, 461, 61))
        font = QtGui.QFont()
        font.setFamily("微软雅黑")
        font.setPointSize(25)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(50, 70, 631, 381))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.lineEdit = QtWidgets.QLineEdit(self.tab)
        self.lineEdit.setGeometry(QtCore.QRect(10, 40, 113, 20))
        self.lineEdit.setObjectName("lineEdit")
        self.erfenbutton = QtWidgets.QPushButton(self.tab)
        self.erfenbutton.setGeometry(QtCore.QRect(150, 40, 75, 23))
        self.erfenbutton.setObjectName("erfenbutton")
        self.label_2 = QtWidgets.QLabel(self.tab)
        self.label_2.setGeometry(QtCore.QRect(10, 10, 201, 16))
        self.label_2.setObjectName("label_2")
        self.newtownbutton = QtWidgets.QPushButton(self.tab)
        self.newtownbutton.setGeometry(QtCore.QRect(240, 40, 75, 23))
        self.newtownbutton.setObjectName("newtownbutton")
        self.label_3 = QtWidgets.QLabel(self.tab)
        self.label_3.setGeometry(QtCore.QRect(10, 80, 54, 12))
        self.label_3.setObjectName("label_3")
        self.lineEdit_2 = QtWidgets.QLineEdit(self.tab)
        self.lineEdit_2.setGeometry(QtCore.QRect(10, 110, 113, 20))
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.label_4 = QtWidgets.QLabel(self.tab_2)
        self.label_4.setGeometry(QtCore.QRect(30, 10, 151, 20))
        self.label_4.setObjectName("label_4")
        self.line = QtWidgets.QFrame(self.tab_2)
        self.line.setGeometry(QtCore.QRect(22, 30, 21, 16))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.line_2 = QtWidgets.QFrame(self.tab_2)
        self.line_2.setGeometry(QtCore.QRect(13, 38, 20, 91))
        self.line_2.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.line_3 = QtWidgets.QFrame(self.tab_2)
        self.line_3.setGeometry(QtCore.QRect(22, 120, 21, 16))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.line_4 = QtWidgets.QFrame(self.tab_2)
        self.line_4.setGeometry(QtCore.QRect(162, 37, 20, 91))
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.line_5 = QtWidgets.QFrame(self.tab_2)
        self.line_5.setGeometry(QtCore.QRect(150, 30, 21, 16))
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.line_6 = QtWidgets.QFrame(self.tab_2)
        self.line_6.setGeometry(QtCore.QRect(150, 118, 21, 16))
        self.line_6.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.lineEdit_3 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_3.setGeometry(QtCore.QRect(60, 50, 21, 20))
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.lineEdit_4 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_4.setGeometry(QtCore.QRect(110, 50, 21, 20))
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.lineEdit_5 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_5.setGeometry(QtCore.QRect(60, 90, 21, 20))
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.lineEdit_6 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_6.setGeometry(QtCore.QRect(110, 90, 21, 20))
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.lineEdit_7 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_7.setGeometry(QtCore.QRect(220, 50, 21, 20))
        self.lineEdit_7.setObjectName("lineEdit_7")
        self.lineEdit_8 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_8.setGeometry(QtCore.QRect(220, 90, 21, 20))
        self.lineEdit_8.setObjectName("lineEdit_8")
        self.label_5 = QtWidgets.QLabel(self.tab_2)
        self.label_5.setGeometry(QtCore.QRect(210, 14, 54, 12))
        self.label_5.setObjectName("label_5")
        self.pushButton = QtWidgets.QPushButton(self.tab_2)
        self.pushButton.setGeometry(QtCore.QRect(30, 150, 75, 23))
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(self.tab_2)
        self.pushButton_2.setGeometry(QtCore.QRect(150, 150, 75, 23))
        self.pushButton_2.setObjectName("pushButton_2")
        self.label_6 = QtWidgets.QLabel(self.tab_2)
        self.label_6.setGeometry(QtCore.QRect(70, 200, 54, 12))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.tab_2)
        self.label_7.setGeometry(QtCore.QRect(210, 200, 54, 12))
        self.label_7.setObjectName("label_7")
        self.lineEdit_9 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_9.setGeometry(QtCore.QRect(10, 230, 113, 20))
        self.lineEdit_9.setObjectName("lineEdit_9")
        self.lineEdit_10 = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_10.setGeometry(QtCore.QRect(160, 230, 113, 20))
        self.lineEdit_10.setObjectName("lineEdit_10")
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.label_8 = QtWidgets.QLabel(self.tab_3)
        self.label_8.setGeometry(QtCore.QRect(20, 10, 391, 16))
        self.label_8.setObjectName("label_8")
        self.lineEdit_11 = QtWidgets.QLineEdit(self.tab_3)
        self.lineEdit_11.setGeometry(QtCore.QRect(20, 40, 113, 20))
        self.lineEdit_11.setObjectName("lineEdit_11")
        self.label_9 = QtWidgets.QLabel(self.tab_3)
        self.label_9.setGeometry(QtCore.QRect(20, 110, 54, 12))
        self.label_9.setObjectName("label_9")
        self.pushButton_3 = QtWidgets.QPushButton(self.tab_3)
        self.pushButton_3.setGeometry(QtCore.QRect(30, 73, 75, 23))
        self.pushButton_3.setObjectName("pushButton_3")
        self.textEdit = QtWidgets.QTextEdit(self.tab_3)
        self.textEdit.setGeometry(QtCore.QRect(20, 140, 291, 161))
        self.textEdit.setObjectName("textEdit")
        self.tabWidget.addTab(self.tab_3, "")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.lineEdit_12 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_12.setGeometry(QtCore.QRect(20, 20, 31, 31))
        self.lineEdit_12.setObjectName("lineEdit_12")
        self.lineEdit_13 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_13.setGeometry(QtCore.QRect(70, 20, 31, 31))
        self.lineEdit_13.setObjectName("lineEdit_13")
        self.lineEdit_14 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_14.setGeometry(QtCore.QRect(120, 20, 31, 31))
        self.lineEdit_14.setObjectName("lineEdit_14")
        self.lineEdit_15 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_15.setGeometry(QtCore.QRect(20, 70, 31, 31))
        self.lineEdit_15.setObjectName("lineEdit_15")
        self.lineEdit_16 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_16.setGeometry(QtCore.QRect(70, 70, 31, 31))
        self.lineEdit_16.setObjectName("lineEdit_16")
        self.lineEdit_17 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_17.setGeometry(QtCore.QRect(120, 70, 31, 31))
        self.lineEdit_17.setObjectName("lineEdit_17")
        self.lineEdit_18 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_18.setGeometry(QtCore.QRect(20, 120, 31, 31))
        self.lineEdit_18.setObjectName("lineEdit_18")
        self.lineEdit_19 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_19.setGeometry(QtCore.QRect(70, 120, 31, 31))
        self.lineEdit_19.setObjectName("lineEdit_19")
        self.lineEdit_20 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_20.setGeometry(QtCore.QRect(120, 120, 31, 31))
        self.lineEdit_20.setObjectName("lineEdit_20")
        self.pushButton_4 = QtWidgets.QPushButton(self.tab_4)
        self.pushButton_4.setGeometry(QtCore.QRect(50, 170, 75, 23))
        self.pushButton_4.setObjectName("pushButton_4")
        self.lineEdit_21 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_21.setGeometry(QtCore.QRect(30, 230, 113, 20))
        self.lineEdit_21.setObjectName("lineEdit_21")
        self.line_7 = QtWidgets.QFrame(self.tab_4)
        self.line_7.setGeometry(QtCore.QRect(10, 150, 21, 16))
        self.line_7.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")
        self.line_8 = QtWidgets.QFrame(self.tab_4)
        self.line_8.setGeometry(QtCore.QRect(1, 8, 20, 151))
        self.line_8.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_8.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_8.setObjectName("line_8")
        self.line_9 = QtWidgets.QFrame(self.tab_4)
        self.line_9.setGeometry(QtCore.QRect(10, 0, 21, 16))
        self.line_9.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_9.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_9.setObjectName("line_9")
        self.line_10 = QtWidgets.QFrame(self.tab_4)
        self.line_10.setGeometry(QtCore.QRect(140, 150, 21, 16))
        self.line_10.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_10.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_10.setObjectName("line_10")
        self.line_11 = QtWidgets.QFrame(self.tab_4)
        self.line_11.setGeometry(QtCore.QRect(150, 8, 20, 151))
        self.line_11.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_11.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_11.setObjectName("line_11")
        self.line_12 = QtWidgets.QFrame(self.tab_4)
        self.line_12.setGeometry(QtCore.QRect(140, 0, 21, 16))
        self.line_12.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_12.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_12.setObjectName("line_12")
        self.label_10 = QtWidgets.QLabel(self.tab_4)
        self.label_10.setGeometry(QtCore.QRect(70, 210, 54, 12))
        self.label_10.setObjectName("label_10")
        self.tabWidget.addTab(self.tab_4, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.lineEdit_22 = QtWidgets.QLineEdit(self.tab_5)
        self.lineEdit_22.setGeometry(QtCore.QRect(20, 40, 113, 20))
        self.lineEdit_22.setObjectName("lineEdit_22")
        self.label_11 = QtWidgets.QLabel(self.tab_5)
        self.label_11.setGeometry(QtCore.QRect(20, 10, 291, 16))
        self.label_11.setObjectName("label_11")
        self.pushButton_5 = QtWidgets.QPushButton(self.tab_5)
        self.pushButton_5.setGeometry(QtCore.QRect(30, 70, 75, 23))
        self.pushButton_5.setObjectName("pushButton_5")
        self.lineEdit_24 = QtWidgets.QLineEdit(self.tab_5)
        self.lineEdit_24.setGeometry(QtCore.QRect(20, 180, 113, 20))
        self.lineEdit_24.setObjectName("lineEdit_24")
        self.label_12 = QtWidgets.QLabel(self.tab_5)
        self.label_12.setGeometry(QtCore.QRect(20, 150, 241, 16))
        self.label_12.setObjectName("label_12")
        self.lineEdit_25 = QtWidgets.QLineEdit(self.tab_5)
        self.lineEdit_25.setGeometry(QtCore.QRect(150, 180, 113, 20))
        self.lineEdit_25.setObjectName("lineEdit_25")
        self.pushButton_6 = QtWidgets.QPushButton(self.tab_5)
        self.pushButton_6.setGeometry(QtCore.QRect(280, 180, 75, 23))
        self.pushButton_6.setObjectName("pushButton_6")
        self.pushButton_7 = QtWidgets.QPushButton(self.tab_5)
        self.pushButton_7.setGeometry(QtCore.QRect(360, 180, 75, 23))
        self.pushButton_7.setObjectName("pushButton_7")
        self.label_13 = QtWidgets.QLabel(self.tab_5)
        self.label_13.setGeometry(QtCore.QRect(20, 110, 54, 12))
        self.label_13.setObjectName("label_13")
        self.label_14 = QtWidgets.QLabel(self.tab_5)
        self.label_14.setGeometry(QtCore.QRect(20, 210, 54, 12))
        self.label_14.setObjectName("label_14")
        self.tabWidget.addTab(self.tab_5, "")
        self.tab_6 = QtWidgets.QWidget()
        self.tab_6.setObjectName("tab_6")
        self.label_15 = QtWidgets.QLabel(self.tab_6)
        self.label_15.setGeometry(QtCore.QRect(20, 10, 131, 16))
        self.label_15.setObjectName("label_15")
        self.lineEdit_23 = QtWidgets.QLineEdit(self.tab_6)
        self.lineEdit_23.setGeometry(QtCore.QRect(20, 40, 113, 20))
        self.lineEdit_23.setObjectName("lineEdit_23")
        self.pushButton_8 = QtWidgets.QPushButton(self.tab_6)
        self.pushButton_8.setGeometry(QtCore.QRect(20, 80, 75, 23))
        self.pushButton_8.setObjectName("pushButton_8")
        self.pushButton_9 = QtWidgets.QPushButton(self.tab_6)
        self.pushButton_9.setGeometry(QtCore.QRect(110, 80, 75, 23))
        self.pushButton_9.setObjectName("pushButton_9")
        self.label_16 = QtWidgets.QLabel(self.tab_6)
        self.label_16.setGeometry(QtCore.QRect(20, 120, 54, 12))
        self.label_16.setObjectName("label_16")
        self.tabWidget.addTab(self.tab_6, "")
        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(483, 460, 201, 20))
        self.label_17.setObjectName("label_17")
        title.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(title)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 713, 23))
        self.menubar.setObjectName("menubar")
        title.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(title)
        self.statusbar.setObjectName("statusbar")
        title.setStatusBar(self.statusbar)

        self.retranslateUi(title)
        self.tabWidget.setCurrentIndex(5)
        self.erfenbutton.clicked.connect(title.erfen)
        self.newtownbutton.clicked.connect(title.newton)
        self.pushButton.clicked.connect(title.sr)
        self.pushButton_2.clicked.connect(title.sdm)
        self.pushButton_3.clicked.connect(title.rk)
        self.pushButton_4.clicked.connect(title.pm)
        self.pushButton_5.clicked.connect(title.bsa)
        self.pushButton_6.clicked.connect(title.che)
        self.pushButton_7.clicked.connect(title.leg)
        self.pushButton_8.clicked.connect(title.bsa2)
        self.pushButton_9.clicked.connect(title.leg2)
        QtCore.QMetaObject.connectSlotsByName(title)

    def retranslateUi(self, title):
        _translate = QtCore.QCoreApplication.translate
        title.setWindowTitle(_translate("title", "MainWindow"))
        self.label.setText(_translate("title", "数值分析算法模拟系统"))
        self.erfenbutton.setText(_translate("title", "二分法"))
        self.label_2.setText(_translate("title", "请输入要计算的值"))
        self.newtownbutton.setText(_translate("title", "牛顿法"))
        self.label_3.setText(_translate("title", "答案"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("title", "迭代法求平方根"))
        self.label_4.setText(_translate("title", "请输入未知数系数矩阵A"))
        self.label_5.setText(_translate("title", "请输入b"))
        self.pushButton.setText(_translate("title", "共轭梯度法"))
        self.pushButton_2.setText(_translate("title", "最速下降法"))
        self.label_6.setText(_translate("title", "x"))
        self.label_7.setText(_translate("title", "y"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("title", "迭代法解线性方程组"))
        self.label_8.setText(_translate("title", "y\'=y-2*x/y 以空格为间隔，输入(a,b)和步长，例如0 1 0.2"))
        self.label_9.setText(_translate("title", "结果"))
        self.pushButton_3.setText(_translate("title", "四阶显式RK"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("title", "常微分方程"))
        self.pushButton_4.setText(_translate("title", "幂法"))
        self.label_10.setText(_translate("title", "答案"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), _translate("title", "矩阵特征值计算"))
        self.label_11.setText(_translate("title", "输入多项式次数"))
        self.pushButton_5.setText(_translate("title", "BSA"))
        self.label_12.setText(_translate("title", "输入多项式次数与辛普森公式间隔大小"))
        self.pushButton_6.setText(_translate("title", "切比雪夫"))
        self.pushButton_7.setText(_translate("title", "勒让德"))
        self.label_13.setText(_translate("title", "结果如图"))
        self.label_14.setText(_translate("title", "结果如图"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), _translate("title", "函数逼近"))
        self.label_15.setText(_translate("title", "输入多项式次数"))
        self.pushButton_8.setText(_translate("title", "BSA"))
        self.pushButton_9.setText(_translate("title", "Leg"))
        self.label_16.setText(_translate("title", "结果如图"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6), _translate("title", "曲线拟合"))
        self.label_17.setText(_translate("title", "Designed by 琚烈 潘恒 凌宇 贾亚光"))

