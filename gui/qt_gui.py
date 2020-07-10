
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QTextEdit, QGraphicsView
from PyQt5 import QtCore, QtGui, QtWidgets, uic
import sys 
import pyqtgraph as pg
import numpy as np
from deerlab import *
import matplotlib.pyplot as plt

def deerlabsim():
    # Simulation
    t = np.linspace(-0.5,5,300)
    r = np.linspace(2,6,300)
    P = dd_gauss2(r,[4, 0.5, 0.4, 5, 0.6, 0.6])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.00)
    return t,V,r,P

class UI(QMainWindow):
    def __init__(self):
        super(UI, self).__init__()
        uic.loadUi("qt_gui.ui", self)
 
        # find the widgets in the xml file
 
        self.graph1 = self.findChild(QGraphicsView, "timeDomain")
        self.graph2 = self.findChild(QGraphicsView, "distanceDomain")
        self.button = self.findChild(QPushButton, "loadButton")
        self.button.clicked.connect(self.clickedBtn)
        self.show()
 
 
 
    def clickedBtn(self):
        t,V,r,P = deerlabsim()
        self.graph1.plot(t,V)
        self.graph2.plot(t,P)
 
 
app = QApplication(sys.argv)
window = UI()
app.exec_()