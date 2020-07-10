from PyQt5.QtWidgets import QLineEdit,QAction, QMainWindow, QApplication, QPushButton, QTextEdit, QGraphicsView, QFileDialog,QTableWidgetItem
from PyQt5 import QtCore, QtGui, QtWidgets, uic
import sys 
import pyqtgraph as pg
import numpy as np
from deerlab import *
from deerlab_default import deerlab_default
import os
dir = os.path.dirname(__file__)

pg.setConfigOption('background','w')
pg.setConfigOptions(antialias=True)

class UI(QMainWindow):
    def __init__(self):
    #===============================================
        super(UI, self).__init__()
        uic.loadUi(os.path.join(dir, 'qt_gui.ui'), self)
 
        # set the title 
        self.setWindowTitle('DeerLab - Standalone GUI prototype') 

        # find the widgets in the xml file
        self.graph1 = self.findChild(QGraphicsView, "timeDomain")
        self.graph2 = self.findChild(QGraphicsView, "distanceDomain")
        self.loadaction = self.findChild(QAction, "actionLoad_file")
        self.exitaction = self.findChild(QAction, "actionExit")
        self.statusbar = self.findChild(QLineEdit, "statusBar")
        self.loadaction.triggered.connect(self.clickedBtn)

        self.show()
    #===============================================

    def clickedBtn(self):
    #===============================================
        # Load file via UI window
        fileName, _ = QFileDialog.getOpenFileNames()
        t, V, _ = deerload(fileName[0])
        # Remove singleton dimensions
        t = np.squeeze(t)
        V = np.real(np.squeeze(V))

        self.statusbar.setText('Processing... please wait')
        self.statusbar.update()
        QApplication.processEvents()
        
        # Run default DeerLab analysis
        V,Vfit,Pfit,r,P95,P50,parfit,parci = deerlab_default(t,V)
        
        self.graph1.clear()
        self.graph2.clear()

        # Update time-domain plot
        self.graph1.plot(t, V, symbolBrush=(80,80,80), symbol='o', symbolSize=3)
        self.graph1.plot(t, Vfit, pen=pg.mkPen(color='b'))
        self.graph1.setLabel('bottom','Time [us]')
        self.graph1.setLabel('left','V(t)')
        self.graph1.showGrid(x=True, y=True)
        self.graph2.setXRange(min(t),max(t),padding=0)
        self.graph2.setYRange(min(V),max(V),padding=0)

        # Update distance-domain plot
        self.graph2.plot(r,Pfit,pen = pg.mkPen(color='b'))
        phigh = pg.PlotCurveItem(r, P95[:,0])           
        plow = pg.PlotCurveItem(r, P95[:,1])   
        pfill = pg.FillBetweenItem(phigh, plow, brush= (38, 64, 171, 100))
        self.graph2.addItem(pfill)
        phigh = pg.PlotCurveItem(r, P50[:,0])           
        plow = pg.PlotCurveItem(r, P50[:,1])   
        pfill = pg.FillBetweenItem(phigh, plow, brush= (38, 64, 171, 100))
        self.graph2.addItem(pfill)
        self.graph2.setLabel('bottom','Distance [nm]')
        self.graph2.setLabel('left','P(r)')
        font=QtGui.QFont()
        font.setPixelSize(40)
        self.graph2.getAxis("bottom").setStyle = font
        self.graph2.setXRange(min(r),max(r),padding=0)
        self.graph2.setYRange(0,max(P95[:,1]),padding=0)
        self.graph2.showGrid(x=True, y=True)

        paramnames = ('Modulation depth','Background decay rate')
        units = ('','us-1')
        self.FitParamTable.setColumnCount(5)  
        self.FitParamTable.setRowCount(len(parfit))  
        for i in range(len(parfit)):
            self.FitParamTable.setItem(i,0,QTableWidgetItem(paramnames[i])) 
            self.FitParamTable.setItem(i,1,QTableWidgetItem(units[i])) 
            self.FitParamTable.setItem(i,2,QTableWidgetItem(np.array2string(parfit[i], formatter={'float_kind':lambda x: "%.2f" % x}))) 
            self.FitParamTable.setItem(i,3,QTableWidgetItem(np.array2string(parci[i,0], formatter={'float_kind':lambda x: "%.2f" % x}))) 
            self.FitParamTable.setItem(i,4,QTableWidgetItem(np.array2string(parci[i,1], formatter={'float_kind':lambda x: "%.2f" % x}))) 
    
        rmsd = np.sqrt(np.sum((V - Vfit)**2))
        Ndof = len(V) - len(paramnames)
        chi2 = rmsd**2/Ndof/np.std(V-Vfit)**2
        self.GoodnessOfFitTable.setItem(0,0,QTableWidgetItem(np.array2string(rmsd, formatter={'float_kind':lambda x: "%.3f" % x}))) 
        self.GoodnessOfFitTable.setItem(1,0,QTableWidgetItem(np.array2string(chi2, formatter={'float_kind':lambda x: "%.3f" % x}))) 

        self.statusbar.setText('Done')
        self.statusbar.update()
    #===============================================


app = QApplication(sys.argv)
window = UI()
app.exec_()