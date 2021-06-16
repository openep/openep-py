'''
GUI code for OpenEp
'''


from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyvistaqt as pyvqt

# import scipy.io as sio
# import pyqtgraph as pg
# import numpy as np
# import trimesh as tm

import sys
sys.path.append('../openep')

from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw



class App(QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'OpenEp Application'
        self.left = 10
        self.top = 10
        self.width = 840
        self.height = 480
        self.initUI()
        self.mesh_obj = None
        self.voltage_new = None


    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Defining the Vertical & Horizontal Box layout
        vbox = QFormLayout(self)
        hbox = QHBoxLayout(self)

        # Main Label
        mainLabel1 = QLabel('OpenEp Tool',self)
        mainLabel2 = QLabel( 'The Open Source solution for electrophysiology data analysis',self)

        vbox.addWidget(mainLabel1)
        vbox.addWidget(mainLabel2)



        # Button 1
        button1 = QPushButton('Load OpenEp Data', self)
        button1.setGeometry(200,150,100,40)
        a = button1.clicked.connect(self.on_click)


        # Button 2
        button2 = QPushButton('Plot 3D Mesh', self)
        button2.setGeometry(200,150,100,40)
        button2.clicked.connect(self.on_click2)

        # Adding buttons to the horizontal layout
        hbox.addWidget(button1)
        hbox.addWidget(button2)

        vbox.addRow(hbox)

        # Plot
        # self.plotter = pyvqt.BackgroundPlotter()
        vbox.addRow(QLabel('3-D Plot'))
        # vbox.addWidget(self.plotter)
        

        # vbox.addRow(pg.PlotWidget())

        # create a mesh object,
        # display the mesh

        # vbox hidden once data is loaded, self.vbox so you have access to it later on.







        # # Plot
        # # self.graphWidget = pg.PlotWidget()
        # self.graph = pg.PlotWidget()
        # # self.setCentralWidget(self.graph)
        # hour = [1,2,3,4,5,6,7,8,9,10]
        # temperature = [30,32,34,32,33,31,29,32,35,45]
        # print('hour\n',hour)
        # print('temperature\n',temperature)
        # self.graph.plot(hour, temperature)
        #

        self.show()

    def on_click(self):
        print('Please wait: Loading Data ... ')

        # Loading file from a Dialog Box
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)

        if dlg.exec_():
            filenames = dlg.selectedFiles()
            self.ep_case = openep_io.load_case(filenames[0])


    def on_click2(self):
        c = draw.DrawMap(self.ep_case,freeboundary_color='black',freeboundary_width=5)
        # import pyvista as pv
        # from pyvistaqt import MultiPlotter
        # plotter = MultiPlotter()
        # _ = plotter[0, 0].add_mesh(c)




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
