import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QApplication, QWidget

class Canvas(FigureCanvas):
    def __init__(self, parent):
        fig, self.ax = plt.subplots(figsize=(5, 4), dpi=200)
        super().__init__(fig)
        self.setParent(parent)

        """ 
        Matplotlib Script
        """
        t = np.arange(0.0, 2.0, 0.01)
        s = 1 + np.sin(2 * np.pi * t)
        
        self.ax.plot(t, s)

        self.ax.set(xlabel='time (s)', ylabel='voltage (mV)',
               title='About as simple as it gets, folks')
        self.ax.grid()

class AppDemo(QWidget)      :
    def __init__(self):
        super().__init__()
        self.resize(1600, 800)

        chart = Canvas(self)

app = QApplication(sys.argv)        
demo = AppDemo()
demo.show()
sys.exit(app.exec_())






# import sys
# from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QPushButton, QAction
# from PyQt5.QtGui import QIcon
# from PyQt5.QtCore import pyqtSlot

# class App(QMainWindow):
    
#     def __init__(self):
#         super().__init__()
#         self.title = 'PyQt5 menu - pythonspot.com'
#         self.left = 10
#         self.top = 10
#         self.width = 640
#         self.height = 400
#         self.initUI()
    
#     def initUI(self):
#         self.setWindowTitle(self.title)
#         self.setGeometry(self.left, self.top, self.width, self.height)
        
#         mainMenu = self.menuBar()
#         fileMenu = mainMenu.addMenu('File')
#         editMenu = mainMenu.addMenu('Edit')
#         viewMenu = mainMenu.addMenu('View')
#         searchMenu = mainMenu.addMenu('Search')
#         toolsMenu = mainMenu.addMenu('Tools')
#         helpMenu = mainMenu.addMenu('Help')
        
#         exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
#         exitButton.setShortcut('Ctrl+Q')
#         exitButton.setStatusTip('Exit application')
#         exitButton.triggered.connect(self.close)
#         fileMenu.addAction(exitButton)
    
#         self.show()
    
# if __name__ == '__main__':
#     app = QApplication(sys.argv)
#     ex = App()
#     sys.exit(app.exec_())

















# from PyQt5.QtCore import *
# from PyQt5.QtWidgets import *
# import sys

# class Window(QWidget):
#     def __init__(self):
#         QWidget.__init__(self)
#         layout = QGridLayout()
#         self.setLayout(layout)

#         # create menu
#         menubar = QMenuBar()
#         layout.addWidget(menubar, 0, 0)
#         actionFile = menubar.addMenu("File")
#         actionFile.addAction("New")
#         actionFile.addAction("Open")
#         actionFile.addAction("Save")
#         actionFile.addSeparator()
#         actionFile.addAction("Quit")
#         menubar.addMenu("Edit")
#         menubar.addMenu("View")
#         menubar.addMenu("Help")

#         # add textbox
#         tbox = QPlainTextEdit()
#         layout.addWidget(tbox, 1, 0)

        
# app = QApplication(sys.argv)
# screen = Window()
# screen.show()
# sys.exit(app.exec_())







# '''
# GUI code for OpenEp
# '''


# from PyQt5.QtGui import *
# from PyQt5.QtWidgets import *
# from PyQt5.QtCore import *
# import vtk
# from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
# import scipy.io as sio
# import pyqtgraph as pg
# import numpy as np
# import trimesh as tm
# from examples.main import visualise_vtk

# import sys
# sys.path.append('../openep')

# class App(QWidget):

#     def __init__(self):
#         super().__init__()
#         self.title = 'OpenEp Application'
#         self.left = 10
#         self.top = 10
#         self.width = 840
#         self.height = 480
#         self.initUI()
#         self.mesh_obj = None
#         self.voltage_new = None


#     def initUI(self):
#         self.setWindowTitle(self.title)
#         self.setGeometry(self.left, self.top, self.width, self.height)

#         # Defining the Vertical & Horizontal Box layout
#         self.vbox = QFormLayout(self)
#         hbox = QHBoxLayout(self)

#         # Main Label
#         mainLabel1 = QLabel('OpenEp Tool',self)
#         mainLabel2 = QLabel( 'The Open Source solution for electrophysiology data analysis',self)

#         self.vbox.addWidget(mainLabel1)
#         self.vbox.addWidget(mainLabel2)


#         # Button 1
#         button1 = QPushButton('Load OpenEp Data', self)
#         button1.setGeometry(200,150,100,40)
#         a = button1.clicked.connect(self.on_click)


#         # Button 2
#         button2 = QPushButton('Plot Voltage Map', self)
#         button2.setGeometry(200,150,100,40)
#         button2.clicked.connect(self.on_click2)

#         # Adding buttons to the horizontal layout
#         hbox.addWidget(button1)
#         hbox.addWidget(button2)

#         self.vbox.addRow(hbox)

#         # Plot
#         # self.vbox.addRow(QLabel('Plot'))
#         # self.vbox.addRow(pg.PlotWidget())

#         # create a mesh object,
#         # display the mesh

#         # vbox hidden once data is loaded, self.vbox so you have access to it later on.







#         # # Plot
#         # # self.graphWidget = pg.PlotWidget()
#         # self.graph = pg.PlotWidget()
#         # # self.setCentralWidget(self.graph)
#         # hour = [1,2,3,4,5,6,7,8,9,10]
#         # temperature = [30,32,34,32,33,31,29,32,35,45]
#         # print('hour\n',hour)
#         # print('temperature\n',temperature)
#         # self.graph.plot(hour, temperature)
#         #

#         self.show()

#     def visualise_surface(self, data_path):
#         frame = QFrame()
#         vtkWidget = QVTKRenderWindowInteractor(frame)
#         self.vbox.addWidget(vtkWidget)
#         vtkWidget.resize(vtkWidget.sizeHint())

#         ren = vtk.vtkRenderer()
#         vtkWidget.GetRenderWindow().AddRenderer(ren)
#         iren = vtkWidget.GetRenderWindow().GetInteractor()
#         style = vtk.vtkInteractorStyleTrackballCamera()
#         iren.SetInteractorStyle(style)

#         visualise_vtk(data_path, renderer=ren)

#         self.show()
#         iren.Initialize()
#         iren.Start()


#     def on_click(self):
#         print('Loading Data ... ')

#         # Loading file from a Dialog Box
#         dlg = QFileDialog()
#         dlg.setFileMode(QFileDialog.AnyFile)

#         if dlg.exec_():
#             filenames = dlg.selectedFiles()
#             self.visualise_surface(filenames[0])








#     def on_click2(self):
#         self.mesh_obj
#         self.voltage_new

#         # have an if condition if its none, do not run visualisation or have an
#         # advisory in the GUI.



#         pass
#         # print('Loading Voltage Map ... ')
#         # hour = [1,2,3,4,5,6,7,8,9,10]
#         # temperature = [30,32,34,32,33,31,29,32,35,45]
#         # print('hour\n',hour)
#         # print('temperature\n',temperature)
#         # self.graph.plot(hour, temperature)

#     # def load_data(filename):
#     #     main_file = sio.loadmat(filenames[0])
#     #     data_tri_X = main_file['userdata']['surface_triRep_X'][0][0]
#     #     t = main_file['userdata']['surface_triRep_Triangulation'][0][0] - 1
#     #     data_act_bip = main_file['userdata']['surface'][0][0]['act_bip'][0][0]
#     #     voltage_data = data_act_bip[:,1]
#     #     x = data_tri_X[:,0]
#     #     y = data_tri_X[:,1]
#     #     z = data_tri_X[:,2]
#     #     print('Nodes\n',data_tri_X)
#     #     print('tri_data\n',t)
#     #     print('Data Successfully Loaded')
#     #     mesh_obj = tm.Trimesh(vertices=data_tri_X,faces=t)
#     #     voltage_new = np.asarray(list(map(lambda x:0 if np.isnan(x) else x, voltage_data)))
#     #     return data_tri_X, data_tri_T, voltage_new



# if __name__ == '__main__':
#     app = QApplication(sys.argv)
#     ex = App()
#     sys.exit(app.exec_())
