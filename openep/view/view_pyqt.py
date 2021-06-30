'''
GUI code for OpenEp
'''

from PyQt5 import QtWidgets as qtw
import pyvista as pv
from pyvistaqt import QtInteractor, MainWindow


import sys
sys.path.append('../openep')

from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw


class OpenEpGUI(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.title = 'OpenEp Application'
        self.left = 10
        self.top = 10
        self.width = 840
        self.height = 480
        self.initUI()
        self.thresholds = False
        self.minval = 0
        self.maxval = 2


    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # LAYOUTS
        self.mainLayout = qtw.QVBoxLayout(self)
        self.labelLayout = qtw.QFormLayout(self)
        self.buttonLayout = qtw.QHBoxLayout(self)

        
        mainLabel1 = qtw.QLabel('OpenEp Tool',self)
        self.labelLayout.addWidget(mainLabel1)
        mainLabel2 = qtw.QLabel( 'The Open Source solution for electrophysiology data analysis',self)
        self.labelLayout.addWidget(mainLabel2)



        # Button 1
        button1 = qtw.QPushButton('Load OpenEp Data', self)
        button1.setGeometry(200,150,100,40)
        button1.clicked.connect(self.on_click)


        # Button 2
        button2 = qtw.QPushButton('Plot Voltage Map', self)
        button2.setGeometry(200,150,100,40)
        button2.clicked.connect(self.on_click2)

        # Adding buttons to the horizontal layout
        self.buttonLayout.addWidget(button1)
        self.buttonLayout.addWidget(button2)

        self.labelLayout.addRow(self.buttonLayout)

        # Plot
        self.frame = qtw.QFrame()
        self.plotLayout = qtw.QHBoxLayout()
        self.plotter = QtInteractor(self.frame)
        self.plotLayout.addWidget(self.plotter.interactor)

        # QDock Widget
        
        self.plotter1 = QtInteractor(self.frame)
        
        self.dock_plot = qtw.QDockWidget("Dockable", self)
        self.dock_plot.setFloating(False)
        self.dock_plot.setWidget(self.plotter1)

        self.plotLayout.addWidget(self.dock_plot)
        
        



        # VoltThresholds limit
        self.limitLayout = qtw.QFormLayout(self)
        self.lowerlimit = qtw.QLineEdit()
        self.upperlimit = qtw.QLineEdit()
        self.lowerlimit.setText(str(0))
        self.upperlimit.setText(str(2))
        self.limitLayout.addRow('Voltage Threshold - Lower', self.lowerlimit)
        self.limitLayout.addRow('Voltage Threshold - Upper', self.upperlimit)

        button3 = qtw.QPushButton('Set Voltage Thresholds', self)
        button3.setGeometry(200,150,100,40)
        button3.clicked.connect(self.on_click3)
        self.limitLayout.addRow(button3)


        # Nesting Layouts and setting the main layout
        self.mainLayout.addLayout(self.labelLayout)
        self.mainLayout.addLayout(self.plotLayout)
        self.mainLayout.addLayout(self.limitLayout)
        self.setLayout(self.mainLayout)







    def on_click(self):
        print('Please wait: Loading Data ... ')
        # Loading file from a Dialog Box
        dlg = qtw.QFileDialog()
        dlg.setFileMode(qtw.QFileDialog.AnyFile)

        if dlg.exec_():
            filenames = dlg.selectedFiles()
            self.ep_case = openep_io.load_case(filenames[0])


    def on_click2(self):
        surf = draw.DrawMap(self.ep_case,
                            cmap='jet_r',
                            freeboundary_color='black',
                            freeboundary_width=5,
                            minval=self.minval,
                            maxval=self.maxval,
                            volt_below_color='brown', 
                            volt_above_color='magenta', 
                            nan_color='gray',
                            plot=False)

        self.mesh = surf['pyvista-mesh']
        self.volt = surf['volt']
        self.nan_color = surf['nan_color']
        self.minval = surf['minval']
        self.maxval = surf['maxval']
        self.cmap = surf['cmap']
        self.below_color = surf['volt_below_color']
        self.above_color = surf['volt_above_color']
        self.freeboundary_points = draw.getAnatomicalStructures(self.ep_case,plot=False)

        self.sargs = dict(interactive=True, 
                          n_labels=2,
                          label_font_size=18,
                          below_label='  ',
                          above_label='  ')

        self.plotter.add_mesh(self.mesh,
                              scalar_bar_args=self.sargs,
                              annotations=False,
                              show_edges=False,
                              smooth_shading=True,
                              scalars=self.volt,
                              nan_color=self.nan_color,
                              clim=[self.minval,self.maxval],
                              cmap=self.cmap,
                              below_color=self.below_color,
                              above_color=self.above_color)
        
        for indx in range(len(self.freeboundary_points['FreeboundaryPoints'])):
            self.plotter.add_lines(self.freeboundary_points['FreeboundaryPoints'][indx],color='black',width=5)
            self.plotter.reset_camera()

        
        self.plotter1.add_mesh(self.mesh,
                              scalar_bar_args=self.sargs,
                              annotations=False,
                              show_edges=False,
                              smooth_shading=True,
                              scalars=self.volt,
                              nan_color=self.nan_color,
                              clim=[self.minval,self.maxval],
                              cmap=self.cmap,
                              below_color=self.below_color,
                              above_color=self.above_color)





    def on_click3(self):
        self.minval = float(self.lowerlimit.text())
        self.maxval = float(self.upperlimit.text())
        # self.plotter.clear()
        self.plotter.add_mesh(self.mesh,
                              scalar_bar_args=self.sargs,
                              show_edges=False,
                              smooth_shading=True,
                              scalars=self.volt,
                              nan_color=self.nan_color,
                              clim=[self.minval,self.maxval],
                              cmap=self.cmap,
                              below_color=self.below_color,
                              above_color=self.above_color)

def main():
    # Create an instance of Qapplication
    app = qtw.QApplication(sys.argv)
    # Create an instance of GUI
    window = OpenEpGUI()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()