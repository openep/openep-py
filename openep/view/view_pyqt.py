'''
GUI code for OpenEp
'''

from logging import setLogRecordFactory
from PyQt5 import QtWidgets as qtw
from matplotlib.figure import Figure
import pyvista as pv
from pyvistaqt import QtInteractor, MainWindow
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar

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
        self.egm_point = 0



    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # # LAYOUTS
        self.mainLayout = qtw.QVBoxLayout(self)
        self.menulayout = qtw.QGridLayout(self)
        self.labelLayout = qtw.QFormLayout(self)
        self.buttonLayout = qtw.QHBoxLayout(self)

        # MenuBar
        menubar=qtw.QMenuBar()
        self.menulayout.addWidget(menubar, 0, 0)
        file_menu = menubar.addMenu("File")
        
        plot_menu =qtw.QMenu('Plots',self)
        voltage_map_plot_act = qtw.QAction('3-D Voltage Map',self)
        self.voltage_map_electroanatomic_plot_act  = qtw.QAction('Electroanatomic Plot',self)
        self.historgram_plot_act = qtw.QAction('Histogram Plot',self)
        self.egm_plot_act = qtw.QAction('EGM Plot',self)


        voltage_map_plot_act.triggered.connect(self.on_click2)
        plot_menu.addAction(voltage_map_plot_act)
        
        self.voltage_map_electroanatomic_plot_act.triggered.connect(self.plot_electroanatomic)
        plot_menu.addAction(self.voltage_map_electroanatomic_plot_act)
        
        self.historgram_plot_act.triggered.connect(self.plot_histogram)
        plot_menu.addAction(self.historgram_plot_act)


        self.egm_plot_act.triggered.connect(self.plot_egm)
        plot_menu.addAction(self.egm_plot_act)

        file_menu.addMenu(plot_menu)



        file_menu.addSeparator()
        file_menu.addAction("Quit")

        menubar.addMenu("Edit")
        menubar.addMenu("View")
        menubar.addMenu("Help")

        # Label
        mainLabel1 = qtw.QLabel('OpenEp Tool',self)
        self.labelLayout.addWidget(mainLabel1)
        mainLabel2 = qtw.QLabel( 'The Open Source solution for electrophysiology data analysis',self)
        self.labelLayout.addWidget(mainLabel2)

        # # Button 1
        button1 = qtw.QPushButton('Load OpenEp Data', self)
        button1.setGeometry(200,150,100,40)
        button1.clicked.connect(self.on_click)


        # # Button 2
        button2 = qtw.QPushButton('Plot Voltage Map', self)
        button2.setGeometry(200,150,100,40)
        button2.clicked.connect(self.on_click2)

        # # Adding buttons to the horizontal layout
        self.buttonLayout.addWidget(button1)
        self.buttonLayout.addWidget(button2)

        self.labelLayout.addRow(self.buttonLayout)

        # # Plot
        self.frame = qtw.QFrame()
        self.plotLayout = qtw.QGridLayout()
        self.plotLayout.setColumnStretch(0,5)
        self.plotLayout.setColumnStretch(1,5)


        self.plotter = QtInteractor(self.frame)
        self.plotLayout.addWidget(self.plotter.interactor,0,0)
        
        
        # # VoltThresholds limit
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

        # selecting the Egm points
        self.egmselect = qtw.QLineEdit()
        self.egmselect.setText(str(0))
        self.limitLayout.addRow('Egm Point',self.egmselect)

        button_egmselect = qtw.QPushButton('Set Egm Point', self)
        button_egmselect.setGeometry(200,150,100,40)
        button_egmselect.clicked.connect(self.plot_egm)
        self.limitLayout.addRow(button_egmselect)



        # Nesting Layouts and setting the main layout
        self.mainLayout.addLayout(self.menulayout)
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
                            volt='bip',
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



    def plot_electroanatomic(self):

        distance_thresh = 10
        # Anatomic descriptions (Mesh) - nodes and indices
        pts = self.ep_case.nodes
        indices = self.ep_case.indices

        # Electric data
        # Locations – Cartesian co-ordinates, projected on to the surface 
        locations = case_routines.get_electrogram_coordinates(self.ep_case,'type','bip')

        i_egm = self.ep_case.electric['egm'].T
        i_vp = case_routines.getMappingPointsWithinWoI(self.ep_case)
        # macthing the shape of ivp with data
        i_vp_egm = np.repeat(i_vp, repeats=i_egm.shape[1], axis=1)
        # macthing the shape of ivp with coords
        i_vp_locations = np.repeat(i_vp, repeats=locations.shape[1],axis=1)

        # Replacing the values outside the window of interest with Nan values
        i_egm[~i_vp_egm] = np.nan
        locations[~i_vp_locations] = np.nan

        # For each mapping point, n, find the voltage amplitude
        max_volt = np.amax(a=i_egm,axis=1).reshape(len(i_egm),1)
        min_volt = np.amin(a=i_egm,axis=1).reshape(len(i_egm),1)

        amplitude_volt = np.subtract(max_volt,min_volt)

        for indx in range(amplitude_volt.shape[1]):
            temp_data = amplitude_volt[:,indx]
            temp_coords = locations
            i_nan = np.isnan(temp_data)
            temp_data=temp_data[~i_nan]
            temp_coords=temp_coords[~i_nan]


            interp = case_routines.OpenEPDataInterpolator(method='rbf',distanceThreshold=distance_thresh,rbfConstant=1)
            vertex_voltage_data = interp.interpolate(x0=temp_coords,d0=temp_data,x1=pts)

        # # QDock Widget
        self.plotter1 = QtInteractor(self.frame)
        self.dock_plot = qtw.QDockWidget("ElectroAnatomic Plot", self)
        self.dock_plot.setFloating(False)
        self.dock_plot.setWidget(self.plotter1)

        self.plotLayout.addWidget(self.dock_plot,0,1)

        self.plotter1.add_mesh(self.mesh,
                        scalar_bar_args=self.sargs,
                        show_edges=False,
                        smooth_shading=True,
                        scalars=vertex_voltage_data,
                        nan_color=self.nan_color,
                        clim=[self.minval,self.maxval],
                        cmap=self.cmap,
                        below_color=self.below_color,
                        above_color=self.above_color)


    def plot_egm(self):
        self.egm_point = int(self.egmselect.text())

        self.fig,self.ax = plt.subplots(ncols=1,nrows=1)
        self.fig.set_facecolor('gray')
        self.canvas = FigureCanvas(self.fig)

        self.egm = case_routines.get_egms_at_points(self.ep_case,"iegm",[self.egm_point])
        self.egm_traces = self.egm['egm_traces']
        self.sample_range = self.egm['sample_range'][0]
        seperation = 7

        for i in range(len(self.egm_traces)):
            y = self.egm_traces[i][0][self.sample_range[0]:self.sample_range[1]]
            t = np.arange(self.sample_range[0],self.sample_range[1],1)
            self.ax.plot(t,y+(seperation*i))
            self.ax.get_yaxis().set_visible(False)
            self.ax.set_xlabel('Samples')
        
        
        # toolbar = NavigationToolbar(self.canvas,self)
        self.dock_plot1 = qtw.QDockWidget("EGM Plot", self)
        self.dock_plot1.setFloating(False)
        self.dock_plot1.setWidget(self.canvas)
        # self.dock_plot1.setWidget(toolbar)

        self.plotLayout.addWidget(self.dock_plot1,1,1)


    def plot_histogram(self):
        self.plotter3 = QtInteractor(self.frame)
        self.dock_plot2 = qtw.QDockWidget("Histogram Plot", self)
        self.dock_plot2.setFloating(False)
        self.dock_plot2.setWidget(self.plotter3)

        self.plotLayout.addWidget(self.dock_plot2,1,0)
        pass

    def on_click3(self):
        self.minval = float(self.lowerlimit.text())
        self.maxval = float(self.upperlimit.text())
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