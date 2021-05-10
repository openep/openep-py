'''
GUI code for OpenEp
'''
import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import scipy.io as sio
# from ui_imagedialog import Ui_ImageDialog

class App(QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'OpenEp Application'
        self.left = 10
        self.top = 10
        self.width = 840
        self.height = 480
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        main_layout = QFormLayout(self)
        mainLabel1 = QLabel('OpenEp Tool')
        main_layout.addWidget(mainLabel1)
        mainLabel2 = QLabel('''The Open Source solution for electrophysiology data analysis''')
        main_layout.addWidget(mainLabel2)

        button1 = QPushButton('Load OpenEp Data', self)
        button1.setGeometry(200,150,100,40)
        button1.clicked.connect(self.on_click)
        main_layout.addWidget(button1)

        # Dialog = QDialog(self)
        # ui = Ui_ImageDialog()
        # ui.setupUi(Dialog)
        # Dialog.show()




        button2 = QPushButton('Plot Voltage Map', self)
        button2.setGeometry(200,150,100,40)
        button2.clicked.connect(self.on_click2)
        main_layout.addWidget(button2)



        # mainLabel.setAlignment(Qt.Qt.AlignCenter)
        # print(dir(mainLabel))






        self.show()

    def on_click(self):
        print('Loading Data ... ')
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        # dlg.setFilter("Text files (*.txt)")
        # filenames = QStringList()

        if dlg.exec_():
            filenames = dlg.selectedFiles()
            print(filenames[0])
            main_file = sio.loadmat(filenames[0])
            data_tri_X = main_file['userdata']['surface_triRep_X'][0][0]
            x = data_tri_X[:,0]
            y = data_tri_X[:,1]
            z = data_tri_X[:,2]
            print('x\n',x)
            print('y\n',y)
            print('z\n',z)



    def on_click2(self):
        print('Loading Voltage Map ... ')




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())




































# from PyQt5.QtWidgets import QApplication,QLabel,QWidget,QPushButton,QVBoxLayout
#
# # Instance of QApplication
# # Pass in the command line arguments within the brackets
# app = QApplication([])
# app.setStyle('Fusion')
#
# # w = 500
# # h = 300
#
#
# # WINDOW
# window = QWidget()
# window.title = 'OpenEp Application'
# window.height = 300
# window.width = 500
#
# # LAYOUT - ARRANGING THE WIDGETS
# layout = QVBoxLayout()
# layout.addWidget(QLabel('OpenEP Tool'))
# layout.addWidget(QPushButton('LoadDataset'))
# layout.addWidget(QPushButton('Voltage Plot'))
#
# window.setLayout(layout)
# window.show()
#
#
# # Run the application
# app.exec()
