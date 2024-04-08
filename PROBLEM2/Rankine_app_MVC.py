#region imports
import sys
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from Rankine_GUI import Ui_Form  # Importing the user interface from Rankine_GUI module
from Rankine_Classes_MVC import rankineController  # Importing a controller class
from UnitConversions import UnitConverter as UC  # Importing UnitConverter class
# These imports are necessary for drawing a Matplotlib graph on the GUI
# No simple widget exists in QT Designer, so it's added in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
#endregion

class MainWindow(qtw.QWidget, Ui_Form):
    """
    Main window class for the application.

    This class inherits from QWidget and Ui_Form, which is the generated UI from Rankine_GUI module.
    """

    def __init__(self):
        """
        Constructor for the MainWindow class.
        """
        super().__init__()  # Call the constructor of the parent class
        self.setupUi(self)  # Set up the UI
        self.AssignSlots()  # Assign slots for signal handling
        self.MakeCanvas()  # Create canvas for plotting

        # Create lists of input and display widgets
        self.input_widgets = [self.rb_SI,self.le_PHigh, self.le_PLow, self.le_TurbineInletCondition, self.rdo_Quality, self.le_TurbineEff, self.cmb_XAxis, self.cmb_YAxis, self.chk_logX, self.chk_logY]
        self.display_widgets=[self.lbl_PHigh, self.lbl_PLow, self.lbl_SatPropLow,self.lbl_SatPropHigh, self.lbl_TurbineInletCondition, self.lbl_H1, self.lbl_H1Units, self.lbl_H2, self.lbl_H2Units, self.lbl_H3, self.lbl_H3Units, self.lbl_H4, self.lbl_H4Units, self.lbl_TurbineWork, self.lbl_TurbineWorkUnits, self.lbl_PumpWork, self.lbl_PumpWorkUnits, self.lbl_HeatAdded, self.lbl_HeatAddedUnits, self.lbl_ThermalEfficiency, self.canvas, self.figure, self.ax]
        
        # Instantiate a rankineController object and pass the lists of widgets to be unpacked into the View
        self.RC = rankineController(self.input_widgets, self.display_widgets)

        self.setNewPHigh()  # Set new PHigh value
        self.setNewPLow()  # Set new PLow value

        self.Calculate()  # Calculate using initial values

        # Variables to store coordinates from the last position on the graph
        self.oldXData = 0.0
        self.oldYData = 0.0

        self.show()  # Show the main window
        
  def AssignSlots(self):
    """
    Setup signals and slots for the program.

    Connects various UI elements to corresponding methods for handling events.

    :return: None
    """
    # Connect UI elements to methods for handling events
    self.btn_Calculate.clicked.connect(self.Calculate)
    self.rdo_Quality.clicked.connect(self.SelectQualityOrTHigh)
    self.rdo_THigh.clicked.connect(self.SelectQualityOrTHigh)
    self.le_PHigh.editingFinished.connect(self.setNewPHigh)
    self.le_PLow.editingFinished.connect(self.setNewPLow)
    self.rb_SI.clicked.connect(self.SetUnits)
    self.rb_English.clicked.connect(self.SetUnits)
    self.cmb_XAxis.currentIndexChanged.connect(self.SetPlotVariables)
    self.cmb_YAxis.currentIndexChanged.connect(self.SetPlotVariables)
    self.chk_logX.toggled.connect(self.SetPlotVariables)
    self.chk_logY.toggled.connect(self.SetPlotVariables)

def MakeCanvas(self):
    """
    Create a canvas for plotting the Rankine cycle graph.

    Steps:
    1. Create a Figure object called self.figure.
    2. Create a FigureCanvasQTAgg object called self.canvas.
    3. Create an axes object for making plot.
    4. Add self.canvas to self.Layout_Plot which is a Vertical box layout.
    5. Attach an event handler for mouse movement on graph.

    :return: None
    """
    # Step 1.
    self.figure = Figure(figsize=(1, 1), tight_layout=True, frameon=True)
    # Step 2.
    self.canvas = FigureCanvasQTAgg(self.figure)
    # Step 3.
    self.ax = self.figure.add_subplot()
    # Step 4.
    self.Layout_Plot.addWidget(NavigationToolbar2QT(self.canvas, self))
    self.Layout_Plot.addWidget(self.canvas)
    # Step 5.
    self.canvas.mpl_connect("motion_notify_event", self.mouseMoveEvent_Canvas)

def mouseMoveEvent_Canvas(self, event):
    """
    Handle mouse movement event on the canvas.

    Displays the current coordinates of the mouse cursor on the graph.

    :param event: Mouse event object
    :return: None
    """
    self.oldXData = event.xdata if event.xdata is not None else self.oldXData
    self.oldYData = event.ydata if event.ydata is not None else self.oldYData
    sUnit = 'kJ/(kg*K)' if self.rb_SI.isChecked() else 'BTU/(lb*R)'
    TUnit = 'C' if self.rb_SI.isChecked() else 'F'
    self.setWindowTitle('s:{:0.2f} {}, T:{:0.2f} {}'.format(self.oldXData, sUnit, self.oldYData, TUnit))

# Remaining methods seem to call corresponding methods of the controller, which handle various aspects of the application.

# Since they merely delegate functionality to the controller, no additional comments are deemed necessary.
#endregion

#if this module is being imported, this won't run. If it is the main module, it will run.
if __name__== '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Rankine calculator')
    sys.exit(app.exec())
