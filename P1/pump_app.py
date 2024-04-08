# region imports
import numpy as np
import PyQt5.QtWidgets as qtw
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from pathlib import Path
import sys
import os

# Custom imports for MVC components and UI generated from Qt Designer
from pump import Ui_Form
from Pump_MVC import Pump_Controller


# endregion

# region class definitions
class PumpCurve_GUI_Class(Ui_Form, qtw.QWidget):
    """
    This class serves as the main window for the Pump Curve application. It integrates the UI
    with the underlying Model-View-Controller (MVC) architecture for processing and presenting
    pump curve data.
    """

    def __init__(self):
        """
        Constructor for PumpCurve_GUI_Class. Sets up the user interface, connects signals,
        and initializes the pump controller.
        """
        super().__init__()  # Initialize the QWidget component
        self.setupUi(self)  # Set up the UI from the inherited Ui_Form
        self.AssignSignals()  # Connect UI signals to the corresponding slots
        self.FilePath = os.getcwd()  # Set the initial file path to the current working directory
        self.FileName = ""  # Initialize the file name as an empty string

        # Matplotlib canvas setup for plotting
        self.canvas = FigureCanvasQTAgg(Figure(figsize=(5, 3), tight_layout=True, frameon=True))
        self.ax = self.canvas.figure.add_subplot()
        self.GL_Output.addWidget(self.canvas, 5, 0, 1, 4)  # Add the canvas to the GridLayout

        # Initialize the Pump_Controller from the MVC components
        self.myPump = Pump_Controller()
        self.setViewWidgets()  # Pass the UI widgets to the controller

        self.show()  # Display the window

    def AssignSignals(self):
        """
        Connects UI buttons to their respective functions (slots).
        """
        self.PB_Exit.clicked.connect(self.Exit)
        self.CMD_Open.clicked.connect(self.ReadAndCalculate)

    def setViewWidgets(self):
        """
        Passes UI widgets to the controller so it can update the view.
        """
        widgets = [
            self.LE_PumpName, self.LE_FlowUnits, self.LE_HeadUnits,
            self.LE_HeadCoefs, self.LE_EffCoefs, self.ax, self.canvas
        ]
        self.myPump.setViewWidgets(widgets)

    def ReadAndCalculate(self):
        """
        Opens a dialog box for selecting a data file and triggers data processing.
        Updates the UI with the processed data.
        :return: A boolean indicating whether the file was opened and processed successfully.
        """
        if self.OpenFile():
            with open(self.FileName, 'r') as file:
                data = file.readlines()
            self.myPump.ImportFromFile(data)
            return True
        return False

    def OpenFile(self):
        """
        Opens a file dialog to select a data file.
        :return: A boolean indicating whether a file was selected.
        """
        fname = qtw.QFileDialog.getOpenFileName(self, 'Open file', self.FilePath)
        fileSelected = len(fname[0]) > 0
        if fileSelected:
            self.FileName = fname[0]
            self.FilePath = str(Path(fname[0]).parents[0]) + '/'
            self.TE_Filename.setText(self.FileName)
        return fileSelected

    def Exit(self):
        """
        Exits the application.
        """
        qapp.exit()


# endregion

# region function definitions
def main():
    """
    The main function that creates an instance of the GUI and starts the application loop.
    """
    PumpCurve_GUI = PumpCurve_GUI_Class()
    qapp.exec_()


# endregion

# region function calls
if __name__ == "__main__":
    qapp = qtw.QApplication(sys.argv)
    main()
# endregion