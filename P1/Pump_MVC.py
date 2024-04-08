# region imports
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import PyQt5.QtWidgets as qtw

# importing from previous work on least squares fit
from LeastSquares import LeastSquaresFit_Class


# endregion

# region class definitions
class Pump_Model():
    """
    This is the pump model. It just stores data.
    """

    def __init__(self):
        """
        Initializes the Pump_Model class.
        """
        # create some class variables for storing information
        self.PumpName = ""
        self.FlowUnits = ""
        self.HeadUnits = ""

        # place to store data from file
        self.FlowData = np.array([])
        self.HeadData = np.array([])
        self.EffData = np.array([])

        # place to store coefficients for cubic fits
        self.HeadCoefficients = np.array([])
        self.EfficiencyCoefficients = np.array([])

        # create two instances (objects) of least squares class
        self.LSFitHead = LeastSquaresFit_Class()
        self.LSFitEff = LeastSquaresFit_Class()


class Pump_Controller():
    def __init__(self):
        """
        Initializes the Pump_Controller class.
        """
        self.Model = Pump_Model()
        self.View = Pump_View()

    # region functions to modify data of the model
    def ImportFromFile(self, data):
        """
        Processes the list of strings in data to build the pump model.

        :param data: List of strings containing pump data.
        """
        self.Model.PumpName = data[0].strip()  # first txt line is pump name
        units_line = data[2].strip()  # third line containing the units
        units = units_line.split()
        L = data[2].split()
        self.Model.FlowUnits = units[0]  # gpm
        self.Model.HeadUnits = units[1]  # ft

        # extracts flow, head and efficiency data and calculates coefficients
        self.SetData(data[3:])
        self.updateView()

    def SetData(self, data):
        '''
        Expects three columns of data in an array of strings with space delimiter
        Parse line and build arrays.
        :param data: Array of strings containing pump data.
        '''
        # erase existing data
        self.Model.FlowData = np.array([])
        self.Model.HeadData = np.array([])
        self.Model.EffData = np.array([])

        # parse new data
        for L in data:
            Cells = L.split()  # split the line into a list of strings
            self.Model.FlowData = np.append(self.Model.FlowData, float(Cells[0]))  # flow data
            self.Model.HeadData = np.append(self.Model.HeadData, float(Cells[1]))  # head data
            self.Model.EffData = np.append(self.Model.EffData, float(Cells[2]))  # efficiency data

        # call least square fit for head and efficiency
        self.LSFit()

    def LSFit(self):
        '''Fit cubic polynomial using Least Squares'''
        self.Model.LSFitHead.x = self.Model.FlowData
        self.Model.LSFitHead.y = self.Model.HeadData
        self.Model.LSFitHead.LeastSquares(3)  # calls LeastSquares function of LSFitHead object

        self.Model.LSFitEff.x = self.Model.FlowData
        self.Model.LSFitEff.y = self.Model.EffData
        self.Model.LSFitEff.LeastSquares(3)  # calls LeastSquares function of LSFitEff object

    # endregion

    # region functions interacting with view
    def setViewWidgets(self, w):
        """
        Sets view widgets for the Pump_Controller.

        :param w: Tuple containing view widgets.
        """
        self.View.setViewWidgets(w)

    def updateView(self):
        """
        Updates the view.
        """
        self.View.updateView(self.Model)
    # endregion


class Pump_View():
    def __init__(self):
        """
        Initializes the Pump_View class.
        """
        self.LE_PumpName = qtw.QLineEdit()
        self.LE_FlowUnits = qtw.QLineEdit()
        self.LE_HeadUnits = qtw.QLineEdit()
        self.LE_HeadCoefs = qtw.QLineEdit()
        self.LE_EffCoefs = qtw.QLineEdit()
        self.ax = None
        self.canvas = None

    def updateView(self, Model):
        """
        Updates the view with model parameters.

        :param Model: Pump_Model object.
        """
        self.LE_PumpName.setText(Model.PumpName)
        self.LE_FlowUnits.setText(Model.FlowUnits)
        self.LE_HeadUnits.setText(Model.HeadUnits)
        self.LE_HeadCoefs.setText(Model.LSFitHead.GetCoeffsString())
        self.LE_EffCoefs.setText(Model.LSFitEff.GetCoeffsString())
        self.DoPlot(Model)

    def DoPlot(self, Model):
        """
        Creates a plot with dual y-axes for Head and Efficiency.

        :param Model: Pump_Model object.
        """
        # Get the fit and plot info for Head and Efficiency
        headx, heady, headRSq = Model.LSFitHead.GetPlotInfo(2, npoints=500)  # Quadric fit for Head
        effx, effy, effRSq = Model.LSFitEff.GetPlotInfo(3, npoints=500)  # Cubic fit for Efficiency

        # Clear any existing plots
        self.ax.clear()
        self.ax.set_xlabel("Flow Rate (gpm)")
        self.ax.set_ylabel("Head (ft)")

        # Plot the Head data and fit on the primary y-axis
        head_data_line, = self.ax.plot(Model.FlowData, Model.HeadData, 'ko', label='Head', linestyle='None')
        head_fit_line, = self.ax.plot(headx, heady, 'k--', label=f'Head (R²={headRSq:.3f})')

        # Create a secondary y-axis for the Efficiency data and fit
        ax2 = self.ax.twinx()
        ax2.set_ylabel("Efficiency (%)")
        eff_data_line, = ax2.plot(Model.FlowData, Model.EffData, 'k^', label='Efficiency', linestyle='None')
        eff_fit_line, = ax2.plot(effx, effy, 'k:', label=f'Efficiency (R²={effRSq:.3f})')

        # Set tick marks to point inwards
        self.ax.tick_params(direction='in')
        ax2.tick_params(direction='in')

        # Add legends separately
        self.ax.legend(handles=[head_fit_line, head_data_line], loc='center left')  # adds head legend  to center left
        ax2.legend(handles=[eff_fit_line, eff_data_line], loc='upper right')  # adds efficiency legend to top right

        # Draw the canvas
        self.canvas.draw()

    def setViewWidgets(self, w):
        """
        Sets view widgets.

        :param w: Tuple containing view widgets.
        """
        self.LE_PumpName, self.LE_FlowUnits, self.LE_HeadUnits, self.LE_HeadCoefs, self.LE_EffCoefs, self.ax, self.canvas = w
# endregion