# region imports
import math  # Provides access to mathematical functions (e.g., sin, cos, sqrt).
from Calc_state import *  # Imports functions and classes for calculating thermodynamic state properties.
from UnitConversions import UnitConverter as UC  # For performing unit conversions, aliased as UC for simplicity.
import numpy as np  # Used for numerical operations, array manipulations.
from matplotlib import pyplot as plt  # For plotting graphs.
from copy import deepcopy as dc  # To create deep copies of objects, ensuring no unintended shared references.

# These imports are required for embedding matplotlib graphs into a Qt-based graphical user interface (GUI).
# Since Qt Designer does not offer a direct widget for matplotlib, it needs to be manually integrated.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # Embeds matplotlib into Qt widgets.
from matplotlib.figure import Figure  # Represents a figure in matplotlib.
#endregion

# region class definitions
class rankineModel:
    def __init__(self):
        """
        Initializes a model for a Rankine power cycle. This class is designed according to the Model component
        of the Model-View-Controller (MVC) architecture, focusing solely on data storage. The Controller component
        is responsible for modifying this model based on user input, whereas the View component displays the data.

        Attributes:
            p_low (float): Low pressure isobar for the cycle, in kPa. Default is None.
            p_high (float): High pressure isobar for the cycle, in kPa. Default is None.
            t_high (float): Optional temperature at turbine inlet (State 1), in degrees C. Default is None.
            name (str): A label for the cycle. Default is None.
            efficiency (float): The overall cycle efficiency. Default is None.
            turbine_eff (float): Isentropic efficiency of the turbine. Default is None.
            turbine_work (float): Work done by the turbine, in the same units as heat_added. Default is None.
            pump_work (float): Work required by the pump, in the same units as heat_added. Default is None.
            heat_added (float): Heat added during the cycle, typically in kJ/kg. Default is None.
            steam (Steam_SI): An object for calculating steam properties. Instantiated by default.
            state1, state2s, state2, state3, state4 (stateProps): Objects representing key thermodynamic states in the cycle.
            SI (bool): True if using SI units, False for English units. Affects pressure and temperature inputs.
            satLiqPlotData, satVapPlotData, upperCurve, lowerCurve (StateDataForPlotting): Containers for data useful in plotting.
        """
        # Initialize all attributes to None or appropriate default values.
        self.p_low = None
        self.p_high = None
        self.t_high = None
        self.name = None
        self.efficiency = None
        self.turbine_eff = None
        self.turbine_work = None
        self.pump_work = None
        self.heat_added = None
        # Instantiate a steam object and initialize state properties for key cycle states.
        self.steam = Steam_SI()
        self.state1 = stateProps()
        self.state2s = stateProps()
        self.state2 = stateProps()
        self.state3 = stateProps()
        self.state4 = stateProps()
        # Set SI unit flag and prepare containers for plotting data.
        self.SI = True
        self.satLiqPlotData = StateDataForPlotting()
        self.satVapPlotData = StateDataForPlotting()
        self.upperCurve = StateDataForPlotting()
        self.lowerCurve = StateDataForPlotting()

class rankineView:
    def __init__(self):
        """
        Constructor for the rankineView class. This class represents the View in the MVC architecture,
        handling the graphical user interface components. The constructor is intentionally left empty
        because the initialization of this class doesn't require any specific setup. The GUI components
        are set up externally and passed to this class through other methods.
        """

    def setWidgets(self, *args):
        """
        Initializes class variables to reference the GUI widgets. This method categorizes widgets into
        input and display widgets for further manipulation.

        :param args: A tuple containing two lists:
                     args[0] - A list of input widget references.
                     args[1] - A list of display widget references.
        """
        # Assign input widget references to class variables for easy access.
        (self.rb_SI, self.le_PHigh, self.le_PLow, self.le_TurbineInletCondition,
         self.rdo_Quality, self.le_TurbineEff, self.cmb_XAxis, self.cmb_YAxis,
         self.chk_logX, self.chk_logY) = args[0]

        # Assign display widget references to class variables for easy access.
        (self.lbl_PHigh, self.lbl_PLow, self.lbl_SatPropLow, self.lbl_SatPropHigh,
         self.lbl_TurbineInletCondition, self.lbl_H1, self.lbl_H1Units, self.lbl_H2,
         self.lbl_H2Units, self.lbl_H3, self.lbl_H3Units, self.lbl_H4, self.lbl_H4Units,
         self.lbl_TurbineWork, self.lbl_TurbineWorkUnits, self.lbl_PumpWork,
         self.lbl_PumpWorkUnits, self.lbl_HeatAdded, self.lbl_HeatAddedUnits,
         self.lbl_ThermalEfficiency, self.canvas, self.figure, self.ax) = args[1]

    def selectQualityOrTHigh(self, Model=None):
        """
        Handles the selection between Quality (x) and high temperature (THigh) for the turbine inlet
        condition. It updates the GUI accordingly and recalculates the state properties if necessary.

        :param Model: The Rankine cycle model instance to query or update based on the selection. This
                      allows for direct interaction with the model to fetch or modify cycle parameters.
        """
        # Check the system of units selected by the user (SI or English units).
        SI = self.rb_SI.isChecked()

        # If Quality is selected, set the Turbine Inlet Condition field to 1.0 and disable it.
        if self.rdo_Quality.isChecked():
            self.le_TurbineInletCondition.setText("1.0")
            self.le_TurbineInletCondition.setEnabled(False)
        else:
            # If THigh is selected, enable the field and update its value based on the saturated temperature at PHigh.
            # Convert the pressure if needed and fetch the corresponding saturated temperature.
            PCF = 1 if SI else UC.psi_to_bar
            satPHigh = Model.steam.getsatProps_p(float(self.le_PHigh.text()) * PCF)
            Tsat = satPHigh.tsat
            Tsat = Tsat if SI else UC.C_to_F(Tsat)
            CurrentT = float(self.le_TurbineInletCondition.text())
            CurrentT = max(CurrentT, Tsat)  # Ensure the temperature is not below saturation.
            self.le_TurbineInletCondition.setText("{:0.2f}".format(CurrentT))
            self.le_TurbineInletCondition.setEnabled(True)

        # Update the label for the Turbine Inlet Condition to reflect the current selection.
        x = self.rdo_Quality.isChecked()
        self.lbl_TurbineInletCondition.setText(
            "Turbine Inlet: {}{} =".format('x' if x else 'THigh', '' if x else ('(C)' if SI else '(F)')))


   def setNewPHigh(self, Model=None):
    """
    Updates the displayed high saturation properties based on the high pressure value entered by the user.
    This function adapts the display to reflect changes in units (SI or English) and recalculates the
    saturation properties for the new pressure.

    :param Model: The Rankine cycle model instance, used to access the steam property calculations.
    """
    # Determine if the system is using SI units based on the GUI checkbox.
    SI = self.rb_SI.isChecked()
    # Define the pressure conversion factor based on the unit system.
    PCF = 1 if SI else UC.psi_to_bar
    # Retrieve the saturation properties for the new high pressure value.
    satProp = Model.steam.getsatProps_p(float(self.le_PHigh.text()) * PCF)
    # Update the GUI label to show the new saturation properties.
    self.lbl_SatPropHigh.setText(satProp.getTextOutput(SI=SI))
    # Trigger any additional GUI updates needed due to the change in high pressure.
    self.selectQualityOrTHigh(Model)

def setNewPLow(self, Model=None):
    """
    Updates the displayed low saturation properties based on the low pressure value entered by the user.
    Similar to setNewPHigh, it adjusts for unit differences and recalculates saturation properties.

    :param Model: The Rankine cycle model instance for accessing steam properties.
    """
    # Check if SI units are being used.
    SI = self.rb_SI.isChecked()
    # Define pressure conversion factor based on unit system.
    PCF = 1 if SI else UC.psi_to_bar
    # Retrieve and display the saturation properties for the new low pressure value.
    satProp = Model.steam.getsatProps_p(float(self.le_PLow.text()) * PCF)
    self.lbl_SatPropLow.setText(satProp.getTextOutput(SI=SI))

def outputToGUI(self, Model=None):
    """
    Outputs the current state of the Rankine cycle model to the GUI. This includes updating
    enthalpy values for key states, work done by the turbine and pump, heat added to the cycle,
    and the overall thermal efficiency. Additionally, updates the plot with the new cycle data.

    :param Model: The Rankine cycle model instance containing the current state and calculations.
    """
    # Early exit if the cycle has not been evaluated yet.
    if Model.state1 is None:
        return
    
    # Determine the enthalpy conversion factor based on the system of units.
    HCF = 1 if Model.SI else UC.kJperkg_to_BTUperlb
    # Update GUI labels with the latest model data, converted to the appropriate units.
    self.lbl_H1.setText("{:0.2f}".format(Model.state1.h * HCF))
    self.lbl_H2.setText("{:0.2f}".format(Model.state2.h * HCF))
    self.lbl_H3.setText("{:0.2f}".format(Model.state3.h * HCF))
    self.lbl_H4.setText("{:0.2f}".format(Model.state4.h * HCF))
    self.lbl_TurbineWork.setText("{:0.2f}".format(Model.turbine_work * HCF))
    self.lbl_PumpWork.setText("{:0.2f}".format(Model.pump_work * HCF))
    self.lbl_HeatAdded.setText("{:0.2f}".format(Model.heat_added * HCF))
    self.lbl_ThermalEfficiency.setText("{:0.2f}".format(Model.efficiency))
    # Also, update saturation properties labels for both high and low pressures.
    satPropsLow = Model.steam.getsatProps_p(Model.p_low)
    satPropsHigh = Model.steam.getsatProps_p(Model.p_high)
    self.lbl_SatPropLow.setText(satPropsLow.getTextOutput(SI=Model.SI))
    self.lbl_SatPropHigh.setText(satPropsHigh.getTextOutput(SI=Model.SI))
    
    # Finally, update the cycle plot with the latest data.
    self.plot_cycle_XY(Model=Model)

   def updateUnits(self, Model=None):
    """
    Updates the units displayed in the GUI based on the current selection between SI and English units.
    This includes converting and displaying the pressures, temperatures, and enthalpies in the chosen units.

    :param Model: An instance of the Rankine cycle model containing current state and configuration.
    """
    # Refresh the GUI to ensure all displayed values are current.
    self.outputToGUI(Model=Model)

    # Determine the conversion factor for pressure based on selected units.
    pCF = 1 if Model.SI else UC.bar_to_psi
    # Update displayed pressures with the correct units and converted values.
    self.le_PHigh.setText("{:0.2f}".format(pCF * Model.p_high))
    self.le_PLow.setText("{:0.2f}".format(pCF * Model.p_low))

    # If turbine inlet condition is based on temperature (THigh), convert and update this value.
    if not self.rdo_Quality.isChecked():
        T = float(self.le_TurbineInletCondition.text())
        T = UC.F_to_C(T) if Model.SI else UC.C_to_F(T)
        TUnits = "C" if Model.SI else "F"
        self.le_TurbineInletCondition.setText("{:0.2f}".format(T))
        self.lbl_TurbineInletCondition.setText(f"Turbine Inlet: THigh ({TUnits}):")

    # Update labels to reflect the current unit system.
    self.lbl_PHigh.setText(f"P High ({'bar' if Model.SI else 'psi'})")
    self.lbl_PLow.setText(f"P Low ({'bar' if Model.SI else 'psi'})")
    HUnits = "kJ/kg" if Model.SI else "BTU/lb"
    # Apply unit updates to all relevant labels.
    for label in [self.lbl_H1Units, self.lbl_H2Units, self.lbl_H3Units, self.lbl_H4Units,
                  self.lbl_TurbineWorkUnits, self.lbl_PumpWorkUnits, self.lbl_HeatAddedUnits]:
        label.setText(HUnits)


    def print_summary(self, Model=None):
    """
    Outputs a summary of the Rankine cycle's performance metrics to the command line interface. This includes
    efficiency, work done by the turbine and pump, and the heat added, along with the state properties at key points.

    :param Model: The Rankine cycle model instance to pull data from.
    """
    # Calculate efficiency if not already done.
    if Model.efficiency is None:
        Model.calc_efficiency()
    
    # Print the cycle summary with formatted values for clarity.
    print(f'Cycle Summary for: {Model.name}')
    print(f'\tEfficiency: {Model.efficiency:.3f}%')
    print(f'\tTurbine Eff:  {Model.turbine_eff:.2f}')
    print(f'\tTurbine Work: {Model.turbine_work:.3f} kJ/kg')
    print(f'\tPump Work: {Model.pump_work:.3f} kJ/kg')
    print(f'\tHeat Added: {Model.heat_added:.3f} kJ/kg')
    # Print state properties for each key state in the cycle.
    for state in [Model.state1, Model.state2, Model.state3, Model.state4]:
        state.print()


    def plot_cycle_TS(self, axObj=None, Model=None):
    """
    Plots the Rankine cycle on T-S coordinates, including the vapor dome and the process paths between key states.
    The plot can either be displayed using pyplot directly or drawn on a pre-existing Matplotlib axes object.

    :param axObj: Optional. A Matplotlib axes object to draw the plot on. If None, a new plot will be created.
    :param Model: The Rankine cycle model instance containing state properties and other relevant data.
    """
    # Configure plot settings based on selected unit system.
    SI = Model.SI
    steam = Model.steam  # Shortcut to the steam property calculations in the model.

    # Load saturation properties to plot the vapor dome.
    ts, ps, hfs, hgs, sfs, sgs, vfs, vgs = np.loadtxt('sat_water_table.txt', skiprows=1, unpack=True)
    ax = plt.subplot() if axObj is None else axObj  # Determine the plotting context.

    # Apply unit conversions to the saturation properties based on the selected unit system.
    unitConversions = [1 if SI else conversion for conversion in (UC.kJperkg_to_BTUperlb, UC.kpa_to_psi, UC.kJperkgK_to_BTUperlbR, UC.kgperm3_to_lbperft3)]
    hfs, hgs, sfs, sgs, vfs, vgs, ps = [prop * conversion for prop, conversion in zip((hfs, hgs, sfs, sgs, vfs, vgs, ps), unitConversions)]
    ts = [t if SI else UC.C_to_F(t) for t in ts]  # Convert temperature if necessary.

    # Plot saturated liquid and vapor lines to represent the vapor dome.
    ax.plot(sfs, ts, color='blue')  # Saturated liquid line.
    ax.plot(sgs, ts, color='red')   # Saturated vapor line.
    # Additional steps would involve plotting the process paths and annotating the plot.

        # endregion

        # step 3:  I'll just make a straight line between state3 and state3p
        st3p = steam.getState(Model.p_high, x=0)  # saturated liquid state at p_high
        svals = np.linspace(Model.state3.s, st3p.s, 20)
        hvals = np.linspace(Model.state3.h, st3p.h, 20)
        pvals = np.linspace(Model.p_low, Model.p_high, 20)
        vvals = np.linspace(Model.state3.v, st3p.v, 20)
        tvals = np.linspace(Model.state3.T, st3p.T, 20)
        line3 = np.column_stack([svals, tvals])

        # step 4:
        sat_pHigh = steam.getState(Model.p_high, x=1.0)
        st1 = Model.state1
        svals2p = np.linspace(st3p.s, sat_pHigh.s, 20)
        hvals2p = np.linspace(st3p.h, sat_pHigh.h, 20)
        pvals2p = [Model.p_high for i in range(20)]
        vvals2p = np.linspace(st3p.v, sat_pHigh.v, 20)
        tvals2p = [st3p.T for i in range(20)]
        line4 = np.column_stack([svals2p, tvals2p])
        if st1.T > sat_pHigh.T:  # need to add data points to state1 for superheated
            svals_sh = np.linspace(sat_pHigh.s, st1.s, 20)
            tvals_sh = np.array([steam.getState(Model.p_high, s=ss).T for ss in svals_sh])
            line4 = np.append(line4, np.column_stack([svals_sh, tvals_sh]), axis=0)
        # plt.plot(line4[:,0], line4[:,1])

        # step 5:
        svals = np.linspace(Model.state1.s, Model.state2.s, 20)
        tvals = np.linspace(Model.state1.T, Model.state2.T, 20)
        line5 = np.array(svals)
        line5 = np.column_stack([line5, tvals])
        # plt.plot(line5[:,0], line5[:,1])

        # step 6:
        svals = np.linspace(Model.state2.s, Model.state3.s, 20)
        tvals = np.array([Model.state2.T for i in range(20)])
        line6 = np.column_stack([svals, tvals])
        # plt.plot(line6[:,0], line6[:,1])

        # step 7:
        topLine = np.append(line3, line4, axis=0)
        topLine = np.append(topLine, line5, axis=0)
        xvals = topLine[:, 0]
        y1 = topLine[:, 1]
        y2 = [Model.state3.T for s in xvals]

        if not SI:
            xvals *= UC.kJperkgK_to_BTUperlbR
            for i in range(len(y1)):
                y1[i] = UC.C_to_F(y1[i])
            for i in range(len(y2)):
                y2[i] = UC.C_to_F(y2[i])

        ax.plot(xvals, y1, color='darkgreen')
        ax.plot(xvals, y2, color='black')
        # ax.fill_between(xvals, y1, y2, color='gray', alpha=0.5)

        if SI:
            ax.plot(Model.state1.s, Model.state1.T, marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state2.s, Model.state2.T, marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state3.s, Model.state3.T, marker='o', markeredgecolor='k', markerfacecolor='w')
        else:
            ax.plot(Model.state1.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state1.T), marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state2.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state2.T), marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state3.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state3.T), marker='o', markeredgecolor='k', markerfacecolor='w')

        tempUnits = r'$\left(^oC\right)$' if SI else r'$\left(^oF\right)$'
        entropyUnits = r'$\left(\frac{kJ}{kg\cdot K}\right)$' if SI else r'$\left(\frac{BTU}{lb\cdot ^oR}\right)$'
        ax.set_xlabel(r's ' + entropyUnits, fontsize=18)  # different than plt
        ax.set_ylabel(r'T ' + tempUnits, fontsize=18)  # different than plt
        ax.set_title(Model.name)  # different than plt
        ax.grid(visible='both', alpha=0.5)
        ax.tick_params(axis='both', direction='in', labelsize=18)

        sMin = min(sfs)
        sMax = max(sgs)
        ax.set_xlim(sMin, sMax)  # different than plt

        tMin = min(ts)
        tMax = max(max(ts), st1.T)
        ax.set_ylim(tMin, tMax * 1.05)  # different than plt

        energyUnits = r'$\frac{kJ}{kg}$' if SI else r'$\frac{BTU}{lb}$'
        energyCF = 1 if SI else UC.kJperkg_to_BTUperlb

        if axObj is None:  # this allows me to show plot if not being displayed on a figure
            plt.show()

    def plot_cycle_XY(self, Model=None):
        """
        I want to plot any two thermodynaimc properties on X and Y
        :param X: letter for which variable to plot on X axis
        :param Y: letter for which variable to plot on Y axis
        :return:
        """
        ax = self.ax
        X = self.cmb_XAxis.currentText()
        Y = self.cmb_YAxis.currentText()
        logx = self.chk_logX.isChecked()
        logy = self.chk_logY.isChecked()
        SI = Model.SI
        if X == Y:
            return
        QTPlotting = True  # assumes we are plotting onto a QT GUI form
        if ax == None:
            ax = plt.subplot()
            QTPlotting = False  # actually, we are just using CLI and showing the plot

        ax.clear()
        ax.set_xscale('log' if logx else 'linear')
        ax.set_yscale('log' if logy else 'linear')
        YF = Model.satLiqPlotData.getDataCol(Y, SI=SI)
        YG = Model.satVapPlotData.getDataCol(Y, SI=SI)
        XF = Model.satLiqPlotData.getDataCol(X, SI=SI)
        XG = Model.satVapPlotData.getDataCol(X, SI=SI)
        # plot the vapor dome
        ax.plot(XF, YF, color='b')
        ax.plot(XG, YG, color='r')
        # plot the upper and lower curves
        ax.plot(Model.lowerCurve.getDataCol(X, SI=SI), Model.lowerCurve.getDataCol(Y, SI=SI), color='k')
        ax.plot(Model.upperCurve.getDataCol(X, SI=SI), Model.upperCurve.getDataCol(Y, SI=SI), color='g')
        # ax.fill_between(Model.upperCurve.getDataCol(X), Model.upperCurve.getDataCol(Y), self.lowerCurve.getDataCol(Y), color='grey', alpha=0.2)

        # add axis labels
        ax.set_ylabel(Model.lowerCurve.getAxisLabel(Y, SI=SI), fontsize='large' if QTPlotting else 'medium')
        ax.set_xlabel(Model.lowerCurve.getAxisLabel(X, SI=SI), fontsize='large' if QTPlotting else 'medium')
        # put a title on the plot
        Model.name = 'Rankine Cycle - ' + Model.state1.region + ' at Turbine Inlet'
        ax.set_title(Model.name, fontsize='large' if QTPlotting else 'medium')

        # modify the tick marks
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True,
                       labelsize='large' if QTPlotting else 'medium')  # format tick marks

        # plot the circles for states 1, 2, 3, and 4
        ax.plot(Model.state1.getVal(X, SI=SI), Model.state1.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        ax.plot(Model.state2.getVal(X, SI=SI), Model.state2.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        ax.plot(Model.state3.getVal(X, SI=SI), Model.state3.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        ax.plot(Model.state4.getVal(X, SI=SI), Model.state4.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        # set limits on x and y
        xmin = min(min(XF), min(XG), min(Model.upperCurve.getDataCol(X, SI=SI)), max(Model.lowerCurve.getDataCol(X, SI=SI)))
        xmax = max(max(XF), max(XG), max(Model.upperCurve.getDataCol(X, SI=SI)), max(Model.lowerCurve.getDataCol(X, SI=SI)))
        ymin = min(min(YF), min(YG), min(Model.upperCurve.getDataCol(Y, SI=SI)), max(Model.lowerCurve.getDataCol(Y, SI=SI)))
        ymax = max(max(YF), max(YG), max(Model.upperCurve.getDataCol(Y, SI=SI)),
                   max(Model.lowerCurve.getDataCol(Y, SI=SI))) * 1.1
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        deltax = xmax - xmin
        deltay = ymax - ymin
        # add the summary text to the plot

        # show the plot
        if QTPlotting == False:
            plt.show()
        else:
            self.canvas.draw()


class rankineController():
    def __init__(self, *args):
        """
        Initializes the controller for the Rankine cycle application, following the Model-View-Controller (MVC) pattern.
        This controller is responsible for creating the model and view components, linking them, and handling the flow of data
        between the user interface and the model based on user interactions.

        :param *args: Expected to receive two tuples/lists as arguments:
                      - The first contains widgets related to user inputs.
                      - The second contains widgets used for displaying outputs.
        """
        self.Model = rankineModel()  # Initialize the Rankine cycle model.
        self.View = rankineView()    # Initialize the view component for the GUI.
        self.IW = args[0]  # Input Widgets: A collection of widgets that receive user input.
        self.DW = args[1]  # Display Widgets: A collection of widgets used to display information to the user.
        self.View.setWidgets(self.IW, self.DW)  # Setup the view with references to the input and display widgets.

        self.buildVaporDomeData()  # Precompute data for the vapor dome, used in T-S diagrams.

    def updateModel(self):
        """
        Updates the Rankine cycle model with new parameters based on the current state of the input widgets.
        This includes reading the values for high and low pressures, the turbine inlet temperature or quality,
        and the isentropic efficiency of the turbine. After updating the model, it recalculates the cycle efficiency
        and updates the view to reflect the new state of the model.
        """
        # Determine the system of units based on the GUI state and set the model accordingly.
        self.Model.SI = self.View.rb_SI.isChecked()

        # Convert input values from the GUI to SI units if necessary before updating the model.
        PCF = 1 if self.Model.SI else UC.psi_to_bar  # Pressure Conversion Factor: Adjusts based on selected units.
        self.Model.p_high = float(self.View.le_PHigh.text()) * PCF  # Update high pressure in the model.
        self.Model.p_low = float(self.View.le_PLow.text()) * PCF    # Update low pressure in the model.
        
        # Update turbine inlet condition: Temperature if specified, or quality (x=1) by default.
        T = float(self.View.le_TurbineInletCondition.text())  # Turbine inlet temperature or quality.
        self.Model.t_high = None if self.View.rdo_Quality.isChecked() else (T if self.Model.SI else UC.F_to_C(T))
        
        # Update turbine isentropic efficiency based on user input.
        self.Model.turbine_eff = float(self.View.le_TurbineEff.text())

        # Recalculate cycle efficiency with new parameters and update the view.
        self.calc_efficiency()
        self.updateView()


    def updateUnits(self):
    """
    Updates the unit system (SI or English) in the model based on the current selection in the view,
    and then instructs the view to refresh all unit-dependent displays. This method ensures that
    changing the units immediately reflects across all relevant GUI components without altering
    the underlying model data.
    """
    # Update the model's unit flag based on the GUI selection.
    self.Model.SI = self.View.rb_SI.isChecked()
    # Call the view's method to update all GUI elements to display in the correct units.
    self.View.updateUnits(Model=self.Model)

def selectQualityOrTHigh(self):
    """
    Triggers the view to update its display based on the selection between specifying the turbine
    inlet condition by quality or temperature. This method facilitates dynamic updates to the
    GUI in response to user interaction.
    """
    # Delegate to the view the task of adjusting GUI components for the quality/temperature selection.
    self.View.selectQualityOrTHigh(self.Model)

def setNewPHigh(self):
    """
    Informs the view to update its display based on a new high pressure value specified by the user.
    This is part of handling dynamic user inputs to adjust the simulation parameters.
    """
    # Trigger an update in the view to reflect the new high pressure setting.
    self.View.setNewPHigh(self.Model)

def setNewPLow(self):
    """
    Instructs the view to update its display to accommodate a new low pressure value input by the user.
    This method allows for real-time adjustments to simulation parameters through the GUI.
    """
    # Notify the view to update its elements to show the new low pressure setting.
    self.View.setNewPLow(self.Model)

def calc_efficiency(self):
    """
    Calculates the efficiency of the Rankine cycle based on current model parameters. This involves
    computing state properties at key points in the cycle (turbine inlet, turbine exit, pump inlet,
    and pump exit) and using these to calculate work done by the turbine, work required by the pump,
    and heat added to the cycle.

    The method uses a single instance of a steam properties object (held within the model) to
    perform these calculations, ensuring consistency and reusability of the steam properties data.

    :return: The calculated efficiency of the Rankine cycle as a percentage.
    """
    steam = self.Model.steam  # Reference to the steam object for property calculations.

    # Determine state properties at key points in the cycle based on model inputs.
    # State 1: Turbine inlet. Could be superheated or saturated vapor.
    self.Model.state1 = steam.getState(P=self.Model.p_high, x=1.0 if self.Model.t_high is None else None, T=self.Model.t_high, name='Turbine Inlet')
    
    # State 2: Turbine exit. Calculated for isentropic and actual conditions based on turbine efficiency.
    self.Model.state2s = steam.getState(P=self.Model.p_low, s=self.Model.state1.s, name="Turbine Exit")
    h2 = self.Model.state1.h - self.Model.turbine_eff * (self.Model.state1.h - self.Model.state2s.h) if self.Model.turbine_eff < 1.0 else None
    self.Model.state2 = steam.getState(P=self.Model.p_low, h=h2, name="Turbine Exit") if h2 is not None else self.Model.state2s
    
    # State 3: Pump inlet. Saturated liquid at low pressure.
    self.Model.state3 = steam.getState(P=self.Model.p_low, x=0, name='Pump Inlet')
    
    # State 4: Pump exit. Calculated as saturated liquid at high pressure, assuming isentropic compression.
    self.Model.state4 = steam.getState(P=self.Model.p_high, s=self.Model.state3.s, name='Pump Exit')

    # Perform efficiency calculations.
    self.Model.turbine_work = self.Model.state1.h - self.Model.state2.h
    self.Model.pump_work = self.Model.state4.h - self.Model.state3.h
    self.Model.heat_added = self.Model.state1.h - self.Model.state4.h
    self.Model.efficiency = 100.0 * (self.Model.turbine_work - self.Model.pump_work) / self.Model.heat_added
    
    return self.Model.efficiency
    def updateView(self):
        """
        This is a pass-through function that calls and identically named function in the View, but passes along the
        Model as an argument.
        :param args: A tuple of Widgets that get unpacked and updated in the view
        :return:
        """
        self.buildDataForPlotting()
        self.View.outputToGUI(Model=self.Model)

    def setRankine(self, p_low=8, p_high=8000, t_high=None, eff_turbine=1.0, name='Rankine Cycle'):
        '''
        Set model values for rankine power cycle.  If t_high is not specified, the State 1
        is assigned x=1 (saturated steam @ p_high).  Otherwise, use t_high to find State 1.
        :param p_low: the low pressure isobar for the cycle in kPa
        :param p_high: the high pressure isobar for the cycle in kPa
        :param t_high: optional temperature for State1 (turbine inlet) in degrees C
        :param eff_turbine: isentropic efficiency of the turbine
        :param name: a convenient name
        '''
        self.Model.p_low = p_low
        self.Model.p_high = p_high
        self.Model.t_high = t_high
        self.Model.name = name
        self.Model.efficiency = None
        self.Model.turbine_eff = eff_turbine
        self.Model.turbine_work = 0
        self.Model.pump_work = 0
        self.Model.heat_added = 0
        self.Model.state1 = None  # entrance to turbine
        self.Model.state2s = None  # entrance to condenser (isentropic turbine)
        self.Model.state2 = None  # entrance to condenser (non-isentropic turbine)
        self.Model.state3 = None  # entrance to pump (saturated liquid at plow)
        self.Model.state4 = None  # entrance to boiler (isentropic)

    def print_summary(self):
        """
        A pass-thrugh method for accessing View and passing Model.
        :return:
        """
        self.View.print_summary(Model=self.Model)

    def buildVaporDomeData(self, nPoints=500):
        """
        I'll build the vapor dome from just above the triple point up to the critical point
        :param nPoints:
        :return:
        """
        steam = self.Model.steam
        tp = triplePt_PT()
        cp = criticalPt_PT()
        steam.state.p = cp.p
        steam.state.t = cp.t
        steam.calcState_1Phase()
        critProps = dc(steam.state)
        P = np.logspace(math.log10(tp.p * 1.001), math.log10(cp.p * 0.99), nPoints)
        for p in P:
            sat = steam.getsatProps_p(p)
            self.Model.satLiqPlotData.addPt((sat.tsat, p, sat.uf, sat.hf, sat.sf, sat.vf))
            self.Model.satVapPlotData.addPt((sat.tsat, p, sat.uf, sat.hg, sat.sg, sat.vg))
        self.Model.satLiqPlotData.addPt((critProps.t, critProps.p, critProps.u, critProps.h, critProps.s, critProps.v))
        self.Model.satVapPlotData.addPt((critProps.t, critProps.p, critProps.u, critProps.h, critProps.s, critProps.v))

    def buildDataForPlotting(self):
        """
        I want to create data for plotting the Rankine cycle.  The key states are:
        State 1.  Entrance to Turbine (either saturated vapor or superheated steam at p_High)
        State 2.  Entrance to Condenser (probably two-phase at p_Low)
        State 3.  Entrance to the pump (saturated liquid at p_Low)
        State 4.  Entrance to the boiler (sub-cooled liquid at p_High, isentropic pump)

        I want to create h, s, v, p, T data between states 1-2, 2-3, 3-4, 4-1
        I'll piece together an upperCurve data set from 3-4 + 4-1 + 1-2
        The lowerCurve data set is 2-3
        :return:
        """
        # clear out any old data
        self.Model.upperCurve.clear()
        self.Model.lowerCurve.clear()

        # get saturated properties at PHigh and PLow
        satPLow = self.Model.steam.getsatProps_p(self.Model.p_low)
        satPHigh = self.Model.steam.getsatProps_p(self.Model.p_high)

        steam = self.Model.steam

        # region build upperCurve
        # region states from 3-4
        nPts = 15
        DeltaP = (satPHigh.psat - satPLow.psat)
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=(satPLow.psat + z * DeltaP), s=satPLow.sf)
            self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
        # endregion

        # region states from 4-1
        # first from T4 to T5 where T5 is the saturated liquid at p_High
        T4 = state.t
        T5 = satPHigh.tsat
        DeltaT = (T5 - T4)
        nPts = 20
        P = satPHigh.psat
        for n in range(nPts - 1):
            z = n * 1.0 / (nPts - 2)
            T = T4 + z * DeltaT
            if T < T5:
                state = steam.getState(P=P, T=T)
                self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(satPHigh.psat, x=z)
            self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
        if self.Model.state1.t > (satPHigh.tsat + 1):
            T6 = satPHigh.tsat
            DeltaT = self.Model.state1.t - T6
            for n in range(0, nPts):
                z = n * 1.0 / (nPts - 1)
                if z > 0:
                    state = steam.getState(satPHigh.psat, T=T6 + z * DeltaT)
                    self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
        # endregion

        # region states between 1 and 2
        # I'm assuming a linear change in Pressure from P1 to P2, along with linear change in s,
        # but not sure of details inside the turbine, so this is just a guess.
        s1 = self.Model.state1.s
        s2 = self.Model.state2.s
        P1 = self.Model.state1.p
        P2 = self.Model.state2.p
        Deltas = s2 - s1
        DeltaP = P2 - P1
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=P1 + z * DeltaP, s=s1 + z * Deltas)
            self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
        # endregion
        # endregion

        # region build lowerCurve between states 2 and 3
        x2 = self.Model.state2.x
        state = self.Model.state2
        # account for possibility that T>TSatPHigh
        if state.t > satPLow.tsat:
            nPts = 20
            DeltaT = (state.t - satPLow.tsat) / nPts
            self.Model.lowerCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
            for n in range(nPts):
                t = self.Model.state2.t - n * DeltaT
                if t > satPLow.tsat:
                    state = steam.getState(P=satPLow.psat, T=t)
                    self.Model.lowerCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

        nPts = len(self.Model.upperCurve.t)
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=satPLow.psat, x=(1.0 - z) * x2)
            self.Model.lowerCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
        # endregion
        pass

    def updatePlot(self):
        self.View.plot_cycle_XY(Model=self.Model)


# endregion

# region function definitions
def main():
    RC=rankineController()
    RC.setRankine(8*UC.kpa_to_bar,8000*UC.kpa_to_bar,t_high=500, eff_turbine=0.9,name='Rankine Cycle - Superheated at turbine inlet')
    #t_high is specified
    #if t_high were not specified, then x_high = 1 is assumed
    eff=RC.calc_efficiency()
    print(eff)
    RC.print_summary()
    RC.plot_cycle_TS()


# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion
