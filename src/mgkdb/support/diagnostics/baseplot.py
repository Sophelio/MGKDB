""" Base module for plotting

Contains:
- Matplotib imports and settings

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from packaging.version import parse as parse_version

if parse_version(mpl.__version__) >= parse_version("1.5"):
    from cycler import cycler


class Plotting:
    """ Base class for plotting routines

    Contains the parameters used for the plots such as
    line colors, and font sizes
    Sets the line and point marker widths/sizes
    Sets the default color cycle (respecting an API change in matplotlib 1.5)
    """

    def __init__(self):
        # Matplotlib 2.x changes the default plot styles
        if parse_version(mpl.__version__) >= parse_version("2.0"):
            # Colormap for quantities >0 (e.g. heat flux)
            self.cmap_unidirect = mpl.cm.viridis
            mpl.rc('font', size=15)
        else:
            self.cmap_unidirect = mpl.cm.viridis
            mpl.rc('font', size=15)
        # for quantities with critical value (e.g. 0 for phi)
        self.cmap_bidirect = mpl.cm.bwr
        self.color_list = plt.cm.Dark2(np.linspace(0, 1.0, 9))  # Plot line colors
        # Set some basic plot properties
        # Those are set globally everytime a Plotting object is instantiated!
        mpl.rc('lines', linewidth=2, markeredgewidth=1.5, markersize=10)
        # This sets the figure frames to transparent both on screen and in calls to savefig
        mpl.rc('figure', facecolor=(1, 1, 1, 0))
        mpl.rc('figure', frameon=False)
        mpl.rc('savefig', facecolor=(1, 1, 1, 0))
        #mpl.rc('savefig', frameon=False)
        # mpl 1.5 introduces a new, more flexible prop_cycle parameter, so different line styles
        # can be defined in the else case as well
        if parse_version(mpl.__version__) <= parse_version("1.4"):
            mpl.rc('axes', color_cycle=self.color_list)
        else:
            mpl.rc('axes', prop_cycle=cycler('color', self.color_list))
        # Dictionary that converts internal names to well-formatted LaTeX for the plots
        self.titles = {'omts': r'$a/L_T$', 'omns': r'$a/L_n$', 'ns': r'$n$', 'Ts': r'$T$',
                       'Gammaturb': r'$\Gamma_{\mathrm{turb}}/\Gamma_{gB}$',
                       'Gammanc': r'$\Gamma_{\mathrm{NC}}/\Gamma_{gB}$',
                       'Qturb': r'$Q_{\mathrm{turb}}/Q_{gB}$', 'Qnc': r'$Q_{\mathrm{NC}}/Q_{gB}$',
                       'Piturb': r'$\Pi_{\mathrm{turb}}/\Pi_{gB}$',
                       'Pinc': r'$\Pi_{\mathrm{turb}}/\Pi_{gB}$', 
                       "jbs": r'$j_\mathrm{BS}$',
                       "phi": r"$\phi$", 
                       "apar": r"$A_\parallel$", 
                       "dens": r"$n$",
                       "tpar": r"$T_\parallel$", 
                       "tperp": r"$T_\perp$",
                       "qpar": r"$q_\parallel + 1.5p_0 u_\parallel$",
                       "qperp": r"$q_\perp + p_0 u_\parallel$", 
                       "upar": r"$u_\parallel$",
                       "densI1": r"$nI1$",
                       "TparI1": r"$TparI1$",
                       "TppI1": r"$TppI1$",
                       "bpar":  r"$B_\parallel$"}

    def plot_contour_quant(self, ax, fig, x, y, data, is_complex, xlabel=None, ylabel=None,
                           title=None, is_bidirect=False, axcb=None):
        """ Plot a 2d quantity as a color map

        :param ax: Axis object to plot into
        :param fig: Figure object ax belongs to (needed for colorbar)
        :param x: x-axis grid
        :param y: y-axis grid
        :param data: 2d-array in x-y order
        :param is_complex: Takes the absolute value of data before plotting (for complex arrays)
        :param xlabel: Label for x-axis, None will leave it unchanged
        :param ylabel: Label for y-axis, None will leave it unchanged
        :param title: Axis title, None will leave it unchanged
        :param is_bidirect: Use a bidirectional colormap for data
        :param axcb: Axes object of the colorbar, needed for time series
        :return: Axes object of the colorbar for reuse
        """
        mapdata = np.abs(np.squeeze(data)).T if is_complex else np.squeeze(data).T
        if is_bidirect:
            posval = np.abs(np.nanmax(mapdata))
            negval = np.abs(np.nanmin(mapdata))
            vmax = posval if posval > negval else negval
            vmin = -vmax
            cmap = self.cmap_bidirect
        else:
            vmax, vmin = None, None
            cmap = self.cmap_unidirect
        cm1 = ax.pcolormesh(x, y, mapdata, cmap=cmap, vmax=vmax, vmin=vmin)
        ax.set_rasterization_zorder(z=-10)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title)
        axcb = fig.colorbar(cm1, cax=axcb).ax
        return axcb


    def plot_1d_quant(self, ax, fig, x, data, is_complex, xlabel=None, title=None, axcb=None):
        """ Plot a 1d quantity

        :param ax: Axis object to plot into
        :param fig: Figure object ax belongs to (needed for colorbar)
        :param x: x-axis grid
        :param data: 1d-arra
        :param is_complex: Takes the absolute value of data before plotting (for complex arrays)
        :param xlabel: Label for x-axis, None will leave it unchanged
        :param title: Axis title, None will leave it unchanged
        """

        ax.plot(x, np.absolute(data.T) if is_complex else data.T, color='b')
        ax.set_rasterization_zorder(z=-10)

        if xlabel:
            ax.set_xlabel(xlabel)
        if title:
            ax.set_title(title)
        axcb = ax.axes
        return axcb

    def plot_contour_quant_timeseries(self, ax, fig, x, y, data, is_complex, times, xlabel=None,
                                      ylabel=None, title=None, is_bidirect=False):
        """ Plot a 2d quantity as a navigatable series of color maps

        :param ax: Axis object to plot into
        :param fig: Figure object ax belongs to (needed for colorbar)
        :param x: x-axis grid
        :param y: y-axis grid
        :param data: 2d-array in x-y order
        :param is_complex: Takes the absolute value of data before plotting (for complex arrays)
        :param times:
        :param xlabel: Label for x-axis, None will leave it unchanged
        :param ylabel: Label for y-axis, None will leave it unchanged
        :param title: Axis title, None will leave it unchanged
        :param is_bidirect: Use a bidirectional colormap for data
        """
        def key_event(event):
            nonlocal icurr_time
            keyresult = self._keyhandler(event, icurr_time, len(times))
            if keyresult is None:
                return
            else:
                icurr_time = keyresult
            # Clear both the axes (leaving labels) and the colorbar object
            for artist in ax.lines+ax.collections:
                artist.remove()
            axcb.clear()
            self.plot_contour_quant(ax, fig, x, y, data[icurr_time], is_complex,
                                    title="{}, t={:.3f}".format(title, times[icurr_time]),
                                    is_bidirect=is_bidirect, axcb=axcb)
            fig.canvas.draw()

        icurr_time = 0
        axcb = self.plot_contour_quant(ax, fig, x, y, data[icurr_time], is_complex, xlabel, ylabel,
                                title="{}, t={:.3f}".format(title, times[icurr_time]),
                                is_bidirect=is_bidirect)
        fig.canvas.mpl_connect('key_press_event', key_event)




    def plot_1d_quant_timeseries(self, ax, fig, x, data, is_complex, times, xlabel=None, title=None):
        """ Plot a 2d quantity as a navigatable series of color maps

        :param ax: Axis object to plot into
        :param fig: Figure object ax belongs to (needed for colorbar)
        :param x: x-axis grid
        :param data: 1d-array
        :param is_complex: Takes the absolute value of data before plotting (for complex arrays)
        :param times:
        :param xlabel: Label for x-axis, None will leave it unchanged
        """
        def key_event(event):
            nonlocal icurr_time
            keyresult = self._keyhandler(event, icurr_time, len(times))
            if keyresult is None:
                return
            else:
                icurr_time = keyresult
            # Clear both the axes (leaving labels)
            for artist in ax.lines+ax.collections:
                artist.remove()
            axcb.clear()
            self.plot_1d_quant(ax, fig, x, data[icurr_time], is_complex, xlabel,
                                    title="{}, t={:.3f}".format(title, times[icurr_time]), axcb=axcb)
            fig.canvas.draw()

        icurr_time = 0
        axcb = self.plot_1d_quant(ax, fig, x, data[icurr_time], is_complex, xlabel,
                           title="{}, t={:.3f}".format(title, times[icurr_time]))
        fig.canvas.mpl_connect('key_press_event', key_event)

    @staticmethod
    def _keyhandler(event, icurr_time, timelen):
        new_ind = icurr_time
        if event.key == "pageup":
            new_ind += 10
        elif event.key == "pagedown":
            new_ind -= 10
        elif event.key == "right":
            new_ind += 1
        elif event.key == "left":
            new_ind -= 1

        elif event.key == "home":
            new_ind = 0
        elif event.key == "end":
            new_ind = timelen-1
        else:
            return None
        return new_ind%timelen

    @staticmethod
    def _minlogaxis(arr):
        """ Calculate lower limit for plotting on a log axis

        lower Integer * 10^minimumexponent

        :param arr: Array to base the limit on
        """
        minim = arr[arr != 0].min()  # ignore 0
        min_exp = np.floor(np.log10(minim))
        min_mant = np.floor(minim*10 ** (-min_exp))
        floormin = min_mant*10 ** min_exp
        return floormin

    @staticmethod
    def _maxlogaxis(arr):
        """ Calculate upper limit for plotting on a log axis

         higher Integer * 10^maximumexponent

        :param arr: Array to base the limit on
        """
        maxim = arr[arr != 0].max()
        max_exp = np.floor(np.log10(maxim))
        max_mant = np.ceil(maxim*10 ** (-max_exp))
        ceilmax = max_mant*10 ** max_exp
        return ceilmax

    @staticmethod
    def gyrobohm_SI(rundata, quantity):
        """ Convert gyroBohm unit to SI unit

        :param rundata: Run object
        :param quantity: Name of the quantity to calculate gB for
        """
        # pylint: disable=invalid-name
        pnt = rundata
        elementary_charge = 1.60217662e-19
        temp_unit = 1e3*elementary_charge  # temp is in keV
        dens_unit = 1e19  # m**(-3)
        mass_unit = 1.6726219e-27  # proton mass
        # Magnetic field and reference length are already in T and m
        nref = pnt.pars["Lref"]*dens_unit
        Tref = pnt.pars["Tref"]*temp_unit
        mref = pnt.pars["mref"]*mass_unit
        cref = np.sqrt(Tref/mref)
        rhoref = mref*cref/elementary_charge/pnt.pars["Bref"]
        Gamma_gB = nref*cref*rhoref ** 2/pnt.pars["Lref"] ** 2
        gBunits = {"Ts": Tref, "ns": nref, "omts": 1, "omns": 1, "Gammanc": Gamma_gB/1e19,
                   "Gammaturb": Gamma_gB/1e19, "Qturb": Gamma_gB*Tref/1e3, "Qnc": Gamma_gB*Tref/1e3,
                   "Pinc": Gamma_gB*cref*mref, "Piturb": Gamma_gB*cref*mref,
                   "jbs": pnt.pars["Bref"]*nref*cref*rhoref/pnt.pars["Lref"]}
        return gBunits[quantity]
