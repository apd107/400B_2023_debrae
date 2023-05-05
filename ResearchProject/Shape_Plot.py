# import modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import scipy.optimize as so

class Shape_Plot:
    def __init__(self):
        # plot the particle density and contour fitting lines
        self.fig, (self.ax1, self.ax2)= plt.subplots(1, 2, figsize=(10, 5))

    def find_confidence_interval(self, x, pdf, confidence_level):
        return pdf[pdf > x].sum() - confidence_level

    def density_contour(self, xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
        """ Create a density contour plot.
        Parameters
        ----------
        xdata : numpy.ndarray
        ydata : numpy.ndarray
        nbins_x : int
            Number of bins along x dimension
        nbins_y : int
            Number of bins along y dimension
        ax : matplotlib.Axes (optional)
            If supplied, plot the contour to this axis. Otherwise, open a new figure
        contour_kwargs : dict
            kwargs to be passed to pyplot.contour()
            
        Example Usage
        -------------
        density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
        e.g.:
        density_contour(xD, yD, 80, 80, ax=ax, 
            colors=['red','orange', 'yellow', 'orange', 'yellow'])

        """

        H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
        x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
        y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

        pdf = (H*(x_bin_sizes*y_bin_sizes))
        
        X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
        Z = pdf.T
        fmt = {}
        
        
        # Contour Levels Definitions
        one_sigma = so.brentq(self.find_confidence_interval, 0., 1., args=(pdf, 0.80))
        two_sigma = so.brentq(self.find_confidence_interval, 0., 1., args=(pdf, 0.99))

        # Array of Contour levels. Adjust according to the above
        levels = [one_sigma, two_sigma][::-1]
        
        # contour level labels  Adjust accoding to the above.
        strs = ['0.80', '1.0'][::-1]


        
        ###### 
        
        if ax == None:
            contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
            for l, s in zip(contour.levels, strs):
                fmt[l] = s
            plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

        else:
            contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
            for l, s in zip(contour.levels, strs):
                fmt[l] = s
            ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
        
        return contour

    def plot(self, x, y, z, x_ellipse, y_ellipse, z_ellipse1, z_ellipse2, title,
             xlim_lo, xlim_hi, ylim_lo, ylim_hi, bins, xadj=0, yadj=0):
        '''
        This function plots the 2D density of dark matter halos, the density
        contour lines, and the ellipse which fits to the density contour lines. 

        INPUT
        -----
        x: 'np.array'
            The x coordinate of the particles.

        y: 'np.array'
            The y coordinate of the particles.
        
        z: 'np.array'
            The z coordinate of the particles.
        
        x_ellipse: 'np.array'
            The x coordinate of the ellipse.

        y_ellipse: 'np.array'
            The y coordinate of the ellipse.
        
        z_ellipse(1 & 2): 'np.array'
            The z coordinate of the ellipse.

        _lim_lo: 'integer'
            I did not include each input. This values simply determines the lower
            limit of the given axis. 
        
        _lim_hi: 'integer'
            I did not include each input. This values simply determines the upper
            limit of the given axis.

        bins: 'integer'
            The bin count for the histogram to plot the density of particles.

        _adj: 'integer'
            I did not include each input. This value adjusts the given axis value 
            according to the center of the ellipse. 
        
        OUTPUT
        ------
        No output, just plots. 
        '''
        self.ax1.set_title(title + ' x vs. z')
        self.ax1.set_xlabel('x (kpc)', fontsize=12)
        self.ax1.set_ylabel('z (kpc)', fontsize=12)
        self.ax1.hist2d(x, z, bins=bins, norm=LogNorm(), cmap='inferno')
        self.density_contour(x, z, 50, 50, ax=self.ax1, colors=['cyan', 'limegreen'])
        self.ax1.plot(x_ellipse+xadj, z_ellipse1+yadj, color='red', linestyle='dashdot', linewidth=3)
        self.ax1.set_xlim(xlim_lo, xlim_hi)
        self.ax1.set_ylim(ylim_lo, ylim_hi)
        #######################
        self.ax2.set_title(title + ' y vs. z')
        self.ax2.set_xlabel('y (kpc)', fontsize=12)
        self.ax2.hist2d(y, z, bins=bins, norm=LogNorm(), cmap='inferno')
        self.density_contour(y, z, 50, 50, ax=self.ax2, colors=['cyan', 'limegreen'])
        self.ax2.plot(y_ellipse+xadj, z_ellipse2+yadj, color='red', linestyle='dashdot', linewidth=3)
        self.ax2.set_xlim(xlim_lo, xlim_hi)
        self.ax2.set_ylim(ylim_lo, ylim_hi)
        label_size = 12
        matplotlib.rcParams['xtick.labelsize'] = label_size 
        matplotlib.rcParams['ytick.labelsize'] = label_size
        plt.savefig('Start_shape_plot')