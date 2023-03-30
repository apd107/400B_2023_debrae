# QUESTIONS #
# - What are the ideal constraints for the axes in my plot especially 
#   considering the halo particles seem to extend very far out?
# - What is the best resolution (i.e. bin size) for analyzing the
#   the halo?
# - Are there normal conventions for the contour line colors? 
# - Is it useful to rotate the halo on a given axis like what we did
#   for the disk in Lab 7? Or does the orientation not matter? 

#%%
# import modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import scipy.optimize as so


# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
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
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.60))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.82))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.90))

    # Array of Contour levels. Adjust according to the above
    levels = [one_sigma, two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    strs = ['0.60', '0.80', '0.90'][::-1]

    
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


def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel

# Create a COM of object for M33
M33_COM = CenterOfMass('M33_VLowRes/M33_000.txt', 1.0)

# Compute COM of M33 using halo particles
M33_COM_pos = M33_COM.COM_P(0.1, 4)
M33_COM_vel = M33_COM.COM_V(M33_COM_pos[0],
                            M33_COM_pos[1],
                            M33_COM_pos[2])

xH = M33_COM.x - M33_COM_pos[0].value 
yH = M33_COM.y - M33_COM_pos[1].value 
zH = M33_COM.z - M33_COM_pos[2].value

r_tot = np.sqrt(xH**2 + yH**2 + zH**2)

vxH = M33_COM.vx - M33_COM_vel[0].value 
vyH = M33_COM.vy - M33_COM_vel[1].value 
vzH = M33_COM.vz - M33_COM_vel[2].value

v_tot = np.sqrt(vxH**2 + vyH**2 + vzH**2)

r = np.array([xH,yH,zH]).T 
v = np.array([vxH,vyH,vzH]).T

rn, vn = RotateFrame(r, v)

fig, ax= plt.subplots(figsize=(10, 10))
plt.xlabel('x', fontsize=22)
plt.ylabel('z', fontsize=22)
plt.hist2d(rn[:,0], rn[:,2], bins=600, norm=LogNorm(), cmap='viridis')
density_contour(rn[:,0], rn[:,2], 80, 80, ax=ax, colors=['black', 'red', 'white'])
plt.xlabel('y', fontsize=22)
plt.ylabel('z', fontsize=22)
plt.hist2d(rn[:,1], rn[:,2], bins=600, norm=LogNorm(), cmap='viridis')
density_contour(rn[:,1], rn[:,2], 80, 80, ax=ax, colors=['black', 'red', 'white'])
plt.ylim(-500,500)
plt.xlim(-500,500)
plt.colorbar()
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size
# %%
