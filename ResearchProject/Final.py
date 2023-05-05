#%%
import numpy as np
from Approach import Approach
from Halo_Shape import Halo_Shape
from Shape_Plot import Shape_Plot
from Shape_Characterization import Shape_Char
from Jacobi import Jacobi
from math import trunc

# Here, I use the Approach class and the orbit files to calculate the starting, 
# apocenter, and pericenter distance and time. This is so that I can see the 
# most dramatic changes to the halo shape as the local group evolves toward the 
# collision event. 
M33_M31_approach = Approach('New_Orbit_M31.txt', 'New_Orbit_M33.txt')
start_dist, peri_t, peri_dist, apo_t, apo_dist = M33_M31_approach.approach()
start_dist = np.around(start_dist)
apo_dist = np.around(apo_dist)
peri_dist = np.around(peri_dist)

# Print the time and position values for the starting, pericenter, and apocenter 
# to plot the proper snapshat file and distance.
print('Start distance (kpc): ' + str(start_dist))
print('Pericenter time snapshot (Gyr): ' + str(peri_t))
print('Pericenter distance (kpc): ' + str(peri_dist))
print('Apocenter time snapshot (Gyr): ' + str(apo_t))
print('Apocenter distance (kpc): ' + str(apo_dist))

# cut down particles outside of jacobi radius (fix contours)
# download new data files


### START ###
# Initialize the particles components using the Halo_Shape class. These are going 
# to be the components used to plot the shape of the halo at a given snapshot.
M33_start_shape = Halo_Shape('M33_000.txt')
start_r, start_v = M33_start_shape.r, M33_start_shape.v
# Rotate the halo to be oriented along the z-axis so that I can effectively 
# charactize the shape of the halo according to the dimensions of an ellipsoid. 
start_rrot, start_vrot = M33_start_shape.RotateFrame(start_r, start_v)

# Start Jacobi Radius
# Here, I use the Jacobi class to calculate the Jacobi radius of the halo at a given
# snapshot because the Jacobi radius encompasses the region which best represents
# the shape of the halo. 
start_system = Jacobi('M33_000.txt', 'M31_000.txt', start_dist)
start_Msat = start_system.M33_mass.value
start_Mhost = start_system.M31_mass.value
start_r = start_system.r
start_Rj = start_system.Jacobi_R(start_Msat, start_Mhost, start_r)
start_Rj = np.around(start_Rj)
limit = start_Rj + 10
print('Initial Jacobi radius (kpc): ' + str(start_Rj))

# SCALE SHOULD BE CHANGED ONCE EVERYTHING IS SET UP AND IS READY TO BE ADJUSTED
# Now, I'll define the exact x, y, z coordinates of the entire halo according to the 
# given snapshot and after the halo has been rotated along the z-axis.
start_x = start_rrot[:,0]
start_y = start_rrot[:,1]
start_z = start_rrot[:,2]
start_r = np.sqrt(start_x**2 + start_y**2 + start_z**2)
# Here, I find the index values of the particles that are withihn the jacobi radius. 
# This not only improves the resolution of the plots but also allows for a more precise
# definiton of the region that represents the shape of the halo. 
start_jacobi_index = np.where(start_r < start_Rj)
start_jacobi_index = start_jacobi_index[0]
start_jx = []
start_jy = []
start_jz = []
# Using the index values that are found from the Jacobi radius I'll now extract the 
# x, y, z coordinates of the particles. I'll also conver these back into np.arrays.
for i in start_jacobi_index:
    start_jx.append(start_x[i])
    start_jy.append(start_y[i])
    start_jz.append(start_z[i])
start_jx = np.array(start_jx)
start_jy = np.array(start_jy)
start_jz = np.array(start_jz)
# Here, I use the ellipse function defined within the Halo_Shape class to define the 
# components of the ellipse that will be used to define the shape of the halo along 
# with the radii of the major, minor, and intermediate axes of the ellipsoid. 
# An important feature of the ellipse function is scaling the ellipse such that it
# fits the contour line of the halo. 
start_x_ell, start_z_ellx, start_x_rad, start_z_radx = M33_start_shape.ellipse(start_jx, start_jz, 0.8, 0.8)
start_y_ell, start_z_elly, start_y_rad, start_z_rady = M33_start_shape.ellipse(start_jy, start_jz, 0.809, 0.8)


# PLOT THE SHAPE
# Now, using all of the above definitions of each component I will plot the 2D projection
# of the density of particles, the corresponding contour line fit, and the ellipse
# shape of the halo. 
start_shape_plot = Shape_Plot()
start_shape_plot.plot(start_jx, start_jy, start_jz, start_x_ell, start_y_ell, start_z_ellx, start_z_elly, 
                      'Start', -limit, limit, -limit, limit, 80, 0, 0)

# CHARACTERIZE SHAPE
# Here, using the Shape_Characterization class I will define the dimensions of the shape
# and the define it as being triaxial, prolate, or oblate
start_shape = Shape_Char(start_x_rad, start_y_rad, start_z_radx)
start_shape.T_calc()
start_shape = Shape_Char(start_x_rad, start_y_rad, start_z_rady)
start_shape.T_calc()
print('The starting shape of the halo is Triaxial.')


# ALL OF THE ABOVE METHODS ARE REPEATED FOR SNAPSHOT 184 AND 794 WHICH CORRESPOND TO THE 
# APOCENTER AND PERICENTER RESPECTIVELY. 

### APOCENTER ###
M33_apo_shape = Halo_Shape('M33_184.txt')
apo_r, apo_v = M33_apo_shape.r, M33_apo_shape.v
apo_rrot, apo_vrot = M33_apo_shape.RotateFrame(apo_r, apo_v)

# Apocenter Jacobi Radius
apo_system = Jacobi('M33_184.txt', 'M31_184.txt', apo_dist)
apo_Msat = apo_system.M33_mass.value
apo_Mhost = apo_system.M31_mass.value
apo_r = apo_system.r
apo_Rj = apo_system.Jacobi_R(apo_Msat, apo_Mhost, apo_r)
apo_Rj = np.around(apo_Rj)
limit = apo_Rj + 5
print('Apocenter Jacobi Radius (kpc): ' + str(apo_Rj))

# SCALE SHOULD BE CHANGED ONCE EVERYTHING IS SET UP AND IS READY TO BE ADJUSTED
apo_x = apo_rrot[:,0]
apo_y = apo_rrot[:,1]
apo_z = apo_rrot[:,2]
apo_r = np.sqrt(apo_x**2 + apo_y**2 + apo_z**2)
apo_jacobi_index = np.where(apo_r < apo_Rj)
apo_jacobi_index = apo_jacobi_index[0]
apo_jx = []
apo_jy = []
apo_jz = []
for i in apo_jacobi_index:
    apo_jx.append(apo_x[i])
    apo_jy.append(apo_y[i])
    apo_jz.append(apo_z[i])
apo_jx = np.array(apo_jx)
apo_jy = np.array(apo_jy)
apo_jz = np.array(apo_jz)

apo_x_ell, apo_z_ellx, apo_x_rad, apo_z_radx = M33_apo_shape.ellipse(apo_jx, apo_jz, 0.75, 0.75)
apo_y_ell, apo_z_elly, apo_y_rad, apo_z_rady = M33_apo_shape.ellipse(apo_jy, apo_jz, 0.75, 0.75)


# PLOT THE SHAPE
apo_shape_plot = Shape_Plot()
apo_shape_plot.plot(apo_jx, apo_jy, apo_jz, apo_x_ell, apo_y_ell, apo_z_ellx, apo_z_elly,
                    'Apocenter', -limit, limit, -limit, limit, 80, 0, 0)


# CHARACTERIZE SHAPE
# Here, using the Shape_Characterization class I will define the dimensions of the shape
# and the define it as being triaxial, prolate, or oblate
apo_xz_shape = Shape_Char(apo_x_rad, apo_y_rad, apo_z_radx)
apo_xz_shape.T_calc()
apo_yz_shape = Shape_Char(apo_x_rad, apo_y_rad, apo_z_rady)
apo_yz_shape.T_calc()
print('The apocenter shape of the halo is Oblate.')



### PERICENTER ###
M33_peri_shape = Halo_Shape('M33_755.txt')
peri_r, peri_v = M33_peri_shape.r, M33_peri_shape.v
peri_rrot, peri_vrot = M33_peri_shape.RotateFrame(peri_r, peri_v)

# Pericenter Jacobi Radius
peri_system = Jacobi('M33_755.txt', 'M31_755.txt', peri_dist)
peri_Msat = peri_system.M33_mass.value
peri_Mhost = peri_system.M31_mass.value
peri_r = peri_system.r
peri_Rj = peri_system.Jacobi_R(peri_Msat, peri_Mhost, peri_r)
peri_Rj = np.around(peri_Rj)
limit = peri_Rj + 20
print('Pericenter Jacobi radius (kpc): ' + str(peri_Rj))

# SCALE SHOULD BE CHANGED ONCE EVERYTHING IS SET UP AND IS READY TO BE ADJUSTED
peri_x = peri_rrot[:,0]
peri_y = peri_rrot[:,1]
peri_z = peri_rrot[:,2]
peri_r = np.sqrt(peri_x**2 + peri_y**2 + peri_z**2)
peri_jacobi_index = np.where(peri_r < peri_Rj + 20)
peri_jacobi_index = peri_jacobi_index[0]
peri_jx = []
peri_jy = []
peri_jz = []
for i in peri_jacobi_index:
    peri_jx.append(peri_x[i])
    peri_jy.append(peri_y[i])
    peri_jz.append(peri_z[i])
peri_jx = np.array(peri_jx)
peri_jy = np.array(peri_jy)
peri_jz = np.array(peri_jz)

peri_x_ell, peri_z_ellx, peri_x_rad, peri_z_radx = M33_peri_shape.ellipse(peri_jx, peri_jz, 0.45, 0.45)
peri_y_ell, peri_z_elly, peri_y_rad, peri_z_rady = M33_peri_shape.ellipse(peri_jy, peri_jz, 0.45, 0.45)

# PLOT THE SHAPE
peri_shape_plot = Shape_Plot()
peri_shape_plot.plot(peri_jx, peri_jy, peri_jz, peri_x_ell, peri_y_ell, peri_z_ellx, peri_z_elly, 
                     'Pericenter', -limit, limit, -limit, limit, 40, 0, 0)


# CHARACTERIZE SHAPE
# Here, using the Shape_Characterization class I will define the dimensions of the shape
# and the define it as being triaxial, prolate, or oblate
peri_xz_shape = Shape_Char(peri_x_rad, peri_y_rad, peri_z_radx)
peri_xz_shape.T_calc()
peri_yz_shape = Shape_Char(peri_x_rad, peri_y_rad, peri_z_rady)
peri_yz_shape.T_calc()
print('The pericenter shape of the halo is Oblate.')

# %%


# centroiding is off (include in analysis)
# dont forget to include that i have not precisely 
# accounted for the change in mass 