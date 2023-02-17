#%%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from MassProfile import MassProfile

### HW ANSWERS ###

# MW
MW = MassProfile('MW', 0)
radii = np.arange(0.1, 30, 0.1)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Rotation Curve for MW')
plt.xlabel('Distance from COM (kpc)', fontsize=13)
plt.ylabel('Circular Velocity (km/s)', fontsize=13)

# Initialize circular velocities for each galaxy component & total
Mhalo = MW.MassEnclosed(1, radii) # Halo mass for Hernquist circular velocity
haloVcirc = MW.CircularVelocity(1, radii)
diskVcirc = MW.CircularVelocity(2, radii)
bulgeVcirc = MW.CircularVelocity(3, radii)
totalVcirc = MW.CircularVelocityTotal(radii)

# Initialize hernquist mass using the halo mass
# Determine the scale factor
scale_factor = 0.01
hernquistVcirc = MW.HernquistVCirc(radii, scale_factor, Mhalo)

ax.semilogy(radii, haloVcirc, color='orange', linewidth=3, 
            label='Halo v$_{circ}$', linestyle='-')
ax.semilogy(radii, diskVcirc, color='blue', linewidth=3, 
            label='Disk v$_{circ}$', linestyle='-')
ax.semilogy(radii, bulgeVcirc, color='green', linewidth=3, 
            label='Bulge v$_{circ}$', linestyle='-')
ax.semilogy(radii, totalVcirc, color='red', linewidth=3, 
            label='Total v$_{circ}$', linestyle='-')
ax.semilogy(radii, hernquistVcirc, color='black', linewidth=4, 
            label='Hernquist v$_{circ}$, a=0.01', linestyle='dotted')
ax.legend()
plt.savefig('MW_rotation_curve.png')

# M31
M31 = MassProfile('M31', 0)
radii = np.arange(0.1, 30, 0.1)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Rotation Curve for M31')
plt.xlabel('Distance from COM (kpc)', fontsize=13)
plt.ylabel('Circular Velocity (km/s)', fontsize=13)

# Initialize circular velocities for each galaxy component & total
Mhalo = M31.MassEnclosed(1, radii) # Halo mass for Hernquist circular velocity
haloVcirc = M31.CircularVelocity(1, radii)
diskVcirc = M31.CircularVelocity(2, radii)
bulgeVcirc = M31.CircularVelocity(3, radii)
totalVcirc = M31.CircularVelocityTotal(radii)

# Initialize hernquist mass using the halo mass
# Determine the scale factor
scale_factor = 0.01
hernquistVcirc = M31.HernquistVCirc(radii, scale_factor, Mhalo)

ax.semilogy(radii, haloVcirc, color='orange', linewidth=3, 
            label='Halo v$_{circ}$', linestyle='-')
ax.semilogy(radii, diskVcirc, color='blue', linewidth=3, 
            label='Disk v$_{circ}$', linestyle='-')
ax.semilogy(radii, bulgeVcirc, color='green', linewidth=3, 
            label='Bulge v$_{circ}$', linestyle='-')
ax.semilogy(radii, totalVcirc, color='red', linewidth=3, 
            label='Total v$_{circ}$', linestyle='-')
ax.semilogy(radii, hernquistVcirc, color='black', linewidth=4, 
            label='Hernquist v$_{circ}$, a=0.01', linestyle='dotted')
ax.legend()
plt.savefig('M31_rotation_curve.png')

# M31
M33 = MassProfile('M33', 0)
radii = np.arange(0.1, 30, 0.1)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Rotation Curve for M33')
plt.xlabel('Distance from COM (kpc)', fontsize=13)
plt.ylabel('Circular Velocity (km/s)', fontsize=13)

# Initialize circular velocities for each galaxy component & total
Mhalo = M33.MassEnclosed(1, radii) # Halo mass for Hernquist circular velocity
haloVcirc = M33.CircularVelocity(1, radii)
diskVcirc = M33.CircularVelocity(2, radii)
totalVcirc = M33.CircularVelocityTotal(radii)

# Initialize hernquist mass using the halo mass
# Determine the scale factor
scale_factor = 0.01
hernquistVcirc = M33.HernquistVCirc(radii, scale_factor, Mhalo)

ax.semilogy(radii, haloVcirc, color='orange', linewidth=3, 
            label='Halo v$_{circ}$', linestyle='-')
ax.semilogy(radii, diskVcirc, color='blue', linewidth=3, 
            label='Disk v$_{circ}$', linestyle='-')
ax.semilogy(radii, totalVcirc, color='red', linewidth=3, 
            label='Total v$_{circ}$', linestyle='-')
ax.semilogy(radii, hernquistVcirc, color='black', linewidth=4, 
            label='Hernquist v$_{circ}$, a=0.01', linestyle='dotted')
ax.legend()
plt.savefig('M33_rotation_curve.png')










# %%
