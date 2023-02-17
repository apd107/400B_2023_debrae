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
plt.title('Mass Profile for MW')
plt.xlabel(u'Radius (kpc)', fontsize=13)
plt.ylabel(u'Log(Mass Enclosed ($M_\u2609$))', fontsize=13)

# Initialize masses for each galaxy component & the total mass
halo_masses = MW.MassEnclosed(1, radii)
disk_masses = MW.MassEnclosed(2, radii)
bulge_masses = MW.MassEnclosed(3, radii)
total_masses = MW.MassEnclosedTotal(radii)
# Initialize hernquist mass using the halo mass
# Determine the scale factor
scale_factor = 0.1
hernquist_masses = MW.HernquistMass(radii, scale_factor, halo_masses)

ax.semilogy(radii, halo_masses, color='orange', linewidth=3, 
            label='Halo Mass', linestyle='-')
ax.semilogy(radii, disk_masses, color='blue', linewidth=3, 
            label='Disk Mass', linestyle='-')
ax.semilogy(radii, bulge_masses, color='green', linewidth=3, 
            label='Bulge Mass', linestyle='-')
ax.semilogy(radii, total_masses, color='red', linewidth=3, 
            label='Total Mass', linestyle='-')
ax.semilogy(radii, hernquist_masses, color='black', linewidth=4, 
            label='Hernquist Mass', linestyle='dotted')
ax.legend()
plt.savefig('MW_mass_profile.png')

# M31
M31 = MassProfile('M31', 0)
radii = np.arange(0.1, 30, 0.1)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Mass Profile for M31')
plt.xlabel(u'Radius (kpc)', fontsize=13)
plt.ylabel(u'Log(Mass Enclosed ($M_\u2609$))', fontsize=13)

# Initialize masses for each galaxy component & the total mass
halo_masses = M31.MassEnclosed(1, radii)
disk_masses = M31.MassEnclosed(2, radii)
bulge_masses = M31.MassEnclosed(3, radii)
total_masses = M31.MassEnclosedTotal(radii)
# Initialize hernquist mass using the halo mass
# Determine the scale factor
scale_factor = 0.1
hernquist_masses = M31.HernquistMass(radii, scale_factor, halo_masses)

ax.semilogy(radii, halo_masses, color='orange', linewidth=3, 
            label='Halo Mass', linestyle='-')
ax.semilogy(radii, disk_masses, color='blue', linewidth=3, 
            label='Disk Mass', linestyle='-')
ax.semilogy(radii, bulge_masses, color='green', linewidth=3, 
            label='Bulge Mass', linestyle='-')
ax.semilogy(radii, total_masses, color='red', linewidth=3, 
            label='Total Mass', linestyle='-')
ax.semilogy(radii, hernquist_masses, color='black', linewidth=4, 
            label='Hernquist Mass', linestyle='dotted')
ax.legend()
plt.savefig('M31_mass_profile.png')

# M33
M33 = MassProfile('M33', 0)
radii = np.arange(0.1, 30, 0.1)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Mass Profile for M33')
plt.xlabel(u'Radius (kpc)', fontsize=13)
plt.ylabel(u'Log(Mass Enclosed ($M_\u2609$))', fontsize=13)

# Initialize masses for each galaxy component & the total mass
halo_masses = M33.MassEnclosed(1, radii)
disk_masses = M33.MassEnclosed(2, radii)
total_masses = M33.MassEnclosedTotal(radii)
# Initialize hernquist mass using the halo mass
# Determine the scale factor
scale_factor = 0.1
hernquist_masses = M33.HernquistMass(radii, scale_factor, halo_masses)

ax.semilogy(radii, halo_masses, color='orange', linewidth=3, 
            label='Halo Mass', linestyle='-')
ax.semilogy(radii, disk_masses, color='blue', linewidth=3, 
            label='Disk Mass', linestyle='-')
ax.semilogy(radii, total_masses, color='red', linewidth=3, 
            label='Total Mass', linestyle='-')
ax.semilogy(radii, hernquist_masses, color='black', linewidth=4, 
            label='Hernquist Mass', linestyle='dotted')
ax.legend()
plt.savefig('M33_mass_profile.png')









# %%
