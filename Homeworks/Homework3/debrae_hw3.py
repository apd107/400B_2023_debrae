# debrae_hw3.py
# Author: Aidan DeBrae
# This program calculates the mass component of each galaxy 
# according to homework 3. 

import numpy as np
from GalaxyMass import ComponentMass

# Using the function 'ComponentMass' each component of 
# MW (halo, disk, bulge) is calculated below
mw_halo_mass = ComponentMass('MW_000.txt', 1.0)
mw_disk_mass = ComponentMass('MW_000.txt', 2.0)
mw_bulge_mass = ComponentMass('MW_000.txt', 3.0)
# To calculate the total mass of the galaxy each component
# must be added together. 
mw_total_mass = mw_halo_mass + mw_disk_mass + mw_bulge_mass
# The baryon fraction is computed by calculating the ratio
# of the disk & bulge to the total mass of the galaxy. 
mw_f_bar = np.around((mw_disk_mass + mw_bulge_mass) / (mw_total_mass), 3)
print('MW halo mass: ' + f'{mw_halo_mass:.3e}')
print('MW disk mass: ' + f'{mw_disk_mass:.3e}')
print('MW bulge mass: ' + f'{mw_bulge_mass:.3e}')
print('MW total mass: ' + f'{mw_total_mass:.3e}')
print('MW fbar mass: ' + f'{mw_f_bar:.3e}')

print()

# Using the function 'ComponentMass' each component of 
# M31 (halo, disk, bulge) is calculated below
m31_halo_mass = ComponentMass('M31_000.txt', 1.0)
m31_disk_mass = ComponentMass('M31_000.txt', 2.0)
m31_bulge_mass = ComponentMass('M31_000.txt', 3.0)
# To calculate the total mass of the galaxy each component
# must be added together. 
m31_total_mass = m31_halo_mass + m31_disk_mass + m31_bulge_mass
# The baryon fraction is computed by calculating the ratio
# of the disk & bulge to the total mass of the galaxy.
m31_f_bar = np.around((m31_disk_mass + m31_bulge_mass) / (m31_total_mass), 3)
print('M31 halo mass: ' + f'{m31_halo_mass:.3e}')
print('M31 disk mass: ' + f'{m31_disk_mass:.3e}')
print('M31 bulge mass: ' + f'{m31_bulge_mass:.3e}')
print('M31 total mass: ' + f'{m31_total_mass:.3e}')
print('M31 fbar mass: ' + f'{m31_f_bar:.3e}')
print()

# Using the function 'ComponentMass' each component of 
# M33 (halo, disk, bulge, fbar) is calculated below
m33_halo_mass = ComponentMass('M33_000.txt', 1.0)
m33_disk_mass = ComponentMass('M33_000.txt', 2.0)
m33_bulge_mass = ComponentMass('M33_000.txt', 3.0)
# To calculate the total mass of the galaxy each component
# must be added together. 
m33_total_mass = m33_halo_mass + m33_disk_mass + m33_bulge_mass
# The baryon fraction is computed by calculating the ratio
# of the disk & bulge to the total mass of the galaxy.
m33_f_bar = np.around((m33_disk_mass + m33_bulge_mass) / (m33_total_mass), 3)
print('M33 halo mass: ' + f'{m33_halo_mass:.3e}')
print('M33 disk mass: ' + f'{m33_disk_mass:.3e}')
print('M33 bulge mass: ' + f'{m33_bulge_mass:.3e}')
print('M33 total mass: ' + f'{m33_total_mass:.3e}')
print('M33 fbar mass: ' + f'{m33_f_bar:.3e}')
print()

# To calculate each mass component of the Local Group each galaxy 
# component is simply added together.
lg_halo_mass = mw_halo_mass + m31_halo_mass + m33_halo_mass
lg_disk_mass = mw_disk_mass + m31_disk_mass + m33_disk_mass
lg_bulge_mass = mw_bulge_mass + m31_bulge_mass + m33_bulge_mass
lg_total_mass = lg_halo_mass + lg_disk_mass + lg_bulge_mass
lg_stellar = (mw_disk_mass + mw_bulge_mass +
              m31_disk_mass + m31_bulge_mass + 
              m33_disk_mass + m33_bulge_mass)
lg_f_bar = (lg_stellar/lg_total_mass)
print('Local Group halo mass: ' + f'{lg_halo_mass:.3e}')
print('Local Group disk mass: ' + f'{lg_disk_mass:.3e}')
print('Local Group bulge mass: ' + f'{lg_bulge_mass:.3e}')
print('Local Group total mass: ' + f'{lg_total_mass:.3e}')
print('M33 fbar mass: ' + f'{lg_f_bar:.3e}')
print()


