#%%
from M33AnalyticOrbit import M33AnalyticOrbit
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u

# Create an instance of the M33AnalyticOrbit class to 
# write a file which contains the position and velocity 
# information over time. 
M33_orbit = M33AnalyticOrbit('M33_orbit.txt')
M33_orbit.OrbitIntegration(0, 0.1, 10)
# Generate the information from the orbit file created by 
# the class. 
M33 = np.genfromtxt('M33_orbit.txt', dtype=None, names=True)
new_time = M33['t']
new_pos = (np.sqrt(M33['x']**2 + M33['y']**2 + M33['z']**2))*u.kpc
new_vel = (np.sqrt(M33['vx']**2 + M33['vy']**2 + M33['vz']**2))*(u.km/u.s)

# Copy HW6 M31 - M33 plotting scheme
def vector_diff(vector1, vector2):
    diff = vector1 - vector2
    magnitude = np.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
    return magnitude


hw6_M31_data = np.genfromtxt('Orbit_M31.txt', dtype=None, names=True)
hw6_M33_data = np.genfromtxt('Orbit_M33.txt', dtype=None, names=True)
hw6_pos = []
hw6_vel = []
for i in range(160):
    M33_x = hw6_M33_data['x'][i]
    M33_y = hw6_M33_data['y'][i]
    M33_z = hw6_M33_data['z'][i]
    M33_vx = hw6_M33_data['vx'][i]
    M33_vy = hw6_M33_data['vy'][i]
    M33_vz = hw6_M33_data['vz'][i]
    M33_pos_vector = np.array([M33_x, M33_y, M33_z])
    M33_vel_vector = np.array([M33_vx, M33_vy, M33_vz])
    M31_x = hw6_M31_data['x'][i]
    M31_y = hw6_M31_data['y'][i]
    M31_z = hw6_M31_data['z'][i]
    M31_vx = hw6_M31_data['vx'][i]
    M31_vy = hw6_M31_data['vy'][i]
    M31_vz = hw6_M31_data['vz'][i]
    M31_pos_vector = np.array([M31_x, M31_y, M31_z])
    M31_vel_vector = np.array([M31_vx, M31_vy, M31_vz]) 
    pos = vector_diff(M31_pos_vector, M33_pos_vector)
    vel = vector_diff(M31_vel_vector, M33_vel_vector)
    hw6_pos.append(pos)
    hw6_vel.append(vel)

hw6_pos = np.array(hw6_pos)*u.kpc
hw6_vel = np.array(hw6_vel)*(u.km/u.s)
hw6_time = hw6_M31_data['t']


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('M31-M33 Separation')
plt.xlabel('time (Gyr)', fontsize=13)
plt.ylabel('Position (kpc)', fontsize=13)
plt.xlim([0, 10.1])
ax.plot(hw6_time, hw6_pos, color='Black', linewidth=3, label='hw6_M31-M33', linestyle='-')
ax.plot(new_time, new_pos, color='Red', linewidth=3, label='M33_orbit', linestyle='-.')
plt.savefig('Orbit_Separation(M33-M31).png')
plt.legend()

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('M31-M33 Relative Velocity')
plt.xlabel('time (Gyr)', fontsize=13)
plt.ylabel('Velocity (km/s)', fontsize=13)
plt.xlim([0, 10.1])
ax.plot(hw6_time, hw6_vel, color='Black', linewidth=3, label='hw6_M31-M33', linestyle='-')
ax.plot(new_time, new_vel, color='Red', linewidth=3, label='M33_orbit', linestyle='-.')
plt.savefig('Orbit_Velocity(M33-M31).png')
plt.legend()

### QUESTIONS ###
print()
print('Q2')
print('The first pericenter approach is nearly identical for both plots.',
      'Afterwards, however, the M33AnalyticOrbit plot exhibits a much, much',
      'longer orbital period than the previous plot from homework 6.',
      'Additionally, the separation plot from M33AnalyticOrbit exhibits',
      'much larger separation values.')
print()
print('Q3')
print('I think there are a few missing physics components from the M33AnalyticOrbit.',
      'First, we are only taking the particle information from snapshot 0 which does',
      'not account for the change in the system over time. Possible affects of this',
      'particular issue are not accurately calculating the changes in mass of the galaxy',
      'or even the exchange of particles within the system. This could affect the orbit',
      'in substantial ways. Secondly, we are treating the system as point masses.',
      'Calculating the orbit from point masses excludes the potential influence from',
      'the extended parts of the halo.')
print()
print('Q4')
print('I believe there are two ways in which the effects of the Milky Way could be included.',
      'First, the full range of the snapshots should be included to capture the changes of',
      'the system over time. The simulation was ran such that the future of the Local Group',
      'is illustrated and thus accounts for the affects of the entirety of the system including',
      'the Milky Way. Additonally, we should calculate the acceleration of M33 from it\'s host',
      'halo M31 as well as MW. We know that the two galaxies will eventually collide and thus',
      'must consider how the Milky Way might influence M33.')
print()

# %%
