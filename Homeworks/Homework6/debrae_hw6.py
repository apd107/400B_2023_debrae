#%%
from Orbits import OrbitCOM
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u

# I already wrote the orbit files with the functions
# below. Because it takes some time, I have commented 
# them out to save time witht the rest of the assignment. 
# If the orbit files need to be rewritten simply un-comment
# the below lines and update the parameters. 
########################################################

# OrbitCOM('MW', 0, 800, 5)
# OrbitCOM('M31', 0, 800, 5)
# OrbitCOM('M33', 0, 800, 5)

def vector_diff(vector1, vector2):
    diff = vector1 - vector2
    magnitude = np.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
    return magnitude


MW_data = np.genfromtxt('Orbit_MW.txt', dtype=None, names=True)
M31_data = np.genfromtxt('Orbit_M31.txt', dtype=None, names=True)
M33_data = np.genfromtxt('Orbit_M33.txt', dtype=None, names=True)

########################################################
# For the sake of ease and clarity I've commented out the
# code used to create the log plots. Uncomment these pieces
# to create a log plot
########################################################

# M31 - MW
pos_magnitudes = []
vel_magnitudes = []
for i in range(160):
    MW_x = MW_data['x'][i]
    MW_y = MW_data['y'][i]
    MW_z = MW_data['z'][i]
    MW_vx = MW_data['vx'][i]
    MW_vy = MW_data['vy'][i]
    MW_vz = MW_data['vz'][i]
    MW_pos_vector = np.array([MW_x, MW_y, MW_z])
    MW_vel_vector = np.array([MW_vx, MW_vy, MW_vz])
    M31_x = M31_data['x'][i]
    M31_y = M31_data['y'][i]
    M31_z = M31_data['z'][i]
    M31_vx = M31_data['vx'][i]
    M31_vy = M31_data['vy'][i]
    M31_vz = M31_data['vz'][i]
    M31_pos_vector = np.array([M31_x, M31_y, M31_z])
    M31_vel_vector = np.array([M31_vx, M31_vy, M31_vz])
    # magnitude 
    pos = vector_diff(M31_pos_vector, MW_pos_vector)
    vel = vector_diff(M31_vel_vector, MW_vel_vector)
    pos_magnitudes.append(pos)
    vel_magnitudes.append(vel)


pos_magnitudes = np.array(pos_magnitudes)*u.kpc
vel_magnitudes = np.array(vel_magnitudes)*(u.km/u.s)
time = MW_data['t']

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('M31-MW Separation')
# plt.title('M31-MW Separation (Logy)')
plt.xlabel('time (Gyr)', fontsize=13)
plt.ylabel('Position (kpc)', fontsize=13)
ax.plot(time, pos_magnitudes, color='Red', linewidth=3, label='M31-MW', linestyle='-')
# ax.semilogy(time, pos_magnitudes, color='Red', linewidth=3, label='M31-MW', linestyle='-')
plt.legend()
plt.savefig('M31-MW_Separation.png')
# plt.savefig('M31-MW_Separation_Logy.png')

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('M31-MW Relative Velocity')
# plt.title('M31-MW Relative Velocity (Logy)')
plt.xlabel('time (Gyr)', fontsize=13)
plt.ylabel('Velocity (km/s)', fontsize=13)
ax.plot(time, vel_magnitudes, color='Red', linewidth=3, label='M31-MW', linestyle='-')
# ax.semilogy(time, vel_magnitudes, color='Red', linewidth=3, label='M31-MW', linestyle='-')
plt.legend()
plt.savefig('M31-MW_Relative_Velocity.png')
# plt.savefig('M31-MW_Relative_Velocity_Logy.png')

########################################################
# For the sake of ease and clarity I've commented out the
# code used to create the log plots. Uncomment these pieces
# to create a log plot
########################################################

# M31 - M33
pos_magnitudes = []
vel_magnitudes = []
for i in range(160):
    M33_x = M33_data['x'][i]
    M33_y = M33_data['y'][i]
    M33_z = M33_data['z'][i]
    M33_vx = M33_data['vx'][i]
    M33_vy = M33_data['vy'][i]
    M33_vz = M33_data['vz'][i]
    M33_pos_vector = np.array([M33_x, M33_y, M33_z])
    M33_vel_vector = np.array([M33_vx, M33_vy, M33_vz])
    M31_x = M31_data['x'][i]
    M31_y = M31_data['y'][i]
    M31_z = M31_data['z'][i]
    M31_vx = M31_data['vx'][i]
    M31_vy = M31_data['vy'][i]
    M31_vz = M31_data['vz'][i]
    M31_pos_vector = np.array([M31_x, M31_y, M31_z])
    M31_vel_vector = np.array([M31_vx, M31_vy, M31_vz])
    # magnitude 
    pos = vector_diff(M31_pos_vector, M33_pos_vector)
    vel = vector_diff(M31_vel_vector, M33_vel_vector)
    pos_magnitudes.append(pos)
    vel_magnitudes.append(vel)

pos_magnitudes = np.array(pos_magnitudes)*u.kpc
vel_magnitudes = np.array(vel_magnitudes)*(u.km/u.s)
time = M31_data['t']
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('M31-M33 Separation')
# plt.title('M31-M33 Separation (Logy)')
plt.xlabel('time (Gyr)', fontsize=13)
plt.ylabel('Position (kpc)', fontsize=13)
ax.plot(time, pos_magnitudes, color='Blue', linewidth=3, label='M31-M33', linestyle='-')
# ax.semilogy(time, pos_magnitudes, color='Blue', linewidth=3, label='M31-M33', linestyle='-')
plt.legend()
plt.savefig('M31-M33_Separation.png')
# plt.savefig('M31-M33_Separation_Logy.png')


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('M31-M33 Relative Velocity')
# plt.title('M31-M33 Relative Velocity (Logy)')
plt.xlabel('time (Gyr)', fontsize=13)
plt.ylabel('Velocity (km/s)', fontsize=13)
ax.plot(time, vel_magnitudes, color='Blue', linewidth=3, label='M31-M33', linestyle='-')
# ax.semilogy(time, vel_magnitudes, color='Blue', linewidth=3, label='M31-M33', linestyle='-')
plt.legend()
plt.savefig('M31-M33_Relative_Velocity.png')
# plt.savefig('M31-M33_Relative_Velocity_Logy.png')



# %%
