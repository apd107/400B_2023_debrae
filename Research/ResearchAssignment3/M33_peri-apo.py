import numpy as np


M31_data = np.genfromtxt('Orbit_M31.txt', dtype=None, names=True)
M33_data = np.genfromtxt('Orbit_M33.txt', dtype=None, names=True)

def vector_diff(vector1, vector2):
    diff = vector1 - vector2
    magnitude = np.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
    return magnitude

peri_t = 0
peri_dist = 100000 
apo_t = 0
apo_dist = 0
for i in range(160):
    M33_x = M33_data['x'][i]
    M33_y = M33_data['y'][i]
    M33_z = M33_data['z'][i]
    M33_pos_vector = np.array([M33_x, M33_y, M33_z])
    M31_x = M31_data['x'][i]
    M31_y = M31_data['y'][i]
    M31_z = M31_data['z'][i]
    M31_pos_vector = np.array([M31_x, M31_y, M31_z])
    pos = vector_diff(M31_pos_vector, M33_pos_vector)
    if pos > apo_dist:
        apo_t = M33_data['t'][i]
        apo_dist = pos
    if pos < apo_dist:
        peri_t = M33_data['t'][i]
        peri_dist = pos

print(peri_t)
print(peri_dist)
print(apo_t)
print(apo_dist)
