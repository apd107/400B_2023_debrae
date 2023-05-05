# import modules
import numpy as np

class Approach:
    def __init__(self, orbit_file1, orbit_file2):
        # generate data from the orbit files written in hw6
        self.M31_data = np.genfromtxt(orbit_file1, dtype=None, names=True)
        self.M33_data = np.genfromtxt(orbit_file2, dtype=None, names=True)

    # create an easy to use vector difference function
    def vector_diff(self, vector1, vector2):
        diff = vector1 - vector2
        magnitude = np.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
        return magnitude

    def approach(self):
        '''
        This function calculates the time and distance of the starting, 
        apocenter, and pericenter events.  

        INPUT
        -----
        None
        
        OUTPUT
        ------
        starting_dist: 'float'
            Distance of M33 to M31 at the starting point of the simulation.
        
        peri_t: 'float'
            Time at which the pericenter approach occurs. 

        peri_dist: 'float'
            Distance at which the pericenter approach occurs.

        apo_t: 'float'
            Time at which the apocenter approach occurs. 

        apo_dist: 'float'
            Distance at which the apocenter approach occurs. 
        '''
        # I'm going to try and calculate the peri/apo-center approaches
        # of the galaxies in an attempt to observe the most dramatic affects.
        start_dist = 0
        peri_t = 0
        peri_dist = 100000 
        apo_t = 0
        apo_dist = 0
        for i in range(160):
            M33_x = self.M33_data['x'][i]
            M33_y = self.M33_data['y'][i]
            M33_z = self.M33_data['z'][i]
            M33_pos_vector = np.array([M33_x, M33_y, M33_z])
            M31_x = self.M31_data['x'][i]
            M31_y = self.M31_data['y'][i]
            M31_z = self.M31_data['z'][i]
            M31_pos_vector = np.array([M31_x, M31_y, M31_z])
            dist = self.vector_diff(M31_pos_vector, M33_pos_vector)
            if i == 0:
                start_dist = dist
            if dist > apo_dist:
                apo_t = self.M33_data['t'][i]
                apo_dist = dist
            if dist < apo_dist:
                peri_t = self.M33_data['t'][i]
                peri_dist = dist
        return start_dist, peri_t, peri_dist, apo_t, apo_dist