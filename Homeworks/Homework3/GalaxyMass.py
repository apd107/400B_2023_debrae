# GalaxyMass.py
# Author: Aidan DeBrae
# This program has been written so that each mass component of a 
# galaxy of interest can be calculated. 

import numpy as np
import astropy.units as u
from ReadFile import Read

def ComponentMass(filename, type):
    '''
    This function will return the total mass of any desired
    Galaxy component. 

        Inputs: 
            filename: 'string'
                Name of the data file containing Galaxy information. 
            type: 'float'
                A float value which represents the type of the particle
                of interest. 

        Outputs:
            mass: 'float'
                Total mass of the galaxy component of interest.
    '''
    time, total_particles, data = Read(filename) # Reading in the information of the Galaxy file 
    particle_index = np.where(data['type'] == type) # Create an array of indexes which correspond to the particle type
    particles = data[particle_index] # Create an array of the particles of interest
    mass_values = np.array(particles['m']) # Create an array of the particle masses 
    total_mass = np.around(sum(mass_values)*1e10*u.Msun, 3) # Add together all masses for a total mass
    return total_mass


