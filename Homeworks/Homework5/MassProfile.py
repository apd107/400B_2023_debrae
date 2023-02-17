#%%
# Author: Aidan DeBrae
# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
import matplotlib
import matplotlib.pyplot as plt
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass


class MassProfile:
# Class to analyze the mass profile of a given galaxy
# at a given snap number

    def __init__(self, galaxy, snap):
        # add a string of the filenumber to the value '000'
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = str(galaxy) + '_' + ilbl + '.txt'
        # initialize data from galaxy file
        self.time, self.total, self.data = Read(self.filename)
        self.m = self.data['m']
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        # galaxy name 
        self.gname = str(galaxy)
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

    def MassEnclosed(self, ptype, radii):
        '''
        This function calculates the mass enclosed a galaxy for a given particle type.

        Inputs:
            ptype: int
                Integer which represents the particle type
            radii: array
                An array which contains the radii of the enclosed particles

        Outputs:
            mass_sums: array
                An array of masses at each radius interval
        '''
        # Calaculate center of mass from the class of previous homework
        galaxy_COM = CenterOfMass(self.filename, ptype)
        COM_pos = galaxy_COM.COM_P(0.1)
        x_COM, y_COM, z_COM = COM_pos[0], COM_pos[1], COM_pos[2]
        

        # Establish index of the particle type
        index = np.where(self.data['type'] == ptype)

        # Calculate the radius & mass of the particles according to the 
        # particle type index
        x = self.x[index]
        y = self.y[index]
        z = self.z[index]
        particle_masses = self.m[index]

        # redefine using the COM
        x_new = x - x_COM
        y_new = y - y_COM
        z_new = z - z_COM
        particle_radii = np.sqrt(x_new**2 + y_new**2 + z_new**2)

        # Initialize mass array according to the length of the radius
        mass_sums = np.zeros(len(radii))


        # This loop iterates over the radii list and computes the enclosed 
        # mass of the particles within the given radius 
        for i in range(len(radii)):
            # Establish the index of the particles within the radius
            index = np.where(particle_radii <= radii[i]*u.kpc)
            # Calulate masses for the enclosed particles according to
            # the new index within the radius
            masses = particle_masses[index]
            # Sum and store the mass of the particles
            mass_sums[i] = np.sum(masses)
        return mass_sums * 1e10 * u.Msun

    def MassEnclosedTotal(self, radii):
        '''
        This function calculates the total mass enclosed a galaxy for a given particle type.

        Inputs:
            radii: array
                An array which contains the radii of the enclosed particles

        Outputs:
            mass_sums: array
                An array of masses at each radius interval
        '''
        # Calculate the mass of each type of particle using the 
        # MassEnclosed function within the class. DO NOT compute
        # the bulge mass for M33 as it doesn't have a bulge
        if self.gname == 'M33':
            halo_mass = self.MassEnclosed(1, radii)
            disk_mass = self.MassEnclosed(2, radii)
            total_mass = halo_mass + disk_mass
        else:
            halo_mass = self.MassEnclosed(1, radii)
            disk_mass = self.MassEnclosed(2, radii)
            bulge_mass = self.MassEnclosed(3, radii)
            total_mass = halo_mass + disk_mass + bulge_mass
        return total_mass

    def HernquistMass(self, r, a, Mhalo):
        '''
        This function calculates the Hernquist Mass.

        Inputs:
            r: array 
                An array which contains the radii of the enclosed particles
            scale factor: float
                A scle factor in kpc
            Mhalo: array
                Halo masses of the galaxy
                
        Outputs:
            mass_sums: array
                An array of masses at each radius interval
        '''
        # density: (M*a/2*pi*r)*(1/(r+a)**3)
        # mass: (Mhalo*r**2)/(r+a)**2
        M = ((Mhalo*(r*u.kpc)**2)/((a+r)*u.kpc)**2)
        return M

    def CircularVelocity(self, ptype, radii):
        '''
        This function calculates circular velocity.

        Inputs:
            ptype: int
                Integer which represents the particle type
            radii: array
                An array which contains the radii of the enclosed particles
                
        Outputs:
            circular_velocities: array
                An array of velocities at each radius interval
        '''
        # Circular velocity: v = (G*M/r)
        masses = self.MassEnclosed(ptype, radii)
        circular_velocities = np.sqrt((self.G*masses)/radii)
        return np.around(circular_velocities, 2)
    
    def CircularVelocityTotal(self, radii):
        '''
        This function calculates the total circular velocity.

        Inputs:
            radii: array
                An array which contains the radii of the enclosed particles
                
        Outputs:
            circular_velocities: array
                An array of velocities at each radius interval
        '''
        # Circular velocity: v = (G*M/r)
        mass_totals = self.MassEnclosedTotal(radii)
        circular_velocities = np.sqrt((self.G*mass_totals)/radii)
        return np.around(circular_velocities, 2)
    
    def HernquistVCirc(self, r, a, Mhalo):
        '''
        This function calculates the Hernquist velocity.

        Inputs:
            r: array 
                An array which contains the radii of the enclosed particles
            scale factor: float
                A scle factor in kpc
            Mhalo: array
                Halo masses of the galaxy
                
        Outputs:
            circular_velocities: array
                An array of masses at each radius interval
        '''
        # Define the mass using the Hernquist mass function
        M = self.HernquistMass(r, a, Mhalo)
        circular_velocities = np.sqrt((self.G*M)/r*u.kpc)
        return np.around(circular_velocities, 2)




# %%
