# import modules
import numpy as np
from GalaxyMass import ComponentMass

class Jacobi:
    def __init__(self, M33_snap_file, M31_snap_file, dist):
        M33_filename = 'HighRes/M33/' + str(M33_snap_file)
        M31_filename = 'HighRes/M31/' + str(M31_snap_file)
        self.r = dist
        
        M33_halo_mass = ComponentMass(M33_filename, 1.0)
        M33_disk_mass = ComponentMass(M33_filename, 2.0)
        M33_bulge_mass = ComponentMass(M33_filename, 3.0)
        self.M33_mass = M33_halo_mass + M33_disk_mass + M33_bulge_mass

        M31_halo_mass = ComponentMass(M31_filename, 1.0)
        M31_disk_mass = ComponentMass(M31_filename, 2.0)
        M31_bulge_mass = ComponentMass(M31_filename, 3.0)
        self.M31_mass = M31_halo_mass + M31_disk_mass + M31_bulge_mass
    
    def Jacobi_R(self, Msat, Mhost, r):
        '''
        This function calculates the jacobi radius of M33.  

        INPUT
        -----
        Msat: 'float'
            Mass of M33 at this given time. 
        
        Mhost: 'float'
            Mass of M31 at this given time. 
        
        r: 'float' 
            Distance of M33 to M31. 
        
        OUTPUT
        ------
        Rj: 'float' 
            The jacobi radius of M33 at the given snapshot time 
            in kpc. 
        '''
        # Rj = r*(Msat/(2*Mhost))**1/3
        # Mhost mass changes make sure to write about it
        Rj = r*(Msat/(2*Mhost))**(1/3)
        return Rj





# get M31/M31_794.txt Desktop/school/Astro_400B/research/HighRes/M31
# get M33/M33_794.txt Desktop/school/Astro_400B/research/HighRes/M33