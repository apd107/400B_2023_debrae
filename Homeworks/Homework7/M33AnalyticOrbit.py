# Author: Aidan DeBrae
# M33AnalyticOrbit program developed using the Homework7_Template.py   

# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from IPython.display import Latex
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass

# M33AnalyticOrbit

class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): 
        '''
        A class which initializes the relative position and velocity 
        between M33 & M31 and calculates the projected orbit of M33.

        INPUT
        -----
        filename: 'string'
            Desired file name for the M33 orbit information. 
        
        OUTPUT
        ------
        A data file containing the time, position, and velocity
        information of M33's orbit. 
        '''
        # get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        # define the output file name according to the input name
        self.filename = filename
        # get the current pos/vel of M33  
        M33_COM = CenterOfMass('M33_VLowRes/M33_000.txt', 2)
        M33_pos_COM = M33_COM.COM_P(0.1, 4)
        M33_vel_COM = M33_COM.COM_V(M33_pos_COM[0],
                                    M33_pos_COM[1],
                                    M33_pos_COM[2])
        
        M31_COM = CenterOfMass('M31_VLowRes/M31_000.txt', 2)
        M31_pos_COM = M31_COM.COM_P(0.1, 2)
        M31_vel_COM = M31_COM.COM_V(M31_pos_COM[0],
                                    M31_pos_COM[1],
                                    M31_pos_COM[2])

        
        # define the relative position and velocity (M33-M31)
        self.r = (M33_pos_COM - M31_pos_COM).value
        self.v = (M33_vel_COM - M31_vel_COM).value
        
        # define the scale factor and mass of each component of M31 

        # disk
        self.rdisk = (5*u.kpc).value
        self.Mdisk = ComponentMass('M31_VLowRes/M31_000.txt', 2.0).value
        # bulge
        self.rbulge = (1*u.kpc).value
        self.Mbulge = ComponentMass('M31_VLowRes/M31_000.txt', 3.0).value
        # Halo
        self.rhalo = 61.58
        self.Mhalo = ComponentMass('M31_VLowRes/M31_000.txt', 1.0).value
    
    
    def HernquistAccel(self, M, r_a, r):  
        '''
        This function calculates the acceleration according to a 
        Hernquist mass profile. 

        INPUT
        -----
        M: 'float'
            Mass of the galaxy component of which the acceleration 
            is calculated. 
        r_a: 'float'
            Scale factor of the galaxy component. 
        r: 'vector' (preferably an np.array)
            The x,y,z position vector of a galaxy. 
        
        OUTPUT
        ------
        a_Hern: 'vector' (np.array)
            An x,y,z acceleration vector according to the Hernquist
            mass profile.   
        '''
        rmag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        a_Hern =  ((-self.G*M)/(rmag*(r_a + rmag)**2))*r  
        
        return a_Hern
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r): 
        '''
        This function calculates the acceleration according to a 
        Miyamoto Nagai 1975 profile. 

        INPUT
        -----
        M: 'float'
            Mass of the galaxy component of which the acceleration 
            is calculated. 
        r_d: 'float'
            Scale factor of the galaxy component. 
        r: 'vector' (preferably an np.array)
            The x,y,z position vector of a galaxy. 
        
        OUTPUT
        ------
        a_MN: 'vector' (np.array)
            An x,y,z acceleration vector according to the Miyamoto
            Nagai 1975 profile.   
        '''
        R = np.sqrt(r[0]**2 + r[1]**2)
        z_d = r_d/5.0
        B = r_d + np.sqrt(r[2]**2 + z_d**2)
        z_stuff = np.array([1, 1, (B/np.sqrt(r[2]**2 + z_d**2))])

        a_MN = ((-self.G*M)/((R**2 + B**2)**1.5))*r*z_stuff

        return a_MN
     
    
    def M31Accel(self, r):
        '''
        This function calculates the total acceleration from each 
        acceleration profile. 

        INPUT
        ----- 
        r: 'vector' (preferably an np.array)
            The x,y,z position vector of a galaxy as it's being updated. 
        
        OUTPUT
        ------
        total_accel: 'vector' (np.array)
            An x,y,z acceleration vector.   
        '''

        # Call the previous acceleration functions for the halo, bulge and disk
        halo_accel = self.HernquistAccel(self.Mhalo, self.rhalo, r)
        bulge_accel = self.HernquistAccel(self.Mbulge, self.rbulge, r)
        disk_accel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r) 

        total_accel = halo_accel + bulge_accel + disk_accel
        return total_accel
    
    
    
    def LeapFrog(self, dt, r, v): 
        '''
        This function is an "integrator function" which updates the 
        position and velocity of the galaxy according to it's most
        current position, velocity, and acceleration

        INPUT
        -----
        dt: 'float'
            Value by which the time of integration is increased.   
        r: 'vector' (preferably an np.array)
            The most current x,y,z position vector of the galaxy. 
        v: 'vector' (preferably an np.array)
            The most current x,y,z velocity vector of the galaxy.
        
        OUTPUT
        ------
        rnew, vnew: 'vectors' (np.array)
            Two vectors which correspond to the updated position 
            and velocity of the galaxy according to the integration.
        '''
        # predict the position at the next half timestep
        rhalf = r + v*(dt/2)
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        ahalf = self.M31Accel(rhalf)
        vnew = v + ahalf*dt
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew*(dt/2)
        
        return rnew, vnew
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        '''
        This function makes use of the integrator function to
        update and store the position and velocity information
        over time. 

        INPUT
        -----
        t0: 'float'
            Starting time value.   
        dt: 'float' 
            Value by which the time of integration is increased. 
        tmax: 'float'
            Max time value a.k.a the stopping time for integration.
        
        OUTPUT
        ------
        N/A
        '''
        # initialize the time to the input starting time
        t = t0
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2, 7))
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r), *tuple(self.v)
        r = self.r
        v = self.v
        # initialize a counter for the orbit.  
        i = 1
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while t < tmax: 
            t += dt
            # store the new time in the first column of the ith row
            orbit[i, 0] = t
            # advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector 
            r, v = self.LeapFrog(dt, r, v)
            # store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i, 1:4] = r
            orbit[i, 4:7] = v
            # update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1

        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
