
# # Homework 7 Template
# 
# Rixin Li & G . Besla
#
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 

# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex
# import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass
# import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# # M33AnalyticOrbit

class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): 
        """ **** ADD COMMENTS """

        # get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        # define the output file name according to the input name
        self.filename = filename
        # get the current pos/vel of M33  
        M33_COM = CenterOfMass('M33_VLowRes/M33_000.txt', 2)
        M33_pos_COM = M33_COM.COM_P(0.1, 4).value # WHAT IS THE PROPER VOLDEC VALUE?
        M33_vel_COM = M33_COM.COM_V(M33_pos_COM[0],
                                    M33_pos_COM[1],
                                    M33_pos_COM[2]).value
        
        M31_COM = CenterOfMass('M31_VLowRes/M31_000.txt', 2)
        M31_pos_COM = M31_COM.COM_P(0.1, 2).value # WHAT IS THE PROPER VOLDEC VALUE?
        M31_vel_COM = M31_COM.COM_V(M31_pos_COM[0],
                                    M31_pos_COM[1],
                                    M31_pos_COM[2]).value
        
        # store the DIFFERENCE between the vectors posM33 - posM31
        self.r = M33_pos_COM - M31_pos_COM
        self.v = M33_vel_COM - M31_vel_COM
        
        # get the mass of each component in M31 
        # disk
        self.rdisk = (5*u.kpc).value
        self.Mdisk = ComponentMass('M31_VLowRes/M31_000.txt', 2.0)
        # bulge
        self.rbulge = (1*u.kpc).value
        self.Mbulge = ComponentMass('M31_VLowRes/M31_000.txt', 3.0)
        # Halo
        self.rhalo = 0.1
        self.Mhalo = ComponentMass('M31_VLowRes/M31_000.txt', 1.0)     
    
    
    def HernquistAccel(self, M, r_a, r):  
        """ **** ADD COMMENTS """

        # Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        # Store the Acceleration VECTOR
        Hern =  ((-self.G*M)/(rmag*(r_a + rmag)**2))*r  
        
        return Hern
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r): 
        """ **** ADD COMMENTS """
        
        R = np.sqrt(r[0]**2 + r[1]**2)
        z_d = r_d/5.0
        B = r_d + np.sqrt(r[2]**2 + z_d**2)
        z_stuff = np.array([1, 1, (B/np.sqrt(r[2]**2 + z_d**2))])

        MN = ((-self.G*M)/((R**2 + B**2)**1.5))*r*z_stuff

        return MN
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """ **** ADD COMMENTS """

        # SHOULD THERE BE AN INPUT FOR R OR ARE WE USING THE INITIALIZED R?

        # Call the previous acceleration functions for the halo, bulge and disk
        halo_accel = self.HernquistAccel(self.Mhalo, self.rhalo, self.r)
        bulge_accel = self.HernquistAccel(self.Mbulge, self.rbulge, self.r)
        disk_accel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, self.r)  
        # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        total_accel = halo_accel + bulge_accel + disk_accel
        return total_accel
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """ **** ADD COMMENTS """
        
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
        """ **** ADD COMMENTS """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r), *tuple(self.v)
        # the above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while t < tmax: 
            t += dt
            # store the new time in the first column of the ith row
            orbit[i, 0] = t
            
            # advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # ???????? rnew, vnew = self.LeapFrog(dt, orbit[i-1][1:4], orbit[i-1][4:7]) ????????
            self.r, self.v = self.LeapFrog(dt, self.r, self.v)
         
    
            # store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i, 1:4] = self.r
            orbit[i, 4:7] = self.v
            # update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


