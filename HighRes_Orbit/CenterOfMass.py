
# Aidan DeBrae (I used the given class template)
# Homework 4
# Center of Mass Position and Velocity

# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read




class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        
        self.m = self.data['m'][self.index]
        self.x = (self.data['x'][self.index])
        self.y = (self.data['y'][self.index])
        self.z = (self.data['z'][self.index])
        self.vx = (self.data['vx'][self.index])
        self.vy = (self.data['vy'][self.index])
        self.vz = (self.data['vz'][self.index])


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''

        # xcomponent Center of mass
        a_com = sum(a*m)/sum(m)
        # ycomponent Center of mass
        b_com = sum(b*m)/sum(m)
        # zcomponent Center of mass
        c_com = sum(c*m)/sum(m)
        
        # return the 3 components separately
        return a_com, b_com, c_com
       
    
    
    def COM_P(self, delta, volDec):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        volDec: float
            determines the amount by which Rmax should be decreased
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        r_COM = (x_COM**2 + y_COM**2 + z_COM**2)**0.5


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = (x_new**2 + y_new**2 + z_new**2)**0.5

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/volDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        while (change > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            array_radii = (x_new**2 + y_new**2 + z_new**2)**0.5
            index2 = np.where(array_radii < r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)

            # compute the new 3D COM position
            r_COM2 = (x_COM2**2 + y_COM2**2 + z_COM2**2)**0.5

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)                                                                                  

            # Before loop continues, reset : r_max, particle separations and COM                                        
            # reduce the volume by a factor of 2 again                                                                 
            r_max /= volDec
                                                                                    

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            r_new = (x_new**2 + y_new**2 + z_new**2)**0.5

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        p_COM = np.around(p_COM*u.kpc, 2)
        return p_COM
        
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of 
            massposition.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        xV = self.x - x_COM.value
        yV = self.y - y_COM.value
        zV = self.z - z_COM.value
        rV = (xV**2 + yV**2 + zV**2)**0.5
        rV = rV*u.kpc
        
        # determine the index for those particles within the max radius
        indexV = np.where(rV < rv_max)
        
        # determine the velocity and mass of those particles within the mas radius
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new =  self.m[indexV]
        
        # compute the center of mass velocity using those particles
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # create an array to store the COM velocity
        v_COM = np.array([vx_COM, vy_COM, vz_COM])

        # return the COM vector
        v_COM = np.around(v_COM*(u.km/u.s), 2)
        return v_COM                                                                                        
     


