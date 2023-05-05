# import modules
import numpy as np
import matplotlib as plt
import scipy.optimize as so
# my modules from assignments
from ReadFile import Read
from CenterOfMass import CenterOfMass

# snapshot 795 for pericenter
# snapshot 185 for apocenter


class Halo_Shape:
    def __init__(self, snap_file):
        # Create a COM of object for M33
        self.M33_COM = CenterOfMass('HighRes/M33/' + str(snap_file), 2.0)
        M33_COM_pos = self.M33_COM.COM_P(0.1, 4)
        M33_COM_vel = self.M33_COM.COM_V(M33_COM_pos[0],
                            M33_COM_pos[1],
                            M33_COM_pos[2])
        

        self.time, self.total, self.data = Read('HighRes/M33/' + str(snap_file))                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == 1.0)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        
        self.m = self.data['m'][self.index]
        self.x = (self.data['x'][self.index])
        self.y = (self.data['y'][self.index])
        self.z = (self.data['z'][self.index])
        self.vx = (self.data['vx'][self.index])
        self.vy = (self.data['vy'][self.index])
        self.vz = (self.data['vz'][self.index])

        # Set up the coordinates to calulate the correct positions
        xH = self.x - M33_COM_pos[0].value 
        yH = self.y - M33_COM_pos[1].value 
        zH = self.z - M33_COM_pos[2].value

        self.r_tot = np.sqrt(xH**2 + yH**2 + zH**2)

        vxH = self.vx - M33_COM_vel[0].value 
        vyH = self.vy - M33_COM_vel[1].value 
        vzH = self.vz - M33_COM_vel[2].value

        self.v_tot = np.sqrt(vxH**2 + vyH**2 + vzH**2)
        self.r = np.array([xH,yH,zH]).T 
        self.v = np.array([vxH,vyH,vzH]).T
    
    # Rotating the frame of the halo toward the z-axis is useful because I can hopefully
    # more easily calaculate the axes of the ellipsoid to understand the shape. 
    def RotateFrame(self, posI, velI):
        """a function that will rotate the position and velocity vectors
        so that the disk angular momentum is aligned with z axis. 
        
        PARAMETERS
        ----------
            posI : `array of floats`
                3D array of positions (x,y,z)
            velI : `array of floats`
                3D array of velocities (vx,vy,vz)
                
        RETURNS
        -------
            pos: `array of floats`
                rotated 3D array of positions (x,y,z) such that disk is in the XY plane
            vel: `array of floats`
                rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
                is in the +z direction 
        """
        
        # compute the angular momentum
        L = np.sum(np.cross(posI,velI), axis=0)
        # normalize the vector
        L_norm = L/np.sqrt(np.sum(L**2))


        # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
        
        # z unit vector
        z_norm = np.array([0, 0, 1])
        
        # cross product between L and z
        vv = np.cross(L_norm, z_norm)
        s = np.sqrt(np.sum(vv**2))
        
        # dot product between L and z 
        c = np.dot(L_norm, z_norm)
        
        # rotation matrix
        I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
        R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

        # Rotate coordinate system
        pos = np.dot(R, posI.T).T
        vel = np.dot(R, velI.T).T
        
        return pos, vel
    
    def ellipse(self, horizontal, vertical, h_scale=1.0, v_scale=1.0):
        """
        The purpose of this function is to create the horizontal
        and vetical radii of an ellipse according to the data points 
        of the halo particles. It will be used to plot and determine
        the shape of the halo.  
        
        PARAMETERS
        ----------
            horizontal : `np.array`
                The horizontal data points of the halo. 
            vertical : `np.array`
                The vertical data points of the halo.
            h_scale: 'float'
                The scale factor by which the horizontal coordinates 
                should be scaled by 
            v_scale: 'float'
                The scale factor by which the vertical coordinates 
                should be scaled by
            n: 'float'
                Size of the ellipse which corresponds to the number
                of desired standard deviations
                
        RETURNS
        -------
            h_coord: `np.array`
                The horizontal coordinates of the ellipse 
            v_coord: 'np.array'
                The vertical coordinates of the ellipse 
        """
        h_mean = sum(horizontal)/len(horizontal)
        v_mean = sum(vertical)/len(vertical)
        N = len(horizontal)
        h_std = np.std(horizontal)
        v_std = np.std(vertical)
        cov = sum((horizontal-h_mean)*(vertical-v_mean))/N

        p = cov/(h_std*v_std)
        h_radii = np.sqrt(1+p)
        v_radii = np.sqrt(1-p)

        h_radii = h_radii*(2*h_scale*h_std)
        v_radii = v_radii*(2*v_scale*v_std)

        angles = np.linspace(0, 2*np.pi, 200)
        h_coord = h_radii*np.cos(angles)
        v_coord = v_radii*np.sin(angles)
        return h_coord, v_coord, h_radii, v_radii
