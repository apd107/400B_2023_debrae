# Author: Aidan DeBrae
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from CenterOfMassMod import CenterOfMass
from ReadFile import Read

def OrbitCOM(galaxy, start, end, n):
    '''
    This function is designed to parse through the particle files 
    of galaxies and create a file which contains the information of 
    their separation and relative velocities. 

    PARAMETERS
    ----------
    galaxy: string
        String of the galaxy name to be analyzed
    start: integer
        The integer representing the snapshot id 
        to start on while parsing through galaxy files
    end: integer
        The integer representing the snapshot id 
        to end on while parsing through galaxy files
    n: integer
        The integer value which defines the interval of
        which the snap ids are parsed through
    
    RETURNS
    -------
    Does not return any values but rather writes a file 
    which contains the separation and relative velocities 
    between two galaxies as a function of time

    '''
    fileout = 'Orbit_' + galaxy + '.txt'
    delta = 0.1
    # Define volDec to be 4 for M33 because it
    # will be severly tidally stripped toward 
    # the end of the simutlation
    if galaxy == 'M33':
        volDec = 4
    else:
        volDec = 2
    snap_ids = np.arange(start, end, n)
    # Check if the snap id array contains info,
    # if it does not, the code should stop
    if snap_ids.size == 0:
        return
    else:
        orbit = np.zeros([len(snap_ids), 7])
        for i, snap_id in enumerate(snap_ids):
            print('\'for loop\' counter: ' + str(i))
            # Define the galaxy file name using the short 
            # algorithm using the snap_id
            ilbl = '000' + str(snap_id)
            ilbl = ilbl[-3:]
            galaxyfile = str(galaxy) + '_' + ilbl + '.txt'
            galaxypath = str(galaxy) + '_VLowRes/' + galaxyfile
            # Calculate the COM of mass for the position and 
            # velocity for the given galaxy
            galaxyCOM = CenterOfMass(galaxypath, 2)
            COM_pos = galaxyCOM.COM_P(delta, volDec)
            x, y, z = COM_pos[0], COM_pos[1], COM_pos[2]
            COM_vel = galaxyCOM.COM_V(x, y, z)
            vx, vy, vz = COM_vel[0], COM_vel[1], COM_vel[2]
            # Define the time of the given snapshot
            time = galaxyCOM.time
            time = (time.value)/1000
            # Iterate through the galaxy files according to the 
            # snapshot and save each component in the array
            orbit[i][0] = time
            orbit[i][1] = x.value
            orbit[i][2] = y.value
            orbit[i][3] = z.value
            orbit[i][4] = vx.value
            orbit[i][5] = vy.value
            orbit[i][6] = vz.value
        # Define the header format for the txt file as it is very long
        # Save the text file of the orbit information
        hf = '{:10s}{:11s}{:11s}{:11s}{:11s}{:11s}{:11s}'.format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz')
        np.savetxt(fileout, orbit, fmt='%11.3f'*7, comments='#', header=hf) 


