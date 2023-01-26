import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename, type, number):
    time, total_particles, data = Read(filename)
    particle_index = np.where(data['type'] == type)
    particles = data[particle_index]

    x, y, z = float(particles['x'][number])*u.kpc, \
              float(particles['y'][number])*u.kpc, \
              float(particles['z'][number])*u.kpc

    vx, vy, vz = float(particles['vx'][number])*(u.km/u.second), \
                 float(particles['vy'][number])*(u.km/u.second), \
                 float(particles['vz'][number])*(u.km/u.second)
    
    distance = (x**2 + y**2 + z**2)**0.5
    velocity = (vx**2 + vy**2 + vz**2)**0.5
    mass = float(particles['m'][number])*u.Msun

    distance = np.around(distance, 3)
    velocity = np.around(velocity, 3)

    return distance, velocity, mass


