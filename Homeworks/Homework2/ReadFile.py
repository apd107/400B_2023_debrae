import numpy as np
import astropy.units as u


def Read(filename):
    file = open(filename, 'r')
    line_1 = file.readline()
    label, value = line_1.split()
    time = float(value)*u.Myr
    line_2 = file.readline()
    label, value = line_2.split()
    total_particles = float(value)
    file.close()
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    return time, total_particles, data







