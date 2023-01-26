import astropy.units as u
from ParticleProperties import ParticleInfo

distance, velocity, mass = ParticleInfo('MW_000.txt', 2.0, 99)
converted_distance = distance.to(u.lyr)

print(distance)
print(converted_distance)
print(velocity)
print(mass)