import gsd.hoomd
import os
import sys
import numpy as np

# Reading input arguments from terminal
# Execute with python3 gsd_gen_bin.py 0.5 8 0.5
# the second argument must be a number with exact cubic root
print('Argument list: %s' % str(sys.argv))
rho = float(sys.argv[1])
natoms = int(sys.argv[2])
X_a = float(sys.argv[3])

# Creating the set of atom coordinates and types
frame = gsd.hoomd.Frame()
frame.particles.N = natoms
frame.particles.types = ['A', 'B']

# Computing the number of atoms per species
N_a = int(round(natoms*X_a))
N_b = natoms - N_a

types = []
for i in range(natoms):
     if i<N_a:
        types.append(0)
     else:
        types.append(1)
        
frame.particles.typeid = types

#print(types)

# Creating particle coordinates
N_tot = natoms
N_x = int(np.cbrt(N_tot))
N_y = N_x
N_z = N_x
a_param = (1.0/rho)**(1.0/3.0)
L_x = N_x * a_param
L_y = N_y * a_param
L_z = N_z * a_param

pos = []
for i in range(N_x):
   for j in range(N_y):
      for k in range(N_z):
         x = i*a_param - (0.5*L_x)
         y = j*a_param - (0.5*L_y)
         z = k*a_param - (0.5*L_z)
         pos.append([x, y, z])

frame.particles.position = pos
frame.configuration.box = [L_x, L_y, L_z, 0, 0, 0]
f = gsd.hoomd.open(name='file_bin.gsd', mode='w')
f.append(frame)

