import hoomd
import numpy as np
import gsd.hoomd


temp = 0.65
eps_SP = 1.00
archivo_gsd = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}_monom_24.gsd"
file_id = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}"
aspect_ratio = 2
mon_cadena = 24
muestreo = 1_000_000



rho = 0.3
natoms = 240_000
num_monomeros = 2

# Creating the set of atom coordinates and types
frame = gsd.hoomd.Frame()
frame.particles.N = natoms
frame.particles.types = ['S', 'P']

frame.bonds.types = ['P-P']

# Computing the number of atoms per species
N_a = num_monomeros
N_b = natoms - num_monomeros

types = []

        

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
bonds = []
enlace = 0.2
part_id = 0

for i in range(N_x):
   for j in range(N_y):
      for k in range(N_z):
         x = i*a_param - (0.5*L_x)
         y = j*a_param - (0.5*L_y)
         z = k*a_param - (0.5*L_z)
         types.append(0)


         monomer_id = part_id
         part_id += 1
         pos.append([x, y, z])

         if len(bonds) < 1:
            if np.random.rand() < 0.10:
                    x_comomer = x + enlace
                    y_comomer = y
                    z_comomer = z

                    pos.append([x_comomer, y_comomer, z_comomer])
                    comomer_id = part_id
                    types.append(1)
                    part_id += 1

                    bonds.append([monomer_id, comomer_id])

frame.bonds.N = len(bonds)
frame.particles.typeid = types
frame.bonds.group[:] = np.array(bonds)
frame.particles.position = pos
frame.configuration.box = [L_x, L_y, L_z, 0, 0, 0]
f = gsd.hoomd.open(name='file_bin.gsd', mode='w')
f.append(frame)


