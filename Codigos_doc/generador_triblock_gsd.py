"""
Este código genera la configuración inicial de un sistema con 2 solventes y un grupo
de polímeros triblock, y la guarda directamente en un archivo .gsd

FABIO NORIEGA HERNANDEZ
"""
import numpy as np
import gsd.hoomd

# PARAMETERS
total_chains = 12_000
a_monomers = 21
b_monomers = 0
c_monomers = 0
total_monomers = a_monomers + b_monomers + c_monomers  # number of monomers per polymer
n_polymers = total_chains * total_monomers              # total number of polymer beads

density = 3.0
bond_length = 0.75
a_mass = b_mass = c_mass = d_mass = e_mass = 1.0

# In here we select if the system will create a subset of particles inside a droplet to help
# a fast thermalization of the system
droplet = 0

if droplet == 0:
    # Create in function of polymer concentration
    polymer_concentration = round(1 / 3, 5)

    total_beads = int(n_polymers * (1.0 / polymer_concentration))
    ceramic_beads = 1
    solvent_beads = total_beads - n_polymers - ceramic_beads

# Configuration of simulation box
aspect_ratio = 4.0
aspect_ratio_squared = pow(aspect_ratio, 2)
volume = total_beads / density
lx3 = volume * aspect_ratio_squared

lx = pow(lx3, (1 / 3))
ly = lz = lx / aspect_ratio

print(f"Total number of polymer beads: {n_polymers}")
print(f"Total ceramic beads: {ceramic_beads}")
print(f"Total solvent beads: {solvent_beads}")
print(f"Total beads: {total_beads}")
print(f"Volume: {volume:.4f}")
print(f"Box dimensions: lx={lx:.4f}, ly={ly:.4f}, lz={lz:.4f}")


def random_three_vector():
    phi = np.random.uniform(0, np.pi * 2)
    cos_theta = np.random.uniform(-1, 1)

    theta = np.arccos(cos_theta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    x *= bond_length
    y *= bond_length
    z *= bond_length

    return (x, y, z)


def positions():
    abc_monomers = a_monomers + b_monomers + c_monomers
    coordinates = np.zeros((abc_monomers, 3))

    if droplet == 0:
        for h in range(abc_monomers):
            x1, y1, z1 = random_three_vector()
            if h == 0:
                x = np.random.uniform((1 / 3), 1) * (lx / 2)
                x = round(x, 4)

                y = np.random.uniform(-1, 1) * (ly / 2)
                y = round(y, 4)

                z = np.random.uniform(-1, 1) * (lz / 2)
                z = round(z, 4)

                coordinates[h, 0] = x
                coordinates[h, 1] = y
                coordinates[h, 2] = z
            else:
                x = x + x1
                x = round(x - lx * int(round(x / lx)), 4)

                y = y + y1
                y = round(y - ly * int(round(y / ly)), 4)

                z = z + z1
                z = round(z - lz * int(round(z / lz)), 4)

                coordinates[h, 0] = x
                coordinates[h, 1] = y
                coordinates[h, 2] = z

        return coordinates


# This function assigns coordinates (x, y, z) to the solvent particles
def solvent_particles_positions_1():
    solvent_1 = np.zeros((solvent_beads, 3))

    # Generate positions over all volume of the box
    if droplet == 0:
        for p in range(solvent_beads):
            x = np.random.uniform(-1, 1 / 3) * (lx / 2)
            x = round(x, 4)

            y = np.random.uniform(-1, 1) * (ly / 2)
            y = round(y, 4)

            z = np.random.uniform(-1, 1) * (lz / 2)
            z = round(z, 4)

            solvent_1[p, 0] = x
            solvent_1[p, 1] = y
            solvent_1[p, 2] = z

        return solvent_1


# This function assigns coordinates (x, y, z) to the ceramic particles
def solvent_particles_positions_2():
    solvent_2 = np.zeros((ceramic_beads, 3))

    if droplet == 0:
        for p in range(ceramic_beads):
            x = np.random.uniform(-1, 1) * (lx / 2)
            x = round(x, 4)

            y = np.random.uniform(-1, 1) * (ly / 2)
            y = round(y, 4)

            z = np.random.uniform(-1, 1) * (lz / 2)
            z = round(z, 4)

            solvent_2[p, 0] = x
            solvent_2[p, 1] = y
            solvent_2[p, 2] = z

        return solvent_2


# Order of types stored in the GSD file (index used for typeid). 1=A, 2=B, 3=C, 4=D (solvent), 5=E (ceramic)
particle_types = ['A', 'B', 'C', 'D', 'E']
type_index = {t: i for i, t in enumerate(particle_types)}


# This function assigns the molecule typeid and mass to each position
def assign_types(total_particles):
    print(f"Arreglo de tamaño: {total_particles}")

    typeid = np.zeros(total_particles, dtype=np.uint32)
    mass_assignments = np.ones(total_particles, dtype=np.float32)
    counter = 0

    for i in range(total_chains):
        xd = total_monomers * counter

        for a in range(xd, a_monomers + xd):
            typeid[a] = type_index['A']
            mass_assignments[a] = a_mass

        for b in range(a_monomers + xd, a_monomers + b_monomers + xd):
            typeid[b] = type_index['B']
            mass_assignments[b] = b_mass

        for c in range(a_monomers + b_monomers + xd, a_monomers + b_monomers + c_monomers + xd):
            typeid[c] = type_index['C']
            mass_assignments[c] = c_mass

        counter += 1

    for k in range(n_polymers, n_polymers + solvent_beads):
        typeid[k] = type_index['D']
        mass_assignments[k] = d_mass

    for j in range(n_polymers + solvent_beads, total_particles):
        typeid[j] = type_index['E']
        mass_assignments[j] = e_mass

    return typeid, mass_assignments


# This function assigns bonds between monomers (returns an array of particle-index pairs)
def bonds_assigning():
    bond_groups = []
    counter = 0

    for i in range(total_chains):
        for j in range(total_monomers - 1):
            bond_groups.append([counter, counter + 1])
            counter += 1
        counter += 1

    return np.array(bond_groups, dtype=np.uint32)


# This function assigns angles between consecutive monomer triplets
def angle_assigning(typeid):
    angle_groups = []
    angle_type_strs = []
    counter = 0

    for i in range(total_chains):
        for j in range(total_monomers - 2):
            t = "{}-{}-{}".format(
                particle_types[typeid[counter]],
                particle_types[typeid[counter + 1]],
                particle_types[typeid[counter + 2]],
            )
            angle_type_strs.append(t)
            angle_groups.append([counter, counter + 1, counter + 2])
            counter += 1
        counter += 2

    return angle_groups, angle_type_strs


# --- Build positions for all chains + solvent + ceramic beads ---
chain_positions = positions()
for i in range(total_chains - 1):
    chain_positions = np.concatenate((chain_positions, positions()))

solvent_positions = solvent_particles_positions_1()
ceramic_positions = solvent_particles_positions_2()

all_positions = np.concatenate((chain_positions, solvent_positions))
all_positions = np.concatenate((all_positions, ceramic_positions))

total_particles = all_positions.shape[0]

typeid, particle_masses = assign_types(total_particles)

bond_groups = bonds_assigning()
angle_groups, angle_type_strs = angle_assigning(typeid)

# angle types are dynamic (depend on a_monomers/b_monomers/c_monomers), so build the
# type list from what actually occurs
angle_types = sorted(set(angle_type_strs)) if angle_type_strs else []
angle_type_index = {t: i for i, t in enumerate(angle_types)}
angle_typeid = np.array([angle_type_index[t] for t in angle_type_strs], dtype=np.uint32)

# --- Build the GSD frame ---
# NOTE: gsd >= 3.0 uses gsd.hoomd.Frame; older versions use gsd.hoomd.Snapshot.
# If Frame doesn't exist in your installed gsd version, use: frame = gsd.hoomd.Snapshot()
frame = gsd.hoomd.Frame()

frame.particles.N = total_particles
frame.particles.types = particle_types
frame.particles.typeid = typeid
frame.particles.position = all_positions.astype(np.float32)
frame.particles.mass = particle_masses
frame.particles.diameter = np.ones(total_particles, dtype=np.float32)
frame.particles.body = -1 * np.ones(total_particles, dtype=np.int32)
frame.particles.charge = np.zeros(total_particles, dtype=np.float32)

frame.bonds.N = bond_groups.shape[0]
frame.bonds.types = ['polymer']
frame.bonds.typeid = np.zeros(bond_groups.shape[0], dtype=np.uint32)
frame.bonds.group = bond_groups

frame.angles.N = len(angle_groups)
frame.angles.types = angle_types
frame.angles.typeid = angle_typeid
frame.angles.group = np.array(angle_groups, dtype=np.uint32) if angle_groups else np.zeros((0, 3), dtype=np.uint32)

frame.configuration.box = [lx, ly, lz, 0, 0, 0]
frame.configuration.step = 0

with gsd.hoomd.open(name='config.gsd', mode='w') as f:
    f.append(frame)

print("Archivo config.gsd generado correctamente.")


