import hoomd
import numpy as np
import os
import sys

# Reading input arguments from terminal
# Execute with python3 LJbin.py 0.8
print('Argument list: %s' % str(sys.argv))
temperature = float(sys.argv[1])

# Setting up the simulation device parameters 
gpu = hoomd.device.GPU()
gpu.notice_level = 3
sim = hoomd.Simulation(device = gpu, seed=2467)

# Parameters for the initial simulation box
# This only works for crystal. If noption = 1 it creates 
# a cubic simulation box, if noption = 0 it reads a gsd file
print('==========================================================')
print('Reading a GSD file')
   
# Restart the simulation from previous state
sim.create_state_from_gsd(filename="config.gsd")

# Neighbor list param
nl = hoomd.md.nlist.Cell(buffer=0.4)

# Setting up the interaction potential
lj = hoomd.md.pair.LJ(nlist=nl, default_r_cut=4.0, mode="shift")
lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
lj.params[('B', 'B')] = dict(epsilon=1.0, sigma=1.0)
lj.params[('A', 'B')] = dict(epsilon=0.3, sigma=1.0)

# Setting up integration method and thermostat for NVT
thermo_mttk = hoomd.md.methods.thermostats.MTTK(kT=temperature, tau=0.2)
nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=thermo_mttk)
integrator = hoomd.md.Integrator(dt=0.002, methods=[nvt], forces=[lj])
remover = hoomd.md.update.ZeroMomentum(hoomd.trigger.Periodic(1_000))
sim.operations.integrator = integrator

# Creating a data logger for thermodynamic properties
thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermo)

output_file = f"simpleLJ.log"

logger = hoomd.logging.Logger(categories=['scalar'])
logger.add(sim, quantities=['timestep'])
logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure', 'pressure_tensor'])

# Output file for thermodynamic properties
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(1_000), logger=logger, output=open(output_file, 'w'))

sim.operations.writers.append(table)

# Send the logger outputs to the terminal
term_log = hoomd.logging.Logger(categories=['scalar', 'string'])
term_log.add(sim, quantities=['timestep', 'tps'])
term_writer = hoomd.write.Table(trigger=hoomd.trigger.Periodic(100), logger=term_log)
sim.operations.writers.append(term_writer)

# System state output files
restart_gsd = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(10_000), filename="restart.gsd", filter=hoomd.filter.All(), truncate=True)

movie_gsd = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(500), filename="movie.gsd", filter=hoomd.filter.All())

sim.operations.writers.append(restart_gsd)
sim.operations.writers.append(movie_gsd)

# Run simulation
sim.run(200_000)

