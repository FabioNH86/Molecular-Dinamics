from hoomd import *
from hoomd import md
from hoomd import deprecated
import math
import numpy

#context.initialize("--mode=cpu --notice-level=0")
context.initialize("--mode=gpu --notice-level=0")

# launch with python3 dpdc.py 1 25 25 1.0 5 5 3
import sys
print('Argument List: %s' % str(sys.argv))
noption = int(sys.argv[1])
aii = float(sys.argv[2])
aij = float(sys.argv[3])
temp = float(sys.argv[4])
monoA = int(sys.argv[5])
monoB = int(sys.argv[6])
n_chains = int(sys.argv[7])

monomer_bond = 0.75;

n_changes = 120;
n_steps = 5000;
nsteps_tot = n_changes*n_steps;
random_steps = 20000;
phase = 0;
period = n_steps/20;
duty_high = 0.5;
duty_low = 0.5;

ramp_up_duration = n_steps*(1.0-duty_high);
ramp_down_duration = n_steps*(1.0-duty_low);

# 0=new configuration, 1=run from previous configuration 
if noption == 1:
   system = deprecated.init.read_xml(filename='step2.xml', restart='restart.xml')
   print(system.box)
   print('===================================')
   print('Running from previous configuration')
   print('===================================')
else:
   rho = 3.0;
   N_poly = n_chains;
   D_poly = monoA + monoB;
   # Calculate the characteristic size
   C_p = 0.1;
   N_t = (1.0*N_poly*D_poly)/C_p;
   Lx_cubic = N_t/rho;
   Lx = Lx_cubic**(1.0/3.0);
   box_side = Lx;
   print(box_side)
   Vol = box_side*box_side*box_side;
   N_solv = N_t - N_poly*D_poly;

   polymer = dict(bond_len=1.0, type=['A']*monoA+['B']*monoB, bond='linear', count=N_poly)
   solvent = dict(bond_len=1.0, type=['S'], bond="linear", count=N_solv)

   molecules = [polymer,solvent]

   system = deprecated.init.create_random_polymers(box=data.boxdim(Lx=box_side, Ly=box_side, Lz=box_side),polymers=molecules, separation=dict(A=0.1, B=0.1, S=0.1));

   # Bonded interactions for freely jointed chain
   harmonic = md.bond.harmonic()
   harmonic.bond_coeff.set('polymer', k=64.0, r0=monomer_bond)

   nl = md.nlist.cell(check_period=1)
   nl.reset_exclusions(exclusions=None)
                           
   dpd = md.pair.dpd(r_cut=1.0,nlist=nl,kT=temp,seed=2467)#2467
   dpd.pair_coeff.set('A', 'A', A=aii, gamma=4.5)
   dpd.pair_coeff.set('S', 'S', A=25, gamma=4.5)
   dpd.pair_coeff.set('B', 'B', A=aii, gamma=4.5)
   dpd.pair_coeff.set('A', 'B', A=aij, gamma=4.5)
   dpd.pair_coeff.set('A', 'S', A=aij, gamma=4.5)
   dpd.pair_coeff.set('B', 'S', A=aij, gamma=4.5)

   md.integrate.mode_standard(dt=0.01)
   nve = md.integrate.nve(group=group.all())
   movie1 = dump.gsd(group=group.all(),filename="step1.gsd", period=period, phase=phase)
   run(random_steps)
   nve.disable()
   movie1.disable()
   dpd.disable()
   harmonic.disable()

if noption == 1:
   filename = "step3"
else:
   filename = "step2"

logger = analyze.log(filename=filename+".log", period=period, phase=phase,
                     quantities=['time','potential_energy','temperature','pressure'], header_prefix='#', overwrite=False)

#### This saves the configuration at specific times in case that a restart is needed  ####
xml = deprecated.dump.xml(group=group.all(),filename=filename+".xml", vis=True, velocity = True, restart=True, period=1000000, phase=phase)
#gsd_restart = dump.gsd(filename="restart.gsd", group=group.all(), truncate=True, period=1000000, phase=phase, dynamic=['momentum'])

#### This saves the trajectory in dcd format ###
dump.dcd(filename=filename+".dcd", period=period, phase=phase)
dump.gsd(group=group.all(),filename=filename+".gsd", period=period, phase=phase)

if noption==1:
   # Angular potential for rigid chains
   ang_value = 3.1416
   angular = md.angle.harmonic()
   angular.angle_coeff.set('angleA', k=100.0, t0=ang_value)
   angular.angle_coeff.set('angleB', k=1.0, t0=ang_value)
   angular.angle_coeff.set('angleAAB', k=10.0, t0=ang_value)
   angular.angle_coeff.set('angleABB', k=10.0, t0=ang_value)

# Bond potential for chains

# Triangular periodic perturbations
t_high = 1.1
t_low = 0.9
m = 0.0
contador = 0
ref_time = 0
t0 = 0.0
T_var = 0.0
L_var = 0.0
duty = 0;
t_target = 0;
max_step = 0;

nl = md.nlist.cell(check_period=1)
nl.reset_exclusions(exclusions=None)
#nl = md.nlist.cell()

harmonic = md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=128.0, r0=monomer_bond)

dpd = md.pair.dpd(r_cut=1.0,nlist=nl,kT=temp, seed=2467)#2467
dpd.pair_coeff.set('A', 'A', A=aii, gamma=4.5)
dpd.pair_coeff.set('S', 'S', A=25, gamma=4.5)
dpd.pair_coeff.set('B', 'B', A=aii, gamma=4.5)
dpd.pair_coeff.set('A', 'B', A=aij, gamma=4.5)
dpd.pair_coeff.set('A', 'S', A=aij, gamma=4.5)
dpd.pair_coeff.set('B', 'S', A=aij, gamma=4.5)  

def prueba(timestep):
   global m
   global contador
   global ref_time
   global t0
   global T_var
   global L_var
   global t_target
   global duty
   global max_step

   if timestep == (random_steps) or (timestep-random_steps)%(n_steps)==0:
      ref_time = timestep
      print('change time')
      if contador%2 == 0:
         #negative slope (small to big lambda)
         m = -1.0
         t0 = t_high
         t_target = t_low;
         duty = 1.0-duty_low;
         max_step = ramp_down_duration;
      else:
         #positive slope (big to small lambda)
         m = 1.0
         t0 = t_low
         t_target = t_high;
         duty = 1.0-duty_high;
         max_step = ramp_up_duration;

      contador = contador+1
     
   #T_var = t0 + m*(timestep-ref_time)/n_steps * (t_high-t_low)
   if timestep<ref_time+max_step:
      T_var = t0 + m*(1.0/duty)*(timestep-ref_time)/n_steps * (t_high-t_low);
   else:
      T_var = t_target;
 
   # Sigmoidal function for size change
   Tc = 1.0
   w = -50.0
   lambda0 = 11.0
   L_var = lambda0*(1.0+T_var) + lambda0/(1.0+numpy.exp(w*(T_var-Tc)))
   #print(timestep, T_var, L_var)

   dpd.set_params(kT = T_var)
   dpd.pair_coeff.set('A', 'B', A=L_var, gamma=4.5)
   dpd.pair_coeff.set('A', 'S', A=L_var, gamma=4.5)
   dpd.pair_coeff.set('B', 'S', A=L_var, gamma=4.5)  
   dpd.update_coeffs()
                          
 
md.integrate.mode_standard(dt=0.01)
md.integrate.nve(group=group.all())
#run(nsteps_tot)
run(nsteps_tot, callback_period=1, callback=prueba)

# write the final configuration
xml.write_restart()
#xml.write(filename="start.xml", time_step=0)
#gsd_restart.write_restart()
