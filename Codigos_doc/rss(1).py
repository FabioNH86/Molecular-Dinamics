import hoomd
import numpy
from hoomd import md
from hoomd import deprecated
#hoomd.context.initialize("--mode=cpu")
hoomd.context.initialize("--mode=gpu --gpu=0")

# 0=new configuration, 1=run from previous configuration 
noption = 1;
nat = 32768;
rho = 0.1;
lado = numpy.sqrt(nat/rho); 

if noption == 1:
   system = hoomd.init.read_gsd(filename='init.gsd', restart='restart.gsd')
   print(system.box)
#   system.box = hoomd.data.boxdim(Lx=35.0, Ly=35.0, Lz=35.0)
   print(system.box)
   print('===========================')
   print('Running from previous configuration')
   print('===========================')
else:
   # create a square box with N particles of A species
   system = deprecated.init.create_random(N=nat, box=hoomd.data.boxdim(Lx=lado, Ly=lado, dimensions=2), name = 'A')
#system=hoomd.init.create_lattice(unitcell=hoomd.lattice.sq(a=2.0),n=[10,10]);
#system=hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=1.2),n=[10,5]);
#system = deprecated.init.create_random(N=1000, phi_p=0.261799387, name = 'A')
   print('===========================')   
   print('New run')
   print('===========================')

# specify RSS soft interactions between particle pairs

def rss(r, rmin, rmax, epsilon, sigma1, sigma2, k0, n):
    V = epsilon*(sigma1/r)**n + epsilon*(1.0/2.0)*(1.0-numpy.tanh(k0*(r-sigma2)));
    F = n*(sigma1/r)**(n-1.0)*(sigma1/r**2.0)+(k0/2.0)*(1.0-numpy.tanh(k0*(r-sigma2))**2.0);
    return (V, F)
nl = md.nlist.cell()
table = md.pair.table(width=1000, nlist=nl);
table.pair_coeff.set('A', 'A', func=rss, rmin=0.001, rmax=3.0, coeff=dict(epsilon=1.0, sigma1=1.0, sigma2=1.35, k0=10.0, n=14.0))

# integrate
all = hoomd.group.all();
md.integrate.mode_standard(dt=0.001)

#luigi=md.integrate.langevin(group=all, kT=1.0, seed=5)
#hoomd.run(1e3)
#luigi.disable()
#md.integrate.nve(group=all)

md.integrate.nvt(group=all, kT=0.18, tau=1.0)
#md.integrate.npt(group=all, kT=0.18, tau=0.005, tauP=0.5, P=5.0)

# write the output files
hoomd.analyze.log(filename='temperature.dat', quantities=['temperature'], period=1000, phase=1000)
hoomd.analyze.log(filename='pressure.dat', quantities=['pressure'], period=1000, phase=1000)
hoomd.analyze.log(filename='u_energy.dat', quantities=['potential_energy'], period=1000, phase=1000)
#hoomd.analyze.log(filename='volume.dat', quantities=['volume'], period=1000, phase=1000)
# for the surface tension 
#hoomd.analyze.log(filename='pressure_xx.dat', quantities=['pressure_xx'], period=100, phase=100, overwrite=False)
#hoomd.analyze.log(filename='pressure_yy.dat', quantities=['pressure_yy'], period=100, phase=100, overwrite=False)
#hoomd.analyze.log(filename='pressure_zz.dat', quantities=['pressure_zz'], period=100, phase=100, overwrite=False)
# write the configuration of the system 
gsd_restart = hoomd.dump.gsd(filename="restart.gsd", group=hoomd.group.all(), truncate=True, period=1000000, phase=0)
#xml_draw = hoomd.deprecated.dump.xml(group=hoomd.group.all(), filename="movie", period=1000000, phase=0)
pelicula=hoomd.dump.gsd(filename="movie.gsd", period=100000, group=hoomd.group.all(), phase=0)
deprecated.dump.pos(filename='xyz.pos', period=100000)
# time steps
#hoomd.run_upto(3e7)
hoomd.run(1e8)
gsd_restart.write_restart()

