import hoomd

temp = 0.65
eps_SP = 1.00
archivo_gsd = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}_monom_24"
file_id = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}"
aspect_ratio = 2
mon_cadena = 24
muestreo = 1_000_000

# Inicializar HOOMD con el snapshot
gpu = hoomd.device.GPU()
sim = hoomd.Simulation(device=gpu, seed=42)
sim.create_state_from_gsd(archivo_gsd, frame=-1)
    
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp)
    
original_box = sim.state.box


# Revisamos si la caja ya viene estirada o no 
if original_box.Lx > 100.0:
    print(f"📦 El GSD ya cuenta con la caja estirada ([{original_box.Lx:.2f}, {original_box.Ly:.2f}, {original_box.Lz:.2f}]). Reanudando directo...")
else:
    print("📦 El GSD tiene la caja chica del equilibrio. Aplicando estiramiento controlado...")
    new_Lx = original_box.Lx * aspect_ratio
    new_Ly = original_box.Ly + 0.1
    new_Lz = original_box.Lz + 0.1

    new_box = hoomd.Box(
        Lx=new_Lx,
        Ly=new_Ly,
        Lz=new_Lz,
        xy=original_box.xy,
        xz=original_box.xz,
        yz=original_box.yz
    )

    sim.state.set_box(box=new_box)
    print(f"Se estriró la caja a: \n{[new_box.Lx, new_box.Ly, new_box.Lz]}")

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    mie = hoomd.md.pair.Mie(nlist=cell, default_r_cut=4.0, mode='shift')

    mie.params[('S', 'S')] = dict(epsilon=1.0, sigma=1.0, n=12, m=6)
    mie.params[('P', 'P')] = dict(epsilon=1.0, sigma=1.0, n=12, m=6)
    mie.params[('S', 'P')] = dict(epsilon=eps_SP, sigma=1.0, n=12, m=6)

    # Fuerza de enlace para mantener la integridad de los polímeros
    armonico = hoomd.md.bond.Harmonic()
    armonico.params['P-P'] = dict(k=10.0, r0=1.0)
    

    # -- Integrador y termostato del ensamble NVT --
    # termostato = hoomd.md.methods.thermostats.Bussi(kT=temp, tau=0.01)
    termostato = hoomd.md.methods.thermostats.MTTK(kT=temp, tau=0.2)

    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)


    integrator = hoomd.md.Integrator(dt=0.005, methods=[nvt], forces=[mie, armonico])
    sim.operations.integrator = integrator

    # -- Loggers (Muestra de la información) --
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # Logger de datos termodinámicos (parecido a todo.dat)
    logger = hoomd.logging.Logger(categories=['scalar'])
    logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure'])    

    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(10000),
                              logger=logger,
                              output=open(f"log_{file_id}_monom_{mon_cadena}.csv", mode='a'))       
    
    sim.operations.writers.append(table)

    # Para poder guardar la primer configuración en el gsd y ver cómo se acomodaron las partículas
    trigger_combinado = hoomd.trigger.Or([hoomd.trigger.On(0),
                                          hoomd.trigger.On(1),
                                          hoomd.trigger.Periodic(5000)])       
    
    gsd_writer = hoomd.write.GSD(trigger=trigger_combinado,
                                 filename=f"{file_id}_monom_{mon_cadena}.gsd",
                                 mode='ab') 
    
    sim.operations.writers.append(gsd_writer)

    # sim.run(equilibracion)

    sim.run(muestreo)