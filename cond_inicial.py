import hoomd
import numpy as np
from scipy.spatial import cKDTree

def run_polymer_hoomd(temp, equilibracion, muestreo, n_monomeros_totales, monomeros_por_polimero, n_solvente, densidad_líquido=0.6,eps_SP=1.0):
    # --- Identificador para archivos ---
    file_id = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}"
    
    print(f'>>> Ejecutando: {file_id}...')
    # -- Uso de GPU -- 
    device = hoomd.device.GPU()
    sim = hoomd.Simulation(device=device, seed=42)


    padding = 5.0  # Espacio mínimo desde las paredes para evitar solapamientos
    distancia_minima = 0.9 # Distancia mínima entre partículas para evitar solapamientos

    n_polimeros = n_monomeros_totales // monomeros_por_polimero
    
    n_total = n_monomeros_totales + n_solvente

    n_enlaces = n_polimeros * (monomeros_por_polimero - 1) 


    # Para gotícula 
    parti_x, parti_y, parti_z = n_total ** (1/3), n_total ** (1/3), n_total ** (1/3)


    aspect_ratio = 4.0
    # Dimensiones caja cuadrada 
    # rho = N / V => V = N / rho => L = (N / rho)^(1/3)
    # l = (n_total / densidad_líquido) ** (1/3)

    parametro_red = (1/densidad_líquido) ** (1/3)

    lx, ly, lz = parti_x * aspect_ratio, parti_y * parametro_red, parti_z * parametro_red

    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        snap.configuration.box = [lx, ly, lz, 0, 0, 0]
        snap.particles.N = n_total
        snap.particles.types = ['S', 'P'] # S para solvente, P para polímero
        snap.particles.mass[:] = [1.0] * n_total

        # Polímero bond
        snap.bonds.N = n_enlaces
        snap.bonds.types = ['P-P']

        type_S_id = snap.particles.types.index('S')
        type_P_id = snap.particles.types.index('P')

        # Construcción de enlaces para los polímeros
        # --- Construcción Distribución de Polímeros (Garantiza cero solapamientos) ---
        bond_counter = 0
        
        # Creamos una red de puntos espaciados exclusivamente para colocar los polímeros
        # Buscamos una distribución tridimensional para las n_polimeros cadenas
        n_p_eje = int(np.ceil(n_polimeros ** (1/3)))
        px_coords = np.linspace(-lx/2 + padding, lx/2 - padding, n_p_eje)
        py_coords = np.linspace(-ly/2 + padding, ly/2 - (monomeros_por_polimero * 0.9) - padding, n_p_eje)
        pz_coords = np.linspace(-lz/2 + padding, lz/2 - padding, n_p_eje)
        
        PX, PY, PZ = np.meshgrid(px_coords, py_coords, pz_coords, indexing='ij')
        orígenes_polimeros = np.vstack([PX.ravel(), PY.ravel(), PZ.ravel()]).T

        # Acomodo de las partículas en cadenas lineales bien separadas
        for i in range(n_polimeros):
            start_idx = i * monomeros_por_polimero
            
            # En lugar de np.random.uniform, tomamos un origen fijo y seguro de la red de polímeros
            x_start, y_start, z_start = orígenes_polimeros[i]
            
            for j in range(monomeros_por_polimero):
                idx = start_idx + j
                snap.particles.position[idx] = [x_start, y_start + j * 0.9, z_start]
                snap.particles.typeid[idx] = type_P_id

                # Enlace con el monómero anterior (excepto el primero de cada cadena)
                if j > 0:
                    snap.bonds.group[bond_counter] = [idx - 1, idx]
                    snap.bonds.typeid[bond_counter] = 0
                    bond_counter += 1

        posiciones_polimeros = np.array([
            snap.particles.position[i] for i in range(n_monomeros_totales)])

        # Solvente
        # --- Solvente en Red Uniforme (Lattice) ---
        # 1. Calculamos la densidad de número del solvente para estimar el espaciado
        volumen_caja = lx * ly * lz
        distancia_nodos = (volumen_caja / n_solvente) ** (1/3)
        
        # 2. Determinamos cuántas partículas caben idealmente en cada eje según las proporciones de tu caja (1000, 500, 500)
        # Esto nos dará una relación aproximada de 2:1:1 en la cantidad de divisiones
        nx = int(np.round(ly / distancia_nodos))
        ny = int(np.round(ly / distancia_nodos))
        nz = int(np.round(lz / distancia_nodos))
        
        # Ajustamos ligeramente para asegurarnos de tener suficientes nodos en la red espacial
        while (nx * ny * nz) < n_solvente:
            nx += 1
            ny += 1
            nz += 1

        # 3. Generamos las rejillas lineales para cada eje (dejando un margen en las paredes)
        x_coords = np.linspace(-(1/8) * lx, (1/8) * lx, nx)
        y_coords = np.linspace(-ly/2 + padding, ly/2 - padding, ny)
        z_coords = np.linspace(-lz/2 + padding, lz/2 - padding, nz)
        
        # 4. Creamos la matriz tridimensional de puntos
        X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
        red_completa = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
        
        # 5. Tomamos exactamente la cantidad de partículas requeridas (n_solvente)
        # Dado que la densidad es baja, la distancia entre ellas será ~ 7.5 unidades (cero riesgo de overlap)
        pos_solvente = red_completa[:n_solvente]
        
        # Asignamos al snapshot
        arbol = cKDTree(posiciones_polimeros)

        pos_solvente_filtrada = []
        for punto in red_completa:
            dist, _ = arbol.query(punto, k=1)
            if dist >= distancia_minima:
                pos_solvente_filtrada.append(punto)
            if len(pos_solvente_filtrada) == n_solvente:
                break

        if len(pos_solvente_filtrada) < n_solvente:
            raise ValueError(f"No hay suficientes puntos libres: solo {len(pos_solvente_filtrada)}")

        snap.particles.position[n_monomeros_totales:] = np.array(pos_solvente_filtrada)

        snap.particles.typeid[n_monomeros_totales:] = type_S_id

    sim.create_state_from_snapshot(snap)
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp)

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    mie = hoomd.md.pair.Mie(nlist=cell, default_r_cut=4.0, mode='shift')

    mie.params[('S', 'S')] = dict(epsilon=1.0, sigma=1.0, n=12, m=6)
    mie.params[('P', 'P')] = dict(epsilon=1.0, sigma=1.0, n=12, m=6)
    mie.params[('S', 'P')] = dict(epsilon=eps_SP, sigma=1.0, n=12, m=6)

    # Fuerza de enlace para mantener la integridad de los polímeros
    armonico = hoomd.md.bond.Harmonic()
    armonico.params['P-P'] = dict(k=10.0, r0=1.0)
    

    # -- Integrador y termostato del ensamble NVT --
    termostato = hoomd.md.methods.thermostats.Bussi(kT=temp, tau=0.01)
    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)


    integrator = hoomd.md.Integrator(dt=0.001, methods=[nvt], forces=[mie, armonico])
    sim.operations.integrator = integrator

    # -- Loggers (Muestra de la información) --
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # Logger de datos termodinámicos (parecido a todo.dat)
    logger = hoomd.logging.Logger(categories=['scalar'])
    logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure'])    

    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(10000),
                              logger=logger,
                              output=open(f"log_{file_id}.csv", 'w'))       
    
    sim.operations.writers.append(table)

    # Para poder guardar la primer configuración en el gsd y ver cómo se acomodaron las partículas
    trigger_combinado = hoomd.trigger.Or([hoomd.trigger.On(1),
                                          hoomd.trigger.Periodic(5000)])       
    
    gsd_writer = hoomd.write.GSD(trigger=trigger_combinado,
                                 filename=f"traj_{file_id}.gsd",
                                 mode='wb') 
    
    sim.operations.writers.append(gsd_writer)
    print(f'Densidad: {n_total / (lx * ly * lz):.4f} | Lx = {lx:.4f} | Ly = {ly:.4f} | Lz = {lz:.4f} | N = {n_total}')

    # print("--- Relajando el sistema con FIRE para eliminar solapamientos ---")
    # # 1. Creamos un método de velocidad nula para la minimización (un "cojín" para que no salgan disparadas)
    # displacement_capped = hoomd.md.methods.DisplacementCapped(
    # filter=hoomd.filter.All(),
    # maximum_displacement=1e-2)
    
    # # 2. Configuramos el optimizador FIRE
    # fire = hoomd.md.minimize.FIRE(dt=0.0001, 
    #                               force_tol=1e-1, 
    #                               angmom_tol=1e-1, 
    #                               energy_tol=1e-4, 
    #                               methods=[displacement_capped], 
    #                               forces=[mie, armonico])
    
    # # Asignamos FIRE temporalmente a la simulación
    # sim.operations.integrator = fire
    
    # # Corremos FIRE hasta que converja o alcance un máximo de 800 pasos
    # while not fire.converged and sim.timestep < 5000:
    #     sim.run(50)
        
    # print(f"Sistema minimizado en el paso: {sim.timestep}. ¿Convergió?: {fire.converged}")

    # 3. Quitamos FIRE y restablecemos tu integrador NVT original para la producción
    sim.operations.integrator = integrator
    
    # Resetear el contador de pasos 

    # --- AHORA SÍ, CORRES TUS ETAPAS ORIGINALES ---
    sim.run(equilibracion)
    sim.run(muestreo)

    print(f'Simulación finalizada ✅\n >> Condiciones: T={temp}, eps_SP={eps_SP}, Monómeros/Polímero={monomeros_por_polimero}, Concentrción de polímeros={n_polimeros / (lx * ly * lz):.4f}\n')