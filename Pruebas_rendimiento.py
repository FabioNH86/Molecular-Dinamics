import hoomd
import hoomd.md
import numpy as np

# =============================================================================
# 1. CONFIGURACIÓN DEL SISTEMA MASIVO
# =============================================================================
# 1 millón de partículas para poner a prueba los 8GB de la RTX 4060
n_particles = 5_800_000  
densidad = 0.2
box_length = (n_particles / densidad) ** (1/3)  # Caja lo suficientemente grande para una densidad estándar

print(f"--- Creando snapshot para {n_particles:,} partículas ---")
snapshot = hoomd.Snapshot()
snapshot.configuration.box = [box_length, box_length, box_length, 0, 0, 0]
snapshot.particles.N = n_particles
snapshot.particles.types = ['A']

# =============================================================================
# GENERACIÓN DE UNA RED CÚBICA SIMPLE (SC)
# =============================================================================
print(f"--- Generando red ordenada para {n_particles:,} partículas ---")

# 1. Calcular cuántas partículas van por lado del cubo (raíz cúbica)
# Usamos int y redondeo para asegurar que quepa el número objetivo
particles_per_side = int(np.round(n_particles ** (1/3)))
actual_particles = particles_per_side ** 3

if actual_particles != n_particles:
    print(f"Nota: Ajustando n_particles de {n_particles:,} a {actual_particles:,} "
          f"para formar un cubo perfecto de {particles_per_side}x{particles_per_side}x{particles_per_side}.")
    n_particles = actual_particles
    snapshot.particles.N = n_particles  # Actualizar el snapshot

# 2. Generar los puntos lineales espaciados uniformemente dentro de la caja
# Dejamos un pequeño margen en los bordes para evitar solapamientos con las fronteras periódicas
grid_points = np.linspace(-box_length/2 + (box_length / (2 * particles_per_side)), 
                          box_length/2 - (box_length / (2 * particles_per_side)), 
                          particles_per_side)

# 3. Crear la malla tridimensional (coordenadas X, Y, Z)
x, y, z = np.meshgrid(grid_points, grid_points, grid_points, indexing='ij')

# 4. Aplanar las matrices para obtener el arreglo de posiciones (N, 3)
positions = np.vstack((x.flatten(), y.flatten(), z.flatten())).T

# 5. Asignar al snapshot de HOOMD
snapshot.particles.position[:] = positions
print(f"--- Red cúbica inicializada correctamente ---")

# Inicializar la simulación en la GPU
sim = hoomd.Simulation(device=hoomd.device.GPU(), seed=1)
sim.create_state_from_snapshot(snapshot)

print("--- Estado del sistema inicializado en la GPU ---")

# =============================================================================
# 2. FUERZAS E INTERACCIONES (El devorador de VRAM)
# =============================================================================
# El 'Neighbor List' en GPU reserva grandes arreglos de memoria para los pares
nlist = hoomd.md.nlist.Cell(buffer=0.4)
lj = hoomd.md.pair.LJ(nlist=nlist)

# Parámetros estándar de Lennard-Jones
lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
lj.r_cut[('A', 'A')] = 4.0  # Radio de corte estándar

# =============================================================================
# 3. INTEGRADOR DE DINÁMICA MOLECULAR
# =============================================================================
# Usamos ConstantVolume (NVE) estándar para la prueba de carga
integrator = hoomd.md.Integrator(dt=0.001)
integrator.forces.append(lj)

thermo_method = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
integrator.methods.append(thermo_method)
sim.operations.integrator = integrator

# =============================================================================
# 4. EJECUCIÓN (Monitorea nvidia-smi aquí)
# =============================================================================
print("\n>>> Iniciando pasos de simulación. Revisa la VRAM ahora de inmediato. <<<")
try:
    # Corremos pasos en bloques para observar el incremento en la terminal
    for i in range(5):
        print(f"Ejecutando bloque {i+1}/5 (100 pasos)...")
        sim.run(100)
    print(f"\n--- Prueba terminada con éxito. ¡Tu VRAM resistió {n_particles} partículas! ---")
except RuntimeError as e:
    if "Out of memory" in str(e) or "an illegal memory access" in str(e).lower():
        print("\n[¡ÉXITO!] Hemos alcanzado el límite de la VRAM. El sistema se quedó sin memoria (OOM).")
    else:
        print(f"\nOcurrió un error durante la ejecución: {e}")
