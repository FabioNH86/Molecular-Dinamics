import hoomd
import hoomd.md
import numpy as np

# =============================================================================
# 1. CONFIGURACIÓN DEL SISTEMA MASIVO
# =============================================================================
# 1 millón de partículas para poner a prueba los 8GB de la RTX 4060
n_particles = 4_000_000  
box_length = 150.0  # Caja lo suficientemente grande para una densidad estándar

print(f"--- Creando snapshot para {n_particles:,} partículas ---")
snapshot = hoomd.Snapshot()
snapshot.configuration.box = [box_length, box_length, box_length, 0, 0, 0]
snapshot.particles.N = n_particles
snapshot.particles.types = ['A']

# Posiciones aleatorias en la caja
np.random.seed(42)
positions = np.random.uniform(-box_length/2, box_length/2, size=(n_particles, 3))
snapshot.particles.position[:] = positions

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
print("\n>>> Iniciando pasos de simulación. Revisa tu VRAM ahora de inmediato. <<<")
try:
    # Corremos pasos en bloques para observar el incremento en la terminal
    for i in range(5):
        print(f"Ejecutando bloque {i+1}/5 (100 pasos)...")
        sim.run(100)
    print("\n--- Prueba terminada con éxito. ¡Tu VRAM resistió el millón de partículas! ---")
except RuntimeError as e:
    if "Out of memory" in str(e) or "an illegal memory access" in str(e).lower():
        print("\n[¡ÉXITO!] Hemos alcanzado el límite de la VRAM. El sistema se quedó sin memoria (OOM).")
    else:
        print(f"\nOcurrió un error durante la ejecución: {e}")