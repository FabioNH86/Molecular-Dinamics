from funciones import crear_primer_frame, correr_simulacion 
import gsd.hoomd
import os

carpeta_guardado = 'Pruebas/'
if not os.path.exists(carpeta_guardado):
        os.makedirs(carpeta_guardado)

num_prueba = 0

# NOTA DE PRUEBA: Para que corra rápido y no lance errores de espacio con 
# n_monomeros=12000 (que genera casi 1.2 millones de solventes al 1%), 
# bajé temporalmente n_monomeros a 400. Sube este número si tienes potencia/tiempo.
snapshot = crear_primer_frame(
    densidad_goticula=0.3, 
    aspect_ratio=4.0, 
    concentracion_porcentual_monomeros=1, 
    monomeros_en_polimero=5,
    n_monomeros=100 
    )

# --- TEST 1: Guardar la configuración inicial a mano vía GSD (Tu script original) ---
frame = gsd.hoomd.Frame()
frame.particles.N = snapshot.particles.N
frame.particles.position = snapshot.particles.position
frame.particles.typeid = snapshot.particles.typeid
frame.particles.types = snapshot.particles.types
frame.configuration.box = snapshot.configuration.box

nombre_archivo_ini = f"configuracion_inicial_{num_prueba}.gsd"
ruta_completa_ini = os.path.join(carpeta_guardado, nombre_archivo_ini)

with gsd.hoomd.open(name=ruta_completa_ini, mode='w') as archivo_gsd:
        archivo_gsd.append(frame)
print(f"¡Archivo inicial {nombre_archivo_ini} creado con éxito!")

    # --- TEST 2: Pasar el snapshot a HOOMD y correr la simulación ---
    # Usamos pasos cortos (200) para validar que el integrador y los potenciales no rompan.
# correr_simulacion(
#         snapshot=snapshot,
#         temp=1.2,
#         equilibracion=20,
#         muestreo=20,
#         eps_SP=1.0
#     )