import os 
from funciones import run_hoomd_simulation

"""
Este código se encarga de automatizar una serie de simulaciones de dinámica molecular
utilizando la librería HOOMD-blue. El script itera sobre una lista de temperaturas
pre-asignadas, configurando el estado inicial, el potencial de Mie y el ensamble NVT 
de forma nativa en Python.

Los resultados (trayectorias en formato .gsd y datos termodinámicos en .csv) se 
almacenan automáticamente en directorios organizados por temperatura:
(Resultados/P{num_prueba}_HOOMD_Mie/T=*).

Fabio Noriega Hernández
Abril 2026
"""

# -- SEÑALA EL NÚMERO DE ENSAYO QUE HARÁS PARA ALMACENAR LOS RESULTADOS EN SU CARPETA CORRESPONDIENTE --
num_prueba = 1

#temperaturas = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temperaturas = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20] 
rho_objetivo = 0.1875
ruta_base = f"Resultados/P{num_prueba}_HOOMD_Mie"



for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    ruta_destino = os.path.join(ruta_base, nombre_carpeta)

    # Si la carpeta no existe aún:
    if not os.path.exists(ruta_destino):
        os.makedirs(ruta_destino)
        print(f"\n📁 Creada carpeta: {ruta_destino}")

    try:
        # Llamamos a la función que usa hoomd
        run_hoomd_simulation(temp=T, rho=rho_objetivo, modo='barrido', ruta_destino=ruta_destino)

    except Exception as e:
        print(f"❌ Error en la simulación T={T:.2f}: {e}")
        break

print("\n--- Todas las simulaciones han terminado ---")
