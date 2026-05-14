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
num_prueba = 6

temperatura = 1.20
lista_partículas = [[50, 20, 20], [60, 30, 30], [70, 30, 30], [70, 40, 40]] 
ruta_base = f"Resultados/P{num_prueba}_HOOMD_Mie"



for n in lista_partículas:
    nombre_config = f"N_{n[0]}_{n[1]}_{n[2]}"
    ruta_destino = os.path.join(ruta_base, nombre_config)

    # Si la carpeta no existe aún:
    if not os.path.exists(ruta_destino):
        os.makedirs(ruta_destino)
        print(f"\n📁 Creada carpeta: {ruta_destino}")

    try:
        print(f"Se realizarán un total de {len(lista_partículas)} simulaciones :)")
        print(f"🧪 Ejecutando simulación para ndiv={n} a T={temperatura}...")
        # Llamamos a la función que usa hoomd
        run_hoomd_simulation(temp=temperatura, 
                             modo='barrido', 
                             ruta_destino=ruta_destino, 
                             length_minibox=30.0, 
                             equilibracion=1e6,
                             muestreo=1e6,
                             ndiv_entrada=n)

    except Exception as e:
        print(f"❌ Error en la configuración {n}: {e}")
        break

print("\n--- Todas las simulaciones han terminado ---")
