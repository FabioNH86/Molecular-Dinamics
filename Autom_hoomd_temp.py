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
num_prueba = 7

temperaturas = [0.65]
lista_partículas = [
    [64, 45, 25] # 77325 partículas 
]

lista_zero_momentum = [100] # Frecuencia de aplicación del zero momentum (en pasos)

longitudes_para_minibox = [55] # Longitud de la caja para la minimización de energía (en unidades reducidas)

ruta_base = f"Resultados/HOOMD/P{num_prueba}_HOOMD_Mie"

# Parámetros de tiempo
pasos_equil = int(5e5)
pasos_muestreo = int(1e6)

print(f"Se realizarán un total de {len(lista_partículas)*len(lista_zero_momentum)*len(temperaturas)*len(longitudes_para_minibox)} simulaciones :)")

for n in lista_partículas:
    nombre_config = f"N_{n[0]}_{n[1]}_{n[2]}"
    ruta_destino = os.path.join(ruta_base, nombre_config)

    # Si la carpeta no existe aún:
    if not os.path.exists(ruta_destino):
        os.makedirs(ruta_destino)
        print(f"\n📁 Creada carpeta: {ruta_destino}")

    try:
        for each in lista_zero_momentum:
            print(f"   - Zero Momentum cada {each} pasos")

            for temperatura in temperaturas:
                for longitud in longitudes_para_minibox:
                    print(f"🧪 Ejecutando simulación para ndiv={n} a T={temperatura}...")
                    # Llamamos a la función que usa hoomd
                    run_hoomd_simulation(temp=temperatura, 
                                        modo='barrido', 
                                        ruta_destino=ruta_destino, 
                                        length_minibox=longitud,                                     
                                        equilibracion=pasos_equil,
                                        muestreo=pasos_muestreo,
                                        ndiv_entrada=n,
                                        periodic_zeromomentum=each)

    except Exception as e:
        print(f"❌ Error en la configuración {n}: {e}")
        break

print("\n--- Todas las simulaciones han terminado ---")
