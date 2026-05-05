import os 
import numpy as np
from funciones import run_hoomd_simulation


"""
Este código se encarga de automatizar simulaciones alterando el archivo in.dat después de cada 
corrida según una lista de temperaturas pre-asignadas. 
La simulación es llevada acabo por el ejecutable `mie`. Los archivos generados se guardan automáticamente
(movie.gro, presiones.dat, resumen.dat, todo.dat y xyz.dat) 
en una carpeta que señala la temperatura usada (Resultados/P{num_prueba}_LV_Mie/T=*).

Fabio Noriega Hernández
"""

# -- SEÑALA EL NÚMERO DE ENSAYO QUE HARÁS PARA ALMACENAR LOS RESULTADOS EN SU CARPETA CORRESPONDIENTE --
num_prueba = 4

temperatura = 1.20
densidades = [0.10, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2047, 0.21, 0.22, 0.23, 0.24, 0.25, 0.30] # Región para aumentar resolución (Ya calculados): [0.01, 0.05, 0.1, 0.3]
# print(densidades)
print(f'Se realizarán un total de: {len(densidades)} simulaciones.')
ruta_base = f"Resultados/Isotermas/Ronda_{num_prueba}/Temp={temperatura:.2f}"
#ejecutable = "./mie_isotermas"


for rho in densidades:
    nombre_carpeta = f"Rho={rho:.4f}"
    ruta_destino = os.path.join(ruta_base, nombre_carpeta)

    # Si la carpeta no existe aún:
    if not os.path.exists(ruta_destino):
        os.makedirs(ruta_destino)
        print(f"\n📁 Creada carpeta: {ruta_destino}")

    # Se actulaiza el archivo in.dat
    #actualizar_entradas(temp=0.7, modo='isoterma', densidad_obj=rho)
    print(f"🌡️ Iniciando T = {temperatura:.2f} | Rho = {rho:.4f}...")

    original_path = os.getcwd()

    try: 
        os.chdir(ruta_destino)

        run_hoomd_simulation(temp=temperatura, 
                             rho=rho, 
                             modo='isoterma', 
                             ruta_destino=ruta_destino, 
                             length_minibox=0,
                             equilibracion=500000,
                             muestreo=1000000)

        print(f"✅ Finalizada Rho={rho:.4f}. Datos guardados en {nombre_carpeta}/")

    except Exception as e:
        print(f"❌ Error crítico en la simulación para Rho={rho:.4f}: {e}")
        os.chdir(original_path)
        break

    finally:
        os.chdir(original_path)

    

print("\n--- Todas las simulaciones han terminado ---")
