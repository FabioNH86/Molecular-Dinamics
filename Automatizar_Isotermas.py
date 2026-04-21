import os 
import shutil
import subprocess
import numpy as np
from funciones import actualizar_entradas


"""
Este código se encargar de automatizar simulaciones alterando el archivo in.dat después de cada 
corrida según una lista de temperaturas pre-asignadas. 
La simulación es llevada acabo por el ejecutable `mie`. Los archivos generados se guardan automáticamente
(movie.gro, presiones.dat, resumen.dat, todo.dat y xyz.dat) 
en una carpeta que señala la temperatura usada (Resultados/P{num_prueba}_LV_Mie/T=*).

Fabio Noriega Hernández
"""
# -- SEÑALA EL NÚMERO DE ENSAYO QUE HARÁS PARA ALMACENAR LOS RESULTADOS EN SU CARPETA CORRESPONDIENTE --
num_prueba = 9


temperatura = 0.7
densidades = [round(x, 3) for x in np.arange(0.005, 0.9, 0.105)]
# print(densidades)
print(f'Se realizarán un total de: {len(densidades)} simulaciones.')
ruta_base = f"Resultados/Isortermas"
ejecutable = "./mie_exec"


for rho in densidades:
    nombre_carpeta = f"Rho={rho:.2f}"
    ruta_destino = os.path.join(ruta_base, nombre_carpeta)

    # Si la carpeta no existe aún:
    if not os.path.exists(ruta_destino):
        os.makedirs(ruta_destino)
        print(f"\n📁 Creada carpeta: {ruta_destino}")

    # Se actulaiza el archivo in.dat
    actualizar_entradas(densidad=rho, temp=temperatura, densidad_constante=False)
    print(f"🌡️ Configurando T = {temperatura:.2f}...")
    print(f"> Configurando Rho = {rho:.2f} para caja centrada")


    # Se ejecuta la simulación
    try:
        print(f"🚀 Ejecutando simulación para {nombre_carpeta}...")
        subprocess.run(ejecutable, check=True)

        # Se mueven los archivos generados
        archivos_a_mover = ["movie.gro", "presiones.dat", "resumen.dat", "todo.dat", "xyz.dat", "gr_dm.dat", "perfil_dm.dat"]

        for archivo in archivos_a_mover:
            if os.path.exists(archivo):
                # Se actuliza el nombre del archivo para evitar confusiones
                nombre_nuevo = archivo.replace(".dat", f"_T{str(temperatura).replace('.', '-')}.dat")
                shutil.move(archivo, os.path.join(ruta_destino, nombre_nuevo))
            else:
                print(f"⚠️ Advertencia: No se encontró {archivo}")

        print(f"✅ Finalizada T={temperatura:.2f}. Archivos guardados. \n")

    except subprocess.CalledProcessError:
        print(f"❌ Error crítico en la ejecución de 'new_mie' para T={temperatura:.2f}")
        break

print("\n--- Todas las simulaciones han terminado ---")
