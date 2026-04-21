import os 
import shutil
import subprocess
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

#temperaturas = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temperaturas = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]
ruta_base = f"Resultados/P{num_prueba}_LV_Mie"
ejecutable = "./mie_exec"


for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    ruta_destino = os.path.join(ruta_base, nombre_carpeta)

    # Si la carpeta no existe aún:
    if not os.path.exists(ruta_destino):
        os.makedirs(ruta_destino)
        print(f"\n📁 Creada carpeta: {ruta_destino}")

    # Se actulaiza el archivo in.dat
    actualizar_entradas(T)
    print(f"🌡️ Configurando T = {T:.2f}...")


    # Se ejecuta la simulación
    try:
        print(f"🚀 Ejecutando simulación para {nombre_carpeta}...")
        subprocess.run(ejecutable, check=True)

        # Se mueven los archivos generados
        archivos_a_mover = ["movie.gro", "presiones.dat", "resumen.dat", "todo.dat", "xyz.dat", "gr_dm.dat", "perfil_dm.dat"]

        for archivo in archivos_a_mover:
            if os.path.exists(archivo):
                # Se actuliza el nombre del archivo para evitar confusiones
                nombre_nuevo = archivo.replace(".dat", f"_T{str(T).replace('.', '-')}.dat")
                shutil.move(archivo, os.path.join(ruta_destino, nombre_nuevo))
            else:
                print(f"⚠️ Advertencia: No se encontró {archivo}")

        print(f"✅ Finalizada T={T:.2f}. Archivos guardados. \n")

    except subprocess.CalledProcessError:
        print(f"❌ Error crítico en la ejecución de 'new_mie' para T={T:.2f}")
        break

print("\n--- Todas las simulaciones han terminado ---")
