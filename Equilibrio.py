from funciones import calcular_promedios_energía_claude_2, calcular_presiones_vapor, calcular_tension_superficial
import glob 
import os 

os.system('clear')

"""
Este programa busca hacer un análisis de datos para replicar los datos de simulación de la NIST mediante MD (Dinámica Molecular).
Dichos resultados se encuentran reportados en: https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-linear-force-shifted

Fabio Noriega Hernández
29/Mar/2026

Asesor: Dr. Luis Padilla
"""

# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
entrada = 11
num_prueba = int(entrada)

temperaturas_originales = [1.20]
temperaturas = [T for T in temperaturas_originales if T != 0.95] # Se omite la temperatura con error


ruta_comun = f'Resultados/P{num_prueba}_LV_Mie'

for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    ruta_todo = os.path.join(ruta_comun, nombre_carpeta, "todo_T*.dat")
    ruta_presiones = os.path.join(ruta_comun, nombre_carpeta, "presiones_T*.dat")
    archivos_todo = glob.glob(ruta_todo)
    archivos_presiones = glob.glob(ruta_presiones)

    if not archivos_todo: 
        continue
    archivo_actual_todo = archivos_todo[0]

    if not archivos_presiones:
        continue
    archivo_actual_presiones = archivos_presiones[0]

    print(f'Procesando T={T}')
    primer_bloque_estable = calcular_promedios_energía_claude_2(archivo=archivo_actual_todo, ancho_bloques=100000, mostrar_progreso=True)

    if primer_bloque_estable is not None:
        presiones, presion_vapor = calcular_presiones_vapor(archivo=archivo_actual_presiones, configuraciones_consideradas=(primer_bloque_estable))
        calcular_tension_superficial(df_presiones=presiones, longitud_partpendicular_interface=40.0)
    else:
        print(f"⚠️ Saltando cálculo de presiones para T={T} por falta de equilibrio.")

    print('='*100)
    print('\n')