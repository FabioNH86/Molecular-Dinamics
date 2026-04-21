from funciones import calcular_promedios_energía_claude_2, generar_dataframes_todo, calcular_capacidad_calorífica
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import os
import glob

os.system('cls' if os.name == 'nt' else 'clear')

"""
Este programa busca hacer un análisis de datos para replicar los datos de simulación de la NIST mediante MD (Dinámica Molecular).
Dichos resultados se encuentran reportados en: https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-linear-force-shifted

Fabio Noriega Hernández
29/Mar/2026

Asesor: Dr. Luis Padilla
"""

# -- SEÑALA EN NÚMERO DE RONDA QUE DESEAS VISUALIZAR --
ronda = 1
ronda = int(ronda)

temperatura = 0.7
densidades = [round(x, 3) for x in np.arange(0.001, 0.1, 0.01)]

ruta_comun = f"Resultados/Isortermas/Ronda_{ronda}"

# Definimos la ruta de guardado
ruta_graficos = os.path.join(ruta_comun, 'Capacidad Calorífica')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

print(f"Analizando datos para T* = {temperatura}...")

cv_calculados = []
rhos_encontradas = []

for rho in densidades:

    nombre_carpeta = f"Rho={rho:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "todo_T*.dat")
    archivo_encontrado = glob.glob(ruta_busqueda)
    print(f'Archivo actual: {archivo_encontrado}')
    archivo_encontrado = archivo_encontrado[0]

    estabilidad = calcular_promedios_energía_claude_2(archivo=archivo_encontrado, mostrar_progreso=False)
    todo = generar_dataframes_todo(archivo=archivo_encontrado)

    cv = calcular_capacidad_calorífica(dataframe=todo, T=temperatura)

    cv_calculados.append(cv)
    rhos_encontradas.append(rho)
    print('\n')
    print('='*100)
    

# --- GENERACIÓN DEL GRÁFICO ---
if cv_calculados:
    plt.figure(figsize=(7, 6))
    
    # Graficamos usando las listas de datos encontrados
    plt.plot(rhos_encontradas, cv_calculados, 
         marker='o',           # Forma del punto
         linestyle='-',        # Estilo de línea
         color='blue',         # Color
         markersize=6,         # Tamaño del punto
         linewidth=1.5,        # Grosor de línea
         label=f'T* = {temperatura}')
    
    plt.autoscale(enable=True, axis='x', tight=True)

    # Estética estilo NIST
    plt.xlabel(r'$\rho^*$ (Densidad reducida)', fontsize=12)
    plt.ylabel(r'$c_v^*$ (Capacidad calorífica)', fontsize=12)
    plt.title(f'Isoterma de Capacidad Calorífica (Ronda {ronda})', fontsize=13)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Ajuste de límites según tu imagen de referencia
    # plt.xlim(0, 1.0)
    plt.ylim(0, max(cv_calculados) * 1.2 if cv_calculados else 3.5)

    # Guardado y visualización
    nombre_img = f'Isoterma_T{temperatura}_Ronda{ronda}.png'
    plt.savefig(os.path.join(ruta_graficos, nombre_img), dpi=300, bbox_inches='tight')
    print(f"\nÉxito: Gráfico guardado como {nombre_img}")
    plt.show()
else:
    print("\nError: No se encontraron datos para graficar.")