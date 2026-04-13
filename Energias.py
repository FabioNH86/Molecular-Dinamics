import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import os
import glob

os.system('clear')

"""
Este programa busca hacer un análisis de datos para replicar los datos de simulación de la NIST mediante MD (Dinámica Molecular).
Dichos resultados se encuentran reportados en: https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-linear-force-shifted

Fabio Noriega Hernández
29/Mar/2026

Asesor: Dr. Luis Padilla
"""

# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
entrada = 4
num_prueba = int(entrada)
num_bines = 50 # Establece el número de fracciones en que se dividirá la caja de simulación (siempre en el eje x)


# -- VALORES REPORTADOS POR LA NIST
nist_rho_v = [8.450e-04, 1.828e-03, 3.508e-03, 6.146e-03, 1.004e-02, 1.553e-02, 2.304e-02, 4.664e-02, 6.489e-02, 9.047e-02, 1.310e-01, 2.047e-01]
nist_rho_l = [8.643e-01, 8.426e-01, 8.203e-01, 7.970e-01, 7.728e-01, 7.474e-01, 7.203e-01, 6.592e-01, 6.233e-01, 5.807e-01, 5.238e-01, 4.367e-01]

    # Valores de error de la NIST
errores_nist_v = []
errores_nist_l = []


# Lista de Temperaturas
temperaturas_originales = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temperaturas = [T for T in temperaturas_originales if T != 0.95] # Se omite la temperatura con error

ruta_comun = f'Resultados/P{num_prueba}_LV_Mie'

# Definimos la ruta de guardado
ruta_graficos = os.path.join(ruta_comun, 'Energías')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

# Lista de densidades de líquido y vapor
densidades_liquido = []
densidades_vapor = []

for T in temperaturas_originales:

    if T == 0.95: # Se omite el procesado de datos a esa temperatura
        continue


    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "todo_T*.dat")
    archivo_encontrado = glob.glob(ruta_busqueda)
    print(f'Archivo actual: {archivo_encontrado}')
    archivo_encontrado = archivo_encontrado[0]

    # Extraemos la información de cada configuración
    todo = pd.read_csv(archivo_encontrado, sep='\s+', header=None)
    print(todo.head())
    # Definición de variables por columna
    config      = todo[0]
    density     = todo[1]
    kinetic_e   = todo[2]
    potential_e = todo[3]
    total_e     = todo[4]
    temperature = todo[5]
    pressure    = todo[6]
    error       = todo[7]
        
    # Parámetro de corte para la fase de equilibración
    salto: int = 300000
    ventana = -50000

    # --- FIGURA 1: ANÁLISIS DE ENERGÍAS Y DENSIDAD (SUBPLOTS) ---
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')

    # Panel 1: Energía Cinética (Fase inicial)
    axs[0, 0].plot(config[ventana:], kinetic_e[ventana:], color='tab:blue')
    axs[0, 0].set_title('Energía Cinética (Equilibración)')
    axs[0, 0].set_xlabel('Configuración')
    axs[0, 0].set_ylabel('Energía')

    # Panel 2: Energía Potencial (Fase inicial)
    axs[0, 1].plot(config[ventana:], potential_e[ventana:], color='tab:red')
    #axs[0, 1].axhline(y=mean_potential_reported, color='black', linestyle='--', linewidth=1.5, label='Referencia Tabla')
    axs[0, 1].set_title('Energía Potencial (Equilibración)')
    axs[0, 1].set_xlabel('Configuración')
    axs[0, 1].set_ylabel('U / ε')
    #axs[0, 1].legend()

    # Panel 3: Energía Total (Fase inicial)
    axs[1, 0].plot(config[ventana:], total_e[ventana:], color='tab:purple', alpha=0.8)
    axs[1, 0].set_title('Energía Total del Sistema')
    axs[1, 0].set_xlabel('Configuración')
    axs[1, 0].set_ylabel('Energía')

    # Panel 4: Variación de la Presión (Toda la simulación)
    axs[1, 1].plot(config[ventana:], pressure[ventana:], color='tab:green')
    axs[1, 1].set_title('Evolución de la Densidad')
    #axs[1, 1].axhline(y=mean_pressure_reported, color='black', linestyle='--', linewidth=1.5, label='Referencia Tabla')
    axs[1, 1].set_xlabel('Configuración')
    plt.show()