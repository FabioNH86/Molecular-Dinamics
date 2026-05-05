import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob
from funciones import graficar_analisis_termo, visualizar_estabilidad_dinamica, encontrar_equilibrio_hoomd

# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
entrada = 5
num_prueba = int(entrada)
num_bines = 50 # Establece el número de fracciones en que se dividirá la caja de simulación (siempre en el eje x)

# Lista de Temperaturas
# temperaturas_originales = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temperaturas_originales = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20] 

temperaturas = [T for T in temperaturas_originales if T != 0.95] # Se omite la temperatura con error

ruta_comun = f'Resultados/P{num_prueba}_HOOMD_Mie'

# Definimos la ruta de guardado
ruta_graficos = os.path.join(ruta_comun, 'Energías')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

for T in temperaturas:
    nombre_archivo = f"todo_T{T:.2f}.csv"
    ruta_archivo = os.path.join(ruta_comun, f"T={T:.2f}", nombre_archivo)
    
    if os.path.exists(ruta_archivo):
        # 1. Calculamos el equilibrio con la función que afinamos antes
        # (Asumiendo que devuelve el paso entero o None)
        # paso_eq = visualizar_estabilidad_dinamica(archivo_csv=ruta_archivo, 
        #                                           pasos_totales=1500000,
        #                                           T=T,
        #                                           ruta_guardado=ruta_graficos,
        #                                           paso_eq=1300000)

        paso_eq = encontrar_equilibrio_hoomd(archivo_csv=ruta_archivo,
                                             pasos_totales=1500000)
        
        # 2. Graficamos
        graficar_analisis_termo(ruta_archivo, ruta_graficos=ruta_graficos, paso_equilibrio=paso_eq, pasos_totales=1500000)