import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import os
from funciones import calcular_densidades
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
entrada = 10
num_prueba = int(entrada)
num_bines = 100 # Establece el número de fracciones en que se dividirá la caja de simulación (siempre en el eje x)
un_tercio = num_bines // 3

# -- VALORES REPORTADOS POR LA NIST
# nist_rho_v = [8.450e-04, 1.828e-03, 3.508e-03, 6.146e-03, 1.004e-02, 1.553e-02, 2.304e-02, *, 4.664e-02, 6.489e-02, 9.047e-02, 1.310e-01, 2.047e-01]
# nist_rho_l = [8.643e-01, 8.426e-01, 8.203e-01, 7.970e-01, 7.728e-01, 7.474e-01, 7.203e-01, *, 6.592e-01, 6.233e-01, 5.807e-01, 5.238e-01, 4.367e-01]

nist_rho_v = [8.450e-04]
nist_rho_l = [8.643e-01]


    # Valores de error de la NIST
errores_nist_v = []
errores_nist_l = []


# Lista de Temperaturas
#temperaturas_originales = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temperaturas_originales = [0.60]
temperaturas = [T for T in temperaturas_originales if T != 0.95] # Se omite la temperatura con error

ruta_comun = f'Resultados/P{num_prueba}_LV_Mie'

# Definimos la ruta de guardado
ruta_graficos = os.path.join(ruta_comun, 'Graficos_2')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

# Lista de densidades de líquido y vapor
densidades_liquido = []
densidades_vapor = []
centro = num_bines // 2
margen = num_bines // 30

for T in temperaturas_originales:

    if T == 0.95: # Se omite el procesado de datos a esa temperatura
        continue


    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "movie.gro")
    archivo_encontrado = glob.glob(ruta_busqueda)
    print(f'Archivo actual: {archivo_encontrado}')

    if not archivo_encontrado:
        continue

    archivo_encontrado = archivo_encontrado[0]

    # Llamamos la función
    x, rho_prom, rho_std = calcular_densidades(filename=archivo_encontrado, start_conf=900000, num_bines=num_bines)



    bines_liquido = rho_prom[centro - margen : centro + margen]
    rho_liquido = np.mean(bines_liquido)
    std_rho_liquido = np.std(bines_liquido)

    bines_vapor = np.concatenate([rho_prom[:un_tercio], rho_prom[2 * un_tercio:]])
    rho_vapor = np.mean(bines_vapor)
    std_rho_vapor = np.std(bines_vapor)

    # Almacenamos estos valores para el diagrama de coexistencia
    densidades_liquido.append(rho_liquido)
    densidades_vapor.append(rho_vapor)

    print(f"T={T:.2f} | rho_L: {rho_liquido:.4f} (std: {std_rho_liquido:.4e}) | rho_V: {rho_vapor:.4f} (std: {std_rho_vapor:.4e})")

    print("="*100)
    print('\n')

    # # --- GRÁFICO SCATTER (PROYECCIÓN XY) ---
    # plt.figure(figsize=(10, 4)) # Más alargado para resaltar la caja de simulación
    # plt.scatter(atomos['X'], atomos['Y'], s=2, alpha=0.4, c='dodgerblue', edgecolors='none')
    
    # # Dibujamos los límites de la caja de simulación
    # plt.axhline(0, color='black', linewidth=1, linestyle='--')
    # plt.axvline(0, color='black', linewidth=1, linestyle='--')
    # plt.axhline(Ly, color='black', linewidth=1, linestyle='--')
    # plt.axvline(Lx, color='black', linewidth=1, linestyle='--')

    # plt.title(f'Distribución Proyección XY | $T^* = {T:.2f}$ | P{num_prueba}', fontsize=12)
    # plt.xlabel('X ($\sigma$)')
    # plt.ylabel('Y ($\sigma$)')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.grid(True, linestyle=':', alpha=0.5)
    
    # nombre_archivo = f"Scatter_T_{T:.2f}.png"
    # plt.savefig(os.path.join(ruta_graficos, nombre_archivo), dpi=150, bbox_inches='tight')
    # plt.close()



plt.figure(figsize=(9, 7))

# Graficamos las ramas unidas por una línea para formar la "campana"
# Graficamos la rama del líquido con barras de error
plt.errorbar(densidades_liquido, temperaturas, 
             xerr=std_rho_liquido, # <--- Aquí van tus desviaciones de la fase líquida
             fmt='o', color='blue', ecolor='lightblue', elinewidth=1.5, capsize=3,
             label='Simulación: Líquido', markersize=5)

plt.errorbar(densidades_vapor, temperaturas, 
             xerr=std_rho_vapor,
             fmt='o', color='red', ecolor='lightsalmon', elinewidth=1.5, capsize=3,
             label='Simulación: Vapor', markersize=5)

# --- DATOS DE LA NIST (Para validación) ---
# Solo comparamos hasta donde tengamos datos de simulación (usando zip)
plt.scatter(nist_rho_v[:len(temperaturas)], temperaturas, marker='x', s=60, color='darkred', label='NIST (Vapor)', zorder=5)
plt.scatter(nist_rho_l[:len(temperaturas)], temperaturas, marker='x', s=60, color='darkblue', label='NIST (Líquido)', zorder=5)

# --- CÁLCULO DE ERRORES ---
def calcular_metricas(mis_datos, nist_datos):
    errores_abs = [abs(m - n) for m, n in zip(mis_datos, nist_datos)]
    errores_rel = [(abs(m - n) / n) * 100 for m, n in zip(mis_datos, nist_datos)]
    return np.mean(errores_abs), np.mean(errores_rel)

mae_l, mre_l = calcular_metricas(densidades_liquido, nist_rho_l)
mae_v, mre_v = calcular_metricas(densidades_vapor, nist_rho_v)

# Estética del gráfico
plt.title(f'Diagrama de Coexistencia $T^*$ vs $\\rho^*$ | Ensayo {num_prueba}', fontsize=14)
plt.xlabel('Densidad Reducida $\\rho^*$', fontsize=12)
plt.ylabel('Temperatura Reducida $T^*$', fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend(frameon=True, loc='best')

ruta_coexistencia = os.path.join(ruta_graficos, f'Coexistencia_P{num_prueba}_vs_NIST.png')
plt.savefig(ruta_coexistencia, dpi=300, bbox_inches='tight')
plt.show()

# --- REPORTE FINAL ---
print('\n' + '='*50)
print(f"{'VALIDACIÓN FINAL VS NIST':^50}")
print('='*50)
print(f"{'Fase':<15} | {'Error Absoluto (MAE)':<20} | {'Error Relativo':<15}")
print('-'*50)
print(f"{'Líquido':<15} | {mae_l:<20.4f} | {mre_l:<15.2f}%")
print(f"{'Vapor':<15} | {mae_v:<20.4f} | {mre_v:<15.2f}%")
print('='*50)