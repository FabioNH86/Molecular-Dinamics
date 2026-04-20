import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import numpy as np 
import pandas as pd 
import os
import glob
from funciones import calcular_densidades

os.system('clear')

"""
Este programa busca hacer un análisis de datos para replicar los datos de simulación de la NIST mediante MD (Dinámica Molecular).
Dichos resultados se encuentran reportados en: https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-linear-force-shifted

Fabio Noriega Hernández
29/Mar/2026

Asesor: Dr. Luis Padilla
"""

# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
entrada = 9
num_prueba = int(entrada)
num_bines = 90 # Establece el número de fracciones en que se dividirá la caja de simulación (siempre en el eje x)
num_conf_buscado = 500000

# -- VALORES REPORTADOS POR LA NIST
nist_rho_v = [8.450e-04, 1.828e-03, 3.508e-03, 6.146e-03, 1.004e-02, 1.553e-02, 2.304e-02, 4.664e-02, 6.489e-02, 9.047e-02, 1.310e-01, 2.047e-01]
nist_rho_l = [8.643e-01, 8.426e-01, 8.203e-01, 7.970e-01, 7.728e-01, 7.474e-01, 7.203e-01, 6.592e-01, 6.233e-01, 5.807e-01, 5.238e-01, 4.367e-01]

# Lx = 63.496      
# Ly = 15.874
# Lz = 15.874

# Lista de Temperaturas
#temperaturas_originales = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temperaturas_originales = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]
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


# -- CONFIGURACIÓN DEL GRÁFICO CONJUNTO --
plt.figure(figsize=(12, 8))
norm = plt.Normalize(min(temperaturas), max(temperaturas))
cmap = cm.viridis
ax = plt.gca()

for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "movie.gro")
    archivos = glob.glob(ruta_busqueda)
    
    if not archivos: 
        continue
    archivo_actual = archivos[0]
    
    print(f'Procesando T={T}')

    # Llamada a tu función
    x, rho_prom, rho_std = calcular_densidades(filename=archivo_actual, start_conf=900000)

    if x is not None:
        color_T = cmap(norm(T))
        
        # Graficamos la línea de esta temperatura específica
        plt.plot(x, rho_prom, label=f'T* = {T:.2f}', color=color_T, linewidth=1.8, alpha=0.9)
        
        # Agregamos la desviación estándar como un área sombreada muy tenue
        plt.fill_between(x, rho_prom - rho_std, rho_prom + rho_std, color=color_T, alpha=0.1)

# --- ESTÉTICA FINAL (FUERA DEL BUCLE) ---
plt.title(f'Evolución de Perfiles de Densidad - Ensayo {num_prueba}', fontsize=14, fontweight='bold')
plt.xlabel('Posición X (unidades reducidas)', fontsize=12)
plt.ylabel('Densidad reducida $\\rho^*$', fontsize=12)

# Añadimos una barra de colores lateral para que se entienda el gradiente de T
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Temperatura Reducida $T^*$')

plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(bbox_to_anchor=(1.25, 1), loc='upper left', fontsize='small', ncol=1)
plt.tight_layout()

# Guardar resultado
ruta_perfil = os.path.join(ruta_graficos, f"Perfiles_Densidad_P{num_prueba}.png")
plt.savefig(ruta_perfil, dpi=300, bbox_inches='tight')
plt.show()