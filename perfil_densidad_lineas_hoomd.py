import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import numpy as np 
import pandas as pd 
import os
import glob
# Asegúrate de que la nueva función esté en funciones.py
from funciones import calcular_perfil_densidad_gsd 

os.system('clear')

"""
Visualización de la evolución de perfiles de densidad (Gradiente de T).
Adaptado para resultados de HOOMD-blue.

Fabio Noriega Hernández
Abril 2026
"""

# -- CONFIGURACIÓN --
num_prueba = 2
num_bines = 100
start_frame = 10  # Frame inicial para el promedio (ajustar según tu simulación)

temperaturas_originales = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]
# Filtramos si alguna temperatura dio error o no se desea graficar
temperaturas = [T for T in temperaturas_originales if T != 0.95]

ruta_comun = f'Resultados/P{num_prueba}_HOOMD_Mie'
ruta_graficos = os.path.join(ruta_comun, 'Graficos_Analisis')

if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos, exist_ok=True)

# -- CONFIGURACIÓN DEL GRÁFICO --
plt.figure(figsize=(12, 7))
norm = plt.Normalize(min(temperaturas), max(temperaturas))
cmap = cm.viridis
ax = plt.gca()

for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, f"trajectory_T{T:.2f}.gsd")
    archivos = glob.glob(ruta_busqueda)
    
    if not archivos: 
        print(f"⚠️ Saltando T={T:.2f}: Archivo GSD no encontrado.")
        continue
    
    archivo_actual = archivos[0]
    print(f'📊 Generando perfil para T={T:.2f}...')

    try:
        # Usamos la función adaptada para GSD
        x, rho_prom, rho_std = calcular_perfil_densidad_gsd(
            gsd_file=archivo_actual, 
            start_frame=start_frame, 
            num_bines=num_bines
        )

        color_T = cmap(norm(T))
        
        # Graficamos el perfil
        plt.plot(x, rho_prom, label=f'T* = {T:.2f}', color=color_T, linewidth=2, alpha=0.9)
        
        # Área de desviación estándar (sombreado tenue)
        plt.fill_between(x, rho_prom - rho_std, rho_prom + rho_std, color=color_T, alpha=0.1)

    except Exception as e:
        print(f"❌ Error en T={T}: {e}")

# --- ESTÉTICA DEL GRÁFICO ---
plt.title(f'Evolución de Perfiles de Densidad | Ensayo {num_prueba}', fontsize=14, fontweight='bold')
plt.xlabel('Posición X (unidades reducidas $\sigma$)', fontsize=12)
plt.ylabel('Densidad reducida $\\rho^*$', fontsize=12)

# Barra de colores lateral (Colorbar)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Temperatura Reducida $T^*$', fontsize=11)

plt.grid(True, linestyle=':', alpha=0.5)
# Ajustamos la leyenda para que no tape los perfiles
plt.legend(bbox_to_anchor=(1.15, 1), loc='upper left', fontsize='small', frameon=False)

plt.tight_layout()

# Guardar y mostrar
ruta_salida = os.path.join(ruta_graficos, f"Evolucion_Perfiles_P{num_prueba}.png")
plt.savefig(ruta_salida, dpi=300, bbox_inches='tight')
print(f"\n✅ Gráfico guardado en: {ruta_salida}")
plt.show()