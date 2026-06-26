import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import os
from funciones import calcular_perfil_densidad_gsd
import glob

os.system('clear')

"""
Análisis de Coexistencia LV para simulaciones HOOMD-blue.
Comparación con datos de referencia NIST (Mie potential).

Fabio Noriega Hernández
Abril 2026
"""

# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
entrada = 1
num_prueba = int(entrada)
num_bines = 100 # Establece el número de fracciones en que se dividirá la caja de simulación (siempre en el eje x)
start_frame = 10  # Ajustar según cuándo empiece el equilibrio en tus GSD

# -- VALORES REPORTADOS POR LA NIST
# nist_rho_l = [8.643e-01, 8.426e-01, 8.203e-01, 7.970e-01, 7.728e-01, 7.474e-01, 7.203e-01, *, 6.592e-01, 6.233e-01, 5.807e-01, 5.238e-01, 4.367e-01]
# nist_rho_v = [8.450e-04, 1.828e-03, 3.508e-03, 6.146e-03, 1.004e-02, 1.553e-02, 2.304e-02, *, 4.664e-02, 6.489e-02, 9.047e-02, 1.310e-01, 2.047e-01]

nist_rho_v = [8.450e-04, 3.508e-03, 1.004e-02, 2.304e-02, 4.664e-02, 9.047e-02, 2.047e-01]
nist_rho_l = [8.643e-01, 8.203e-01, 7.728e-01, 7.203e-01, 6.592e-01, 5.807e-01, 4.367e-01]


    # Valores de error de la NIST
errores_nist_v = []
errores_nist_l = []


# Lista de Temperaturas
#temperaturas_originales = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temperaturas_originales = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]
temperaturas = [T for T in temperaturas_originales if T != 0.95] # Se omite la temperatura con error

ruta_comun = f'Resultados/P{num_prueba}_HOOMD_Mie'

# Definimos la ruta de guardado
ruta_graficos = os.path.join(ruta_comun, 'Graficos_Análisis')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

# Lista de densidades de líquido y vapor
densidades_liquido = []
densidades_vapor = []
std_liquido = []
std_vapor = []
temps_ejecutadas = []

centro = num_bines // 2
margen = 5

for T in temperaturas_originales:

    if T == 0.95: # Se omite el procesado de datos a esa temperatura
        continue


    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "movie.gro")
    archivo_encontrado = glob.glob(ruta_busqueda)
    print(f'Archivo actual: {archivo_encontrado}')

    if not archivo_encontrado:
        print(f"⚠️ No se encontró archivo para T={T:.2f}")
        continue

    archivo_gsd = archivo_encontrado[0]
    print(f'Procesando GSD: {archivo_gsd}')

    # --- LLAMADA A LA NUEVA FUNCIÓN ---
    try:
        x, rho_prom, rho_std = calcular_perfil_densidad_gsd(
            gsd_file=archivo_gsd, 
            start_frame=start_frame, 
            num_bines=num_bines
        )

        # Cálculo de densidades de fase (Líquido al centro, Vapor a los extremos)
        # Nota: HOOMD suele centrar la fase densa si usaste la lógica de centrado
        bines_liquido = rho_prom[centro - margen : centro + margen]
        rho_l = bines_liquido.mean()
        std_l = rho_std[centro - margen : centro + margen].mean()

        bines_vapor = np.concatenate([rho_prom[:10], rho_prom[-10:]])
        rho_v = bines_vapor.mean()
        std_v = rho_std[:10].mean()

        densidades_liquido.append(rho_l)
        std_liquido.append(std_l)
        densidades_vapor.append(rho_v)
        std_vapor.append(std_v)
        temps_ejecutadas.append(T)

        print(f"T={T:.2f} | rho_L: {rho_l:.4f} | rho_V: {rho_v:.4f}")
        print("-" * 60)

    except Exception as e:
        print(f"❌ Error procesando T={T}: {e}")



# --- GRÁFICA DE COEXISTENCIA ---
plt.figure(figsize=(8, 6))

# Datos Simulación
plt.errorbar(densidades_liquido, temps_ejecutadas, xerr=std_liquido, 
             fmt='o', color='royalblue', label='HOOMD: Líquido', capsize=3)
plt.errorbar(densidades_vapor, temps_ejecutadas, xerr=std_vapor, 
             fmt='o', color='crimson', label='HOOMD: Vapor', capsize=3)

# Datos NIST (Ajustamos el slice para que coincida con lo procesado)
n_puntos = len(temps_ejecutadas)
plt.scatter(nist_rho_l[:n_puntos], temperaturas_originales[:n_puntos], 
            marker='x', color='navy', label='NIST: Líquido', zorder=5)
plt.scatter(nist_rho_v[:n_puntos], temperaturas_originales[:n_puntos], 
            marker='x', color='darkred', label='NIST: Vapor', zorder=5)

plt.title(f'Envolvente de Coexistencia - Ensayo {num_prueba}')
plt.xlabel('Densidad Reducida $\\rho^*$')
plt.ylabel('Temperatura Reducida $T^*$')
plt.legend()
plt.grid(alpha=0.3)

plt.savefig(os.path.join(ruta_graficos, 'campana_coexistencia.png'), dpi=300)
plt.show()

# --- MÉTRICAS DE ERROR ---
def calcular_error(sim, ref):
    sim, ref = np.array(sim), np.array(ref[:len(sim)])
    mae = np.mean(np.abs(sim - ref))
    mre = np.mean(np.abs((sim - ref) / ref)) * 100
    return mae, mre

mae_l, mre_l = calcular_error(densidades_liquido, nist_rho_l)
mae_v, mre_v = calcular_error(densidades_vapor, nist_rho_v)

print(f"\n{' RESULTADOS VS NIST ':^40}")
print("-" * 40)
print(f"Líquido - MAE: {mae_l:.4f} | MRE: {mre_l:.2f}%")
print(f"Vapor   - MAE: {mae_v:.4f} | MRE: {mre_v:.2f}%")