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
ronda = 4
ronda = int(ronda)

temperatura = 0.70
densidades = [0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2]

ruta_comun = f"Resultados/Isortermas/Ronda_{ronda}/Temp={temperatura:.2f}"

# Definimos la ruta de guardado
ruta_graficos = os.path.join(ruta_comun, 'Capacidad Calorífica')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

print(f"Analizando datos para T* = {temperatura} (HOOMD-blue v6.1.1)")

cv_calculados = []
rhos_encontradas = []

for rho in densidades:

    nombre_carpeta = f"Rho={rho:.4f}"
    patron = os.path.join(ruta_comun, nombre_carpeta, "*.csv")
    archivos = glob.glob(patron)

    if not archivos:
        print(f"⚠️ No se encontró archivo para Rho={rho:.4f}")
        continue

    archivo_actual = archivos[0]
    print(f"📄 Procesando: {os.path.basename(archivo_actual)}")

    df = pd.read_csv(archivo_actual, sep=r'\s+')
    df.columns = [c.split('.')[-1] for c in df.columns]
    print(f"Columnas detectadas: {df.columns.tolist()}")

    df_produccion = df.iloc[50:].copy()

    print(df_produccion.head())

    N = 5508  # Tu número de partículas
    var_u = df_produccion['potential_energy'].var()
    
    # Térmico (gas ideal) + Fluctuación (configuracional)
    cv = 1.5 + (var_u / (N * (temperatura**2)))

    cv_calculados.append(cv)
    rhos_encontradas.append(rho)
    print(f"   -> Cv calculado: {cv:.4f}")

# --- GENERACIÓN DEL GRÁFICO ---
if cv_calculados:
    plt.figure(figsize=(8, 6))
    plt.plot(rhos_encontradas, cv_calculados, 'o-', color='#2c3e50', label=f'HOOMD (T*={temperatura})')
    
    plt.xlabel(r'Densidad $\rho^*$', fontsize=12)
    plt.ylabel(r'Capacidad Calorífica $c_v^*$', fontsize=12)
    plt.title(f'Capacidad Calorífica vs Densidad - Ronda {ronda}', fontsize=13)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    
    nombre_img = f'Cv_Isoterma_T{temperatura}.png'
    plt.savefig(os.path.join(ruta_graficos, nombre_img), dpi=300)
    print(f"\n✅ Gráfico guardado en: {ruta_graficos}")
    plt.show()