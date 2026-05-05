from funciones import calcular_cv_hoomd
import pandas as pd 
import os
import glob
import matplotlib.pyplot as plt

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

temperatura = 0.80
densidades = [0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.017, 0.02, 0.03, 0.035, 0.04, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.055, 0.06, 0.07, 0.08, 0.09, 0.1]

ruta_comun = f"Resultados/Isotermas/Ronda_{ronda}/Temp={temperatura:.2f}"

# Definimos la ruta de guardado
ruta_graficos = os.path.join(ruta_comun, 'Capacidad Calorífica')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

print(f"Analizando datos para T* = {temperatura}")

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
    
    cv = calcular_cv_hoomd(df=df, T=temperatura)

    cv_calculados.append(cv)
    rhos_encontradas.append(rho)
    print(f"   -> Cv calculado: {cv:.4f}")

# --- GENERACIÓN DEL GRÁFICO ---
if cv_calculados:
    plt.figure(figsize=(10, 7))
    plt.plot(rhos_encontradas, cv_calculados, 'o-', color='#2c3e50', label=f'HOOMD (T*={temperatura})')
    
    for x, y in zip(rhos_encontradas, cv_calculados):
        label = f"{x:.3f}"
        plt.annotate(label,
                     (x, y),
                     textcoords="offset points",
                     xytext=(0, 10),
                     ha='center',
                     va='bottom',
                     rotation=90,
                     fontsize=9,
                     color='#34495e')

    plt.xlabel(r'Densidad $\rho^*$', fontsize=12)
    plt.ylabel(r'Capacidad Calorífica $c_v^*$', fontsize=12)
    plt.title(f'Capacidad Calorífica vs Densidad - Ronda {ronda}', fontsize=13)

    plt.ylim(min(cv_calculados)*0.9, max(cv_calculados)*1.1)
    plt.grid(True, linestyle=':', alpha=0.2)
    plt.margins(y=2.0)
    plt.legend()
    
    nombre_img = f'Cv_Isoterma_T{temperatura}.png'
    plt.savefig(os.path.join(ruta_graficos, nombre_img), dpi=300)
    print(f"\n✅ Gráfico guardado en: {ruta_graficos}")
    plt.show()