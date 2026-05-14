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

# 1. Definimos los límites "razonables"
rho_min_liquido = 0.70  # Umbral para considerar que hay una fase líquida clara
rho_max_vapor = 0.10    # Umbral para la fase vapor
dx_usado = 0.50

# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
entrada = 3
num_prueba = int(entrada)

# Lista de Temperaturas
temperaturas = [1.20]
ruta_comun = f'Resultados/P{num_prueba}_LV_Mie'
todos_los_datos = {}
propiedades_referencia = None
datos_campana = []

# Variables globales para graficar después del bucle
Lx_global = None 
volumen_bin_global = None
dx_usado_global = None

def obtener_propiedades_sistema(ruta): 
    with open(ruta, 'r') as f:
        lineas = f.readlines()
        num_particulas = int(lineas[0].strip().split()[0])
        cordenadas_xyz = lineas[1].strip().split()
        dimensiones_caja = [float(cordenadas_xyz[0]), float(cordenadas_xyz[1]), float(cordenadas_xyz[2])] 
        return [num_particulas, dimensiones_caja]
    
def obtener_dimensiones_bin(dx, propiedades_caja):
    Lx, Ly, Lz = propiedades_caja[1]    
    dimensiones_bin = [dx, Ly, Lz]
    volumen_bin = dx * Ly * Lz 
    return [volumen_bin, dimensiones_bin]
    
# Iniciar figura principal para los perfiles
plt.figure(figsize=(10, 6))

# Se recorren todas las carpetas de los Resultados
for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "xyz_T*.dat")
    archivos_encontrados = glob.glob(ruta_busqueda)

    if archivos_encontrados:
        ruta_archivo = archivos_encontrados[0]
        
        try:
            propiedades_actuales = obtener_propiedades_sistema(ruta=ruta_archivo)

            if propiedades_referencia is None: 
                propiedades_referencia = propiedades_actuales

            else:
                if propiedades_actuales[0] != propiedades_referencia[0]:
                    print(f'⚠️ ¡El Número de partículas no coincide en T={T:.2f}!')
                elif propiedades_actuales[1] != propiedades_referencia[1]:
                    print(f'⚠️ ¡Las dimensiones de la caja no coinciden en T={T:.2f}!')
                
            volumen_bin, dimensiones_bin = obtener_dimensiones_bin(dx=dx_usado, propiedades_caja=propiedades_referencia)
            num_particulas, dimensiones_caja = propiedades_referencia
            Lx, Ly, Lz = dimensiones_caja
            
            # Guardamos para usar fuera del bucle
            Lx_global = Lx 
            volumen_bin_global = volumen_bin
            dx_usado_global = dx_usado

            # Leer datos y limpiar
            df_temp = pd.read_csv(ruta_archivo, sep=r'\s+', skiprows=2, header=None)
            df_temp = df_temp.iloc[:, :-3] 
            coordenadas_x = df_temp[0].values

            # Centrado
            cuentas_pre, bordes_pre = np.histogram(coordenadas_x, bins=300, range=(0, Lx))
            idx_max = np.argmax(cuentas_pre)
            pos_liquido = (bordes_pre[idx_max] + bordes_pre[idx_max+1]) / 2
            desplazamiento = (Lx / 2) - pos_liquido
            coordenadas_centradas = (coordenadas_x + desplazamiento) % Lx

            # Cálculo de densidad final
            num_bines = int(Lx / dx_usado)
            bordes_bines = np.linspace(0, Lx, num_bines+1)
            particulas_bin, _ = np.histogram(coordenadas_centradas, bins=bordes_bines)
            densidad_perfil = particulas_bin / volumen_bin
            centros_bines = (bordes_bines[:-1] + bordes_bines[1:]) / 2

            # Guardar resultados
            todos_los_datos[T] = {
                'centros': centros_bines,
                'densidad': densidad_perfil
            }
            print(f"✅ Perfil calculado para T={T:.2f}")

        except Exception as e:
            print(f"❌ Error al procesar {ruta_archivo}: {e}")
            continue # Saltamos a la siguiente temperatura si hay error
    else:
        print(f"⚠️ No se encontró archivo en: {nombre_carpeta}")
        continue # Saltamos si no hay archivo

    # --- ANÁLISIS DE LA CURVA (Asegurado dentro de la iteración con datos válidos) ---
    x = todos_los_datos[T]['centros']
    y = todos_los_datos[T]['densidad']
    
    # Redujimos el rango para el líquido para evitar las interfases (0.45 a 0.55)
    mascara_liquido = (x > Lx * 0.45) & (x < Lx * 0.55)
    rho_L = np.mean(y[mascara_liquido])

    mascara_vapor = (x < Lx * 0.1) | (x > Lx * 0.9)
    rho_V = np.mean(y[mascara_vapor])

    # Cambio Por CHAT-GPT
    #rho_L = np.mean(np.sort(y)[-10:])   # top densidades
    #rho_V = np.mean(np.sort(y)[:10])    # bottom densidades
    
    max_rho = np.max(y)
    min_rho = np.min(y)
    
    if max_rho >= rho_min_liquido and min_rho <= rho_max_vapor:
        color_curva = 'dodgerblue'
        label_prefix = "Coexistencia L-V"
        alpha = 1.0
        linewidth = 2
        datos_campana.append([T, rho_L, rho_V])
        print(f"   -> Coexistencia: rho_L = {rho_L:.4f}, rho_V = {rho_V:.4f}")
    else:
        color_curva = 'gray'
        label_prefix = "No coexistente"
        alpha = 0.3
        linewidth = 1
        print(f"   -> Homogénea o crítica.")

    plt.plot(x, y, label=f'{label_prefix} T={T:.2f}', color=color_curva, alpha=alpha, lw=linewidth)


# --- IMPRESIÓN DE RESUMEN ---
if propiedades_referencia:
    print(f'\n--- Resumen del Sistema ---')
    print(f'Dimensiones de la caja: {propiedades_referencia[1]}')
    print(f'Número de partículas N: {propiedades_referencia[0]}')
    print(f'Volumen del bin usado: {volumen_bin_global:.2f}, ancho dx: {dx_usado_global:.2f}')   

# --- CONFIGURACIÓN GRÁFICA DE PERFILES ---
if Lx_global is not None:
    plt.axhline(y=rho_min_liquido, color='red', linestyle='--', alpha=0.6, label=f'Umbral Liquido (rho >= {rho_min_liquido})')
    plt.axhline(y=rho_max_vapor, color='green', linestyle='--', alpha=0.6, label=f'Umbral Vapor (rho <= {rho_max_vapor})')
    plt.fill_between([0, Lx_global], rho_min_liquido, 1.0, color='blue', alpha=0.05, label='Region Liquida')
    plt.fill_between([0, Lx_global], 0, rho_max_vapor, color='green', alpha=0.05, label='Region Vapor')

plt.xlabel('Posicion en X (sigma)', fontsize=12)
plt.ylabel('Densidad local (rho*)', fontsize=12)
plt.title(f'Analisis de Fases - Prueba {num_prueba}', fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
plt.grid(alpha=0.2)
plt.tight_layout()
plt.show()

# --- GRÁFICA DE LA CAMPANA ---
if datos_campana:
    campana = np.array(datos_campana)
    # Ordenar por temperatura por si acaso, ayuda a visualizar mejor las líneas
    campana = campana[campana[:, 0].argsort()]
    
    plt.figure(figsize=(8, 6))
    plt.plot(campana[:, 1], campana[:, 0], 'o-', label='Liquido', color='tab:blue')
    plt.plot(campana[:, 2], campana[:, 0], 'o-', label='Vapor', color='tab:orange')
    plt.xlabel('Densidad (rho*)', fontsize=12)
    plt.ylabel('Temperatura (T*)', fontsize=12)
    plt.title('Curva de Coexistencia L-V (Diagrama de Fases)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()