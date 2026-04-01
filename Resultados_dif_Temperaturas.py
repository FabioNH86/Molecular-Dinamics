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
ancho_bin = 0.50


# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
#entrada = input("Selecciona el número de prueba que deseas graficar: \n >> ")
entrada = 2
num_prueba = int(entrada)


# Lista de Temperaturas
temperaturas = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
temp_colors = plt.cm.viridis(np.linspace(0, 1, len(temperaturas)))
ruta_comun = f'Resultados/P{num_prueba}_LV_Mie'
todos_los_datos = {}
propiedades_referencia = None
datos_campana = []

# Cálculo de densidades 
def obtener_propiedades_sistema(ruta): # Se obtienen las propiedades directamente de cada archivo xyz.dat a analizar 
    with open(ruta, 'r') as f:
        lineas = f.readlines()
        num_particulas = int(lineas[0].strip().split()[0])
        cordenadas_xyz = lineas[1].strip().split()
        dimensiones_caja = [float(cordenadas_xyz[0]), float(cordenadas_xyz[1]), float(cordenadas_xyz[2])] 
       
        return [num_particulas, dimensiones_caja]
    
def obtener_dimensiones_bin(dx, propiedades_caja):
    # El bin tendrá una fracción dx de la longitud total Lx de la caja de simulación.
    Lx, Ly, Lz = propiedades_caja[1]    
    dimensiones_bin = [dx, Ly, Lz]

    volumen_bin = dx * Ly * Lz 
    
    return [volumen_bin, dimensiones_bin]
    

# Se recorren todas las carpetas de los Resultados
for i, T in enumerate(temperaturas):
    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "xyz_T*.dat")
    
    # glob.glob nos devuelve una lista de todos los archivos que coincidan
    archivos_encontrados = glob.glob(ruta_busqueda)

    if archivos_encontrados:
        # Tomamos el primero que encuentre
        ruta_archivo = archivos_encontrados[0]
        
        try:
            propiedades_actuales = obtener_propiedades_sistema(ruta=ruta_archivo)

            if propiedades_referencia is None: 
                propiedades_referencia = propiedades_actuales

            else:
                if propiedades_actuales[0] != propiedades_referencia[0]:
                    print(f'⚠️ El Número de partículas no coinciden!')
                    print(f'El problema se encuentra en T={T:.2f}')
                    print(f"Esperado: {propiedades_referencia[0]} | Encontrado: {propiedades_actuales[0]}")
                    
                elif propiedades_actuales[1] != propiedades_referencia[1]:
                    print(f'⚠️ Las dimensiones de la caja de simulación no coinciden!')
                    print(f'El problema se encuentra en T={T:.2f}')
                    print(f"Esperado: {propiedades_referencia[1]} | Encontrado: {propiedades_actuales[1]}")
                
            # Se define el volumen del bin
            volumen_bin, dimensiones_bin = obtener_dimensiones_bin(dx=ancho_bin, propiedades_caja=propiedades_referencia)
            dx_usado = dimensiones_bin[0]

            num_particulas, dimensiones_caja = propiedades_referencia
            Lx, Ly, Lz = dimensiones_caja
            volumen_caja = Lx * Ly * Lz

            # Leer datos saltando las 2 filas de cabecera
            df_temp = pd.read_csv(ruta_archivo, sep=r'\s+', skiprows=2, header=None)
            
            # Limpiar: Eliminar las últimas dos columnas
            df_temp = df_temp.iloc[:, :-3] # Nos quedamos solo con las posiciones de las partículas

            coordenadas_x = df_temp[0].values

            # 1. Histogramas rápidos para encontrar dónde está el líquido en este archivo
            cuentas_pre, bordes_pre = np.histogram(coordenadas_x, bins=50, range=(0, Lx))
            # 2. Encontramos el bin con la densidad MÁXIMA (el corazón del líquido)
            idx_max = np.argmax(cuentas_pre)
            pos_liquido = (bordes_pre[idx_max] + bordes_pre[idx_max+1]) / 2

            centro_masa_x = np.mean(coordenadas_x)

            # 3. Calculamos cuánto hay que moverlo para que ese máximo esté en Lx/2
            # Usamos el módulo para que la rotación sea fluida
            desplazamiento = (Lx / 2) - pos_liquido
            coordenadas_centradas = (coordenadas_x + desplazamiento) % Lx

            # Guardar en el diccionario
            todos_los_datos[T] = df_temp
            print(f"Cargado con éxito: {os.path.basename(ruta_archivo)} para T={T:.2f}")

            # Dividimos el volumen de la caja entre el volumen del bin 
            num_bines = int(Lx / dx_usado)
            bordes_bines = np.linspace(0, Lx, num_bines+1)

            particulas_bin, _ = np.histogram(coordenadas_centradas, bins=bordes_bines)

            densidad_perfil = particulas_bin / volumen_bin
            centros_bines = (bordes_bines[:-1] + bordes_bines[1:]) / 2

            todos_los_datos[T] = {
                'centros': centros_bines,
                'densidad': densidad_perfil
            }

            
        except Exception as e:
            print(f"Error al leer {ruta_archivo}: {e}")
        
        print(f"✅ Perfil calculado para T={T:.2f} ({num_bines} bines)")


    else:
        print(f"No se encontró ningún archivo xyz_T en: {nombre_carpeta}")


    # Una vez centradas:
    # Rango del líquido (Centro): de Lx*0.4 a Lx*0.6
    mascara_liquido = (centros_bines > Lx*0.4) & (centros_bines < Lx*0.6)
    rho_L = np.mean(densidad_perfil[mascara_liquido])

    # Rango del vapor (Extremos): de 0 a Lx*0.1 y de Lx*0.9 a Lx
    mascara_vapor = (centros_bines < Lx*0.1) | (centros_bines > Lx*0.9)
    rho_V = np.mean(densidad_perfil[mascara_vapor])
    
    print(f"Para T={T:.2f}: rho_L = {rho_L:.4f}, rho_V = {rho_V:.4f}")


    if T in todos_los_datos:
        x = todos_los_datos[T]['centros']
        y = todos_los_datos[T]['densidad']
        
        # 2. Criterio: ¿El máximo de esta curva supera el umbral de líquido?
        max_rho = np.max(y)
        min_rho = np.min(y)
        
        # Solo se pinta de azul si cumple AMBAS condiciones
        if max_rho >= rho_min_liquido and min_rho <= rho_max_vapor:
            color_curva = 'dodgerblue'
            label_prefix = "Coexistencia L-V"
            alpha = 1.0
            linewidth = 2
            datos_campana.append([T, rho_L, rho_V])
        else:
            # Si falla alguna, se va a gris (fase homogénea o cerca del punto crítico)
            color_curva = 'gray'
            label_prefix = "No coexistente"
            alpha = 0.3
            linewidth = 1

        plt.plot(x, y, label=f'{label_prefix} T={T:.2f}', color=color_curva, 
                 alpha=alpha, lw=linewidth)


    



print(f'\nLas dimensiones de las cajas fueron: {propiedades_referencia[1]}')
print(f'El número de partículas N fue: {propiedades_referencia[0]}')
print(f'El volumen del bin usado es: {volumen_bin}, con un ancho de: {dx_usado}')   





# 3. Dibujar líneas horizontales de delimitación (Texto Plano)
plt.axhline(y=rho_min_liquido, color='red', linestyle='--', alpha=0.6, 
            label=f'Umbral Liquido (rho >= {rho_min_liquido})')
plt.axhline(y=rho_max_vapor, color='green', linestyle='--', alpha=0.6, 
            label=f'Umbral Vapor (rho <= {rho_max_vapor})')

# 4. Sombreado de zonas
plt.fill_between([0, Lx], rho_min_liquido, 1.0, color='blue', alpha=0.05, label='Region Liquida')
plt.fill_between([0, Lx], 0, rho_max_vapor, color='green', alpha=0.05, label='Region Vapor')

# Configuración estética 
plt.xlabel('Posicion en X (sigma)', fontsize=12)
plt.ylabel('Densidad local (rho*)', fontsize=12)
plt.title(f'Analisis de Fases - Prueba {num_prueba}', fontsize=14)

# Ajuste de la leyenda
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
plt.grid(alpha=0.2)

plt.tight_layout() # Importante para que la leyenda no se corte
plt.show()

if datos_campana:
    campana = np.array(datos_campana)
    plt.figure()
    plt.plot(campana[:, 1], campana[:, 0], 'o-', label='Liquido')
    plt.plot(campana[:, 2], campana[:, 0], 'o-', label='Vapor')
    plt.xlabel('Densidad (rho*)')
    plt.ylabel('Temperatura (T*)')
    plt.title('Curva de Coexistencia L-V')
    plt.legend()
    plt.show()