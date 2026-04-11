import matplotlib.pyplot as plt 
import matplotlib.cm as cm
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
num_bines = 90 # Establece el número de fracciones en que se dividirá la caja de simulación (siempre en el eje x)


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
ruta_graficos = os.path.join(ruta_comun, 'Graficos_2')

# Si la carpeta no existe, la creamos automáticamente
if not os.path.exists(ruta_graficos):
    os.makedirs(ruta_graficos)
    print(f"Carpeta creada: {ruta_graficos}")

# Lista de densidades de líquido y vapor
densidades_liquido = []
densidades_vapor = []


# -- GENERACIÓN DE UN SOLO PERFIL --
plt.figure(figsize=(16, 10))
centros_bines = np.linspace(0, 63.496, num_bines)
norm = plt.Normalize(min(temperaturas), max(temperaturas))
cmap = cm.viridis # 'viridis' es excelente para escalas térmicas
perfiles_por_temperatura = {} # Diccionario para guardar {T: [densidades]}
ax = plt.gca()

for T in temperaturas_originales:

    if T == 0.95: # Se omite el procesado de datos a esa temperatura
        continue


    nombre_carpeta = f"T={T:.2f}"
    ruta_busqueda = os.path.join(ruta_comun, nombre_carpeta, "xyz_T*.dat")
    archivo_encontrado = glob.glob(ruta_busqueda)
    print(f'Archivo actual: {archivo_encontrado}')
    archivo_encontrado = archivo_encontrado[0]

    # Extraemos las dimensiones de la caja 
    with open(archivo_encontrado, 'r') as f:
        f.readline()
        segunda_fila = f.readline()

        Lx, Ly, Lz = map(float, segunda_fila.split())

    volumen_caja = Lx * Ly * Lz
    print(f"Dimensiones extraídas: Lx={Lx:.4f}, Ly={Ly:.4f}, Lz={Lz:.4f}")
    print(f"Volumen de la caja: {volumen_caja:.4f}")

    # Hacemos las dataframes de las coordenadas
    posiciones_velocidades = pd.read_csv(archivo_encontrado, sep=r'\s+', skiprows=2, header=None)
    coordenadas = posiciones_velocidades.iloc[:, :-3] 
    coordenadas.columns = ['x', 'y', 'z']
    #print(coordenadas.head())

    #  -- INICIA EL ACOMODO DE LAS PARTÍCULAS --
    # Se buscan los valores más negativos 
    x_min = coordenadas['x'].min()
    y_min = coordenadas['y'].min()
    z_min = coordenadas['z'].min()
    print(f'Coordenadas orinales: {x_min:.4f}, {y_min:.4f}, {z_min:.4f}')

    # Se acomodan las partículas con el origen en una esquina de la caja de simulación
    new_x = coordenadas['x'] - x_min
    new_y = coordenadas['y'] - y_min
    new_z = coordenadas['z'] - z_min


    coordenadas_norm = pd.DataFrame({
    'x': new_x,
    'y': new_y,
    'z': new_z
    })
    print('Nuevas coordenadas')
    print(coordenadas_norm.head())

    # --- GENERACIÓN DEL GRÁFICO POR TEMPERATURA (DENTRO DEL CICLO) ---
    print(f'Generando gráfico para T = {T:.2f}...')

    cortes_x = np.linspace(0, Lx, num_bines + 1)

    #plt.figure(figsize=(12, 6))
    

    num_particulas_global = len(coordenadas_norm)
    densidad_caja = num_particulas_global/volumen_caja

    print(f'La densidad de partículas en la caja es: {densidad_caja:.4f} con {num_particulas_global} partículas totales')

        # -- MEDICIÓN DE DENSIDADES POR BINES (SECCIONES DE LA CAJA DE SIMULACIÓN) --
    particulas_por_bin = []
    densidades_locales = []
    dx = Lx / num_bines

    volumen_bin = dx * Ly * Lz

    for i in range(num_bines):
        # Definimos los límites del bin
        x_inicial = cortes_x[i]
        x_final = cortes_x[i+1]

        # DEBUG: Vemos los límites del bin en 2D
        #plt.axvline(x=x_final, color='red', linewidth=1.0, alpha=0.3)

        bin_actual = coordenadas_norm[(coordenadas_norm['x'] >= x_inicial) & (coordenadas_norm['x'] < x_final)]

        num_particulas_bin = len(bin_actual)

            # Calculamos las densidades de los bines
        densidad_local = num_particulas_bin / volumen_bin

            # Almacenamos los datos
        densidades_locales.append(densidad_local)
        particulas_por_bin.append(num_particulas_bin)

    print(f"Conteo completado. Densidad más alta registrada: {max(densidades_locales):.4f}")
    print(f"Conteo completado. Densidad más baja registrada: {min(densidades_locales):.4f}")

    # -- COMIENZA LA SELECCIÓN DE DENSIDADES LÍQUIDAS Y GASEOSAS --
        # Se convierten en serie las densidades de cada bin
    densidades_series = pd.Series(densidades_locales)



        # Tomamos las 10 densidades más altas y las 10 más bajas para líquido y vapor respectivamente
    top_bines = densidades_series.nlargest(5)
    rho_liquido = top_bines.mean()
    std_rho_liquido = top_bines.std() 

    rho_vapor = densidades_series.nsmallest(5).mean()
    std_rho_vapor = rho_vapor.std()

    # Almacenamos estos valores para graficar 
    densidades_liquido.append(rho_liquido)
    densidades_vapor.append(rho_vapor)

    print(f"Densidad Líquida estimada: {rho_liquido:.4f}, con una desviación estándar de: {std_rho_liquido}")
    print(f"Densidad Vapor estimada: {rho_vapor:.4f}, con una desviación estándar de: {std_rho_vapor}")

    print('')

    
    # --- GENERACIÓN DEL PERFIL DE DENSIDAD (Línea continua) ---
    #plt.figure(figsize=(10, 5))
    
    # Calculamos el centro de cada bin para que el eje X sea la posición exacta
    # cortes_x tiene los bordes, así que promediamos para obtener el centro.
    # centros_bines = (cortes_x[:-1] + cortes_x[1:]) / 2

    # # Graficamos la línea continua de densidad local
    # plt.plot(centros_bines, densidades_locales, color='teal', linewidth=2, label='Perfil de densidad')
    
    # # Opcional: Rellenar el área bajo la curva para mejor visualización
    # plt.fill_between(centros_bines, densidades_locales, color='teal', alpha=0.2)

    # # --- Verificación visual de los promedios calculados ---
    # # Dibujamos líneas horizontales para mostrar qué valores tomó el script como rho_L y rho_V
    # if not np.isnan(rho_liquido):
    #     plt.axhline(rho_liquido, color='blue', linestyle='--', alpha=0.7, label=f'Promedio Líquido: {rho_liquido:.3f}')
    # plt.axhline(rho_vapor, color='red', linestyle='--', alpha=0.7, label=f'Promedio Vapor: {rho_vapor:.3f}')

    # # Configuración de etiquetas y título
    # plt.title(f'Perfil de Densidad Local\nTemperatura = {T:.2f} (Prueba P{num_prueba})', fontsize=14, fontweight='bold')
    # plt.xlabel('Posición en el eje X de la caja ($L_x$)', fontsize=12)
    # plt.ylabel('Densidad Reducida $\\rho^*$', fontsize=12)
    
    # # Ajustar límites
    # plt.xlim(0, Lx)
    # plt.ylim(0, max(densidades_locales) * 1.1) # Un 10% más arriba del máximo para que respire el gráfico
    
    # plt.grid(True, linestyle=':', alpha=0.6)
    # plt.legend()
    # plt.tight_layout()
    
    # # Cambiamos el nombre del archivo para reflejar que es un perfil de densidad
    # nombre_archivo = f"Perfil_Densidad_T_{T:.2f}.png"
    # ruta_final_archivo = os.path.join(ruta_graficos, nombre_archivo)
    
    # # Guardamos la imagen
    # plt.savefig(ruta_final_archivo, dpi=300, bbox_inches='tight')
    # print(f"Perfil de densidad guardado en: {ruta_final_archivo}")
    
    # # IMPORTANTE: Cerrar la figura
    # plt.close()
    
    # Aplicamos el centrado que discutimos antes
    idx_max = np.argmax(densidades_locales)
    shift = (num_bines // 2) - idx_max
    densidades_centradas = np.roll(densidades_locales, shift)
    
    # Guardamos el perfil asociado a su temperatura
    perfiles_por_temperatura[T] = densidades_centradas

for T, perfil in perfiles_por_temperatura.items():
    if T > 1.10:
        color = 'darkgray'
    else:
        color = cmap(norm(T))
    plt.plot(centros_bines, perfil, label=f'T* = {T:.2f}', color=color, alpha=0.8, linewidth=1.0)

# Configuración estética del gráfico conjunto
plt.title(f'Evolución de los Perfiles de Densidad con la Temperatura (Ensayo {num_prueba})', fontsize=14, fontweight='bold')
plt.xlabel('Posición en el eje X', fontsize=12)
plt.ylabel('Densidad Reducida $\\rho^*$', fontsize=12)

# Añadimos una barra de colores para indicar la temperatura
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Temperatura Reducida $T^*$')

plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(bbox_to_anchor=(1.25, 1), loc='upper left', fontsize='small', ncol=1) # Leyenda fuera del gráfico
#plt.tight_layout()

# Guardamos el gráfico final que contiene todas las líneas
ruta_perfil_total = os.path.join(ruta_graficos, "Perfil_Densidad_Conjunto_Inestables_grises.png")
plt.savefig(ruta_perfil_total, dpi=300, bbox_inches='tight')
plt.show()
