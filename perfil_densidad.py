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

# -- SEÑALA EN NÚMERO DE PRUEBA/ENSAYO QUE DESEAS VISUALIZAR --
entrada = 1
num_prueba = int(entrada)
num_bines = 50 # Establece el número de fracciones en que se dividirá la caja de simulación (siempre en el eje x)


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

    plt.figure(figsize=(12, 6))
    

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
        plt.axvline(x=x_final, color='red', linewidth=1.0, alpha=0.3)

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

    
    # Graficamos las nuevas coordenadas x e y (Scatter Plot)
    plt.scatter(new_x, new_y, s=8, alpha=0.5, c='dodgerblue', edgecolors='none')
    
    # --- Verificación visual de los ejes ---
    # Dibujamos los ejes principales en negro para verificar el origen (0,0)
    plt.axhline(0, color='black', linewidth=1.5, linestyle='--')
    plt.axvline(0, color='black', linewidth=1.5, linestyle='--')
    plt.axhline(Ly, color='black', linewidth=1.5, linestyle='--')
    plt.axvline(Lx, color='black', linewidth=1.5, linestyle='--')

    
    # Configuración de etiquetas y título
    plt.title(f'Distribución de Partículas (Proyección XY)\nTemperatura = {T:.2f} (Prueba P{num_prueba})', fontsize=14, fontweight='bold')
    plt.xlabel('Coordenada X (Trasladada al origen)', fontsize=12)
    plt.ylabel('Coordenada Y (Trasladada al origen)', fontsize=12)
    
    # Ajustar límites para asegurar que se vea el origen (0,0) claramente
    L_max = max(new_x.max(), new_y.max()) # Aproximación del tamaño de la caja
    #plt.xlim(left=-0.05 * L_max, right=1.05 * L_max)
    #plt.ylim(bottom=-0.05 * L_max, top=1.05 * L_max)
    
    # Aspecto 1:1 para que la caja no se vea distorsionada (si la caja es cúbica)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    
    # Mostrar gráfico y pausar (esperar a que el usuario cierre la ventana antes de la siguiente T)
    # Si quieres que se muestren todas a la vez (no recomendado si son muchas T), quita el plt.show()
    # y usa plt.savefig() para guardarlas en archivos.
    #plt.show() 
    
    # Limpiar la figura actual de la memoria antes de pasar a la siguiente T
    #plt.close()

    # Definimos el nombre del archivo según la temperatura
    nombre_archivo = f"Scatter_T_{T:.2f}.png"
    ruta_final_archivo = os.path.join(ruta_graficos, nombre_archivo)
    
    # Guardamos la imagen con buena resolución (300 dpi)
    plt.savefig(ruta_final_archivo, dpi=300, bbox_inches='tight')
    print(f"Gráfico guardado en: {ruta_final_archivo}")
    
    # IMPORTANTE: Cerrar la figura para liberar memoria RAM
    plt.close()


# --- GRÁFICO DE COEXISTENCIA LÍQUIDO-VAPOR ---
plt.figure(figsize=(8, 6))

# Graficamos la rama del líquido
plt.plot(densidades_liquido, temperaturas, 'o', color='blue', label='Líquido Saturado')

# Graficamos la rama del vapor
plt.plot(densidades_vapor, temperaturas, 'o', color='red', label='Vapor Saturado')


# --- CÁLCULO DE ERRORES ABSOLUTOS MEDIOS (SIMULACIÓN vs NIST) ---

errores_l = []
errores_l_relativos = []
for mi_rho_l, nist_rho_l_val in zip(densidades_liquido, nist_rho_l):
    error = (mi_rho_l - nist_rho_l_val) 
    error_relativo = error/nist_rho_l_val * 100
    errores_l.append(abs(error))
    errores_l_relativos.append(abs(error_relativo))

error_liquido_promedio = sum(errores_l) / len(errores_l)
error_l_rel_promedio = sum(errores_l_relativos)/len(errores_l_relativos)


errores_v = []
errores_v_relativos = []
for mi_rho_v, nist_rho_v_val in zip(densidades_vapor, nist_rho_v):
    error = (mi_rho_v - nist_rho_v_val) 
    error_relativo = error/nist_rho_v_val * 100
    errores_v.append(abs(error))
    errores_v_relativos.append(abs(error_relativo))

error_vapor_promedio = sum(errores_v) / len(errores_v)
error_v_rel_promedio = sum(errores_v_relativos)/len(errores_v_relativos)


# --- DATOS DE LA NIST ---
plt.scatter(nist_rho_v, temperaturas, marker='x', color='darkred', label='NIST (Vapor)', zorder=5)
plt.scatter(nist_rho_l, temperaturas, marker='x', color='darkblue', label='NIST (Líquido)', zorder=5)

plt.title(f'Diagrama de Coexistencia $T$ vs $\\rho$ (Ensayo {num_prueba})')
plt.xlabel('Densidad Reducida $\\rho^*$')
plt.ylabel('Temperatura Reducida $T^*$')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

# --- Guardado de la Curva de Coexistencia ---
ruta_coexistencia = os.path.join(ruta_graficos, f'Coexistencia_P{num_prueba}_vs_NIST.png')
plt.savefig(ruta_coexistencia, dpi=300, bbox_inches='tight')
print(f"\n--- Diagrama de Coexistencia guardado en: {ruta_coexistencia} ---")


plt.show()

print('='*40)
print("--- VALIDACIÓN FINAL (Error Absoluto Medio) ---")
print(f"Error promedio en Fase Líquida: {error_liquido_promedio:.2f}")
print(f"Error promedio en Fase Vapor: {error_vapor_promedio:.2f}")
print('='*40)
print("--- VALIDACIÓN FINAL (Error Relativo Medio) ---")
print(f"Error relativo promedio en Líquido: {error_l_rel_promedio:.2f}%")
print(f"Error relativo promedio en Vapor: {error_v_rel_promedio:.2f}%")
print('='*40)

print("\n--- Procesamiento de todas las temperaturas finalizado ---")