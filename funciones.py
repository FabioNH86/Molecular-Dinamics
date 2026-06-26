import pandas as pd
import numpy as np
from numpy.polynomial import polynomial as P
import os
import matplotlib.pyplot as plt
import math
import hoomd
import gsd.hoomd
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist 


# from hoomd import deprecated


def calcular_dimensiones_por_densidad(ndiv, densidad_objetivo):
    """
    Calcula el lado L necesario para que un sistema con 
    N partículas (nx*ny*nz) tenga la densidad deseada.
    """
    n_particulas = ndiv[0] * ndiv[1] * ndiv[2]
    
    if n_particulas <= 0:
        raise ValueError("El número de partículas (ndiv) debe ser mayor a 0.")
    
    # rho = N / V  =>  V = N / rho
    volumen_necesario = n_particulas / densidad_objetivo
    
    # Para una caja cúbica: L = V^(1/3)
    # Para tu caja alargada (donde Ly=Lz y Lx = 2*Ly):
    # V = Lx * Ly * Lz = (2L) * L * L = 2L^3  => L = (V/2)^(1/3)
    
    # Asumamos caja cúbica por defecto para cálculos de capacidad calorífica
    l_lado = math.pow(volumen_necesario, 1/3)
    
    return l_lado, n_particulas


def actualizar_entradas(temp, densidad_obj=0.5, modo='barrido'):
    """
    modo 'barrido': N y V constantes (Geometría fija de coexistencia).
    modo 'isoterma': N constante, V variable (Ajusta L para la densidad).
    """
    
    if modo == 'isoterma':
        # Definimos un N estándar para capacidad calorífica (Caja Cúbica)
        # 12x12x12 es un buen equilibrio entre estadística y velocidad
        ndiv = [18, 17, 17] 
        
        # Calculamos el lado L necesario para la densidad_obj
        L, n_total = calcular_dimensiones_por_densidad(ndiv, densidad_obj)
        lx = ly = lz = L
        
    else: # modo == 'barrido'
        # Geometría fija para coexistencia (L y N constantes)
        lx, ly, lz = 48.92, 24.46, 24.46
        ndiv = [28, 14, 14]
        n_total = ndiv[0] * ndiv[1] * ndiv[2]

    # Formateo preciso para el in.dat
    ndiv_str = f"{ndiv[0]:<4} {ndiv[1]:<4} {ndiv[2]:<4}"
    l_str = f"{lx:<8.4f} {ly:<8.4f} {lz:<8.4f}"
    
    contenido = f"""2         ncolectivo_1_NVE_2_NVT_3_NPT
1         nopcion_1_inicializacion_2_continuacion
3         dofx_numero_de_dimensiones
12.0      expn_exponente_parte_repulsiva_12_para_LJ
6.0       expm_exponente_parte_atractiva_6_en_general
0.001     time_tiempo_de_integracion
{l_str}  Lx_Ly_Lz_unidades_reducidas
{ndiv_str}  ndivx_ndivy_ndivz_numero_de_atomos_por_lado
{temp:<10.2f}   temp_en_unidades_reducidas_para_asignar_velocidades
0.01      taut_para_termostato_berendsen
0.0       presion_en_unidades_reducidas_para_presostato
0.0       taup_para_presostato_berendsen
4.0       rcut_r_de_corte_del_potencial_menor_a_la_mitad_de_la_caja
1000000    nconfequi_numero_de_configuraciones_para_equilibrar
1500000   nconf_numero_de_configuraciones
50000     nperfil_frecuencia_para_calcular_distribuciones
0.1       deltar_ancho_del_intervalo_para_perfil_de_densidad
10000     nmovie_frecuencia_para_tomar_fotos
10000     nprint_frecuencia_para_imprimir_en_pantalla"""

    with open("in.dat", "w") as f:
        f.write(contenido.strip())

    print(f"--- Configuración {modo.upper()} lista ---")
    print(f"T={temp} | Rho={densidad_obj:.4f} | L={lx:.4f} | N={n_total}")


def calcular_densidades(filename, start_conf=500000, end_conf=1500000, num_atom=5488, num_bines=100):
    with open(filename, 'r') as f:
        # 1. Leer cabecera para obtener Lx y parámetros
        linea1 = f.readline() 
        num_atom_file = int(f.readline())

        for _ in range(num_atom_file): f.readline() # Se saltan las líneas de esos átomos 
        Lx, Ly, Lz = [float(v) / 0.3405 for v in f.readline().split()]

        print("Las dimensiones de la caja son:")                    
        print(f'Lx: {Lx:.4f}\nLy: {Ly:.4f}\nLz: {Lz:.4f} \n')

        volumen_total = Lx * Ly * Lz 
        densidad_global = num_atom / volumen_total
        print(f'El volumen total de la caja es: {volumen_total:.4f} y se dividirá en {num_bines} bines \n')
        print(f"La densidad global es: {densidad_global:.4f}")

        volumen_bin = (Lx / num_bines) * Ly * Lz
        print(f'El volumen de cada bin será: {volumen_bin:.4f} y su ancho será de: {(Lx / num_bines):.4f} \n')
        
        # Cortes a lo largo de Lx
        cortes_x = np.linspace(0, Lx, num_bines + 1)
        # Calcular centros de los bines (Eje X)
        centros_x = (cortes_x[:-1] + cortes_x[1:]) / 2

        # -- HASTA AQUÍ SE ESTABLECEN LAS CONDICIONES DE INICIO --
        # (Cuando la caja de simulación está acomodadita)

        # -- SALTO DE CONFIGURACIONES -- 
        print(f"Saltando hasta la configuración {start_conf}...")
        f.seek(0)
        n_conf = 0
        while n_conf < start_conf:
            linea = f.readline()
            if not linea: break # Fin si el archivo tiene un salto de línea inesperado
            if "CONFIGURACION" in linea:
                n_conf = int(linea.split()[-1])

                if n_conf < start_conf:
                    f.readline() # Saltamos num_atom
                    for _ in range(num_atom): f.readline()
                    f.readline() # Línea de posiciones

        # -- AQUI INICIA EL PROCESADO DE LAS DEMÁS POSICIONES --
        print(f"Iniciando procesado desde {n_conf} hasta {end_conf}...")
        
        todas_las_densidades = []
        while n_conf < end_conf:
            linea_conf = f.readline()
            if not linea_conf: break
            n_conf = int(linea_conf.split()[-1])

            f.readline()


            posiciones_x = []

            for _ in range(num_atom_file):
                linea_atomo = f.readline()

                posiciones_x.append(float(linea_atomo.split()[-3]))

            f.readline() # Saltamos las dimensiones de la caja 

            # Convertimos a unidades reducidas
            x_arr = np.array(posiciones_x) / 0.3405
            centro_actual = x_arr.mean()
            desplazamiento = (Lx / 2) - centro_actual
            x_centrado = (x_arr + desplazamiento) % Lx


            # CALCULAMOS LA DENSIDAD CON NUMPY 
            counts, _ = np.histogram(x_centrado, bins=cortes_x)
            densidades_por_bin = counts / volumen_bin
            todas_las_densidades.append(densidades_por_bin)


            if n_conf % 100000 == 0:
                print(f'Procesado: {n_conf}')

            
        # Conjunto de densidades
        df_densidades = pd.DataFrame(todas_las_densidades)

        #perfil_promedio = densidades_por_configuracion.mean(axis=0) 
        perfil_promedio = df_densidades.mean()
        desviacion_estandar = df_densidades.std()
        
    return centros_x, perfil_promedio, desviacion_estandar


def calcular_promedios_energía(archivo, ancho_bloques=100000, variacion_permitida=0.03):
    # Leemos el archivo 
    try:
        todo = pd.read_csv(archivo, header=None, sep=r'\s+')
        print(f"Trabajando en: {archivo}")
    except FileNotFoundError:
        print(f"No se encontró el archivo en: {archivo} \n")
        return None


    # Eliminamos la columna de la densidad
    energias = todo.drop(columns=[1])
    energias.columns = ['iconf', 'eki', 'epi', 'etot', 'tempi', 'presi', 'error']
    #print(energias.head())

    ultimo_etot_promedio = None 
    ultimo_eki_promedio = None 
    ultimo_epi_promedio = None


    print(f"{'Configuración':<13} | {'E_total':<15} | {'E_Cinetica':<15} | {'E_Potencial':<15}")
    print("-" * 65)

    try: # Se calcula el promedio cada ancho_bloques
        for i in range(0, len(energias), ancho_bloques):
            bloque = energias.iloc[i: i + ancho_bloques]

            # Si por alguna razón no se completa el último bloque
            if len(bloque) < ancho_bloques:
                break # Se evita su cálculo

            resumen_bloque = bloque.drop(columns=['iconf']).mean()


            # Se extrae el valor de la energía total    
            promedio_etot_actual = resumen_bloque['etot']
            promedio_eki_actual = resumen_bloque['eki']
            promedio_epi_actual = resumen_bloque['epi']
            configuracion_inicial = bloque['iconf'].iloc[0]

            if ultimo_etot_promedio is not None and ultimo_eki_promedio is not None and ultimo_epi_promedio is not None: 
                # Se calula la variación relativa
                variacion_etot = abs((promedio_etot_actual - ultimo_etot_promedio) / ultimo_etot_promedio)
                variacion_eki = abs((promedio_eki_actual - ultimo_eki_promedio) / ultimo_eki_promedio)
                variacion_epi = abs((promedio_epi_actual - ultimo_epi_promedio) / ultimo_epi_promedio)

                # Imprimimos el progreso
                print(f"{configuracion_inicial:<13.2f} | {promedio_etot_actual:<15.2f} | {promedio_eki_actual:<15.2f} | {promedio_epi_actual:<15.2f}")


                # if variacion_etot <= variacion_permitida and variacion_eki <= variacion_permitida and variacion_epi <= variacion_permitida:
                #     print(f"\n✅ ¡Estabilidad encontrada!")
                #     print(f"La configuración '{configuracion_inicial}' es la primera con variación < {variacion_permitida}")
                #     return configuracion_inicial # Retorna la iconf donde se cumple
                if (bloque_es_estable(bloque['etot']) and
                    bloque_es_estable(bloque['eki']) and
                    bloque_es_estable(bloque['epi']) and
                    variacion_etot <= variacion_permitida and 
                    variacion_eki <= variacion_permitida and 
                    variacion_epi <= variacion_permitida):
                    print(f"✅ Estabilidad encontrada en configuración {configuracion_inicial}")
                    return configuracion_inicial
        
            ultimo_etot_promedio = promedio_etot_actual
            ultimo_eki_promedio = promedio_eki_actual
            ultimo_epi_promedio = promedio_epi_actual


        print("\n❌ No se encontró ninguna configuración con esa estabilidad.")
        return None
        
    except ValueError:
        print("Error en el cálculo!")


def bloque_es_estable(valores, umbral_pendiente=0.0000001):
    x = np.arange(len(valores))

    # Ajuste lineal: coef[1] es la pendiente normalizada
    coef = np.polyfit(x, valores, 1)
    pendiente_normalizada = abs(coef[0] / np.mean(valores))
    return pendiente_normalizada < umbral_pendiente


def calcular_promedios_energía_claude(archivo, ancho_bloques=1000, variacion_permitida=0.03, bloques_consecutivos=5):
    try:
        todo = pd.read_csv(archivo, header=None, sep=r'\s+')
        print(f"Trabajando en: {archivo}")
    except FileNotFoundError:
        print(f"No se encontró el archivo en: {archivo}")
        return None

    energias = todo.drop(columns=[1])
    energias.columns = ['iconf', 'eki', 'epi', 'etot', 'tempi', 'presi', 'error']

    # Calcular promedios por bloque
    bloques = []
    for i in range(0, len(energias), ancho_bloques):
        bloque = energias.iloc[i: i + ancho_bloques]
        if len(bloque) < ancho_bloques:
            break
        resumen = bloque.drop(columns=['iconf']).mean()
        resumen['iconf_ini'] = bloque['iconf'].iloc[0]
        bloques.append(resumen)

    bloques_df = pd.DataFrame(bloques).reset_index(drop=True)

    print(f"{'Configuración':<13} | {'E_total':<15} | {'E_Cinetica':<15} | {'E_Potencial':<15}")
    print("-" * 65)

    # Buscar la primera configuración desde donde TODOS los bloques
    # siguientes son estables (variación < umbral entre bloques consecutivos)
    primer_bloque_estable = None

    for i in range(1, len(bloques_df)):
        prev = bloques_df.iloc[i - 1]
        curr = bloques_df.iloc[i]

        var_etot = abs((curr['etot'] - prev['etot']) / prev['etot']) if prev['etot'] != 0 else float('inf')
        var_eki  = abs((curr['eki']  - prev['eki'])  / prev['eki'])  if prev['eki']  != 0 else float('inf')
        var_epi  = abs((curr['epi']  - prev['epi'])  / prev['epi'])  if prev['epi']  != 0 else float('inf')

        estable = var_etot <= variacion_permitida and var_eki <= variacion_permitida and var_epi <= variacion_permitida

        marca = "✅" if estable else "❌"
        print(f"{curr['iconf_ini']:<13.0f} | {curr['etot']:<15.4f} | {curr['eki']:<15.4f} | {curr['epi']:<15.4f} {marca}")

    # Buscar desde qué índice TODOS los bloques restantes son estables
    for i in range(1, len(bloques_df) - bloques_consecutivos + 1):
        ventana = []
        for j in range(i, i + bloques_consecutivos):
            prev = bloques_df.iloc[j - 1]
            curr = bloques_df.iloc[j]

            var_etot = abs((curr['etot'] - prev['etot']) / prev['etot']) if prev['etot'] != 0 else float('inf')
            var_eki  = abs((curr['eki']  - prev['eki'])  / prev['eki'])  if prev['eki']  != 0 else float('inf')
            var_epi  = abs((curr['epi']  - prev['epi'])  / prev['epi'])  if prev['epi']  != 0 else float('inf')

            ventana.append(var_etot <= variacion_permitida and
                           var_eki  <= variacion_permitida and
                           var_epi  <= variacion_permitida)

        # Solo cuenta si además los bloques FINALES también son estables
        bloques_finales_estables = all(ventana) and all(
            abs((bloques_df.iloc[j]['etot'] - bloques_df.iloc[j-1]['etot']) / bloques_df.iloc[j-1]['etot']) <= variacion_permitida
            for j in range(i, len(bloques_df))
            if bloques_df.iloc[j-1]['etot'] != 0
        )

        if bloques_finales_estables:
            primer_bloque_estable = bloques_df.iloc[i]['iconf_ini']
            break

    if primer_bloque_estable:
        print(f"\n✅ Estabilidad sostenida desde la configuración: {primer_bloque_estable:.0f}")
        return primer_bloque_estable
    else:
        print("\n❌ No se encontró estabilidad sostenida hasta el final.")
        return None


def calcular_promedios_energía_claude_2(archivo, ancho_bloques=1000, variacion_permitida=0.03, fraccion_cola=0.2, mostrar_progreso=True):

    try:
        todo = pd.read_csv(archivo, header=None, sep=r'\s+')
    except FileNotFoundError:
        print(f"No se encontró el archivo en: {archivo}")
        return None

    energias = todo.drop(columns=[1])
    energias.columns = ['iconf', 'eki', 'epi', 'etot', 'tempi', 'presi', 'error']

    # Calcular promedios por bloque
    bloques = []
    for i in range(0, len(energias), ancho_bloques):
        bloque = energias.iloc[i: i + ancho_bloques]
        if len(bloque) < ancho_bloques:
            break
        resumen = bloque.drop(columns=['iconf']).mean()
        resumen['iconf_ini'] = bloque['iconf'].iloc[0]
        bloques.append(resumen)

    bloques_df = pd.DataFrame(bloques).reset_index(drop=True)

    # ✅ Referencia: promedio del último 20% de los datos (la "cola" equilibrada)
    n_cola = max(1, int(len(bloques_df) * fraccion_cola))
    referencia = bloques_df.tail(n_cola)[['etot', 'eki', 'epi']].mean()

    print(f"\nReferencia (promedio del último {fraccion_cola*100:.0f}%):")
    print(f"  E_total={referencia['etot']:.4f} | E_cin={referencia['eki']:.4f} | E_pot={referencia['epi']:.4f}")

    if mostrar_progreso:
        print(f"\n{'Configuración':<13} | {'E_total':<15} | {'E_Cinetica':<15} | {'E_Potencial':<15}")
        print("-" * 65)

    primer_bloque_estable = None

    for i, row in bloques_df.iterrows():
        # Variación respecto al promedio de la cola final
        var_etot = abs((row['etot'] - referencia['etot']) / referencia['etot'])
        var_eki  = abs((row['eki']  - referencia['eki'])  / referencia['eki'])
        var_epi  = abs((row['epi']  - referencia['epi'])  / referencia['epi'])

        estable = var_etot <= variacion_permitida and var_eki <= variacion_permitida and var_epi <= variacion_permitida
        marca = "✅" if estable else "❌"

        if mostrar_progreso:
            print(f"{row['iconf_ini']:<13.0f} | {row['etot']:<15.4f} | {row['eki']:<15.4f} | {row['epi']:<15.4f} {marca}")

        # Primera configuración estable desde la que TODOS los siguientes también lo son
        if estable and primer_bloque_estable is None:
            # Verificar que todos los bloques restantes también sean estables
            resto = bloques_df.iloc[i:]
            todos_estables = all(
                abs((r['etot'] - referencia['etot']) / referencia['etot']) <= variacion_permitida and
                abs((r['eki']  - referencia['eki'])  / referencia['eki'])  <= variacion_permitida and
                abs((r['epi']  - referencia['epi'])  / referencia['epi'])  <= variacion_permitida
                for _, r in resto.iterrows()
            )
            if todos_estables:
                primer_bloque_estable = row['iconf_ini']

    if primer_bloque_estable:
        print(f"\n✅ Sistema estabilizado desde la configuración: {primer_bloque_estable:.0f}")
    else:
        print(f"\n❌ No se encontró estabilidad sostenida. Considera aumentar variacion_permitida o ancho_bloques.")

    if primer_bloque_estable is not None:
        return int(primer_bloque_estable)


def calcular_presiones_vapor(archivo, configuraciones_consideradas=500000, eje_perpendicular_interface='x'): # configuraciones_consideradas: Se refieren a aquellas que ocurren después de alcanzado el equilibrio
    try:                                                              # Para conocer la configuración a la que se alcanzó el equilibrio, usar función calcular_promedios_energía
        presiones = pd.read_csv(archivo, header=None, sep=r'\s+')
        print(f"Trabajando presiones en: {archivo}")
    except FileNotFoundError:
        print(f"No se encontró el archivo en: {archivo} \n")
        return None
    
    presiones = presiones.drop(columns=[0])
    presiones.columns = ['x', 'y', 'z'] # Se asignan nómbres a los ejes

    # Filtrado 
    presiones = presiones.tail(configuraciones_consideradas)

    presión_vapor = presiones[eje_perpendicular_interface].mean()
    p_desv_estand = presiones[eje_perpendicular_interface].std()

    print(f'La presión de vapor promedio es: {presión_vapor:.4f}')
    print(f'Con una desviación estándar de: {p_desv_estand:.4f}')

    return presiones, presión_vapor, p_desv_estand


def calcular_presion_vapor_hoomd(archivo, Lx, eje_normal='x', configuraciones_consideradas=None):
    try:
        df = pd.read_csv(archivo)
        print(f"Trabajando en: {archivo}")
        print(f'Frames totales: {len(df)}')

    except FileNotFoundError:
        print(f"No se encontró el archivo: {archivo}")
        return None, None
    
    df = df.dropna() # Sacanos los frames que hayan tenido NAN

    if configuraciones_consideradas is not None:
        df = df.tail(configuraciones_consideradas) # Nos desplazamos a las configuraciones que vamos a muestrear

    print(f"-> Frames para promediar: {len(df)}")

    # --- Presiones normales y tangenciales según el eje normal a la interfaz ---
    if eje_normal == 'x':
        P_normal     = df['Pxx']
        P_tangencial = (df['Pyy'] + df['Pzz']) / 2.0
        Lx_col       = None  # Lx no está en el CSV, se pasa aparte si se necesita
    elif eje_normal == 'y':
        P_normal     = df['Pyy']
        P_tangencial = (df['Pxx'] + df['Pzz']) / 2.0
    elif eje_normal == 'z':
        P_normal     = df['Pzz']
        P_tangencial = (df['Pxx'] + df['Pyy']) / 2.0
    else:
        raise ValueError(f"eje_normal debe ser 'x', 'y' o 'z', no '{eje_normal}'")

    presion_vapor = P_normal.mean()
    presion_vapor_std = P_normal.std()

    # Calculamos la presión tangencial promedio
    presion_tang = P_tangencial.mean()
    presion_tang_std = P_tangencial.std()

    presion_total = (df['Pxx'] + df['Pyy'] + df['Pzz']).mean() / 3.0

    print(f"\n  P_normal  ({eje_normal})   : {presion_vapor:.6f} ± {presion_vapor_std:.6f}")
    print(f"  P_tangencial             : {presion_tang:.6f}  ± {presion_tang_std:.6f}")
    print(f"  P_isotropica (traza/3)   : {presion_total:.6f}")

    resultados = {
        'P_normal_mean' : presion_vapor,
        'P_normal_std'  : presion_vapor_std,
        'P_tangencial_mean': presion_tang,
        'P_tangencial_std': presion_tang_std,
        'P_isotropca'   : presion_total,
        'n_frames'      : len(df)
    }

    # Calculamos de una vez la tension superficial
    tension_serie = (Lx / 2.0) * (P_normal - P_tangencial)

    tension_superficial_media = tension_serie.mean()
    tension_superficial_std = tension_serie.std()

    return df, resultados, tension_superficial_media, tension_superficial_std




def graficar_evolucion_presion(df_presiones):
    # Creamos una figura con 3 subgráficos verticales
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    ejes = ['x', 'y', 'z']
    colores = ['#e74c3c', '#2ecc71', '#3498db'] # Rojo, Verde, Azul
    
    for i, eje in enumerate(ejes):
        axes[i].plot(df_presiones[eje], color=colores[i], linewidth=0.5, alpha=0.8)
        axes[i].set_ylabel(f'Presión {eje.upper()}')
        axes[i].grid(True, linestyle='--', alpha=0.6)
        axes[i].legend([f'Eje {eje}'], loc='upper right')

    axes[2].set_xlabel('Configuraciones')
    fig.suptitle('Evolución Temporal de la Presión por Eje', fontsize=16)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


def calcular_tension_superficial(df_presiones, longitud_partpendicular_interface):
    print(f'Se tienen {len(df_presiones)} frames')
    # 1. Calculamos la tensión instantánea completa
    P_tangencial_serie = 0.5 * (df_presiones['y'] + df_presiones['z'])
    tension_instantanea = (longitud_partpendicular_interface / 2) * (df_presiones['x'] - P_tangencial_serie)
    
    # 2. Recortamos la serie para ignorar el inicio usando iloc
    tension_estabilizada = tension_instantanea.iloc[:]
    
    # 3. Extraemos los estadísticos solo de la zona estable
    tension_superficial_media = tension_estabilizada.mean()
    tension_superficial_std   = tension_estabilizada.std()
    
    print(f"Frames usados para el promedio: {len(tension_estabilizada)}")
    print(f"La tensión superficial es: {tension_superficial_media:.4f} ± {tension_superficial_std:.4f}")
    
    return tension_superficial_media, tension_superficial_std


def generar_dataframes_todo(archivo):
    # Extraemos la información de cada configuración
    todo = pd.read_csv(archivo, sep=r'\s+', header=None)
    
    todo.columns = ['iconf','rho', 'eki', 'epi', 'etot', 'tempi', 'presi', 'error']
    
    return todo

    
def calcular_capacidad_calorífica(dataframe, T=0.7):
    kB = 1.0
    
    # Cv = (<U²> - <U>²) / kB * T
    eki = dataframe['eki']
    
    #  <U²>
    promedio_cuadrado = np.mean(np.square(eki))

    #  <U>²
    cuadrado_del_promedio = np.square(eki.mean())

    Cv_total = (promedio_cuadrado - cuadrado_del_promedio) / (kB * T**2)

    print(f'La capacidad calorifica es: {Cv_total}')
    return Cv_total


def calcular_cv_hoomd(df, T=0.7):
    kB = 1.0
    df.columns = [c.split('.')[-1] for c in df.columns]
    # print(f"Columnas detectadas: {df.columns.tolist()}")

    df_produccion = df.iloc[50:].copy()

    # print(df_produccion.head())

    # Cv = (<U²> - <U>²) / kB * T
    total_energy = df_produccion['kinetic_energy'] + df_produccion['potential_energy']
    
    # La parte de arriba del cálculo de la Cv es equivalente a calcular la variación de la energía interna
    Cv_total = total_energy.var() / (kB * T**2)

    print(f'La capacidad calorifica es: {Cv_total}')
    return Cv_total


def run_hoomd_simulation(temp, ruta_destino, length_minibox, equilibracion, muestreo, periodic_zeromomentum, modo='isoterma', rho=0.5, ndiv_entrada=[]):
    print(f'Iniciando simulación a T={temp:.2f}')
    # -- Uso de GPU -- 
    device = hoomd.device.GPU()
    sim = hoomd.Simulation(device=device, seed=42)

    # -- Definición de la geometría del sistema -- 
    if modo == 'isoterma':
        ndiv = [18, 18, 17] 
        n_total = ndiv[0] * ndiv[1] * ndiv[2]

        # L = (N / rho)^(1/3)
        L = (n_total / rho)**(1/3)
        lx = ly = lz = L

    else:
        lx, ly, lz = 48.92, 24.46, 24.46
        ndiv = ndiv_entrada
        n_total = ndiv[0] * ndiv[1] * ndiv[2]

        # Cálculo de espaciado #== Se centran las partículas en un rectángulo interior. ==
        ly_minibox = 24
        lz_minibox = 24

        dx = length_minibox / ndiv[0]
        dy = ly_minibox / ndiv[1]
        dz = lz_minibox / ndiv[2]
        
        # Agregando un espacio entre las partículas y las paredes de la caja
        offset_x = -length_minibox / 2
        offset_y = -ly_minibox / 2
        offset_z = -lz_minibox / 2


    # Se crea un estado incial
    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0: # Verifica que se esté en el primer proceso
        snap.configuration.box = [lx, ly, lz, 0, 0, 0] # Se asignan las dimensiones e inclinación a los lados de la caja de simulación
        snap.particles.N = n_total # Se asigna un número de partículas
        snap.particles.types = ['A']

        if modo == 'isoterma':
            # Acomodo de las partículas
            x = np.linspace(-lx/2, lx/2, ndiv[0], endpoint=False)
            y = np.linspace(-ly/2, ly/2, ndiv[1], endpoint=False)
            z = np.linspace(-lz/2, lz/2, ndiv[2], endpoint=False)
        
        else:
            pos = []
            for i in range(ndiv[0]):
                for j in range(ndiv[1]):
                    for k in range(ndiv[2]):
                        x = i * dx + offset_x + (dx / 2)
                        y = j * dy + offset_y + (dy / 2)
                        z = k * dz + offset_z + (dz / 2)
                        pos.append([x, y, z])

        # xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        # pos = np.stack((xx.flatten(), yy.flatten(), zz.flatten()), axis=-1)
            snap.particles.position[:] = pos # Se copian las coordenadas generadas en un snap de hoomd
            

    sim.create_state_from_snapshot(snap) # Se crea la simulación a partir de las posiciones del snap generado
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp) # Se asignan velocidades iniciales



    cell = hoomd.md.nlist.Cell(buffer=0.4) # Se espera que la partícula se mueva 0.4 unidades 
    # cell: Divide la caja en rejillas de interacción 
    lj = hoomd.md.pair.LJ(nlist=cell, default_r_cut=4.0, mode="shift") # Se define el potecial mie como en el código en C con r_cut = 4.0
    lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)

    # -- Integrador y termostato del ensamble NVT --
    # dt=0.001 y taut=0.01
    termostato = hoomd.md.methods.thermostats.Bussi(kT=temp, tau=0.01)
    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)
    integrator = hoomd.md.Integrator(dt=0.001, methods=[nvt], forces=[lj])
    sim.operations.integrator = integrator

    # -- Loggers (Muestra de la información) --
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # Logger de datos termodinámicos (parecido a todo.dat)
    logger = hoomd.logging.Logger(categories=['scalar'])

    logger.add(sim, quantities=['timestep'])
    logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure', 'pressure_tensor'])

    # El archivo de salida es una tabla
    if modo == 'isoterma':
        log_filename = f"todo_T{temp:.2f}_rho{rho:.4f}.csv"
        gsd_filename = f"trajectory_T{temp:.2f}_rho{rho:.4f}.gsd"
    else: 
        log_filename = os.path.join(ruta_destino, f"todo_T{temp:.2f}_long{length_minibox}.csv")
        gsd_filename = os.path.join(ruta_destino, f"trajectory_T{temp:.2f}_long{length_minibox}.gsd")


    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(10000),
                              logger=logger,
                              output=open(log_filename, 'w'))

    sim.operations.writers.append(table)


    # Para poder guardar la primer configuración en el gsd y ver cómo se acomodaron las partículas
    trigger_combinado = hoomd.trigger.Or([hoomd.trigger.On(1),
                                          hoomd.trigger.Periodic(50000)])
    
    gsd_writer = hoomd.write.GSD(trigger=trigger_combinado,
                                 filename=gsd_filename,
                                 mode='wb')
    
    sim.operations.writers.append(gsd_writer)

    terminal_log = hoomd.logging.Logger(categories=['scalar', 'string'])
    terminal_log.add(sim, quantities=['timestep', 'tps'])
    terminal_writter = hoomd.write.Table(trigger=hoomd.trigger.Periodic(5_000), logger=terminal_log)
    sim.operations.writers.append(terminal_writter)
    
    zero_momentum = hoomd.md.update.ZeroMomentum(trigger=hoomd.trigger.Periodic(periodic_zeromomentum))
    
    sim.operations.updaters.append(zero_momentum)
    

    # Fase de suavizado inicial
    integrator.dt = 0.00001
    sim.run(5000)

    integrator.dt = 0.0001
    sim.run(5000)

    integrator.dt = 0.001 

    sim.run(equilibracion)
    sim.run(muestreo)

    print(f'Simulación finalizada: Rho = {n_total / (lx * ly * lz):.4f} | L = {lx:.4f} | N = {n_total}')


def calcular_perfil_densidad_gsd(gsd_file, start_frame=0, num_bines=100):
    """
    Calcula el perfil de densidad a lo largo del eje X usando el archivo GSD de HOOMD.
    
    Parameters:
    -----------
    gsd_file : str
        Ruta al archivo .gsd
    start_frame : int
        Índice del frame inicial (0 es el primero grabado, no el paso de simulación).
    num_bines : int
        Número de divisiones en el eje X.
    """
    with gsd.hoomd.open(name=gsd_file, mode='r') as trayecto:
        # 1. Obtener datos del primer frame para inicializar
        snap_ref = trayecto[0]
        lx = snap_ref.configuration.box[0]
        ly = snap_ref.configuration.box[1]
        lz = snap_ref.configuration.box[2]
        num_atom = snap_ref.particles.N
        
        print(f"Dimensiones de la caja: Lx={lx:.4f}, Ly={ly:.4f}, Lz={lz:.4f}")
        print(f'Con un total de {num_atom} partículas')
        
        volumen_bin = (lx / num_bines) * ly * lz
        cortes_x = np.linspace(-lx/2, lx/2, num_bines + 1)
        centros_x = (cortes_x[:-1] + cortes_x[1:]) / 2
        
        todas_las_densidades = []
        
        # 2. Iterar sobre los frames (HOOMD usa índices de 0 a N_frames)
        total_frames = len(trayecto)
        print(f"Procesando {total_frames - start_frame} frames de {total_frames} totales...")

        for i in range(start_frame, total_frames):
            snap = trayecto[i]
            
            # HOOMD centra la caja en (0,0,0), las posiciones van de -L/2 a L/2
            posiciones_x = snap.particles.position[:, 0]
            
            # Centrado opcional: En simulaciones LV, la interfase puede moverse.
            # Este paso centra la masa en el origen de la caja.
            centro_masa = np.mean(posiciones_x)
            x_centrado = ((posiciones_x - centro_masa + lx/2) % lx) - lx/2

            # Histograma para obtener densidades
            counts, _ = np.histogram(x_centrado, bins=cortes_x)
            densidades_por_bin = counts / volumen_bin
            todas_las_densidades.append(densidades_por_bin)

            if i % 10 == 0:
                print(f"Frame procesado: {i}/{total_frames}", end="\r")

        # 3. Procesamiento estadístico
        df_densidades = pd.DataFrame(todas_las_densidades)
        perfil_promedio = df_densidades.mean()
        desviacion_estandar = df_densidades.std()
        
        print("\n✅ Cálculo de perfil completado.")
        
    return centros_x, perfil_promedio, desviacion_estandar


def calcular_perfil_densidad_multi_especie(gsd_file, tipos_interes, start_frame=0, num_bines=50):
    """
    Calcula el perfil de densidad a lo largo del eje X para múltiples especies.
    
    Parameters:
    -----------
    gsd_file : str
        Ruta al archivo .gsd
    tipos_interes : list of str
        Lista con los nombres de los tipos de partículas a separar (ej. ['A', 'B'] o ['solvente', 'polimero']).
    start_frame : int
        Índice del frame inicial.
    num_bines : int
        Número de divisiones en el eje X.
    """
    with gsd.hoomd.open(name=gsd_file, mode='r') as trayecto:
        snap_ref = trayecto[0]
        segundo_frame = trayecto[1]
        print(f'Cada frame tiene: {segundo_frame.configuration.step} pasos')        
        lx, ly, lz = snap_ref.configuration.box[0:3]
        print(f'La longitud de la caja es: {lx}')
        
        # Mapear los nombres de los tipos ('solvente', 'polimero') a sus IDs numéricos (0, 1, etc.)
        nombres_tipos = snap_ref.particles.types
        mapa_tipos = {nombre: i for i, nombre in enumerate(nombres_tipos)}
        
        # Validar que los tipos pedidos existan en el GSD
        for tipo in tipos_interes:
            if tipo not in mapa_tipos:
                raise ValueError(f"El tipo '{tipo}' no se encuentra en el archivo GSD. Tipos disponibles: {nombres_tipos}")
        
        volumen_bin = (lx / num_bines) * ly * lz
        cortes_x = np.linspace(-lx/2, lx/2, num_bines + 1)
        centros_x = (cortes_x[:-1] + cortes_x[1:]) / 2
        
        # Inicializar diccionarios para guardar los datos de cada frame por especie
        historial_densidades = {tipo: [] for tipo in tipos_interes}
        
        total_frames = len(trayecto)
        print(f"Procesando {total_frames - start_frame} frames para las especies: {tipos_interes}...")

        for i in range(start_frame, total_frames):
            snap = trayecto[i]
            
            # 1. Centrado usando el centro de masa de TODAS las partículas (para mantener el mismo sistema de referencia)
            pos_x_todas = snap.particles.position[:, 0]
            centro_masa = np.mean(pos_x_todas)
            x_centrado_todas = ((pos_x_todas - centro_masa + lx/2) % lx) - lx/2
            
            # 2. Separar por especie usando el typeid
            typeids = snap.particles.typeid
            
            for tipo in tipos_interes:
                id_numerico = mapa_tipos[tipo]
                # Filtro: solo las posiciones X de la especie actual
                x_especie = x_centrado_todas[typeids == id_numerico]
                
                # Histograma para esta especie
                counts, _ = np.histogram(x_especie, bins=cortes_x)
                densidades_por_bin = counts / volumen_bin
                historial_densidades[tipo].append(densidades_por_bin)

            if i % 10 == 0:
                print(f"Frame procesado: {i}/{total_frames}", end="\r")

        # 3. Procesamiento estadístico por especie
        resultados = {}
        for tipo in tipos_interes:
            df_densidades = pd.DataFrame(historial_densidades[tipo])
            resultados[tipo] = {
                'promedio': df_densidades.mean(),
                'desviacion': df_densidades.std()
            }
        
        print("\n✅ Cálculo completado con éxito.")
        # print("RRESULTADOS")
        # print("-" * 40)
        # for tipo in tipos_interes:
        #     # Calculamos la densidad media total promediando todos los bines
        #     densidad_liquido = resultados[tipo]['promedio'].max()
        #     densidad_vapor = resultados[tipo]['promedio'].min()
        #     print(f"  -> Especie '{tipo}': Densidad de líquido = {densidad_liquido:.5f}")
        #     print(f"  -> Especie '{tipo}': Densidad de vapor = {densidad_vapor:.5f}")
        # print("-" * 40)
    
    return centros_x, resultados


def encontrar_equilibrio_hoomd(archivo_csv, pasos_totales, ancho_bloques=10, variacion_permitida=0.03, fraccion_cola=0.2, mostrar_progreso=True):
    """
    Analiza el archivo CSV de HOOMD para encontrar el paso donde las energías se estabilizan.
    """
    try:
        # 1. Leer CSV. HOOMD a veces usa espacios como delimitadores o comas.
        # skipinitialspace=True ayuda con los espacios después de las comas.
        df = pd.read_csv(archivo_csv, skipinitialspace=True)
        
        # 2. LIMPIEZA DE NOTACIÓN CIENTÍFICA:
        # Eliminamos espacios en los nombres de las columnas
        df.columns = df.columns.str.strip()
        
        # Convertimos todo a numérico. 
        # 'coerce' transformará cualquier cosa que no entienda en NaN.
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Eliminamos filas con NaN (posibles encabezados repetidos o basura)
        df = df.dropna().reset_index(drop=True)
        
    except Exception as e:
        print(f"❌ Error al procesar el archivo: {e}")
        return None

    try:
        # HOOMD genera CSV con encabezados, así que es más fácil de leer
        df = pd.read_csv(archivo_csv)
    except Exception as e:
        print(f"❌ Error al leer {archivo_csv}: {e}")
        return None

    cols_step_found = [c for c in df.columns if 'step' in c.lower() or 'timestep' in c.lower()]

    if cols_step_found:
        col_step = cols_step_found[0]
        df['step_axis'] = df[col_step]
    
    else:
        print(f"⚠️ Columna 'step' no encontrada. Generando eje basado en {pasos_totales} pasos.")
        df['step_axis'] = np.linspace(0, pasos_totales, len(df))

    try:
        col_pot = [c for c in df.columns if 'potential_energy' in c.lower()][0]
        col_kin = [c for c in df.columns if 'kinetic_energy' in c.lower()][0]
        
    except IndexError:
        print("❌ No se encontraron columnas de energía potencial o cinética.")
        return None
    
    # Creamos la columna de Energía Total si no existe
    df['etot'] = df[col_pot] + df[col_kin]
    
    # Renombramos para compatibilidad
    energias = df[['step_axis', col_kin, col_pot, 'etot']].copy()
    energias.columns = ['step', 'eki', 'epi', 'etot']

    # Calcular promedios por bloque
    bloques = []
    for i in range(0, len(energias), ancho_bloques):
        bloque = energias.iloc[i : i + ancho_bloques]
        if len(bloque) < ancho_bloques:
            break
        bloques.append(bloque.mean())

    bloques_df = pd.DataFrame(bloques)

    # ✅ Referencia: promedio del último % de los datos
    n_cola = max(1, int(len(bloques_df) * fraccion_cola))
    referencia = bloques_df.tail(n_cola)[['etot', 'eki', 'epi']].mean()

    print(f"\n--- Análisis de Equilibrio: {os.path.basename(archivo_csv)} ---")
    print(f"Referencia final: Etot={referencia['etot']:.4f} | Ekin={referencia['eki']:.4f}")

    if mostrar_progreso:
        print(f"{'Paso':<12} | {'E_total':<12} | {'E_Pot':<12} | Estabilidad")
        print("-" * 55)

    primer_paso_estable = None

    for i, row in bloques_df.iterrows():
        # Cálculo de variaciones
        var_etot = abs((row['etot'] - referencia['etot']) / referencia['etot'])
        var_eki  = abs((row['eki']  - referencia['eki'])  / referencia['eki'])
        
        # Criterio de estabilidad
        estable = var_etot <= variacion_permitida and var_eki <= variacion_permitida
        
        if mostrar_progreso and i % 5 == 0: # Mostrar cada 5 bloques para no saturar
            marca = "✅" if estable else "❌"
            print(f"{row['step']:<12.0f} | {row['etot']:<12.4f} | {row['epi']:<12.4f} | {marca}")

        # Verificación de estabilidad sostenida
        if estable and primer_paso_estable is None:
            resto = bloques_df.iloc[i:]
            # Todos los bloques siguientes deben cumplir el criterio
            todos_estables = all(
                abs((r['etot'] - referencia['etot']) / referencia['etot']) <= variacion_permitida 
                for _, r in resto.iterrows()
            )
            if todos_estables:
                primer_paso_estable = row['step']

    if primer_paso_estable:
        print(f"✅ Sistema estabilizado en el paso: {primer_paso_estable:.0f}")
        return int(primer_paso_estable)
    else:
        print("❌ No se detectó estabilidad sostenida.")
        return None


def leer_csv_seguro(archivo_csv):
    try:
        # sep=None con engine='python' detecta automáticamente si es coma, espacio o tab
        # comment='#' ignora las líneas de metadatos de HOOMD
        df = pd.read_csv(archivo_csv, sep=r'\s+', engine='python')
        df.columns.str.strip()

        # Limpiar nombres de columnas (quitar espacios y posibles '#' que queden)
        #df.columns = [c.strip().replace('#', '').strip() for c in df.columns]
        
        # Convertir a numérico solo lo que sea necesario
        df = df.apply(pd.to_numeric, errors='coerce')
        df = df.dropna(how='all').reset_index(drop=True) # Solo borra si TODA la fila es NaN
        
        return df
    except Exception as e:
        print(f"❌ Error real en lectura: {e}")
        return None
    

def encontrar_equilibrio_hoomd(archivo_csv, pasos_totales=1000000, ancho_bloques=500, variacion_permitida=0.03, fraccion_cola=0.2, mostrar_progreso=True):
    try:
        # 1. Leer con manejo de espacios
        df = pd.read_csv(archivo_csv, sep=r'\s+', engine='python')
        df.columns = df.columns.str.strip() # Limpiar nombres de columnas
        
        # 2. Conversión forzada a numérico de todo el DataFrame
        df = df.apply(pd.to_numeric, errors='coerce')
        df = df.dropna().reset_index(drop=True)
        
        if df.empty:
            print(f"❌ El archivo {archivo_csv} quedó vacío tras limpiar datos no numéricos.")
            return None

    except Exception as e:
        print(f"❌ Error crítico leyendo el archivo: {e}")
        return None

    # --- Búsqueda de Columnas ---
    cols_step = [c for c in df.columns if 'step' in c.lower() or 'timestep' in c.lower()]
    cols_pot  = [c for c in df.columns if 'potential_energy' in c.lower() or 'pe' == c.lower()]
    cols_kin  = [c for c in df.columns if 'kinetic_energy' in c.lower() or 'ke' == c.lower()]

    if not cols_pot or not cols_kin:
        print(f"❌ No se encontraron energías. Columnas detectadas: {list(df.columns)}")
        return None

    # Asignar valores
    step_vals = df[cols_step[0]] if cols_step else np.linspace(0, pasos_totales, len(df))
    epi_vals = df[cols_pot[0]]
    eki_vals = df[cols_kin[0]]
    etot_vals = epi_vals + eki_vals

    # Crear DataFrame de trabajo limpio y garantizado
    energias = pd.DataFrame({
        'step': step_vals,
        'eki': eki_vals,
        'epi': epi_vals,
        'etot': etot_vals
    })

    # --- Cálculo de Bloques (Forma robusta para evitar pérdida de nombres) ---
    indices_bloques = np.arange(len(energias)) // ancho_bloques
    bloques_df = energias.groupby(indices_bloques).mean()
    
    # Filtrar el último bloque si quedó incompleto
    if len(energias) % ancho_bloques != 0:
        bloques_df = bloques_df.iloc[:-1]

    if bloques_df.empty:
        print("❌ Datos insuficientes para crear bloques de promedio.")
        return None

    # --- Referencia y Estabilidad ---
    n_cola = max(1, int(len(bloques_df) * fraccion_cola))
    # Forzamos que sea una serie con los nombres correctos
    referencia = bloques_df.tail(n_cola).mean()

    print(f"\n--- Análisis: {os.path.basename(archivo_csv)} ---")
    # Usamos .get() o acceso directo ahora que garantizamos la creación arriba
    print(f"Referencia final: Etot={referencia['etot']:.4e} | Ekin={referencia['eki']:.4e}")

    if mostrar_progreso:
        print(f"{'Paso':<12} | {'E_total':<12} | {'E_Pot':<12} | Estabilidad")
        print("-" * 60)

    primer_paso_estable = None
    
    # Valores de referencia para evitar KeyErrors y divisiones por cero
    ref_etot = referencia['etot'] if referencia['etot'] != 0 else 1e-9
    ref_eki  = referencia['eki']  if referencia['eki']  != 0 else 1e-9

    for i, row in bloques_df.iterrows():
        var_etot = abs((row['etot'] - ref_etot) / ref_etot)
        var_eki  = abs((row['eki']  - ref_eki)  / ref_eki)
        
        estable = var_etot <= variacion_permitida and var_eki <= variacion_permitida
        
        if mostrar_progreso and i % 5 == 0:
            marca = "✅" if estable else "❌"
            print(f"{row['step']:<12.0f} | {row['etot']:<12.4e} | {row['epi']:<12.4e} | {marca}")

        if estable and primer_paso_estable is None:
            # Comprobar si el resto de la simulación sigue siendo estable
            resto = bloques_df.iloc[i:]
            if all((abs((r['etot'] - ref_etot) / ref_etot) <= variacion_permitida) for _, r in resto.iterrows()):
                primer_paso_estable = row['step']

    if primer_paso_estable is not None:
        print(f"✅ Equilibrio detectado en el paso: {primer_paso_estable:.0f}")
        return int(primer_paso_estable)
    
    print("❌ No se detectó una región de estabilidad sostenida.")
    return None


def graficar_analisis_termo(archivo_csv, ruta_graficos, paso_equilibrio, pasos_totales=1500000):
    """
    Visualiza las energías y presión, marcando el punto donde el sistema se estabilizó.
    """
    try:
        # 1. Carga con el separador correcto corregido (\s+)
        df = pd.read_csv(archivo_csv, sep=r'\s+', engine='python')
        df.columns = df.columns.str.strip()
        df = df.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)
        
        # 2. Identificación flexible e inmune a los nombres largos de HOOMD
        cols_step = [c for c in df.columns if 'step' in c.lower()]
        
        # Si no existe la columna timestep en el CSV, calculamos el paso real 
        # multiplicando el índice de la fila por el intervalo de guardado (10000)
        if cols_step:
            step = df[cols_step[0]]
        else:
            print("⚠️ Columna 'timestep' no encontrada en el CSV. Reconstruyendo eje X cada 10k pasos.")
            step = df.index * 10000  # <--- Esto mapea el eje X con tus datos reales de 10k en 10k
        
        # Buscar las columnas por coincidencia parcial de texto de manera limpia
        col_pe = [c for c in df.columns if 'potential' in c.lower()][0]
        col_ke = [c for c in df.columns if 'kinetic_energy' in c.lower() or 'kinetic_' in c.lower() and 'temperature' not in c.lower()][0]
        col_pres = [c for c in df.columns if 'pressure' in c.lower()]
        
        pe = df[col_pe]
        ke = df[col_ke]
        etot = pe + ke
        press = df[col_pres[0]] if col_pres else np.zeros(len(df))

    except Exception as e:
        print(f"❌ Error al procesar datos para gráfica: {e}")
        return

    # --- CONFIGURACIÓN DE LA FIGURA ---
    fig, axs = plt.subplots(2, 2, figsize=(14, 9), layout='constrained')
    t_titulo = os.path.basename(archivo_csv).replace('.csv', '')
    fig.suptitle(f'Análisis de Termodinámica: {t_titulo}', fontsize=16, fontweight='bold')

    # Diccionario para iterar los paneles fácilmente
    datasets = [
        (ke, 'Energía Cinética', 'tab:blue', axs[0, 0]),
        (pe, 'Energía Potencial', 'tab:red', axs[0, 1]),
        (etot, 'Energía Total', 'tab:purple', axs[1, 0]),
        (press, 'Presión', 'tab:green', axs[1, 1])
    ]

    for data, title, color, ax in datasets:
        # Graficar toda la serie
        ax.plot(step, data, color=color, alpha=0.3, label='Equilibración')
        
        # Si hay un paso de equilibrio, resaltar la zona de "producción"
        if paso_equilibrio:
            mask_estable = step >= paso_equilibrio
            ax.plot(step[mask_estable], data[mask_estable], color=color, linewidth=1.5, label='Producción (Estable)')
            ax.axvline(x=paso_equilibrio, color='black', linestyle='--', alpha=0.6, label='Punto Eq.')
        
        ax.set_title(title, fontweight='semibold')
        ax.set_xlabel('Timestep')
        ax.grid(True, alpha=0.3)
        if title == 'Energía Total': ax.legend(loc='best', fontsize='small')

    # Guardar resultado
    # (Asegúrate de que ruta_graficos esté definida en tu script principal)
    plt.savefig(os.path.join(ruta_graficos, f"analisis_{t_titulo}.png"), dpi=300)
    plt.show()


def visualizar_estabilidad_dinamica(archivo_csv, paso_eq, T, ruta_guardado, pasos_totales=1500000):
    """
    Replica tu estilo de visualización pero con detección de equilibrio automática.
    """
    df = leer_csv_seguro(archivo_csv)
    if df is None or df.empty:
        print(f"⚠️ No se pudo graficar {archivo_csv}")
        return

    # Identificar columnas (ajusta según tus nombres reales en el CSV)
    c_step = [c for c in df.columns if 'step' in c.lower()][0] if any('step' in c.lower() for c in df.columns) else None
    c_pot  = [c for c in df.columns if 'potential' in c.lower() or 'pe' == c.lower()][0]
    c_kin  = [c for c in df.columns if 'kinetic' in c.lower() or 'ke' == c.lower()][0]
    c_pres = [c for c in df.columns if 'pressure' in c.lower()]
    
    # Eje X
    x = df[c_step] if c_step else np.linspace(0, pasos_totales, len(df))
    df['etot'] = df[c_pot] + df[c_kin]

    # --- FIGURA ---
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    fig.suptitle(f'Análisis de Estabilidad (Temperatura T = {T})', fontsize=16)

    # Configuración de los 4 paneles
    config_paneles = [
        (axs[0, 0], df[c_kin], 'Energía Cinética', 'tab:blue'),
        (axs[0, 1], df[c_pot], 'Energía Potencial', 'tab:red'),
        (axs[1, 0], df['etot'], 'Energía Total', 'tab:purple'),
        (axs[1, 1], df[c_pres[0]] if c_pres else np.zeros(len(df)), 'Presión', 'tab:green')
    ]

    for ax, data, titulo, color in config_paneles:
        # Graficar toda la simulación en sombra
        ax.plot(x, data, color=color, alpha=0.3, label='Equilibración')
        
        if paso_eq is not None:
            # Resaltar la zona estable
            mask = x >= paso_eq
            ax.plot(x[mask], data[mask], color=color, alpha=1, label='Producción')
            ax.axvline(x=paso_eq, color='black', linestyle='--', linewidth=1)
            
        ax.set_title(titulo)
        ax.set_xlabel('Paso (Step)')
        ax.grid(True, linestyle=':', alpha=0.6)

    # Añadir leyenda solo al primer panel para no amontonar
    axs[0, 0].legend(fontsize='small')

    # Guardar y mostrar
    nombre_fig = f"Estabilidad_T{T:.2f}.png"
    plt.savefig(os.path.join(ruta_guardado, nombre_fig), dpi=300)
    plt.show()


def run_sim_binary_sistem(temp, equilibracion, muestreo, eps_AB=1.0, sist_homegeno=True):
    # --- Identificador para archivos ---
    suffix = "Homog" if sist_homegeno else "Separado"
    # Incluimos eps_AB en el nombre para no sobreescribir archivos
    file_id = f"{suffix}_T{temp:.2f}_epsAB{eps_AB:.2f}"
    
    print(f'>>> Ejecutando: {file_id}...')
    # -- Uso de GPU -- 
    device = hoomd.device.GPU()
    sim = hoomd.Simulation(device=device, seed=42)


    lx, ly, lz = 100.0, 25.0, 25.0
    ndiv = [64, 32, 25]
    n_total = ndiv[0] * ndiv[1] * ndiv[2]

    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        snap.configuration.box = [lx, ly, lz, 0, 0, 0]
        snap.particles.N = n_total
        snap.particles.types = ['A', 'B']
        snap.particles.mass[:] = [1.0] * n_total

    # Acomodo de las partículas a los extremos de la caja 
    x = np.linspace(-lx/2 + 0.5, lx/2 - 0.5, ndiv[0])
    y = np.linspace(-ly/2 + 0.5, ly/2 - 0.5, ndiv[1])
    z = np.linspace(-lz/2 + 0.5, lz/2 - 0.5, ndiv[2])

    mesh = np.array(np.meshgrid(x, y, z, indexing='ij')).reshape(3, -1).T
    snap.particles.position[:] = mesh

    if sist_homegeno:
        type_ids = np.zeros(n_total, dtype=int)
        type_ids[n_total // 2:] = 1
        np.random.shuffle(type_ids)
        snap.particles.typeid[:] = type_ids
    
    else:
        # CREAR TIPOS A PARTIR DE LAS POSICIONES
        snap.particles.typeid[:] = (mesh[:, 0] >= 0).astype(int)

    sim.create_state_from_snapshot(snap)
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp)

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    mie = hoomd.md.pair.Mie(nlist=cell, default_r_cut=4.0)
    mie.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0, n=12, m=6)
    mie.params[('B', 'B')] = dict(epsilon=1.0, sigma=1.0, n=12, m=6)

    # Modificación de la interacción cruzada
    mie.params[('A', 'B')] = dict(epsilon=eps_AB, sigma=1.0, n=12, m=6)

    # -- Integrador y termostato del ensamble NVT --
    # dt=0.001 y taut=0.01
    termostato = hoomd.md.methods.thermostats.Bussi(kT=temp, tau=0.01)
    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)
    integrator = hoomd.md.Integrator(dt=0.001, methods=[nvt], forces=[mie])
    sim.operations.integrator = integrator

    # -- Loggers (Muestra de la información) --
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # Logger de datos termodinámicos (parecido a todo.dat)
    logger = hoomd.logging.Logger(categories=['scalar'])
    logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure'])

    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(10000),
                              logger=logger,
                              output=open(f"log_{file_id}.csv", 'w'))
    
    sim.operations.writers.append(table)

    # Para poder guardar la primer configuración en el gsd y ver cómo se acomodaron las partículas
    trigger_combinado = hoomd.trigger.Or([hoomd.trigger.On(1),
                                          hoomd.trigger.Periodic(50000)])
    
    gsd_writer = hoomd.write.GSD(trigger=trigger_combinado,
                                 filename=f"traj_{file_id}.gsd",
                                 mode='wb')
    
 
    sim.operations.writers.append(gsd_writer)

    sim.run(equilibracion)
    sim.run(muestreo)

    print(f'Simulación finalizada ✅\n')


def crear_primer_frame(densidad_goticula, aspect_ratio, concentracion_porcentual_monomeros, grado_polimerizacion, tipos_solventes=2, n_monomeros=12_000):   
    snap = hoomd.Snapshot()
    snap.configuration.step = 1

    # Calculamos cuántos monómeros y polímeros necesitamos
    n_solvente_aprox      = int(((100/concentracion_porcentual_monomeros) - 1) * n_monomeros)
    n_total               = n_solvente_aprox + n_monomeros  
    n_polimeros           = n_monomeros // grado_polimerizacion
    n_enlaces             = n_polimeros * (grado_polimerizacion - 1) 
    n_monomeros_reales    = n_polimeros * grado_polimerizacion
    n_solvente_total      = n_total - n_monomeros_reales

    # Agregamos algún tipo de solvente diferente
    n_solvente_1          = int(n_solvente_total / tipos_solventes)
    n_solvente_2          = n_solvente_total - n_solvente_1

    print(f"Cantidad monomeros: {n_monomeros}")
    print(f"Cantidad part. solvente: {n_solvente_total}")
    print(f"({int(tipos_solventes)} solventes, en proporciones iguales)")
    print(f"Concentración {(n_monomeros / n_total * 100):.2f}")
   
    # Comenzamos con una caja cúbica en el centro del sistema
    volumen_gotícula = n_total / densidad_goticula

    print(f"El volumen de la caja es: {volumen_gotícula:.2f}")

    # Parámetro de red 
    # El espaciado entre partículas es el parámetro de red
    parametro_red = (1 / densidad_goticula) ** (1/3) 

    L = volumen_gotícula ** (1/3) # L: Longitude la minibox de gotícula

    print(f"El parametro de red es: {parametro_red:.2f}")
    
    # Usamos el propio parámetro de red para separar un poco la minibox de la caja de simulación 
    lx = aspect_ratio * L + 2 * parametro_red
    ly = L + 2 * parametro_red
    lz = L + 2 * parametro_red

    # Configuramos la caja
    snap.configuration.box = [lx, ly, lz, 0, 0, 0] # Se genera la caja alargada
    print(f"Lx={lx:.2f}, Ly={ly:.2f}, Lz={lz:.2f}")

    snap.particles.N = n_total
    snap.particles.types = ['S1', 'S2', 'A', 'B', 'C'] # Solvente (S1, S2), Polímero (A, B, C)
    snap.particles.mass[:] = [1.0] * n_total
    # Polímero bond
    snap.bonds.N = n_enlaces
    snap.bonds.types = ['A-A', 'A-B', 'B-B', 'B-C', 'C-C']

    # Cuántas partículas caben en cada eje?
    # Tenemos un n_total
    n_p_eje = int(np.ceil(n_total ** (1/3))) # Al ser un cúbo debe de haber los mismos en cada eje 

    print(f"El total de partículas es: {n_total}")
    print(f"Caben {n_p_eje} partículas por eje")
   
    # Considerando que todo está centrado en el origen
    coord_x = np.linspace(start=-L/2, stop=L/2, num=(n_p_eje + 4)) # Tratamos de poner las partículas perdidas elongando el eje x
    coord_yz = np.linspace(start=-L/2, stop=L/2, num=n_p_eje)
    
    Px, Py, Pz = np.meshgrid(coord_x, coord_yz, coord_yz, indexing='xy')
    # Aplanamos las posiciones y transponemos las posiciones
    posiciones = np.vstack([Px.ravel(), Py.ravel(), Pz.ravel()]).T

    # Checamos que las posiciones estén del lado izquierdo o derecho
    lado_izquierdo = posiciones[:, 0] < 0 
    lado_derecho = posiciones[:, 0] >= 0

    # Asignamos cada grupo a un lado 
    posiciones_izq = posiciones[lado_izquierdo]
    posiciones_der = posiciones[lado_derecho]

    # Metemos aleatoriedad
    np.random.seed(42)
    np.random.shuffle(posiciones_izq) # Tomamos diferentes posiciones al azar
    np.random.shuffle(posiciones_der)

    # Divimos la cantidad de partículas entre la cantidad de solvente
    n_lado = n_total // 2

    particulas_izq = posiciones_izq[:n_lado]
    particulas_der = posiciones_der[:n_lado]

    if n_total % 2 != 0:
        particulas_der = posiciones_der[:n_lado + 1]


    # Jusntamos ambos lados
    posiciones = np.vstack([particulas_izq, particulas_der])

    print(f"Las posiciones totales generadas son: {len(posiciones)} (Izq: {len(particulas_izq)}, Der: {len(particulas_der)})")
    
    # Asignamos las posiciones al snap
    snap.particles.position[:] = posiciones

    # AHORA REVISARÉMOS LAS QUE HAYAN QUEDADO CERCA PARA ENLAZARLAS Y CONSTRUIR POLÍMEROS 
    arbo_espacial = KDTree(posiciones) # Generamos el arbol

    # Generamos los identificadores de cada tipo de partícula
    id_S1 = snap.particles.types.index('S1') # Solvente tipo 1
    id_S2 = snap.particles.types.index('S2') # Solvente tipo 2
    id_A  = snap.particles.types.index('A')
    id_B  = snap.particles.types.index('B')
    id_C  = snap.particles.types.index('C')

    # Creamos un diccionario con los tipos de enlaces para los monomeros del polimero
    mapa_enlaces = {
        ('A', 'A'): snap.bonds.types.index('A-A'), 
        ('A', 'B'): snap.bonds.types.index('A-B'),
        ('B', 'B'): snap.bonds.types.index('B-B'), 
        ('B', 'C'): snap.bonds.types.index('B-C'),
        ('C', 'C'): snap.bonds.types.index('C-C')  
    }

    tipos_iniciales = np.empty(n_total)

    # Como 'posiciones' se construyó apilando [particulas_izq, particulas_der]:
    tipos_iniciales[:n_lado] = id_S1  # La mitad izquierda nace siendo Solvente 1
    tipos_iniciales[n_lado:] = id_S2  # La mitad derecha nace siendo Solvente 2

    snap.particles.typeid[:] = tipos_iniciales # Asignamos todas por defecto a Solvente tipo 1

    # Necesitamos llevar control de partículas ya enlazadas
    particulas_usadas = set()
    enlaces = []
    tipos_enlaces = []

    # Debemos dividir el polímero en el número de bloques de cada tipo de monómero
    bloque_monom = grado_polimerizacion // 3
    primer_seccion = bloque_monom
    segunda_seccion = 2*bloque_monom
    # tercer_seccion = 3*bloque_monom

    # Ahora recorremos todas las partículas hasta alcanzar la cantidad de polímeros que deberíamos de tener
    for _ in range(n_polimeros):
        # Buscamos una semilla inicia (partícula que inicia el polímero)
        semilla = False 
        for particula in range(n_total): # Recorremos todas las partículas
            if particula not in particulas_usadas:
                semilla = particula # Encontramos una partícula para convertila en semilla 
                break
            
        if semilla is None:
            break # En caso de no encontrar partícula

        tipo_actual = 'A' # El polímero siempre comienza con el tipo A
        # Convertimos esta semilla en la primer partícula para el polímero 
        snap.particles.typeid[semilla] = id_A
        particulas_usadas.add(semilla) # La agregamos las que ya están usadas 

        # Comenzamos a elongar el polímero 
        particula_actual = semilla
        for i in range(grado_polimerizacion - 1):
            # Recorremos el polímero asignando un tipo de monómero dependiendo de su posición
            # dependiento de su unicación dentro del polímero
            if i < primer_seccion:
                tipo_siguiente = 'A'
            elif i < segunda_seccion:
                tipo_siguiente = 'B'
            else:
                tipo_siguiente = 'C'


            pos_actual = posiciones[particula_actual] # Tomamos la posicion de la partícula 

            # Buscamos las partículas más cercanas a esta y sacamos su distancia
            distancia, indices_vecinos = arbo_espacial.query(pos_actual, k=15) # Buscamos 15 partículas cercanas

            particula_siguiente = None # Comenzaremos a construir los enlaces 
            for vecina in indices_vecinos:
                if vecina != particula_actual and vecina not in particulas_usadas:
                    particula_siguiente = vecina # Verificamos que la vecina no sea la misma partícula no se encuentre usada
                    break

            if particula_siguiente is None:
                break # Evitamos romper el código si se encierra espacialmente

            # Agregamos el enlace 
            enlaces.append([particula_actual, particula_siguiente])

            id_enlace = mapa_enlaces.get((tipo_actual, tipo_siguiente), 0)
            tipos_enlaces.append(id_enlace)

            # Convertimos las vecinos a Polímero 
            # Agregamos ya con los tipos asignados
            snap.particles.typeid[particula_siguiente] = snap.particles.types.index(tipo_siguiente)
            particulas_usadas.add(particula_siguiente)

            # Avanzamos a la siguiente partícula 
            particula_actual = particula_siguiente
            tipo_actual = tipo_siguiente

    # Guardamos los enlaces generados 
    snap.bonds.group[:] = enlaces
    snap.bonds.typeid[:] = tipos_enlaces

    # Impresión de control para verificar que las proporciones sean correctas
    print(f"\nPolímeros generados: {n_polimeros} de tamaño {grado_polimerizacion}")
    print("--- Conteo Final de Partículas en el Sistema ---")
    for t in snap.particles.types:
        conteo = list(snap.particles.typeid).count(snap.particles.types.index(t))
        print(f"  Tipo {t}: {conteo}")
  
    return snap


def correr_simulacion(snapshot, temp, equilibracion, muestreo, mon_cadena, aspect_ratio, eps_SP=1.0):
    # --- Identificador para archivos ---
    file_id = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}"

    # Inicializar HOOMD con el snapshot
    gpu = hoomd.device.CPU()
    sim = hoomd.Simulation(device=gpu, seed=42)
    sim.create_state_from_snapshot(snapshot)
    print(f"La conf. inicial de la caja es: \n{snapshot.configuration.box}")
    
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp)

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    mie = hoomd.md.pair.Mie(nlist=cell, default_r_cut=4.0, mode='shift')

    mie.params[('S1', 'S1')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)
    mie.params[('S2', 'S2')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)
    mie.params[('S1', 'S2')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)

    mie.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)
    mie.params[('B', 'B')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)
    mie.params[('C', 'C')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)

    mie.params[('A', 'B')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)
    mie.params[('B', 'C')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)
    mie.params[('A', 'C')] = dict(epsilon=1.0, sigma=1.0, n=12.0, m=6.0)

    mie.params[('S1', 'A')] = dict(epsilon=eps_SP, sigma=1.0, n=12.0, m=6.0)
    mie.params[('S1', 'B')] = dict(epsilon=eps_SP, sigma=1.0, n=12.0, m=6.0)
    mie.params[('S1', 'C')] = dict(epsilon=eps_SP, sigma=1.0, n=12.0, m=6.0)

    mie.params[('S2', 'A')] = dict(epsilon=eps_SP, sigma=1.0, n=12.0, m=6.0)
    mie.params[('S2', 'B')] = dict(epsilon=eps_SP, sigma=1.0, n=12.0, m=6.0)
    mie.params[('S2', 'C')] = dict(epsilon=eps_SP, sigma=1.0, n=12.0, m=6.0)

    # Fuerza de enlace para mantener la integridad de los polímeros
    armonico = hoomd.md.bond.Harmonic()
    armonico.params['A-A', 'A-B', 'B-B', 'B-C', 'C-C'] = dict(k=10.0, r0=1.0)
    

    # -- Integrador y termostato del ensamble NVT --
    termostato = hoomd.md.methods.thermostats.Bussi(kT=temp, tau=0.01)
    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)


    integrator = hoomd.md.Integrator(dt=0.005, methods=[nvt], forces=[mie, armonico])
    sim.operations.integrator = integrator

    # -- Loggers (Muestra de la información) --
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # Logger de datos termodinámicos (parecido a todo.dat)
    logger = hoomd.logging.Logger(categories=['scalar'])
    logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure'])  



    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(1000),
                              logger=logger,
                              output=open(f"log_{file_id}_monom_{mon_cadena}.dat", 'w'))       
    
    sim.operations.writers.append(table)

    # Para poder guardar la primer configuración en el gsd y ver cómo se acomodaron las partículas
    trigger_combinado = hoomd.trigger.Or([hoomd.trigger.On(0),
                                          hoomd.trigger.On(1),
                                          hoomd.trigger.On(2),
                                          hoomd.trigger.Periodic(2)])       
    
    gsd_writer = hoomd.write.GSD(trigger=trigger_combinado,
                                 filename=f"{file_id}_monom_{mon_cadena}.gsd",
                                 mode='wb') 
    
    sim.operations.writers.append(gsd_writer)

    sim.run(equilibracion)
    box = snapshot.configuration.box # Accedo a la configuración de la caja actual

    new_box = snapshot.configuration.box = [aspect_ratio * box[0], box[1], box[2], 0.0, 0.0, 0.0] # Altero la configuración de la caja en el muestreo

    sim.state.set_box(box=new_box)
    print(f"Se estriró la caja a: \n{snapshot.configuration.box}")

    # energy_operation = hoomd.update.CustomUpdater(trigger=2)
    # sim.operations += energy_operation



    sim.run(muestreo)


def crear_primer_frame_homopolimero(densidad_goticula, aspect_ratio, concentracion_porcentual_monomeros, monomeros_en_polimero, n_monomeros=12_000):   
    snap = hoomd.Snapshot()
    snap.configuration.step = 1

    # Calculamos cuántos monómeros y polímeros necesitamos
    n_solvente_aprox      = int(((100/concentracion_porcentual_monomeros) - 1) * n_monomeros)
    n_total               = n_solvente_aprox + n_monomeros  
    n_polimeros           = n_monomeros // monomeros_en_polimero
    n_enlaces             = n_polimeros * (monomeros_en_polimero - 1) 
    n_monomeros_reales    = n_polimeros * monomeros_en_polimero
    n_solvente            = n_total - n_monomeros_reales

    print(f"Cantidad monomeros: {n_monomeros}")
    print(f"Cantidad part. solvente: {n_solvente}")
    print(f"Concentración {(n_monomeros / n_total * 100):.2f}")
   
    # Comenzamos con una caja cúbica en el centro del sistema
    volumen_gotícula = n_total / densidad_goticula

    print(f"El volumen de la caja es: {volumen_gotícula:.2f}")

    # Parámetro de red 
    # El espaciado entre partículas es el parámetro de red
    parametro_red = (1 / densidad_goticula) ** (1/3) 

    L = volumen_gotícula ** (1/3) # L: Longitude la minibox de gotícula

    print(f"El parametro de red es: {parametro_red:.2f}")
    
    # Usamos el propio parámetro de red para separar un poco la minibox de la caja de simulación 
    lx = aspect_ratio * L + 2 * parametro_red
    ly = L + 2 * parametro_red
    lz = L + 2 * parametro_red

    # Configuramos la caja
    snap.configuration.box = [lx, ly, lz, 0, 0, 0] # Se genera la caja alargada
    print(f"Lx={lx:.2f}, Ly={ly:.2f}, Lz={lz:.2f}")

    snap.particles.N = n_total
    snap.particles.types = ['S', 'P'] # Solvente, Polímero
    snap.particles.mass[:] = [1.0] * n_total
    # Polímero bond
    snap.bonds.N = n_enlaces
    snap.bonds.types = ['P-P']

    # Cuántas partículas caben en cada eje?
    # Tenemos un n_total
    n_p_eje = int(np.ceil(n_total ** (1/3))) # Al ser un cúbo debe de haber los mismos en cada eje 

    print(f"El total de partículas es: {n_total}")
    print(f"Caben {n_p_eje} partículas por eje")
   
    # Considerando que todo está centrado en el origen
    coord_x = np.linspace(start=-L/2, stop=L/2, num=(n_p_eje)) # Tratamos de poner las partículas perdidas elongando el eje x
    coord_yz = np.linspace(start=-L/2, stop=L/2, num=n_p_eje)
    
    Px, Py, Pz = np.meshgrid(coord_x, coord_yz, coord_yz, indexing='xy')
    # Aplanamos las posiciones y transponemos las posiciones
    posiciones = np.vstack([Px.ravel(), Py.ravel(), Pz.ravel()]).T

    # Metemos aleatoriedad
    np.random.seed(42)
    np.random.shuffle(posiciones) # Tomamos diferentes posiciones al azar

    # Recortamos el exceso de posiciones para tener exactamente n_total
    posiciones = posiciones[:n_total]

    print(f"Las posiciones son: {len(posiciones)}")

    # Asignamos las posiciones al snap
    snap.particles.position[:] = posiciones

    # AHORA REVISARÉMOS LAS QUE HAYAN QUEDADO CERCA PARA ENLAZARLAS Y CONSTRUIR POLÍMEROS 
    arbo_espacial = KDTree(posiciones) # Generamos el arbol

    # Generamos los identificadores de cada tipo de partícula
    id_S = snap.particles.types.index('S') # Solvente
    id_P = snap.particles.types.index('P') # Polímero (Monómero del polímero técnicamente)
    snap.particles.typeid[:] = id_S # Asignamos todas por defecto a Solvente

    # Necesitamos llevar control de partículas ya enlazadas
    particulas_usadas = set()
    enlaces = []

    # Ahora recorremos todas las partículas hasta alcanzar la cantidad de polímeros que deberíamos de tener
    for _ in range(n_polimeros):
        # Buscamos una semilla inicia (partícula que inicia el polímero)
        semilla = False 
        for particula in range(n_total): # Recorremos todas las partículas
            if particula not in particulas_usadas:
                semilla = particula # Encontramos una partícula para convertila en semilla 
                break
            
        if semilla is None:
            break # En caso de no encontrar partícula

        # Convertimos esta semilla en la primer partícula para el polímero 
        snap.particles.typeid[semilla] = id_P
        particulas_usadas.add(semilla) # La agregamos las que ya están usadas 

        # Comenzamos a elongar el polímero 
        particula_actual = semilla
        for _ in range(monomeros_en_polimero - 1):
            pos_actual = posiciones[particula_actual] # Tomamos la posicion de la partícula 

            # Buscamos las partículas más cercanas a esta y sacamos su distancia
            distancia, indices_vecinos = arbo_espacial.query(pos_actual, k=15) # Buscamos 15 partículas cercanas

            particula_siguiente = None # Comenzaremos a construir los enlaces 
            for vecina in indices_vecinos:
                if vecina != particula_actual and vecina not in particulas_usadas:
                    particula_siguiente = vecina # Verificamos que la vecina no sea la misma partícula no se encuentre usada
                    break

            # Agregamos el enlace 
            enlaces.append([particula_actual, particula_siguiente])

            # Convertimos las vecinos a Polímero 
            snap.particles.typeid[particula_siguiente] = id_P
            particulas_usadas.add(particula_siguiente)

            # Avanzamos a la siguiente partícula 
            particula_actual = particula_siguiente

    # Guardamos los enlaces gee=nerados 
    snap.bonds.group[:] = enlaces
    snap.bonds.typeid[:] = [0] * len(enlaces) # Todos pertenecen al tipo P-P

    print(f"Polímeros generados: {n_polimeros} de tamaño {monomeros_en_polimero}")
    print(f"Total de monómeros 'P' asignados: {list(snap.particles.typeid).count(id_P)}")
    print(f"Total de solventes 'S' restantes: {list(snap.particles.typeid).count(id_S)}")
  
    return snap


def correr_simulacion_homoplimero(snapshot, temp, equilibracion, muestreo, mon_cadena, aspect_ratio, eps_SP=1.0):
    # --- Identificador para archivos ---
    file_id = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}"

    # Inicializar HOOMD con el snapshot
    gpu = hoomd.device.GPU()
    sim = hoomd.Simulation(device=gpu, seed=42)
    sim.create_state_from_snapshot(snapshot)
    
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp)

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    lj = hoomd.md.pair.LJ(nlist=cell, default_r_cut=4.0, mode='shift')

    lj.params[('S', 'S')] = dict(epsilon=1.0, sigma=1.0)
    lj.params[('P', 'P')] = dict(epsilon=1.0, sigma=1.0)
    lj.params[('S', 'P')] = dict(epsilon=eps_SP, sigma=1.0)


    # Fuerza de enlace para mantener la integridad de los polímeros
    armonico = hoomd.md.bond.Harmonic()
    armonico.params['P-P'] = dict(k=10.0, r0=1.0)
    

    # -- Integrador y termostato del ensamble NVT --
    termostato = hoomd.md.methods.thermostats.MTTK(kT=temp, tau=0.01)
    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)


    integrator = hoomd.md.Integrator(dt=0.005, methods=[nvt], forces=[lj, armonico])
    sim.operations.integrator = integrator

    # -- Loggers (Muestra de la información) --
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # Logger de datos termodinámicos (parecido a todo.dat)
    logger = hoomd.logging.Logger(categories=['scalar'])
    logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure'])    

    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(10000),
                              logger=logger,
                              output=open(f"log_{file_id}_monom_{mon_cadena}.csv", 'w'))       
    
    sim.operations.writers.append(table)

    # Para poder guardar la primer configuración en el gsd y ver cómo se acomodaron las partículas
    trigger_combinado = hoomd.trigger.Or([hoomd.trigger.On(0),
                                          hoomd.trigger.On(1),
                                          hoomd.trigger.Periodic(5000)])       
    
    gsd_writer = hoomd.write.GSD(trigger=trigger_combinado,
                                 filename=f"{file_id}_monom_{mon_cadena}.gsd",
                                 mode='wb') 
    

    term_log = hoomd.logging.Logger(categories=['scalar', 'string'])
    term_log.add(sim, quantities=['timestep', 'tps'])
    term_writter = hoomd.write.Table(trigger=hoomd.trigger.Periodic(5000), logger=term_log)
    sim.operations.writers.append(term_writter)

    
    sim.operations.writers.append(gsd_writer)

    sim.run(equilibracion)

    box = snapshot.configuration.box # Accedo a la configuración de la caja actual

    new_box = snapshot.configuration.box = [aspect_ratio * box[0], box[1], box[2], 0.0, 0.0, 0.0] # Altero la configuración de la caja en el muestreo

    sim.state.set_box(box=new_box)
    print(f"Se estriró la caja a: \n{snapshot.configuration.box}")

    sim.run(muestreo)


def continue_sim_from_gsd(archivo_gsd, muestreo, temp, eps_SP, mon_cadena, aspect_ratio):
    file_id = f"Poly-Solv_T{temp:.2f}_epsSP{eps_SP:.2f}"

    # Inicializar HOOMD con el snapshot
    gpu = hoomd.device.GPU()
    sim = hoomd.Simulation(device=gpu, seed=42)
    sim.create_state_from_gsd(archivo_gsd, frame=-1)
    
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp)
    
    original_box = sim.state.box


# Revisamos si la caja ya viene estirada o no 
    if original_box.Lx > 100.0:
        print(f"📦 El GSD ya cuenta con la caja estirada ([{original_box.Lx:.2f}, {original_box.Ly:.2f}, {original_box.Lz:.2f}]). Reanudando directo...")
    else:
        print("📦 El GSD tiene la caja chica del equilibrio. Aplicando estiramiento controlado...")
        new_Lx = original_box.Lx * aspect_ratio
        new_Ly = original_box.Ly + 0.1
        new_Lz = original_box.Lz + 0.1

        new_box = hoomd.Box(
            Lx=new_Lx,
            Ly=new_Ly,
            Lz=new_Lz,
            xy=original_box.xy,
            xz=original_box.xz,
            yz=original_box.yz
        )

        sim.state.set_box(box=new_box)
        print(f"Se estriró la caja a: \n{[new_box.Lx, new_box.Ly, new_box.Lz]}")

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    lj = hoomd.md.pair.LJ(nlist=cell, default_r_cut=4.0, mode='shift')

    lj.params[('S', 'S')] = dict(epsilon=1.0, sigma=1.0)
    lj.params[('P', 'P')] = dict(epsilon=1.0, sigma=1.0)
    lj.params[('S', 'P')] = dict(epsilon=eps_SP, sigma=1.0)

    # Fuerza de enlace para mantener la integridad de los polímeros
    armonico = hoomd.md.bond.Harmonic()
    armonico.params['P-P'] = dict(k=10.0, r0=1.0)
    

    # -- Integrador y termostato del ensamble NVT --
    # termostato = hoomd.md.methods.thermostats.Bussi(kT=temp, tau=0.01)
    termostato = hoomd.md.methods.thermostats.MTTK(kT=temp, tau=0.2)

    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)


    integrator = hoomd.md.Integrator(dt=0.005, methods=[nvt], forces=[lj, armonico])
    sim.operations.integrator = integrator

    # -- Loggers (Muestra de la información) --
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # Logger de datos termodinámicos (parecido a todo.dat)
    logger = hoomd.logging.Logger(categories=['scalar'])
    logger.add(thermo, quantities=['potential_energy', 'kinetic_energy', 'kinetic_temperature', 'pressure'])    

    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(5000),
                              logger=logger,
                              output=open(f"log_{file_id}_monom_{mon_cadena}.csv", mode='a'))       
    
    sim.operations.writers.append(table)

    # Send the logger outputs to the terminal
    term_log = hoomd.logging.Logger(categories=['scalar', 'string'])
    term_log.add(sim, quantities=['timestep', 'tps'])
    term_writer = hoomd.write.Table(trigger=hoomd.trigger.Periodic(5000), logger=term_log)
    sim.operations.writers.append(term_writer)


    # Para poder guardar la primer configuración en el gsd y ver cómo se acomodaron las partículas
    trigger_combinado = hoomd.trigger.Or([hoomd.trigger.On(0),
                                          hoomd.trigger.On(1),
                                          hoomd.trigger.Periodic(5000)])       
    
    gsd_writer = hoomd.write.GSD(trigger=trigger_combinado,
                                 filename=f"{file_id}_monom_{mon_cadena}.gsd",
                                 mode='ab') 
    
    sim.operations.writers.append(gsd_writer)

    # sim.run(equilibracion)

    # sim.run(muestreo)
    pasos_restantes = muestreo - sim.timestep
    
    if pasos_restantes > 0:
        sim.run(pasos_restantes)
    else:
        print(f"Aviso: La simulación ya está en el paso {sim.timestep}, no se avanzó.")


def calcular_radio_giro_promedio(trayectoria, longitud_cadenas, dimensiones_caja):
    """
    Calcula el radio de giro promedio (Rg) usando la fórmula de distancias por pares.

    """
    num_frames, num_monomeros, _ = trayectoria.shape
    # Determainamos cuántos polímeros hay dentro de la simulación
    num_cadenas = num_monomeros // longitud_cadenas 
    
    rg_por_frame = []

    for f in range(num_frames):
        # Obtenemos las coordenadas de los monómeros en cada frame
        coordenadas_frame = trayectoria[f]

        # Extraemos las dimensiones de la caja de simulación
        L = dimensiones_caja[f]

        # Calculamos el radio de giro para cada cadena
        for c in range(num_cadenas):
            inicio_cadena = c * longitud_cadenas
            fin_cadena = inicio_cadena + longitud_cadenas
            coordenadas_cadena = coordenadas_frame[inicio_cadena:fin_cadena]

            suma_distancias_cuadrado = 0.0
            for i in range(longitud_cadenas):
                for j in range(i+1, longitud_cadenas):
                    # (Ri - Rj)
                    dr = coordenadas_cadena[i] - coordenadas_cadena[j]

                    # Corregimos la diferencia de distancias si la cadena cruza la caja de simulación
                    dr = dr - L * np.round(dr/L)
                    suma_distancias_cuadrado += np.sum(dr**2)

            rg_cuadrado = suma_distancias_cuadrado / (longitud_cadenas**2)
            # Guardamos los Rg en una lista
            rg_por_frame.append(np.sqrt(rg_cuadrado))
   
    return np.mean(rg_por_frame), np.std(rg_por_frame), rg_por_frame

def procesar_datos_termo(archivo_csv):
    """
    Lee un archivo CSV de HOOMD de forma robusta y extrae el timestep y la energía potencial.
    Devuelve un diccionario con los vectores listos para graficar y el nombre corto del archivo.
    """
    try:
        # 1. Carga con el separador correcto
        df = pd.read_csv(archivo_csv, sep=r'\s+', engine='python')
        df.columns = df.columns.str.strip()
        df = df.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)
        
        # 2. Identificación del Timestep
        cols_step = [c for c in df.columns if 'step' in c.lower()]
        if cols_step:
            step = df[cols_step[0]].to_numpy()
        else:
            # Reconstrucción si no existe la columna
            step = (df.index * 10000).to_numpy()
        
        # 3. Identificación de la Energía Potencial
        col_pe = [c for c in df.columns if 'potential' in c.lower()][0]
        pe = df[col_pe].to_numpy()
        
        nombre_corto = os.path.basename(archivo_csv).replace('.csv', '')
        
        return {
            'nombre': nombre_corto,
            'step': step,
            'pe': pe
        }
        
    except Exception as e:
        print(f"❌ Error al procesar {os.path.basename(archivo_csv)}: {e}")
        return None
    
def crear_primer_frame_solvente(densidad, aspect_ratio, n_particulas=240_000):
    snap = hoomd.Snapshot()
    snap.configuration.step = 1

    # Volumen y parámetro de red
    volumen = n_particulas / densidad
    parametro_red = (1 / densidad) ** (1/3)
    L = volumen ** (1/3)

    print(f"Cantidad de partículas solvente: {n_particulas}")
    print(f"Volumen de la caja: {volumen:.2f}")
    print(f"Parámetro de red: {parametro_red:.2f}")

    lx = aspect_ratio * L + 2 * parametro_red
    ly = L + 2 * parametro_red
    lz = L + 2 * parametro_red

    snap.configuration.box = [lx, ly, lz, 0, 0, 0]
    print(f"Lx={lx:.2f}, Ly={ly:.2f}, Lz={lz:.2f}")

    snap.particles.N = n_particulas
    snap.particles.types = ['S']
    snap.particles.mass[:] = [1.0] * n_particulas
    snap.bonds.N = 0
    snap.bonds.types = []

    # Grid cúbico de posiciones
    n_p_eje = int(np.ceil(n_particulas ** (1/3)))
    print(f"Partículas por eje: {n_p_eje}")

    coord_x  = np.linspace(-L/2, L/2, n_p_eje)
    coord_yz = np.linspace(-L/2, L/2, n_p_eje)

    Px, Py, Pz = np.meshgrid(coord_x, coord_yz, coord_yz, indexing='xy')
    posiciones = np.vstack([Px.ravel(), Py.ravel(), Pz.ravel()]).T

    np.random.seed(42)
    np.random.shuffle(posiciones)
    posiciones = posiciones[:n_particulas]

    snap.particles.position[:] = posiciones
    snap.particles.typeid[:] = [0] * n_particulas  # Todas son 'S'

    print(f"Posiciones asignadas: {len(posiciones)}")
    return snap


def correr_simulacion_solvente(snapshot, temp, equilibracion, muestreo, aspect_ratio):
    file_id = f"Solv_T{temp:.2f}"

    gpu = hoomd.device.GPU()
    sim = hoomd.Simulation(device=gpu, seed=42)
    sim.create_state_from_snapshot(snapshot)

    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp)

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    lj = hoomd.md.pair.LJ(nlist=cell, default_r_cut=4.0, mode='shift')
    lj.params[('S', 'S')] = dict(epsilon=1.0, sigma=1.0)

    termostato = hoomd.md.methods.thermostats.MTTK(kT=temp, tau=0.01)
    nvt = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(), thermostat=termostato)

    integrator = hoomd.md.Integrator(dt=0.005, methods=[nvt], forces=[lj])
    sim.operations.integrator = integrator

    # --- Cómputos termodinámicos ---
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    # --- Logger: escalares + tensor de presión ---
    # El tensor de presión es categoría 'sequence', necesita logger separado
    logger_scalar = hoomd.logging.Logger(categories=['scalar'])
    logger_scalar.add(thermo, quantities=[
        'potential_energy',
        'kinetic_energy',
        'kinetic_temperature',
        'pressure',
    ])

    logger_tensor = hoomd.logging.Logger(categories=['sequence'])
    logger_tensor.add(thermo, quantities=['pressure_tensor'])

    # Escritura de escalares en CSV
    table_scalar = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(5000),
        logger=logger_scalar,
        output=open(f"log_{file_id}.csv", 'w'),
    )
    sim.operations.writers.append(table_scalar)

    # Escritura del tensor de presión en archivo separado
    table_tensor = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(10_000),
        logger=logger_tensor,
        output=open(f"log_{file_id}_pressure_tensor.csv", 'w'),
    )
    sim.operations.writers.append(table_tensor)

    # --- GSD ---
    trigger_combinado = hoomd.trigger.Or([
        hoomd.trigger.On(0),
        hoomd.trigger.On(1),
        hoomd.trigger.Periodic(5000),
    ])
    gsd_writer = hoomd.write.GSD(
        trigger=trigger_combinado,
        filename=f"{file_id}.gsd",
        mode='wb',
    )
    sim.operations.writers.append(gsd_writer)

    # --- Logger de progreso en terminal ---
    term_log = hoomd.logging.Logger(categories=['scalar', 'string'])
    term_log.add(sim, quantities=['timestep', 'tps'])
    term_writer = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(5000),
        logger=term_log,
    )
    sim.operations.writers.append(term_writer)

    # --- Equilibración ---
    sim.run(equilibracion)

    # --- Estiramiento de caja para etapa de muestreo ---
    box = sim.state.box  # Usamos sim.state.box, no el snapshot (ya puede estar desactualizado)
    new_box = hoomd.Box(
        Lx=aspect_ratio * box.Lx,
        Ly=box.Ly,
        Lz=box.Lz,
    )
    sim.state.set_box(box=new_box)
    print(f"Caja estirada a: Lx={new_box.Lx:.2f}, Ly={new_box.Ly:.2f}, Lz={new_box.Lz:.2f}")

    # --- Muestreo ---
    sim.run(muestreo)