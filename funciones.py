import pandas as pd
import numpy as np
from numpy.polynomial import polynomial as P
import io
import matplotlib.pyplot as plt
import math



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
        ndiv = [18, 18, 17]
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
500000    nconfequi_numero_de_configuraciones_para_equilibrar
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

    return presiones, presión_vapor


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
    presion_x = df_presiones['x'].mean()
    presion_y = df_presiones['y'].mean()
    presion_z = df_presiones['z'].mean()

    factor_presiones = presion_x - 0.5 * (presion_y + presion_y)
    tension_superficial = (longitud_partpendicular_interface / 2) * factor_presiones

    print(f'La tensión superficial es: {tension_superficial:.4f}')


def generar_dataframes_todo(archivo):
    # Extraemos la información de cada configuración
    todo = pd.read_csv(archivo, sep='\s+', header=None)
    
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
