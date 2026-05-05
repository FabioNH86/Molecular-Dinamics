import pandas as pd
import numpy as np
from numpy.polynomial import polynomial as P
import os
import matplotlib.pyplot as plt
import math
import hoomd
import gsd.hoomd



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
    print(f"Columnas detectadas: {df.columns.tolist()}")

    df_produccion = df.iloc[50:].copy()

    print(df_produccion.head())

    # Cv = (<U²> - <U>²) / kB * T
    total_energy = df_produccion['kinetic_energy'] + df_produccion['potential_energy']
    
    # La parte de arriba del cálculo de la Cv es equivalente a calcular la variación de la energía interna
    Cv_total = total_energy.var() / (kB * T**2)

    print(f'La capacidad calorifica es: {Cv_total}')
    return Cv_total




def run_hoomd_simulation(temp, ruta_destino, length_minibox, equilibracion, muestreo, modo='isoterma', rho=0.5):
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
        lx, ly, lz = 100.0, 25.0, 25.0
        ndiv = [40, 20, 20]
        n_total = ndiv[0] * ndiv[1] * ndiv[2]

        # Cálculo de espaciado #== Se centran las partículas en un rectángulo interior. ==
        dx = length_minibox / ndiv[0]
        dy = ly / ndiv[1]
        dz = lz / ndiv[2]
        

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
            x = np.linspace(-length_minibox/2 + dx/2, length_minibox/2 - dx/2, ndiv[0])
            y = np.linspace(-ly/2 + dy/2, ly/2 - dy/2, ndiv[1])
            z = np.linspace(-lz/2 + dz/2, lz/2 - dz/2, ndiv[2])

        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        pos = np.stack((xx.flatten(), yy.flatten(), zz.flatten()), axis=-1)
        snap.particles.position[:] = pos # Se copian las coordenadas generadas en un snap de hoomd
            

    sim.create_state_from_snapshot(snap) # Se crea la simulación a partir de las posiciones del snap generado
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=temp) # Se asignan velocidades iniciales

    # Configuración del potencial
    n_exp, m_exp = 12.0, 6.0
    epsilon, sigma = 1.0, 1.0

    # Factor de corrección de largo alcance para Mie
    coef = (n_exp / (n_exp - m_exp)) * (n_exp / m_exp) ** (m_exp / (n_exp - m_exp))

    cell = hoomd.md.nlist.Cell(buffer=0.4) # Se espera que la partícula se mueva 0.4 unidades 
    # cell: Divide la caja en rejillas de interacción 
    mie = hoomd.md.pair.Mie(nlist=cell, default_r_cut=4.0) # Se define el potecial mie como en el código en C con r_cut = 4.0
    mie.params[('A', 'A')] = dict(epsilon=epsilon, sigma=sigma, n=n_exp, m=m_exp) # Se define la manera de interacción entre las partículas 

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

    # El archivo de salida es una tabla
    if modo == 'isoterma':
        log_filename = f"todo_T{temp:.2f}_rho{rho:.4f}.csv"
        gsd_filename = f"trajectory_T{temp:.2f}_rho{rho:.4f}.gsd"
    else: 
        log_filename = os.path.join(ruta_destino, f"todo_T{temp:.2f}.csv")
        gsd_filename = os.path.join(ruta_destino, f"trajectory_T{temp:.2f}.gsd")


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

    sim.run(equilibracion)
    sim.run(muestreo)

    print(f'Simulación finalizada: Rho = {rho}')


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
    

def encontrar_equilibrio_hoomd(archivo_csv, pasos_totales=1000000, ancho_bloques=100, variacion_permitida=0.03, fraccion_cola=0.2, mostrar_progreso=True):
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
        # 1. Carga y limpieza robusta (como la función anterior)
        df = pd.read_csv(archivo_csv, sep=r's\+', engine='python')
        df.columns = df.columns.str.strip()
        df = df.apply(pd.to_numeric, errors='coerce').dropna().reset_index(drop=True)
        
        # 2. Reconstrucción de ejes (basado en lo que ya sabemos de tus archivos)
        cols_step = [c for c in df.columns if 'step' in c.lower() or 'timestep' in c.lower()]
        step = df[cols_step[0]] if cols_step else np.linspace(0, pasos_totales, len(df))
        
        # Identificación flexible de columnas
        pe = df[[c for c in df.columns if 'potential' in c.lower() or 'pe' == c.lower()][0]]
        ke = df[[c for c in df.columns if 'kinetic' in c.lower() or 'ke' == c.lower()][0]]
        etot = pe + ke
        # La presión no siempre está, la buscamos
        col_pres = [c for c in df.columns if 'pressure' in c.lower()]
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