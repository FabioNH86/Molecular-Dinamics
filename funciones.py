import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt

def calcular_densidades(filename, start_conf=500000, end_conf=1000000, num_atom=3000, num_bines=100):
    with open(filename, 'r') as f:
        configuración_1 = f.readline() # Linea 1: CONFIGURACION X
        partes = configuración_1.split()
        n_conf = int(partes[-1]) # Se extrae la primer configuración

        # DEBUG:
        #print(configuración_1)

        num_atom = int(f.readline())
        #DEBUG: 
        #print(num_atom)

        lineas_atomos = [f.readline() for _ in range(num_atom)] # Se saltan todas las posiciones del primer frame
        
        dimensiondes_caja = f.readline()
        #print(dimensiondes_caja)
        Lx, Ly, Lz = map(float, dimensiondes_caja.split()) 
        Lx /= 0.3405 # Se divide por un factor de conversión usado a la hora de escribir la movie.gro
        Ly /= 0.3405
        Lz /= 0.3405
        print("Las dimensiones de la caja son:")                    
        print(f'Lx: {Lx:.4f}\nLy: {Ly:.4f}\nLz: {Lz:.4f} \n')

        volumen_total = Lx * Ly * Lz 
        densidad_global = num_atom / volumen_total
        print(f'El volumen total de la caja es: {volumen_total:.4f} y se dividirá en {num_bines} bines \n')
        print(f"La densidad global es: {densidad_global:.4f}")

        dx = Lx / num_bines
        volumen_bin = dx * Ly * Lz
        print(f'El volumen de cada bin será: {volumen_bin:.4f} y su ancho será de: {dx:.4f} \n')
        
        # Cortes a lo largo de Lx
        cortes_x = np.linspace(0, Lx, num_bines + 1)
        # Calcular centros de los bines (Eje X)
        centros_x = (cortes_x[:-1] + cortes_x[1:]) / 2

        # -- HASTA AQUÍ SE ESTABLECEN LAS CONDICIONES DE INICIO --
        # (Cuando la caja de simulación está acomodadita)

        # -- AQUI INICIA EL PROCESADO DE LAS DEMÁS POSICIONES --
        n_conf = 1
        contador_frames = 0
        todas_las_densidades = []

        while n_conf < end_conf:
            configuracion = f.readline()
            partes = configuracion.split()
            n_conf = int(partes[-1])
            contador_frames += 1

            # Imprimimos cada 100 mil configuraciones
            if n_conf % 100000 == 0:
                print(f"Procesando configuración: {n_conf}")

            new_num_atom = int(f.readline())
            #print(new_num_atom)

            if new_num_atom != num_atom:
                raise ValueError("Cambió el número de átomos")
            
            lineas_atomos = [f.readline() for _ in range(num_atom)]

            atomos = pd.read_csv(
                io.StringIO("".join(lineas_atomos)),
                sep=r'\s+',
                header=None
            )
            
            # Seleccionar las últimas 4 columnas sin importar cuántas haya
            atomos = atomos.iloc[:, -4:].copy()
            atomos.columns = ['ID', 'X', 'Y', 'Z']
            #print(atomos.tail())
            f.readline() # Nos saltamos una línea
            #print("")

            # También se ajustan las coordenadas 
            atomos['X'] /= 0.3405
            atomos['Y'] /= 0.3405
            atomos['Z'] /= 0.3405
            #print(atomos.head())
            
            # -- NORMALIZAMOS LAS POSICIONES --
            atomos['X'] -= atomos['X'].min()
            atomos['Y'] -= atomos['Y'].min()
            atomos['Z'] -= atomos['Z'].min()


            # print(atomos.tail())
            # print("")

        # -- Ahora sí, comenzamos con el cálculo de densidades en cada configuración
            densidades_por_bin = []
            bin = 0
            for i in range(num_bines):
                bin += 1
                # Definimos los límites del bin 
                x_inicial = cortes_x[i]
                x_final = cortes_x[i+1]

                # Revisamos qué átomos tienen coordenadas que entran en el rango del bin actual
                bin_actual = atomos[(atomos['X'] >= x_inicial) & (atomos['X'] < x_final)] 
                # Contamos dichos átomos
                num_atom_bin = len(bin_actual)
                # Calculamos la densidad 
                densidad_bin = num_atom_bin / volumen_bin

                # Almacenamos cada densidad
                densidades_por_bin.append(densidad_bin)

            # Guardamos las listas de densidades por bin en un dataframe como fila, cada fila debe ser una configuración
            # Por lo tanto el dataframe debe ser (numero de bines x configuraciones)
            todas_las_densidades.append(densidades_por_bin)

        
        

        densidades_por_configuracion = pd.DataFrame(todas_las_densidades)

        ultimas_densidades = densidades_por_configuracion.iloc[-500000:]

        #perfil_promedio = densidades_por_configuracion.mean(axis=0) 
        perfil_promedio = ultimas_densidades.mean()
        desviacion_estandar = ultimas_densidades.std()
        
    return centros_x, perfil_promedio, desviacion_estandar





# centros, perfil_prom, desv = calcular_densidades('movie.gro')
# plt.plot(centros, perfil_prom)
# plt.show()


            
                

