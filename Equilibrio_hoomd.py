import glob 
import os 
# Asegúrate de que estas funciones estén actualizadas en funciones.py
from funciones import encontrar_equilibrio_hoomd, calcular_presiones_vapor, calcular_tension_superficial

os.system('clear')

"""
Análisis de Propiedades Termodinámicas (Presión de Vapor y Tensión Superficial)
para simulaciones HOOMD-blue.

Fabio Noriega Hernández
Abril 2026
"""

# -- CONFIGURACIÓN --
num_prueba = 5
temperaturas_originales = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20] 
temperaturas = [T for T in temperaturas_originales if T != 0.95]

# La ruta ahora apunta a la carpeta de HOOMD
ruta_comun = f'Resultados/P{num_prueba}_HOOMD_Mie'

for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    # HOOMD guarda todo en el CSV que definimos como todo_T{temp}.csv
    ruta_csv = os.path.join(ruta_comun, nombre_carpeta, f"todo_T{T:.2f}.csv")
    archivos_todo = glob.glob(ruta_csv)

    if not archivos_todo: 
        print(f"⚠️ No se encontró el archivo CSV para T={T:.2f}")
        continue
        
    archivo_actual = archivos_todo[0]

    print(f'🔍 Procesando Termodinámica para T={T}')
    print('='*100)

    # 1. Encontrar el paso de equilibrio (usando la nueva función adaptada a HOOMD)
    # ancho_bloques=100 suele ser bueno si guardas datos cada 1000 pasos
    paso_estable = encontrar_equilibrio_hoomd(
        archivo_csv=archivo_actual, 
        pasos_totales=1500000,
        ancho_bloques=10, 
        mostrar_progreso=True,
        variacion_permitida=0.1
    )

    # if paso_estable is not None:
    #     # 2. Calcular Presión de Vapor
    #     # En HOOMD, el CSV ya contiene la columna 'pressure' (P_total)
    #     # Si tu función 'calcular_presiones_vapor' espera el archivo original, pásale el CSV.
    #     presiones, presion_vapor = calcular_presiones_vapor(
    #         archivo=archivo_actual, 
    #         paso_inicial=paso_estable
    #     )

    #     # 3. Calcular Tensión Superficial
    #     # Nota: Asegúrate de que la longitud_perpendicular coincida con la lx de tu simulación
    #     # Para el modo 'barrido' usamos 48.92 según tu función anterior.
    #     calcular_tension_superficial(
    #         df_presiones=presiones, 
    #         longitud_perpendicular_interface=48.92 
    #     )
    # else:
    #     print(f"⚠️ Saltando cálculos para T={T}: El sistema no alcanzó estabilidad sostenida.")

    print('='*100)
    print('\n')

print("--- Proceso de análisis finalizado ---")