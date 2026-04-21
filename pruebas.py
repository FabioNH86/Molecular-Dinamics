import os
from funciones import calcular_presiones_vapor, graficar_evolucion_presion, generar_dimensiones_partículas
import numpy as np
os.system('cls' if os.name == 'nt' else 'clear')

# ruta_archivo = 'Resultados/P9_LV_Mie/T=0.90/presiones_T0-9.dat'

# #calcular_promedios_energía_claude(archivo=ruta_archivo, ancho_bloques=100000)

# presiones, presion_vapor = calcular_presiones_vapor(ruta_archivo, configuraciones_consideradas=1000000)

# graficar_evolucion_presion(presiones) # Para conocer la configuración a la que se alcanzó el equilibrio, usar función calcular_promedios_energía


densidades = [round(x, 3) for x in np.arange(0.005, 0.9, 0.105)]
print(densidades)
print(len(densidades))

for each_rho in densidades:
    lados_particulas = generar_dimensiones_partículas(densidad=each_rho)

    print(lados_particulas)
    # Verificación rápida del producto
    n_final = lados_particulas[0] * lados_particulas[1] * lados_particulas[2]
    print(f"Rho: {each_rho:.3f} -> Dimensiones: {lados_particulas} (N aprox: {n_final})")
    print(f"Densidad obtenida: {n_final/16000}")