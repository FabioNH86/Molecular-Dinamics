import os
from funciones import calcular_presiones_vapor, graficar_evolucion_presion
os.system('cls' if os.name == 'nt' else 'clear')

ruta_archivo = 'Resultados/P9_LV_Mie/T=0.90/presiones_T0-9.dat'

#calcular_promedios_energía_claude(archivo=ruta_archivo, ancho_bloques=100000)

presiones, presion_vapor = calcular_presiones_vapor(ruta_archivo, configuraciones_consideradas=1000000)

graficar_evolucion_presion(presiones) # Para conocer la configuración a la que se alcanzó el equilibrio, usar función calcular_promedios_energía