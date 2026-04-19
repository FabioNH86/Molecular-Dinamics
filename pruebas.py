import os
from funciones import calcular_promedios_energía, calcular_promedios_energía_claude, calcular_promedios_energía_claude_2
os.system('cls' if os.name == 'nt' else 'clear')

ruta_archivo = 'Resultados/P7_LV_Mie/T=0.60/todo_T0-6.dat'

calcular_promedios_energía_claude(archivo=ruta_archivo, ancho_bloques=100000)