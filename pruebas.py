import pandas as pd
import numpy as np
import io
from funciones import calcular_promedios_energía

ruta_archivo = 'Resultados/P7_LV_Mie/T=0.70/todo_T0-7.dat'

calcular_promedios_energía(archivo=ruta_archivo)