from sim_simple import run_hoomd_simulation
import os

carpeta_actual = os.getcwd()

run_hoomd_simulation(temp=0.7, 
                     modo='barrido', 
                     ruta_destino=carpeta_actual, 
                     length_minibox=50, equilibracion=100_000,
                     muestreo=1_000_000,
                     ndiv_entrada=[50, 50, 50],
                     periodic_zeromomentum=100)