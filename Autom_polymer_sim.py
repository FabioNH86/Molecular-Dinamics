import os 
from funciones import crear_primer_frame, correr_simulacion

"""
Este código se encarga de automatizar una serie de simulaciones de dinámica molecular
utilizando la librería HOOMD-blue. El script itera sobre configuraciones de polímeros
con solvente puntual, variando la temperatura y la afinidad del solvente (eps_SP).

Los resultados (trayectorias en formato .gsd y datos termodinámicos en .csv) se 
almacenan automáticamente en directorios organizados por configuración.

Fabio Noriega Hernández
Mayo 2026
"""

# -- SEÑALA EL NÚMERO DE ENSAYO QUE HARÁS PARA ALMACENAR LOS RESULTADOS EN SU CARPETA CORRESPONDIENTE --
num_prueba = 8

# Parámetros del sistema a explorar
temperaturas = [0.7, 0.9, 1.1]
lista_eps_SP = [0.1, 1.0]  # Diferentes afinidades solvente-polímero

# Configuración del polímero y solvente
monomeros_por_polimero = [8, 16, 24]  # Número de monómeros por cadena polimérica
n_monomeros_totales = 100  # Total de monómeros en el sistema (ajustar según tu simulación)



ruta_base = f"Resultados/HOOMD/P{num_prueba}_Polimero_Solvente"

# Parámetros de tiempo de simulación
pasos_equil = int(5e5)
pasos_muestreo = int(1e6)

total_simulaciones = len(temperaturas) * len(lista_eps_SP) * len(monomeros_por_polimero)
print(f"Se realizarán un total de {total_simulaciones} simulaciones :)")

# Creamos una carpeta base para esta prueba si no existe
if not os.path.exists(ruta_base):
    os.makedirs(ruta_base)
    print(f"\n📁 Creada carpeta base de resultados: {ruta_base}")

# Cambiamos temporalmente el directorio de trabajo al de destino para que los 
# archivos log_*.csv y traj_*.gsd generados por HOOMD se guarden directamente ahí
directorio_original = os.getcwd()

try:
    for eps in lista_eps_SP:
        print(f"\n💧 Evaluando afinidad Solvente-Polímero eps_SP = {eps}")
        
        for temperatura in temperaturas:
            print(f"  🧪 Ejecutando simulación a T = {temperatura:.2f}...")
            

            for n_monomeros in monomeros_por_polimero:
                n_cadenas = n_monomeros_totales // n_monomeros
                print(f"    - Polímero con {n_monomeros} monómeros por cadena")
            
                # Nos movemos a la ruta base para que HOOMD escriba allí
                os.chdir(os.path.join(directorio_original, ruta_base))
                try:
                    # Llamamos a las funciones 
                    snapshot = crear_primer_frame(densidad_goticula=0.3,
                                        aspect_ratio=4,
                                        concentracion_porcentual_monomeros=1.0,
                                        n_monomeros=n_monomeros_totales,
                                        monomeros_en_polimero=n_monomeros)
                    
                    correr_simulacion(snapshot=snapshot,  
                                      temp=temperatura,
                                      equilibracion=pasos_equil, 
                                      muestreo=pasos_muestreo,
                                      eps_SP=eps,
                                      mon_cadena=n_monomeros)
                
                    print('=========== Nueva simulación ============= \n')
                    

                except Exception as e:
                    print(f"  ❌ Error ejecutando T={temperatura}, eps_SP={eps}: {e}")
                finally:
                    # Siempre regresamos al directorio original por seguridad
                    os.chdir(directorio_original)

except Exception as e:
    print(f"❌ Error general en el bucle de simulación: {e}")
finally:
    # Asegurar que el script termine en el directorio correcto si algo falla radicalmente
    os.chdir(directorio_original)

print("\n--- Todas las simulaciones han terminado ---")