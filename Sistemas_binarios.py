import os 
from funciones import run_sim_binary_sistem

"""
Este script de Python automatiza la ejecución de simulaciones dinámico-moleculares en HOOMD-blue, 
realizando un barrido sistemático de parámetros sobre una lista de temperaturas y potenciales 
de interacción. El código configura el estado inicial del sistema (mezcla homogénea o fase 
separada) y resuelve la integración de las ecuaciones de movimiento utilizando aceleración 
por hardware (GPU).

Los resultados de cada corrida, incluyendo trayectorias y propiedades termodinámicas, se exportan 
y organizan automáticamente en archivos con nombres descriptivos (log_*.csv y traj_*.gsd), 
etiquetados según la temperatura y el valor de epsilon ($\\epsilon_{AB}$) empleados. Esto permite 
un análisis fluido de fenómenos de repulsión y segregación de fases en sistemas binarios bajo 
el potencial Mie.

Fabio Noriega Hernández
"""
# Validación de seguridad: ¿Está el SSD montado?
if not os.path.exists("/media/fabio-noriega/Almacen_Fabio"):
    print("❌ ERROR: El SSD no parece estar montado en la ruta especificada.")
    exit()


num_prueba = 2

# --- LISTAS DE EXPLORACIÓN ---
lista_temperaturas = [0.5, 1.0, 1.5]  # De frío a caliente
lista_repulsiones = [0.1, 1.0]        # 0.1: Muy repulsivo (fomenta desmezcla), 1.0: Ideal

# Parámetros de tiempo
pasos_equil = int(5e5)
pasos_muestreo = int(1e6)

total_sims = len(lista_temperaturas) * len(lista_repulsiones) * 2 # Homogeneo y separado
print(f'Se realizarán un total de: {total_sims} simulaciones.')

ruta_base = f"/media/fabio-noriega/Almacen_Fabio/Simulaciones/Resultados/Sistemas_Binarios/Ronda_{num_prueba}"
original_path = os.getcwd()


# --- BUCLE DE EXPLORACIÓN ---
for t in lista_temperaturas:
    for e_ab in lista_repulsiones:
        
        nombre_carpeta = f"T_{t:.2f}_Eps_{e_ab:.2f}"
        ruta_completa = os.path.join(ruta_base, nombre_carpeta)

        os.makedirs(ruta_completa, exist_ok=True)

        try:
            os.chdir(ruta_completa)
            run_sim_binary_sistem(temp=t, 
                                  equilibracion=pasos_equil, 
                                  muestreo=pasos_muestreo, 
                                  eps_AB=e_ab, 
                                  sist_homegeno=True)
        except Exception as e: 
            print(f"Error en simulación Homogénea T={t}, Eps={e_ab}: {e}")
        finally:
            os.chdir(original_path)

        try:
            os.chdir(ruta_completa)
            run_sim_binary_sistem(temp=t, equilibracion=pasos_equil, muestreo=pasos_muestreo, 
                                  eps_AB=e_ab, sist_homegeno=False)
        except Exception as e:
            print(f"Error en simulación Separada T={t}, Eps={e_ab}: {e}")
        finally:
            os.chdir(original_path)


print('\n (: -- TODAS LAS SIMULACIONES HAN TERMINADO -- :)\n')