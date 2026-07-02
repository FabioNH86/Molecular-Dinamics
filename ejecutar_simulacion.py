import os
import sys
from funciones import crear_primer_frame_homopolimero, correr_simulacion_homoplimero, continue_sim_from_gsd

"""
Script preparado para recibir parámetros desde Bash.
Uso: python3 ejecutar_simulacion.py <num_prueba> <eps> <temperatura> <n_monomeros>
"""

if len(sys.argv) < 5:
    print("❌ Error: Faltan argumentos. Uso: python3 ejecutar_simulacion.py <num_prueba> <eps> <temperatura> <n_monomeros>")
    sys.exit(1)

# Capturar argumentos de la línea de comandos
num_prueba = int(sys.argv[1])
eps = float(sys.argv[2])
temperatura = float(sys.argv[3])
n_monomeros = int(sys.argv[4])
n_monomeros_totales = 2400

ruta_base = f"Resultados/HOOMD/P{num_prueba}_Polimero_Solvente"
directorio_original = os.getcwd()

# Crear directorio si no existe
if not os.path.exists(ruta_base):
    os.makedirs(ruta_base, exist_ok=True)

print(f"\n🧪 [PROCESO AISLADO] Iniciando: T={temperatura}, eps_SP={eps}, Monómeros/Cadena={n_monomeros}")

try: 
    # Moverse a la ruta de resultados antes de inicializar HOOMD
    os.chdir(os.path.join(directorio_original, ruta_base))
        
    # Parámetros de tiempo
    pasos_equil = int(2.5e6)
    pasos_muestreo = int(500_000)

    # Agregamos los nombres de los archivos que ya deberían estar creados
    file_id = f"Poly-Solv_T{temperatura:.2f}_epsSP{eps:.2f}"
    fn_gsd = f"{file_id}_monom_{n_monomeros}.gsd"
    print(f"Buscando: {fn_gsd}")

    # Revisamos si ya están creados los archivos 
    if os.path.exists(fn_gsd):
        print(f"🔄 Archivo previo '{fn_gsd}' detectado. Intentando reanudar...")


        # Llamamos la función que continúa la simulación 
        continue_sim_from_gsd(
            archivo_gsd=fn_gsd,
            temp=temperatura,
            pasos_muestreo=int(3e6),
            eps_SP=eps,
            mon_cadena=n_monomeros,
            aspect_ratio=2.0,
            pasos_equilibracion=pasos_equil
        )
        print(f"✅ ¡Simulación REANUDADA y finalizada con éxito para T={temperatura}, N={n_monomeros}!")

    else:
        print("🆕 No se encontró GSD previo. Iniciando simulación desde cero...")
        # Generar Snapshot
        snapshot = crear_primer_frame_homopolimero(
            densidad_goticula=0.6,
            aspect_ratio=1,
            concentracion_porcentual_monomeros=1.0,
            n_monomeros=n_monomeros_totales,
            monomeros_en_polimero=n_monomeros
        )
        
        # Ejecutar en HOOMD
        correr_simulacion_homoplimero(
            snapshot=snapshot,  
            temp=temperatura,
            equilibracion=pasos_equil, 
            muestreo=pasos_muestreo,
            eps_SP=eps,
            mon_cadena=n_monomeros,
            aspect_ratio=2.0
        )
        print(f"✅ ¡Simulación finalizada con éxito para T={temperatura}, N={n_monomeros}!")

except Exception as e:
    print(f"❌ Error crítico en este proceso: {e}")
    sys.exit(1) # Reportar fallo al script de Bash
finally:
    os.chdir(directorio_original)