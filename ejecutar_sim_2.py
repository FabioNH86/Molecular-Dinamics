import os
import sys
from funciones import crear_primer_frame_solvente, correr_simulacion_solvente

"""
Script para simulación de solvente puro con registro del tensor de presión.
Uso: python3 ejecutar_simulacion_solvente.py <num_prueba> <temperatura> <n_particulas>
"""

if len(sys.argv) < 4:
    print("❌ Error: Faltan argumentos. Uso: python3 ejecutar_simulacion_solvente.py <num_prueba> <temperatura> <n_particulas>")
    sys.exit(1)

# Capturar argumentos
num_prueba   = int(sys.argv[1])
temperatura  = float(sys.argv[2])
n_particulas = int(sys.argv[3])

ruta_base          = f"Resultados/HOOMD/P{num_prueba}_Solvente"
directorio_original = os.getcwd()

if not os.path.exists(ruta_base):
    os.makedirs(ruta_base, exist_ok=True)

print(f"\n🧪 [PROCESO AISLADO] Iniciando: T={temperatura}, N={n_particulas}")

try:
    os.chdir(os.path.join(directorio_original, ruta_base))

    pasos_equil    = int(0)
    pasos_muestreo = int(3e6)

    file_id = f"Solv_T{temperatura:.2f}"
    fn_gsd  = f"{file_id}.gsd"

    print(f"Buscando: {fn_gsd}")

    if os.path.exists(fn_gsd):
        print(f"⚠️  Archivo '{fn_gsd}' ya existe. Esta simulación no admite reanudación por ahora.")
        print("    Elimina el archivo o elige otro num_prueba para correr de nuevo.")
        sys.exit(0)
    else:
        print("🆕 No se encontró GSD previo. Iniciando simulación desde cero...")

        snapshot = crear_primer_frame_solvente(
            densidad    = 0.6,
            aspect_ratio= 1,
            n_particulas= n_particulas,
        )

        correr_simulacion_solvente(
            snapshot     = snapshot,
            temp         = temperatura,
            equilibracion= pasos_equil,
            muestreo     = pasos_muestreo,
            aspect_ratio = 2.0,
        )

        print(f"✅ ¡Simulación finalizada con éxito para T={temperatura}, N={n_particulas}!")

except Exception as e:
    print(f"❌ Error crítico en este proceso: {e}")
    sys.exit(1)

finally:
    os.chdir(directorio_original)