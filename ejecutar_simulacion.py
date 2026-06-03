import os
import sys
from funciones import crear_primer_frame, correr_simulacion

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
n_monomeros_totales = 12000

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
    pasos_equil = int(5e5)
    pasos_muestreo = int(1e6)
    
    # Generar Snapshot
    snapshot = crear_primer_frame(
        densidad_goticula=0.3,
        aspect_ratio=4,
        concentracion_porcentual_monomeros=1.0,
        n_monomeros=n_monomeros_totales,
        monomeros_en_polimero=n_monomeros
    )
    
    # Ejecutar en HOOMD
    correr_simulacion(
        snapshot=snapshot,  
        temp=temperatura,
        equilibracion=pasos_equil, 
        muestreo=pasos_muestreo,
        eps_SP=eps,
        mon_cadena=n_monomeros
    )
    print(f"✅ ¡Simulación finalizada con éxito para T={temperatura}, N={n_monomeros}!")

except Exception as e:
    print(f"❌ Error crítico en este proceso: {e}")
    sys.exit(1) # Reportar fallo al script de Bash
finally:
    os.chdir(directorio_original)