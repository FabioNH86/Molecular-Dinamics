import os 
import shutil
import subprocess
import numpy as np

"""
Este código se encargar de automatizar simulaciones alterando el archivo in.dat después de cada 
corrida según una lista de temperaturas pre-asignadas. 
La simulación es llevada acabo por el ejecutable `mie`. Los archivos generados se guardan automáticamente
(movie.gro, presiones.dat, resumen.dat, todo.dat y xyz.dat) 
en una carpeta que señala la temperatura usada (Resultados/P{num_prueba}_LV_Mie/T=*).

Fabio Noriega Hernández
"""
# -- SEÑALA EL NÚMERO DE ENSAYO QUE HARÁS PARA ALMACENAR LOS RESULTADOS EN SU CARPETA CORRESPONDIENTE --
num_prueba = 9

densidades = [round(x, 2) for x in np.arange(0.7, 0.81, 0.01)]
ruta_base = f"Resultados/P"
ejecutable = "./mie_exec"

def actualizar_entradas(temp):
    """Modifica el archivo in.dat con una temperatura seleccionada."""

    contenido = f"""2         ncolectivo_1_NVE_2_NVT_3_NPT
1         nopcion_1_inicializacion_2_continuacion
3         dofx_numero_de_dimensiones
12.0      expn_exponente_para_la_parte_repulsiva_12_para_LJ
6.0       expm_exponente_para_la_parte_atractiva_6_en_general
0.001     time_tiempo_de_integracion
40 20 20  Lx_Ly_Lz_en_unidades_reducidas
30 10 10  ndivx_ndivy_ndivz_numero_de_atomos_por_lado
{temp:<10.2f}   temp_en_unidades_reducidas_para_asignar_velocidades
0.01      taut_para_termostato_berendsen
0.0       presion_en_unidades_reducidas_para_presostato
0.0       taup_para_presostato_berendsen
4.0       rcut_r_de_corte_del_potencial_menor_a_la_mitad_de_la_caja
500000    nconfequi_numero_de_configuraciones_para_equilibrar
1500000   nconf_numero_de_configuraciones
50000     nperfil_frecuencia_para_calcular_distribuciones
0.1       deltar_ancho_del_intervalo_para_perfil_de_densidad
10000     nmovie_frecuencia_para_tomar_fotos
10000     nprint_frecuencia_para_imprimir_en_pantalla"""

    with open("in.dat", "w") as f:
        f.write(contenido)

print("--- Iniciando batería de simulaciones ---")

for T in temperaturas:
    nombre_carpeta = f"T={T:.2f}"
    ruta_destino = os.path.join(ruta_base, nombre_carpeta)

    # Si la carpeta no existe aún:
    if not os.path.exists(ruta_destino):
        os.makedirs(ruta_destino)
        print(f"\n📁 Creada carpeta: {ruta_destino}")

    # Se actulaiza el archivo in.dat
    actualizar_entradas(T)
    print(f"🌡️ Configurando T = {T:.2f}...")


    # Se ejecuta la simulación
    try:
        print(f"🚀 Ejecutando simulación para {nombre_carpeta}...")
        subprocess.run(ejecutable, check=True)

        # Se mueven los archivos generados
        archivos_a_mover = ["movie.gro", "presiones.dat", "resumen.dat", "todo.dat", "xyz.dat", "gr_dm.dat", "perfil_dm.dat"]

        for archivo in archivos_a_mover:
            if os.path.exists(archivo):
                # Se actuliza el nombre del archivo para evitar confusiones
                nombre_nuevo = archivo.replace(".dat", f"_T{str(T).replace('.', '-')}.dat")
                shutil.move(archivo, os.path.join(ruta_destino, nombre_nuevo))
            else:
                print(f"⚠️ Advertencia: No se encontró {archivo}")

        print(f"✅ Finalizada T={T:.2f}. Archivos guardados. \n")

    except subprocess.CalledProcessError:
        print(f"❌ Error crítico en la ejecución de 'new_mie' para T={T:.2f}")
        break

print("\n--- Todas las simulaciones han terminado ---")
