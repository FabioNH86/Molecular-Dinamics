import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


# -- SELECCIONA EL NÚMERO DE ENSAYO A ANALIZAR --
# num_prueba = 5

# Limpieza de terminal
os.system('clear')

# ruta_comun = f'Dinámicas_moleculares/Resultados/P{num_prueba}_LV_Mie'

# # Definimos la ruta de guardado
# ruta_graficos = os.path.join(ruta_comun, 'Graficos_2')

# # Si la carpeta no existe, la creamos automáticamente
# if not os.path.exists(ruta_graficos):
#     os.makedirs(ruta_graficos)
#     print(f"Carpeta creada: {ruta_graficos}")

# Carga de datos
todo = pd.read_csv('Dinámicas_moleculares/Resultados/P5_LV_Mie/T=1.00/todo_T1-0.dat', sep=r'\s+', header=None)

# Definición de variables por columna
config      = todo[0]
density     = todo[1]
kinetic_e   = todo[2]
potential_e = todo[3]
total_e     = todo[4]
temperature = todo[5]
pressure    = todo[6]
error       = todo[7]

n: int = 3 # <----- SELECCIONA LA PROPIEDAD A ANALIZAR

# mean_potential_reported = -5.44
# mean_pressure_reported = 4.254

# Parámetro de corte para la fase de equilibración
salto: int = 100

# --- FIGURA 1: ANÁLISIS DE ENERGÍAS Y DENSIDAD (SUBPLOTS) ---
fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')

# Panel 1: Energía Cinética (Fase inicial)
axs[0, 0].plot(config[:salto], kinetic_e[:salto], color='tab:blue')
axs[0, 0].set_title('Energía Cinética (Equilibración)')
axs[0, 0].set_xlabel('Configuración')
axs[0, 0].set_ylabel('Energía')

# Panel 2: Energía Potencial (Fase inicial)
axs[0, 1].plot(config[:salto], potential_e[:salto], color='tab:red')
#axs[0, 1].axhline(y=mean_potential_reported, color='black', linestyle='--', linewidth=1.5, label='Referencia Tabla')
axs[0, 1].set_title('Energía Potencial (Equilibración)')
axs[0, 1].set_xlabel('Configuración')
axs[0, 1].set_ylabel('U / ε')
axs[0, 1].legend()

# Panel 3: Energía Total (Fase inicial)
axs[1, 0].plot(config[:salto], total_e[:salto], color='tab:purple', alpha=0.8)
axs[1, 0].set_title('Energía Total del Sistema')
axs[1, 0].set_xlabel('Configuración')
axs[1, 0].set_ylabel('Energía')

# Panel 4: Variación de la Presión (Toda la simulación)
axs[1, 1].plot(config[:salto], pressure[:salto], color='tab:green')
axs[1, 1].set_title('Evolución de la Densidad')
#axs[1, 1].axhline(y=mean_pressure_reported, color='black', linestyle='--', linewidth=1.5, label='Referencia Tabla')
axs[1, 1].set_xlabel('Configuración')

# --- FIGURA 2: HISTOGRAMA DE EQUILIBRIO DE ENERGÍA POTENCIAL ---
# Calculamos estadística solo de la fase de producción (después del salto)
fase_produccion = potential_e[salto:]
sigma = fase_produccion.std()
mean = fase_produccion.mean()

plt.figure(figsize=(8, 5))
plt.hist(fase_produccion, bins=40, color='skyblue', edgecolor='black', alpha=0.7)

# Se añaden líneas de referencia al histograma
plt.axvline(mean, color='red', linestyle='dashed', linewidth=1.5, label=f'Promedio: {mean:.4f}')
plt.axvline(mean + sigma, color='gray', linestyle='dotted', label=f'1 $\sigma$: {sigma:.4f}')
plt.axvline(mean - sigma, color='gray', linestyle='dotted')
#plt.axvline(mean_potential_reported, color='black', linestyle='dotted', label=f'Valor reportado en NIST: {mean_potential_reported}')

plt.xlabel('Energía potencial')
plt.ylabel('Frecuencia')
plt.title('Histograma de Energía Potencial')
plt.legend()
plt.grid(axis='y', alpha=0.3)


# --- FIGURA 3: HISTOGRAMA DE EQUILIBRIO DE PRESIÓN ---
# Calculamos estadística solo de la fase de producción (después del salto)
fase_produccion_p = pressure[salto:]
sigma_p = fase_produccion_p.std()
mean_p = fase_produccion_p.mean()

plt.figure(figsize=(10, 6))
plt.hist(fase_produccion_p, bins=40, color='lightgreen', edgecolor='black', alpha=0.7)

# Se añaden líneas de referencia al histograma
plt.axvline(mean_p, color='red', linestyle='dashed', linewidth=1.5, label=f'Promedio: {mean_p:.4f}')
plt.axvline(mean_p + sigma_p, color='gray', linestyle='dotted', label=f'1 $\sigma$: {sigma_p:.4f}')
plt.axvline(mean_p - sigma_p, color='gray', linestyle='dotted')
#plt.axvline(mean_pressure_reported, color='black', linestyle='dotted', label=f'Valor reportado en NIST: {mean_pressure_reported}')

plt.xlabel('Presión')
plt.ylabel('Frecuencia')
plt.title('Histograma de Presión')
plt.legend()
plt.grid(axis='y', alpha=0.3)

# Mostrar estadísticas en terminal
total_datos = len(todo)
print(f"{' ESTADÍSTICAS DE ENERGÍA POTENCIAL ':=^40}")
print(f"Promedio: {mean:>20.6f}")
print(f"Desviación Estándar: {sigma:>10.6f}")
print(f"{'='*40}")
print(f"El total de configuraciones es: {len(todo[0])}")

# Se divide cada propiedad por bloques para promediar sus avances
tam_bloque: int = 50000
n_bloque = 1
# Usamos solo la fase de producción para que los bloques reflejen el equilibrio
datos_bloques = todo[n][salto:] 
total_puntos = len(datos_bloques)

# Listas de resultados:
promedios = []
desviaciones = []
indices_bloques = []
etiquetas_bloques =  []

print(f"\n{' ANÁLISIS DE ESTABILIDAD: COLUMNA {n} ':=^40}")
for i in range(0, total_datos, tam_bloque):
    bloque_actual = todo[n].iloc[i: i + tam_bloque]

    # Ignoramos el último bloque si es muy pequeño para no sesgar la estadística
    if len(bloque_actual) < (tam_bloque * 0.5): break

    m = bloque_actual.mean()
    s = bloque_actual.std()
    # Mostramos los datos
    print(f"BLOQUE: {n_bloque}  |   PROMEDIO: {m:.6f}   |   DESVIACIÓN EST.: {s:.6f}")
    promedios.append(m)
    desviaciones.append(s)
    indices_bloques.append(f"{i}-{i + tam_bloque}")
    etiquetas_bloques.append(f"B{n_bloque}")
    n_bloque += 1

# Dataframe con la información de los bloques: 
df_bloques = pd.DataFrame({
    'Rango': indices_bloques,
    'Promedio': promedios,
    'Desviación': desviaciones
})

plt.figure(figsize=(10, 6))
plt.errorbar(etiquetas_bloques,
             promedios,
             yerr=desviaciones,
             fmt='o',
             color='tab:red',
             ecolor='gray',
             elinewidth=1.5,
             capsize=4,
             capthick=2,
             markersize=7,
             label='Promedio Bloque +- $\sigma$'
             )

# Línea de Promedio Global de la simulación
promedio_sim = np.mean(promedios)
plt.axhline(y=promedio_sim, color='red', linestyle='--', alpha=0.7, label=f'Promedio Global: {promedio_sim:.4f}')

# Línea de Referencia NIST (Asegúrate que sea el valor de Energía: -0.08455)
#plt.axhline(y=mean_potential_reported, color='black', linestyle=':', linewidth=2, label=f'Referencia NIST: {mean_potential_reported}')

# Ajuste de escala para que se vean los bloques de cerca
# Si el valor de NIST está muy lejos, esto permitirá ver tus datos claramente
margin = (max(promedios) - min(promedios)) * 2 if len(promedios) > 1 else 0.05
plt.ylim(min(promedios) - margin, max(promedios) + margin)

plt.title('Análisis de Estabilidad por Bloques', fontsize=13)
plt.xlabel(f'Número de bloque (Bloques de {tam_bloque} datos)', fontsize=11)
plt.ylabel('Energía Potencial ($U/\epsilon$)', fontsize=11)
plt.grid(axis='y', linestyle='--', alpha=0.3)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()

plt.show()