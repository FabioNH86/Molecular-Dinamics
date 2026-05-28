from funciones import crear_snapshot 
import gsd.hoomd
import os

num_prueba = 3

carpeta_guardado = 'Pruebas/'

snapshot = crear_snapshot(
    densidad_goticula=0.3, 
    aspect_ratio=4.0, 
    concentracion_porcentual_monomeros=1, 
    monomeros_en_polimero=8)

frame = gsd.hoomd.Frame()
frame.particles.N = snapshot.particles.N
frame.particles.position = snapshot.particles.position
frame.particles.typeid = snapshot.particles.typeid
frame.particles.types = snapshot.particles.types
frame.configuration.box = snapshot.configuration.box

# 2. Guardas el snapshot en un archivo GSD usando gsd.hoomd
nombre_archivo = f"configuracion_inicial_{num_prueba}.gsd"

ruta_completa = carpeta_guardado + nombre_archivo

with gsd.hoomd.open(name=ruta_completa, mode='x') as archivo_gsd:
    archivo_gsd.append(frame)

print(f"¡Archivo {nombre_archivo} creado con éxito!")

