## NOTAS ÚTILES

*Crear un ejecutable de c:* `gcc -O3 -fopenmp mie.c -o mie_exec -lm -march=native`

Para Intel: `gcc -O3 -fopenmp -march=native -mtune=native -flto mie.c -o mie_exec -lm`

Para caja de simulación: `vmd > pbc box -center origin`