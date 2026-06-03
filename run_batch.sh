#!/bin/bash

# Automatización de simulaciones HOOMD-blue mediante aislamiento de procesos
# Fabio Noriega Hernández - Mayo 2026

NUM_PRUEBA=8

# Configuración de arreglos (espacios como separadores)
LISTA_EPS=(1.0)
TEMPERATURAS=(0.7 0.9 1.1)
MONOMEROS_POR_POLIMERO=(8 16 24)

# Calcular total de simulaciones
TOTAL_SIM=$(( ${#LISTA_EPS[@]} * ${#TEMPERATURAS[@]} * ${#MONOMEROS_POR_POLIMERO[@]} ))
echo "=========================================================="
echo " Se realizarán un total de $TOTAL_SIM simulaciones en serie :)"
echo "=========================================================="

CONTADOR=1

# Bucles anidados en Bash
for eps in "${LISTA_EPS[@]}"; do
    echo -e "\n💧 [Bash] Evaluando afinidad Solvente-Polímero eps_SP = $eps"
    
    for temp in "${TEMPERATURAS[@]}"; do
        for n_mon in "${MONOMEROS_POR_POLIMERO[@]}"; do
            
            echo -e "\n----------------------------------------------------------"
            echo "🚀 [Simulación $CONTADOR/$TOTAL_SIM]"
            echo "----------------------------------------------------------"
            
            # Ejecutar el script de Python pasando las variables como argumentos
            python3 ejecutar_simulacion.py "$NUM_PRUEBA" "$eps" "$temp" "$n_mon"
            
            # Capturar el código de salida de Python
            STATUS=$?
            if [ $STATUS -eq 0 ]; then
                echo "⭐ Progreso registrado correctamente."
            else
                echo "⚠️ El proceso de Python falló, pero Bash continuará con la siguiente."
            fi
            
            CONTADOR=$((CONTADOR + 1))
        done
    done
done

echo -e "\n=========================================================="
echo "--- 🎉 Todas las simulaciones del script han terminado ---"
echo "=========================================================="