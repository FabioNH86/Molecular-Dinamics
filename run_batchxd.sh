#!/bin/bash
# Automatización de simulaciones de solvente puro con HOOMD-blue
# Fabio Noriega Hernández - Junio 2026

NUM_PRUEBA=15

# Configuración — solo temperatura y número de partículas
TEMPERATURAS=(0.60 0.65)
N_PARTICULAS=(240000)

# Total de simulaciones
TOTAL_SIM=$(( ${#TEMPERATURAS[@]} * ${#N_PARTICULAS[@]} ))

echo "=========================================================="
echo " Simulaciones de solvente puro: $TOTAL_SIM en total :)"
echo "=========================================================="

CONTADOR=1

for temp in "${TEMPERATURAS[@]}"; do
    echo -e "\n🌡️  [Bash] Evaluando temperatura T = $temp"

    for n_part in "${N_PARTICULAS[@]}"; do
        echo -e "\n----------------------------------------------------------"
        echo "🚀 [Simulación $CONTADOR/$TOTAL_SIM] T=$temp, N=$n_part"
        echo "----------------------------------------------------------"

        python3 ejecutar_sim_2.py "$NUM_PRUEBA" "$temp" "$n_part"

        STATUS=$?
        if [ $STATUS -eq 0 ]; then
            echo "⭐ Progreso registrado correctamente."
        else
            echo "⚠️  El proceso de Python falló, pero Bash continuará con la siguiente."
        fi

        CONTADOR=$(( CONTADOR + 1 ))
    done
done

echo -e "\n=========================================================="
echo "--- 🎉 Todas las simulaciones del script han terminado ---"
echo "=========================================================="