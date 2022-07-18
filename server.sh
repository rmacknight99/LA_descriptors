#!/bin/bash
#SBATCH -A CHE220011p
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH -n 4
#SBATCH --job-name LA_complexes 

module load anaconda3
conda activate eda

redis-server redis.config &

funsies wait
funsies clean
funsies worker &

smiles=`cat LAs/*smi`

for smi in $smiles;
do
    python workflow.py --smiles $smi --ntasks "4"
done

# Shutdown funsies
funsies shutdown

# Save result and quit
redis-cli --rdb results.rdb
redis-cli shutdown
