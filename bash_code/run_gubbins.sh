#!/bin/bash --login

#SBATCH --job-name=gubbins
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=4
#SBATCH --mem=64GB
#SBATCH --error=gubbins.err
#SBATCH --output=gubbins.out

conda activate gubbins

baps=$(awk -F '[\t]' '{print $0}' ./st_list | awk "NR == $SLURM_ARRAY_TASK_ID")

cd ./${baps}

run_gubbins.py  -c 32 --prefix ${baps} ${baps}.aln
