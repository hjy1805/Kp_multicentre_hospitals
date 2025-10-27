#!/bin/bash --login

#SBATCH --job-name=call
#SBATCH --time=60:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=4
#SBATCH --mem=64GB
#SBATCH --error=call.err
#SBATCH --output=call.out

conda activate snpsites

acc=$(awk -F '[\t]' '{print $0}' ./st_list | awk "NR == $SLURM_ARRAY_TASK_ID")

cd ${acc}
snp-sites -m -v -p -o ${acc} ${acc}.aln
