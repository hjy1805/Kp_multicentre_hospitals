#!/bin/bash --login

#SBATCH --job-name=kleborate
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-1469
#SBATCH --mem=8GB
#SBATCH --error=typing_kleborate.err
#SBATCH --output=typing_kleborate.out

conda activate kleborate

run_acc=$(awk '{print $0}' ./genome_list | awk "NR == $SLURM_ARRAY_TASK_ID")

kleborate -a ./assemblies/${run_acc}.fasta --all -o ./kleborate_report/${run_acc}_kleborate.txt
