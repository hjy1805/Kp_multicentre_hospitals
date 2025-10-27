#!/bin/bash --login

#SBATCH --job-name=annotation
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-1469
#SBATCH --mem=16GB
#SBATCH --error=annotation.err
#SBATCH --output=annotation.out

module load prokka

run_acc=$(awk '{print $0}' ./genome_list | awk "NR == $SLURM_ARRAY_TASK_ID")

prokka ./assemblies/${run_acc}.fasta --genus Klebsiella --prefix ${run_acc} --force --outdir ./prokka_annotation
