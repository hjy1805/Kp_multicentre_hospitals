#!/bin/bash --login

#SBATCH --job-name=qc
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-211
#SBATCH --mem=8GB
#SBATCH --error=assembly_qc.err
#SBATCH --output=assembly_qc.out

module load quast

run_acc=$(awk -F '[\t]' '{print $0}' ./tag_list | awk "NR == $SLURM_ARRAY_TASK_ID")

mkdir ./${run_acc}/assembly_qc
quast.py ./${run_acc}/hybrid_assembly/assembly.fasta -o ./${run_acc}/assembly_qc
