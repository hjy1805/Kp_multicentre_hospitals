#!/bin/bash --login

#SBATCH --job-name=unicycler
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --array=86
#SBATCH --mem=256GB
#SBATCH --error=assembly_unicycler.err
#SBATCH --output=assembly_unicycler.out

conda activate unicycler

run_acc=$(awk -F '[\t]' '{print $0}' ./tag_list | awk "NR == $SLURM_ARRAY_TASK_ID")

mkdir ${run_acc}
mkdir ./${run_acc}/hybrid_assembly
unicycler -1 ../fastqs/dir_${run_acc}/${run_acc}_1.fq.gz -2 ../fastqs/dir_${run_acc}/${run_acc}_2.fq.gz -l ./${run_acc}.fastq.gz -o ./${run_acc}/hybrid_assembly -t 64 --mode conservative
