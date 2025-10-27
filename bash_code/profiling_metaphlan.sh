#!/bin/bash

#SBATCH --job-name=profiling
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-141
#SBATCH --mem=64GB
#SBATCH --error=profiling.err
#SBATCH --output=profiling.out

module load metaphlan4/4.1.0

run_acc=$(awk '{print $0}' ./samples_b1 | awk "NR == $SLURM_ARRAY_TASK_ID")

mkdir ./${run_acc}/profiling_metaphlan
metaphlan ./${run_acc}/${run_acc}_1.fq.gz,./${run_acc}/${run_acc}_2.fq.gz --bowtie2db /ibex/project/c2205/databases/metaphlan_db --bowtie2out ./${run_acc}/profiling_metaphlan/metagenome.bowtie2.bz2 --nproc 16 --input_type fastq -o ./${run_acc}/profiling_metaphlan/profiled_metagenome.txt
