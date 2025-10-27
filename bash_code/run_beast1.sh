#!/bin/bash --login

#SBATCH --job-name=beast
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=12
#SBATCH --constrain=intel
#SBATCH --mem=128GB
#SBATCH --error=beast_%j.err
#SBATCH --output=beast_%j.out

module load beast

cluster=$(awk -F '[\t]' '{print $1}' ./st_list | awk "NR == $SLURM_ARRAY_TASK_ID")

java -jar $BEAST_JAR/beast.jar -overwrite -working -threads 32 ./${cluster}/beast/${cluster}_added.xml
