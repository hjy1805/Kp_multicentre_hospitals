#!/bin/bash --login

#SBATCH --job-name=mapping
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=3-748
#SBATCH --mem=16GB
#SBATCH --error=mapping.err
#SBATCH --output=mapping.out

conda activate snippy

run_acc=$(awk -F '[\t]' '{print $1}' ./internal_meta_select.tsv | awk "NR == $SLURM_ARRAY_TASK_ID")
st=$(awk -F '[\t]' '{print $7}' ./internal_meta_select.tsv | awk "NR == $SLURM_ARRAY_TASK_ID")
ref=$(awk -F '[\t]' '{print $8}' ./internal_meta_select.tsv | awk "NR == $SLURM_ARRAY_TASK_ID")

if [[ "$run_acc" =~ ^(ERR|Kp) ]]; then
    mkdir ./${st}/${run_acc}_mapping_local_st_cluster
    snippy --cpus 16 --force --outdir ./${st}/${run_acc}_mapping_local_st_cluster --ref ./mapping_ref/${ref}_ref.fasta --R1 /ibex/project/c2205/jiayi/Kp_ml/sharif_samples/${run_acc}/${run_acc}_1.fastq.gz --R2 /ibex/project/c2205/jiayi/Kp_ml/sharif_samples/${run_acc}/${run_acc}_2.fastq.gz
    rm -rf ./${st}/${run_acc}/mapping_local_st_cluster/snps.bam
else
    mkdir ./${st}/${run_acc}_mapping_local_st_cluster
    snippy --cpus 16 --force --outdir ./${st}/${run_acc}_mapping_local_st_cluster --ref ./mapping_ref/${ref}_ref.fasta --R1 /ibex/project/c2205/Malaikah/Kleb_800_project/${run_acc}/${run_acc}_1.fq.gz --R2 /ibex/project/c2205/Malaikah/Kleb_800_project/${run_acc}/${run_acc}_2.fq.gz
    rm -rf ./${st}/${run_acc}_mapping_local_st_cluster/snps.bam
fi
