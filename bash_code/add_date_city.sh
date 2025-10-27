#!/bin/bash

#SBATCH --job-name=add
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=2-17
#SBATCH --mem=4GB
#SBATCH --error=add.err
#SBATCH --output=add.out

baps=$(awk -F '[\t]' '{print $0}' ./st_list | awk "NR == $SLURM_ARRAY_TASK_ID")

fasta_file=./${baps}/${baps}.filtered_polymorphic_sites.fasta
tsv_file=./all_beast_meta.tsv
output_file=./${baps}/${baps}_added.fasta
temp_file=./${baps}/temp.fasta

# Loop through each line in the TSV file
awk 'NR > 1' "${tsv_file}" | while IFS=$'\t' read -r fasta_tag year month day location rest; do
    # Check if the FASTA tag exists in the FASTA file
    if grep -q -w -m 1 "^>${fasta_tag}$" "${fasta_file}"; then
        # Extract the sequence from the FASTA file for the current tag
        sequence=$(awk -v tag=">${fasta_tag}" '$0 ~ tag {flag=1; next} /^>/ {flag=0} flag' "${fasta_file}")

        # Append additional information to the header
        new_header=">${fasta_tag}|${year}-${month}-${day}|${location}"

        # Write the modified header and sequence to the temporary file
        echo -e "${new_header}\n${sequence}" >> "${temp_file}"
    fi
done < "${tsv_file}"

# Replace the original FASTA file with the temporary file
mv "${temp_file}" "${output_file}"

echo "Finished. Modified FASTA file saved to ${output_file}"
