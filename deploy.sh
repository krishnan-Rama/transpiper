#!/bin/bash

# Source config script
source config.parameters_all

# Prompt the user for the HPC partition name
read -p "Enter your preferred HPC partition name: " HPC_partition

# Function to replace <HPC_partition> in a given file
replace_partition() {
    sed -i "s/<HPC_partition>/${HPC_partition}/g" "$1"
}

# Replace in all relevant scripts
for script in "${moduledir}"/*.sh; do
    if grep -q "<HPC_partition>" "$script"; then
        replace_partition "$script"
    fi
done

# Prompt for genome-guided vs de novo assembly
read -p "Do you have a reference genome for genome-guided assembly? (yes/no): " use_reference

if [[ "$use_reference" == "yes" ]]; then
    read -p "Please provide the full path to the reference genome (FASTA file): " REFERENCE_GENOME
    if [[ ! -f "$REFERENCE_GENOME" ]]; then
        echo "ERROR: Reference genome file not found at $REFERENCE_GENOME"
        exit 1
    fi
    export ASSEMBLY_MODE="genome_guided"
    export REFERENCE_GENOME="$REFERENCE_GENOME"
else
    export ASSEMBLY_MODE="de_novo"
fi

# Prompt for species/project identifier
read -p "Please enter the species/project identifier name (e.g., Hsap_120624, Hsap for humans): " species_identifier

# Replace <pipeline> placeholders
replace_pipeline() {
    sed -i "s/<pipeline>/${species_identifier}/g" "$1"
}

for script in "${moduledir}"/*.sh; do
    if grep -q "<pipeline>" "$script"; then
        replace_pipeline "$script"
    fi
done

export SPECIES_IDENTIFIER="$species_identifier"
assembly="${SPECIES_IDENTIFIER}"
export assembly


# Step 0: Data configuration
sbatch -d singleton --error="${log}/rawqc_%J.err" --output="${log}/rawqc_%J.out" "${moduledir}/0-pre.sh"

# Step 1A: FastQC
sbatch -d singleton --error="${log}/rawqc_%J.err" --output="${log}/rawqc_%J.out" "${moduledir}/1A-fastqc_array.sh"

# Step 2A: Fastp trimming
sbatch -d singleton --error="${log}/fastp_%J.err" --output="${log}/fastp_%J.out" "${moduledir}/2A-fastp_array.sh"  

# Step 2B: Kraken2
sbatch -d singleton --error="${log}/kraken2_%J.err" --output="${log}/kraken2__%J.out" "${moduledir}/kraken.sh"
sbatch -d singleton --error="${log}/kraken2_%J.err" --output="${log}/kraken2__%J.out" "${moduledir}/2B-kraken2.sh"

# Step 2C: rcorrector
rcorrector_jobid=$(sbatch --parsable -d singleton --error="${log}/rcor_%J.err" --output="${log}/rcor__%J.out" "${moduledir}/2C-rcorrector.sh")

# Assembly step logic
if [[ "$ASSEMBLY_MODE" == "genome_guided" ]]; then
    echo "Preparing genome-guided assembly..."

    # Submit STAR job after rcorrector
    align_jobid=$(sbatch --parsable --dependency=afterok:${rcorrector_jobid} "${moduledir}/align_with_star.sh")

    if [[ -z "$align_jobid" ]]; then
        echo "❌ ERROR: STAR alignment job did not submit correctly!"
        exit 1
    fi

    echo "Submitting Trinity (genome-guided) with dependency on STAR..."
    sbatch --dependency=afterok:${align_jobid} --error="${log}/assembly_%J.err" --output="${log}/assembly_%J.out" "${moduledir}/3-trinity_assembly.sh"

else
    echo "Submitting Trinity (de novo) with dependency on rcorrector..."
    sbatch --dependency=afterok:${rcorrector_jobid} --error="${log}/assembly_%J.err" --output="${log}/assembly_%J.out" "${moduledir}/3-trinity_assembly.sh"
fi

# Step 1B: FastQC post-processed
sbatch -d singleton --error="${log}/rawqc_2_%J.err" --output="${log}/rawqc_2_%J.out" "${moduledir}/1B-fastqc_array.sh"

# Step 4: Evigene
sbatch -d singleton --error="${log}/evigene_%J.err" --output="${log}/evigene_%J.out" "${moduledir}/4-evigene.sh"
sbatch -d singleton --error="${log}/trinitystats_%J.err" --output="${log}/trinitystats_%J.out" "${moduledir}/4b-trinitystats.sh"

# Step 5: BUSCO
sbatch -d singleton --error="${log}/busco_%J.err" --output="${log}/busco_%J.out" "${moduledir}/5-busco_singularity.sh"

# Step 6: Trinity mapping
sbatch -d singleton --error="${log}/rsem_%J.err" --output="${log}/rsem_%J.out" --array=0-${max_idx} "${moduledir}/6-trinity-mapping.sh"

# Step 7: Diff expression
sbatch -d singleton --error="${log}/deg_%J.err" --output="${log}/deg_%J.out" "${moduledir}/7-rsem-post-reassemble.sh"

# Step 8: MultiQC
sbatch -d singleton --error="${log}/multiqc_%J.err" --output="${log}/multiqc_%J.out" "${moduledir}/8-multiqc.sh"

# Step 9: Blastdb
sbatch -d singleton --error="${log}/blastdb_%J.err" --output="${log}/blastdb_%J.out" "${moduledir}/9-blastdb.sh"

# Step 10: Multispecies blast
sbatch -d singleton --error="${log}/blastp_%J.err" --output="${log}/blastp_%J.out" --array="0-5" "${moduledir}/10-blast.sh"

# Step 11: Uniprot Annotation Extraction
sbatch -d singleton --error="${log}/upimapi_%J.err" --output="${log}/upimapi_%J.out" "${moduledir}/11-upimapi.sh"

# Step 12: Annotation Merge
sbatch -d singleton --error="${log}/merge_%J.err" --output="${log}/merge_%J.out" "${moduledir}/12-datamerge.sh"

# Step 13: Convert to DB
sbatch -d singleton --error="${log}/database_%J.err" --output="${log}/database_%J.out" "${moduledir}/13-csv2db.sh"

# Step 14: Clean-up (optional)
# sbatch -d singleton --error="${log}/cleanup_%J.err" --output="${log}/cleanup_%J.out" "${moduledir}/cleanup.sh"

