#!/bin/bash

#SBATCH --job-name=A_aura_test
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30000
#SBATCH --error=${log}/trinity_%J.err
#SBATCH --output=${log}/trinity_%J.out

echo "Environment Info:"
echo "=================="
echo "Hostname: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE} MB"

# Load Singularity
module load singularity/3.8.7

# Trinity container settings
IMAGE_NAME=trinityrnaseq/trinityrnaseq:latest
SINGULARITY_IMAGE_NAME=trinityrnaseq:latest
SINGIMAGEDIR=${pipedir}/singularities
SINGIMAGENAME=${SINGULARITY_IMAGE_NAME}
WORKINGDIR=${pipedir}

# Pull image if not found
if [ ! -f ${SINGIMAGEDIR}/${SINGULARITY_IMAGE_NAME} ]; then
    echo "Pulling Trinity Singularity image..."
    singularity pull ${SINGIMAGEDIR}/${SINGULARITY_IMAGE_NAME} docker://$IMAGE_NAME
else
    echo "Singularity image already exists."
fi

# RAM in GB
TOTAL_RAM=$(expr ${SLURM_MEM_PER_NODE} / 1024)

# Bind directories
export BINDS="${BINDS},${WORKINGDIR}:${WORKINGDIR}"

# Create command script
CMD_SCRIPT="${log}/trinity_assembly_commands_${SLURM_JOB_ID}.sh"
echo "Writing Trinity command script to $CMD_SCRIPT"

cat > "$CMD_SCRIPT" <<EOF

set -e

if [[ "\$ASSEMBLY_MODE" == "genome_guided" ]]; then
    echo "Running Trinity in GENOME-GUIDED mode..."
    
    Trinity --genome_guided_bam "${rcordir}/aligned.sorted.bam" \\
            --genome_guided_max_intron 10000 \\
            --CPU ${SLURM_CPUS_PER_TASK} \\
            --max_memory ${TOTAL_RAM}G \\
            --output "${assemblydir}/" \\
            --full_cleanup

else
    echo "Running Trinity in DE NOVO mode..."

    # Merge corrected reads if needed
    cat ${rcordir}/*_1.cor.fq.gz > ${rcordir}/all_1.fastq.gz
    cat ${rcordir}/*_2.cor.fq.gz > ${rcordir}/all_2.fastq.gz

    Trinity --seqType fq \\
            --left "${rcordir}/all_1.fastq.gz" \\
            --right "${rcordir}/all_2.fastq.gz" \\
            --max_memory ${TOTAL_RAM}G \\
            --CPU ${SLURM_CPUS_PER_TASK} \\
            --output "${assemblydir}/" \\
            --full_cleanup
fi

# Move and rename outputs
mv "${assemblydir}.Trinity.fasta" "${assemblydir}/${assembly}.fasta"
mv "${assemblydir}.Trinity.fasta.gene_trans_map" "${assemblydir}/${assembly}.gene_trans_map"

# Rename contigs for consistency
sed -i "s/TRINITY_DN/${SPECIES_IDENTIFIER}_/g" "${assemblydir}/${assembly}.fasta"

# Trinity Stats
/usr/local/bin/util/TrinityStats.pl "${assemblydir}/${assembly}.fasta" > "${assemblydir}/${assembly}_stats.txt"

# Copy outputs to outdir
mkdir -p "${outdir}/raw_assembly"
cp "${assemblydir}/${assembly}.fasta" "${outdir}/raw_assembly/"
cp "${assemblydir}/${assembly}_stats.txt" "${outdir}/raw_assembly/"
cp "${assemblydir}/${assembly}.gene_trans_map" "${outdir}/raw_assembly/"

echo "Assembly completed successfully."
EOF

# Run the script inside the Singularity container
singularity exec --contain --bind ${BINDS} --pwd ${WORKINGDIR} ${SINGIMAGEDIR}/${SINGIMAGENAME} bash "$CMD_SCRIPT"

