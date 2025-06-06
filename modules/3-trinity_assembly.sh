#!/bin/bash

#SBATCH --job-name=<pipeline>
#SBATCH --partition=<HPC_partition>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30000

echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE}

module load singularity/3.8.7

IMAGE_NAME=trinityrnaseq/trinityrnaseq:latest
SINGULARITY_IMAGE_NAME=trinityrnaseq:latest

if [ ! -f ${pipedir}/singularities/${SINGULARITY_IMAGE_NAME} ]; then
    echo "Pulling Trinity image..."
    singularity pull ${pipedir}/singularities/${SINGULARITY_IMAGE_NAME} docker://$IMAGE_NAME
else
    echo "✅ Singularity image already exists."
fi

SINGIMAGEDIR=${pipedir}/singularities
SINGIMAGENAME=${SINGULARITY_IMAGE_NAME}
WORKINGDIR=${pipedir}
export BINDS="${BINDS},${WORKINGDIR}:${WORKINGDIR},${rcordir}:${rcordir}"

# Convert SLURM memory to GB
TOTAL_RAM=$(expr ${SLURM_MEM_PER_NODE} / 1024)
CMD_SCRIPT="${log}/trinity_assembly_commands_${SLURM_JOB_ID}.sh"

cat > "$CMD_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail

echo "🧬 Running Trinity with ${TOTAL_RAM}G RAM and ${SLURM_CPUS_PER_TASK} CPUs"
ls -lh "${rcordir}/aligned.sorted.bam" "${rcordir}/aligned.sorted.bam.bai"

# Run Trinity
Trinity --genome_guided_bam "${rcordir}/aligned.sorted.bam" \\
        --genome_guided_max_intron 10000 \\
        --CPU ${SLURM_CPUS_PER_TASK} \\
        --max_memory ${TOTAL_RAM}G \\
        --output "${assemblydir}/" \\
        --full_cleanup

# Standardize output name for consistency
GG_FASTA="\${assemblydir}.Trinity-GG.fasta"
if [ -f "\$GG_FASTA" ]; then
    cp "\$GG_FASTA" "\${assemblydir}/Trinity.fasta"
    cp "\$GG_FASTA.gene_trans_map" "\${assemblydir}/Trinity.fasta.gene_trans_map"
else
    echo "❌ ERROR: Trinity genome-guided output not found."
    exit 1
fi

# Rename with species identifier
mv "\${assemblydir}/Trinity.fasta" "\${assemblydir}/${assembly}.fasta"
mv "\${assemblydir}/Trinity.fasta.gene_trans_map" "\${assemblydir}/${assembly}.gene_trans_map"

# Tag contigs with species name
sed -i "s/TRINITY_DN/${SPECIES_IDENTIFIER}_/g" "\${assemblydir}/${assembly}.fasta"

# Generate assembly stats
/usr/local/bin/util/TrinityStats.pl "\${assemblydir}/${assembly}.fasta" > "\${assemblydir}/${assembly}_stats.txt"

# Copy to output directory
mkdir -p "${outdir}/raw_assembly"
cp "\${assemblydir}/${assembly}.fasta" "\${outdir}/raw_assembly/"
cp "\${assemblydir}/${assembly}_stats.txt" "\${outdir}/raw_assembly/"
cp "\${assemblydir}/${assembly}.gene_trans_map" "\${outdir}/raw_assembly/"
EOF

# Run it inside the container
singularity exec --contain --bind ${BINDS} --pwd ${WORKINGDIR} ${SINGIMAGEDIR}/${SINGIMAGENAME} bash "$CMD_SCRIPT"

