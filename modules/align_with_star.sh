#!/bin/bash

#SBATCH --job-name=<pipeline>
#SBATCH --partition=<HPC_partition>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24000

set -euo pipefail

# Load STAR module instead of using Singularity
module load STAR/2.7.3a

# Settings
WORKINGDIR=${pipedir}
STAR_INDEX="${workdir}/star_index"

# Create STAR index directory
mkdir -p "$STAR_INDEX"

# Decompress genome if needed
if [[ "${REFERENCE_GENOME}" == *.gz ]]; then
    echo "Unzipping reference genome..."
    gunzip -c "$REFERENCE_GENOME" > "${STAR_INDEX}/reference.fa"
    GENOME="${STAR_INDEX}/reference.fa"
else
    GENOME="${REFERENCE_GENOME}"
fi

# Generate STAR index if not already present
if [ ! -f "${STAR_INDEX}/SA" ]; then
    echo "Generating STAR genome index..."
    STAR \
        --runThreadN ${SLURM_CPUS_PER_TASK} \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME"
else
    echo "Using existing STAR genome index at ${STAR_INDEX}"
fi

# Concatenate corrected reads
echo "Concatenating corrected reads..."
if compgen -G "${rcordir}/*_1.cor.fq.gz" > /dev/null && compgen -G "${rcordir}/*_2.cor.fq.gz" > /dev/null; then
    cat ${rcordir}/*_1.cor.fq.gz > ${rcordir}/all_1.fastq.gz
    cat ${rcordir}/*_2.cor.fq.gz > ${rcordir}/all_2.fastq.gz
else
    echo "ERROR: Corrected FASTQ files not found in ${rcordir}"
    exit 1
fi

# Align reads and output sorted BAM
echo "Aligning reads with STAR..."
STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn ${rcordir}/all_1.fastq.gz ${rcordir}/all_2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix "${rcordir}/aligned_" \
    --outSAMtype BAM SortedByCoordinate

# Rename the final BAM for Trinity
OUTBAM="${rcordir}/aligned_Aligned.sortedByCoord.out.bam"
if [ -f "$OUTBAM" ]; then
    mv "$OUTBAM" "${rcordir}/aligned.sorted.bam"
else
    echo "ERROR: STAR did not produce expected BAM output: $OUTBAM"
    exit 1
fi

