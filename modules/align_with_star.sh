#!/bin/bash

#SBATCH --job-name=<pipeline>
#SBATCH --partition=<HPC_partition>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000

set -euo pipefail

# Load STAR module
module load STAR/2.7.3a

# Settings
WORKINGDIR="${pipedir}"
STAR_INDEX="${workdir}/star_index_clean"

# Create STAR index directory
mkdir -p "$STAR_INDEX"

# Clean and decompress reference genome if needed
if [[ "${REFERENCE_GENOME}" == *.gz ]]; then
    echo "🧼 Cleaning and decompressing reference genome..."
    zcat "${REFERENCE_GENOME}" | sed 's/^\(>[0-9A-Za-z]\+\) .*/\1/' > "${STAR_INDEX}/reference_clean.fa"
    GENOME="${STAR_INDEX}/reference_clean.fa"
else
    echo "🧼 Cleaning uncompressed reference genome..."
    sed 's/^\(>[0-9A-Za-z]\+\) .*/\1/' "${REFERENCE_GENOME}" > "${STAR_INDEX}/reference_clean.fa"
    GENOME="${STAR_INDEX}/reference_clean.fa"
fi

# Decompress GTF if provided and needed
if [[ -n "${GTF_FILE:-}" && -f "${GTF_FILE}" ]]; then
    if [[ "${GTF_FILE}" == *.gz ]]; then
        echo "🗜️    Decompressing GTF annotation..."
        gunzip -c "${GTF_FILE}" > "${STAR_INDEX}/annotation.gtf"
        GTF="${STAR_INDEX}/annotation.gtf"
    else
        GTF="${GTF_FILE}"
    fi
    SJDB_OPTS="--sjdbGTFfile ${GTF} --sjdbOverhang 99"
else
    SJDB_OPTS=""
fi

# Generate STAR index if not already present
if [ ! -f "${STAR_INDEX}/SA" ]; then
    echo "🧬 Generating STAR genome index..."
    STAR \
        --runThreadN ${SLURM_CPUS_PER_TASK} \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME" \
        $SJDB_OPTS
else
    echo "✅ Using existing STAR genome index at ${STAR_INDEX}"
fi

# Concatenate corrected reads
echo "📦 Concatenating corrected reads..."
if compgen -G "${rcordir}/*_1.cor.fq.gz" > /dev/null && compgen -G "${rcordir}/*_2.cor.fq.gz" > /dev/null; then
    cat ${rcordir}/*_1.cor.fq.gz > ${rcordir}/all_1.fastq.gz
    cat ${rcordir}/*_2.cor.fq.gz > ${rcordir}/all_2.fastq.gz
else
    echo "❌ ERROR: Corrected FASTQ files not found in ${rcordir}"
    exit 1
fi
# Run STAR alignment
echo "🚀 Aligning reads with STAR..."

STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "${rcordir}/all_1.fastq.gz" "${rcordir}/all_2.fastq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${rcordir}/aligned_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM MD \
    --outSAMprimaryFlag AllBestScore \
    --outSAMattrRGline ID:${SPECIES_IDENTIFIER} SM:${SPECIES_IDENTIFIER}

# Run STAR alignment
echo "🚀 Aligning reads with STAR..."

echo "🔧 Indexing BAM file..."
module load samtools/1.13
samtools index "${rcordir}/aligned.sorted.bam"

# Rename STAR output BAM for Trinity
OUTBAM="${rcordir}/aligned_Aligned.sortedByCoord.out.bam"
if [ -f "$OUTBAM" ]; then
    mv "$OUTBAM" "${rcordir}/aligned.sorted.bam"
    echo "✅ STAR alignment complete. Output: ${rcordir}/aligned.sorted.bam"
else
    echo "❌ ERROR: STAR did not produce the expected BAM output: $OUTBAM"
    exit 1
fi



