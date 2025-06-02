#!/bin/bash

#SBATCH --job-name=<pipeline>
#SBATCH --partition=<HPC_partition>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=64000
#SBATCH --error=${log}/kraken2_%J.err
#SBATCH --output=${log}/kraken2_%J.out

echo "Some Usable Environment Variables:"
echo "=================================="
echo "hostname=$(hostname)"
echo "\$SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "\$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"

# Write jobscript to output file (for reproducibility)
cat $0

module load singularity/3.8.7

# Set variables
WORKINGDIR=${pipedir}
export BINDS="${BINDS},${WORKINGDIR}:${WORKINGDIR}"

KRAKEN_IMAGE=kraken2:2.1.3--pl5321hdcf5f25_0
KRONA_IMAGE=biocontainers/krona:2.8.1--pl5321hdfd78af_1

# Pull Singularity images if missing
if [ ! -f ${pipedir}/singularities/${KRAKEN_IMAGE} ]; then
    wget -O ${pipedir}/singularities/${KRAKEN_IMAGE} https://depot.galaxyproject.org/singularity/${KRAKEN_IMAGE}
fi

if [ ! -f ${pipedir}/singularities/krona_2.8.1.sif ]; then
    singularity pull ${pipedir}/singularities/krona_2.8.1.sif docker://${KRONA_IMAGE}
fi

SINGIMAGEDIR=${pipedir}/singularities
KRAKEN_IMG=${SINGIMAGEDIR}/${KRAKEN_IMAGE}
KRONA_IMG=${SINGIMAGEDIR}/krona_2.8.1.sif

# Create Kraken DB dir if needed
mkdir -p "${pipedir}/kraken_standard/"
if [ ! -f ${pipedir}/kraken_standard/k2_standard_20230605 ]; then
    wget -O ${pipedir}/kraken_standard/kraken_standard.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz
    tar -xzvf ${pipedir}/kraken_standard/kraken_standard.tar.gz -C ${pipedir}/kraken_standard/
fi

# List of sample bases
bases=$(ls ${trimdir}/*_trim_1.fastq.gz | xargs -n 1 basename | sed 's/_trim_1.fastq.gz//' | sort | uniq)

for base in $bases; do
    echo "Processing sample: $base"

    # Kraken2 command
    cat >${log}/kraken_taxa_commands_${SLURM_JOB_ID}_${base}.sh <<EOF
kraken2 --paired --db ${pipedir}/kraken_standard \\
  --output ${krakendir}/${base}_kraken2_output \\
  --report ${krakendir}/${base}_kraken2_report \\
  --classified-out ${krakendir}/${base}_#.classified.fastq \\
  --unclassified-out ${krakendir}/${base}_#.unclassified.fastq \\
  --threads ${SLURM_CPUS_PER_TASK} \\
  ${trimdir}/${base}_trim_1.fastq.gz ${trimdir}/${base}_trim_2.fastq.gz
EOF

    echo "Running Kraken2 for $base..."
    singularity exec --contain --bind ${BINDS} --pwd ${WORKINGDIR} ${KRAKEN_IMG} bash ${log}/kraken_taxa_commands_${SLURM_JOB_ID}_${base}.sh

    # Krona visualization
    echo "Generating Krona HTML for $base..."
    singularity exec --contain --bind ${BINDS} --pwd ${WORKINGDIR} ${KRONA_IMG} ktImportText ${krakendir}/${base}_kraken2_report -o ${krakendir}/${base}_krona.html
done

