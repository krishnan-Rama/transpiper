#!/bin/bash

# Source configs
source config.parameters_all

# Enhanced Usage function
usage() {
    cat <<EOF
Usage: $0 -p <HPC_partition> -n <project_name> [-r <reference_genome.fasta>] [-g <annotation.gtf>] [--trinity-only] [--help]

Options:
  -p <HPC_partition>    Specify the HPC partition to use (required)
  -n <project_name>     Specify the project name (species identifier) (required)
  -r <reference.fasta>  Specify the reference genome FASTA file (for genome-guided assembly)
  -g <annotation.gtf>   Specify the annotation GTF file (for genome-guided assembly)
  --trinity-only        Run only the Trinity assembly step (genome-guided) and exit (for debugging)
  --help, -h            Show this help message and exit
EOF
}

# Initialize flags
TRINITY_ONLY=0
HELP=0

# Process long flags and help option
new_args=()
for arg in "$@"; do
    case "$arg" in
        --trinity-only)
            TRINITY_ONLY=1
            ;;
        --help)
            HELP=1
            ;;
        *)
            new_args+=("$arg")
            ;;
    esac
done
set -- "${new_args[@]}"

# Parse short flags (including -h for help)
while getopts ":p:n:r:g:h" opt; do
    case ${opt} in
        p ) HPC_partition="$OPTARG" ;;
        n ) species_identifier="$OPTARG" ;;
        r ) REFERENCE_GENOME="$OPTARG" ;;
        g ) GTF_FILE="$OPTARG" ;;
        h ) HELP=1 ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

# Handle help request
if [[ $HELP -eq 1 ]]; then
    usage
    exit 0
fi

# Validate required inputs
if [[ -z "${HPC_partition:-}" || -z "${species_identifier:-}" ]]; then
    echo "❌ ERROR: Both -p (HPC partition) and -n (project name) are required."
    usage
    exit 1
fi


# Set genome-guided mode only if both files exist
if [[ -n "${REFERENCE_GENOME:-}" && -f "$REFERENCE_GENOME" && -n "${GTF_FILE:-}" && -f "$GTF_FILE" ]]; then
    ASSEMBLY_MODE="genome_guided"
else
    echo "⚠️  Missing valid reference genome or GTF — defaulting to de novo assembly."
    ASSEMBLY_MODE="de_novo"
fi

# Replace placeholders in module scripts
for script in "${moduledir}"/*.sh; do
    sed -i "s|<HPC_partition>|${HPC_partition}|g" "$script"
    sed -i "s|<pipeline>|${species_identifier}|g" "$script"
done

# Export shared variables
export ASSEMBLY_MODE
export REFERENCE_GENOME
export GTF_FILE
export SPECIES_IDENTIFIER="$species_identifier"
export assembly="$species_identifier"

#-----------------Trinity-Only---------------
if [[ "$TRINITY_ONLY" -eq 1 ]]; then
  echo "🔍 Debug mode: Submitting Trinity genome-guided only..."
  export ASSEMBLY_MODE=genome_guided
  sbatch --error="${log}/trinity_debug.err" \
         --output="${log}/trinity_debug.out" \
         --export=ALL,ASSEMBLY_MODE=genome_guided,REFERENCE_GENOME="${REFERENCE_GENOME}",GTF_FILE="${GTF_FILE}",pipedir="${pipedir}",log="${log}",rcordir="${rcordir}",assemblydir="${assemblydir}",SPECIES_IDENTIFIER="${species_identifier}" \
         "${moduledir}/3-trinity_assembly.sh"
  exit 0
fi
#-----------------Trinity-Only---------------


#-------------- Initial Jobs--------------------
sbatch -d singleton --error="${log}/rawqc_%J.err" --output="${log}/rawqc_%J.out" "${moduledir}/0-pre.sh"
sbatch -d singleton --error="${log}/rawqc_%J.err" --output="${log}/rawqc_%J.out" "${moduledir}/1A-fastqc_array.sh"
sbatch -d singleton --error="${log}/fastp_%J.err" --output="${log}/fastp_%J.out" "${moduledir}/2A-fastp_array.sh"
sbatch -d singleton --error="${log}/kraken2_%J.err" --output="${log}/kraken2_%J.out" "${moduledir}/kraken.sh"
sbatch -d singleton --error="${log}/kraken2_%J.err" --output="${log}/kraken2_%J.out" "${moduledir}/2B-kraken2.sh"

# Run rcorrector
rcorrector_script="${moduledir}/2C-rcorrector.sh"
cor_files_found=$(find "${rcordir}" -maxdepth 1 -type f -name "*_1.cor.fq.gz" | wc -l)

if (( cor_files_found > 0 )); then
    echo "✅ Detected existing .cor.fq.gz files — assuming rcorrector already completed."
    rcorrector_jobid="manual"
else
    if [[ ! -f "$rcorrector_script" ]]; then
        echo "❌ ERROR: rcorrector script not found at $rcorrector_script"
        exit 1
    fi
    rcorrector_jobid=$(sbatch --parsable -d singleton \
        --error="${log}/rcor_%J.err" \
        --output="${log}/rcor__%J.out" \
        --export=ALL,pipedir="${pipedir}",log="${log}",rawdir="${rawdir}",rcordir="${rcordir}",krakendir="${krakendir}" \
        "$rcorrector_script")
fi

# Assembly logic
if [[ "$ASSEMBLY_MODE" == "genome_guided" ]]; then
    echo "🧬 Running genome-guided assembly..."

    if [[ "$rcorrector_jobid" == "manual" ]]; then
        echo "📌 Submitting STAR alignment without dependency..."
        align_jobid=$(sbatch --parsable \
            --error="${log}/star_align_%J.err" \
            --output="${log}/star_align_%J.out" \
            --export=ALL,REFERENCE_GENOME="${REFERENCE_GENOME}",GTF_FILE="${GTF_FILE}",pipedir="${pipedir}",log="${log}",rcordir="${rcordir}",workdir="${workdir}" \
            "${moduledir}/align_with_star.sh")
    else
        echo "📌 Submitting STAR alignment with dependency on rcorrector: $rcorrector_jobid"
        align_jobid=$(sbatch --parsable \
            --dependency=afterok:${rcorrector_jobid} \
            --error="${log}/star_align_%J.err" \
            --output="${log}/star_align_%J.out" \
            --export=ALL,REFERENCE_GENOME="${REFERENCE_GENOME}",GTF_FILE="${GTF_FILE}",pipedir="${pipedir}",log="${log}",rcordir="${rcordir}",workdir="${workdir}" \
            "${moduledir}/align_with_star.sh")
    fi

    if [[ -n "$align_jobid" ]]; then
        echo "🚀 Submitting Trinity (genome-guided) with dependency on STAR..."
        sbatch --dependency=afterok:${align_jobid} \
            --error="${log}/assembly_%J.err" \
            --output="${log}/assembly_%J.out" \
            --export=ALL,pipedir="${pipedir}",log="${log}",rcordir="${rcordir}",assemblydir="${assemblydir}" \
            "${moduledir}/3-trinity_assembly.sh"
    else
        echo "❌ ERROR: STAR alignment failed to return a job ID."
        exit 1
    fi

else
    echo "🧬 Running de novo assembly..."
    if [[ "$rcorrector_jobid" == "manual" ]]; then
        sbatch \
            --error="${log}/assembly_%J.err" \
            --output="${log}/assembly_%J.out" \
            --export=ALL,pipedir="${pipedir}",log="${log}",rcordir="${rcordir}",assemblydir="${assemblydir}" \
            "${moduledir}/3-trinity_assembly.sh"
    else
        sbatch \
            --dependency=afterok:${rcorrector_jobid} \
            --error="${log}/assembly_%J.err" \
            --output="${log}/assembly_%J.out" \
            --export=ALL,pipedir="${pipedir}",log="${log}",rcordir="${rcordir}",assemblydir="${assemblydir}" \
            "${moduledir}/3-trinity_assembly.sh"
    fi
fi


# Downstream jobs
sbatch -d singleton --error="${log}/rawqc_2_%J.err" --output="${log}/rawqc_2_%J.out" "${moduledir}/1B-fastqc_array.sh"
sbatch -d singleton --error="${log}/evigene_%J.err" --output="${log}/evigene_%J.out" "${moduledir}/4-evigene.sh"
sbatch -d singleton --error="${log}/trinitystats_%J.err" --output="${log}/trinitystats_%J.out" "${moduledir}/4b-trinitystats.sh"
sbatch -d singleton --error="${log}/busco_%J.err" --output="${log}/busco_%J.out" "${moduledir}/5-busco_singularity.sh"
sbatch -d singleton --error="${log}/rsem_%J.err" --output="${log}/rsem_%J.out" --array=0-${max_idx} "${moduledir}/6-trinity-mapping.sh"
sbatch -d singleton --error="${log}/deg_%J.err" --output="${log}/deg_%J.out" "${moduledir}/7-rsem-post-reassemble.sh"
sbatch -d singleton --error="${log}/multiqc_%J.err" --output="${log}/multiqc_%J.out" "${moduledir}/8-multiqc.sh"
sbatch -d singleton --error="${log}/blastdb_%J.err" --output="${log}/blastdb_%J.out" "${moduledir}/9-blastdb.sh"
sbatch -d singleton --error="${log}/blastp_%J.err" --output="${log}/blastp_%J.out" --array="0-5" "${moduledir}/10-blast.sh"
sbatch -d singleton --error="${log}/upimapi_%J.err" --output="${log}/upimapi_%J.out" "${moduledir}/11-upimapi.sh"
sbatch -d singleton --error="${log}/merge_%J.err" --output="${log}/merge_%J.out" "${moduledir}/12-datamerge.sh"
sbatch -d singleton --error="${log}/database_%J.err" --output="${log}/database_%J.out" "${moduledir}/13-csv2db.sh"
sbatch -d singleton --error="${log}/GEA_%J.err" --output="${log}/GEA_%J.out" "${moduledir}/14-go-enrichment.sh"
