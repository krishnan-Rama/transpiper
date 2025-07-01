#!/bin/bash

#SBATCH --job-name=GO_enrichment
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16000

module load singularity/3.8.7

# === USER-DEFINED PATHS ===
imgdir=${pipedir}/singularities
singimg=henryc101_clusterprofiler_4.10.0.sif
dockerimg=henryc101/clusterprofiler:4.10.0

# === MOUNT PATHS ===
workdir=${pipedir}/workdir
mergedir=${workdir}/mergedir/${assembly}_combined_final.csv
degfile=${workdir}/rsem/edgeR_results/${assembly}_RSEM.gene.counts.matrix.Nc_vs_Wb.edgeR.DE_results
outdir=${workdir}/GEA
script=${pipedir}/modules/go_enrichment_analysis.R

mkdir -p ${outdir}
mkdir -p ${imgdir}

# === Pull image if not exists ===
if [ ! -f ${imgdir}/${singimg} ]; then
    echo "Pulling Docker image and converting to Singularity..."
    singularity pull ${imgdir}/${singimg} docker://${dockerimg}
else
    echo "Singularity image already exists."
fi

# === RUN THE CONTAINER ===
singularity exec --contain \
  --bind ${workdir}:${workdir} \
  --bind ${script}:${script} \
  ${imgdir}/${singimg} \
  Rscript ${script} ${mergedir} ${degfile} ${outdir}
