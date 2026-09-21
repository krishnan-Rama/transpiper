# transpiper

**Transpipeline v4** — Now supports both **de novo** and **genome-guided** transcriptome assembly, with differential gene expression and enrichment analysiS on SLURM-based HPC clusters.

For additional post-annotation data query and analysis, refer to:
👉 [https://github.com/krishnan-Rama/transpipeline\_containerised.git](https://github.com/krishnan-Rama/transpipeline_containerised.git)

---

## ✅ Installation

1. Clone the repository into your working directory on the HPC system:

```bash
git clone https://github.com/krishnan-Rama/transpiper.git
```
```bash
cd transpiper
```

2. Place your raw reads in the `raw_data/` folder.

3. Run the pipeline:

```bash
./deploy.sh -p <HPC_partition> -n <project_name> -r /path/to/reference.fasta -g /path/to/annotation.gtf
```

#### Run Options

#### Required

* `-p`: SLURM partition to submit jobs (e.g., `short`, `cpu`)
* `-n`: Project name or species identifier (e.g., `Hsap`, `mouse2025`)

#### Optional

* `-r`: Reference genome file (`.fasta`)
* `-g`: Genome annotation file (`.gtf`)

> If both `-r` and `-g` are provided and valid, the pipeline will run in **genome-guided** mode.
> Otherwise, it will default to **de novo** assembly.

---

## 🔁 Reusability

You can run multiple independent projects by cloning this repo into different directories and providing separate `raw_data/` and identifiers:

```bash
git clone https://github.com/krishnan-Rama/transpiper.git my_project_A
cd my_project_A
./deploy_pipeline.sh -p <partition> -n my_project_A
```

---

## 🗂 Output Structure

| Folder     | Description                                            |
| ---------- | ------------------------------------------------------ |
| `log/`     | SLURM `.out` and `.err` job logs                       |
| `workdir/` | Intermediate files for each processing stage           |
| `outdir/`  | Final pipeline outputs (assemblies, annotations, etc.) |
| `modules/` | All job submission scripts used by the pipeline        |

---

## ⚙️ Customizing SLURM Settings

Each pipeline step is modularized in `modules/*.sh`.
To change resources like memory, CPUs, or time limits, edit the respective scripts directly.

---

## 🧬 Workflow Overview

1. Quality control: `FastQC`, `fastp`, `Kraken2`
2. Error correction: `Rcorrector`
3. Transcriptome assembly: `Trinity` (de novo or genome-guided via `STAR`)
4. Evaluation: `BUSCO`, `TrinityStats`, `EviGene`
5. Quantification: `RSEM`, DE analysis
6. Functional annotation: `BLAST`, `UPIMAPI`
7. Summary report: `MultiQC`
8. Final merge and SQLite database creation
9. Gene Enrichment Analysis

---

## transpiper Architecture

```mermaid
flowchart TD

subgraph group_orchestration["Orchestration"]
  node_deploy["Pipeline Launcher<br/>[deploy.sh]"]
  node_slurm["SLURM Cluster"]
end

subgraph group_preprocessing["Read Processing"]
  node_raw_reads["Raw Reads"]
  node_quality_control["Quality Control<br/>[1A-fastqc_array.sh]"]
  node_contamination["Taxonomic Screening<br/>[2B-kraken2.sh]"]
  node_error_correction["Error Correction<br/>[2C-rcorrector.sh]"]
end

subgraph group_assembly["Assembly Analysis"]
  node_star_guidance["Genome Guidance<br/>[align_with_star.sh]"]
  node_trinity["Trinity Assembly"]
  node_evaluation["Assembly Evaluation"]
  node_quantification["Expression Quantification"]
  node_differential_expression["Differential Expression<br/>[process_rna_seq.R]"]
end

subgraph group_annotation["Annotation Results"]
  node_blast["BLAST Annotation<br/>[10-blast.sh]"]
  node_upimapi["UPIMAPI Annotation<br/>[11-upimapi.sh]"]
  node_reporting["MultiQC Reporting<br/>[8-multiqc.sh]"]
  node_data_merge["Result Merge<br/>[12-datamerge.sh]"]
end

subgraph group_query["Data Query"]
  node_csv_database["CSV Database Load<br/>[13-csv2db.sh]"]
  node_gene_database[("Gene SQLite Store")]
  node_gene_query["Gene Query CLI<br/>[query_gene_data.py]"]
  node_enrichment["GO Enrichment"]
  node_merged_csv["Merged Gene CSV"]
end

node_researcher(("Researcher"))

node_researcher -->|"submits project"| node_deploy
node_deploy -->|"submits jobs"| node_slurm
node_raw_reads -->|"processes reads"| node_quality_control
node_quality_control -->|"screens reads"| node_contamination
node_quality_control -->|"passes reads"| node_error_correction
node_contamination -->|"filters reads"| node_error_correction
node_error_correction -->|"corrected reads"| node_trinity
node_star_guidance -.->|"guides assembly"| node_trinity
node_trinity -->|"evaluates assembly"| node_evaluation
node_trinity -->|"quantifies transcripts"| node_quantification
node_quantification -->|"analyzes expression"| node_differential_expression
node_trinity -->|"annotates transcripts"| node_blast
node_trinity -->|"annotates proteins"| node_upimapi
node_quality_control -->|"summarizes QC"| node_reporting
node_evaluation -->|"summarizes metrics"| node_reporting
node_data_merge -->|"writes results"| node_merged_csv
node_blast -->|"adds annotations"| node_data_merge
node_upimapi -->|"adds annotations"| node_data_merge
node_differential_expression -->|"supplies results"| node_enrichment
node_merged_csv -->|"loads CSV"| node_csv_database
node_csv_database -->|"writes database"| node_gene_database
node_researcher -->|"queries genes"| node_gene_query
node_gene_query -->|"reads filters"| node_gene_database
node_gene_query -.->|"creates store"| node_gene_database
node_merged_csv -.->|"creates database"| node_gene_query

click node_deploy "https://github.com/krishnan-rama/transpiper/blob/main/deploy.sh"
click node_raw_reads "https://github.com/krishnan-rama/transpiper/tree/main/raw_data"
click node_quality_control "https://github.com/krishnan-rama/transpiper/blob/main/modules/1A-fastqc_array.sh"
click node_contamination "https://github.com/krishnan-rama/transpiper/blob/main/modules/2B-kraken2.sh"
click node_error_correction "https://github.com/krishnan-rama/transpiper/blob/main/modules/2C-rcorrector.sh"
click node_star_guidance "https://github.com/krishnan-rama/transpiper/blob/main/modules/align_with_star.sh"
click node_trinity "https://github.com/krishnan-rama/transpiper/blob/main/modules/3-trinity_assembly.sh"
click node_evaluation "https://github.com/krishnan-rama/transpiper/blob/main/modules/5-busco_singularity.sh"
click node_quantification "https://github.com/krishnan-rama/transpiper/blob/main/modules/7-rsem-post-reassemble.sh"
click node_differential_expression "https://github.com/krishnan-rama/transpiper/blob/main/modules/process_rna_seq.R"
click node_blast "https://github.com/krishnan-rama/transpiper/blob/main/modules/10-blast.sh"
click node_upimapi "https://github.com/krishnan-rama/transpiper/blob/main/modules/11-upimapi.sh"
click node_reporting "https://github.com/krishnan-rama/transpiper/blob/main/modules/8-multiqc.sh"
click node_data_merge "https://github.com/krishnan-rama/transpiper/blob/main/modules/12-datamerge.sh"
click node_csv_database "https://github.com/krishnan-rama/transpiper/blob/main/modules/13-csv2db.sh"
click node_gene_query "https://github.com/krishnan-rama/transpiper/blob/main/modules/query_gene_data.py"
click node_enrichment "https://github.com/krishnan-rama/transpiper/blob/main/modules/go_enrichment_analysis.R"

classDef toneNeutral fill:#f8fafc,stroke:#334155,stroke-width:1.5px,color:#0f172a
classDef toneBlue fill:#dbeafe,stroke:#2563eb,stroke-width:1.5px,color:#172554
classDef toneAmber fill:#fef3c7,stroke:#d97706,stroke-width:1.5px,color:#78350f
classDef toneMint fill:#dcfce7,stroke:#16a34a,stroke-width:1.5px,color:#14532d
classDef toneRose fill:#ffe4e6,stroke:#e11d48,stroke-width:1.5px,color:#881337
classDef toneIndigo fill:#e0e7ff,stroke:#4f46e5,stroke-width:1.5px,color:#312e81
classDef toneTeal fill:#ccfbf1,stroke:#0f766e,stroke-width:1.5px,color:#134e4a
class node_deploy,node_slurm toneBlue
class node_raw_reads,node_quality_control,node_contamination,node_error_correction toneAmber
class node_star_guidance,node_trinity,node_evaluation,node_quantification,node_differential_expression toneMint
class node_blast,node_upimapi,node_reporting,node_data_merge toneRose
class node_csv_database,node_gene_database,node_gene_query,node_enrichment,node_merged_csv,node_researcher toneIndigo
```
