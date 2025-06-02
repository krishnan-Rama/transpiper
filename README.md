# transpiper

Transpipeline v4: Now supports denovo and genome guided transcriptome assembly 

For installation and deployment, follow: https://github.com/krishnan-Rama/transpipeline_containerised.git

## Installation

1. Install the transpipeline resources into your HPC cluster directory in which you will be performing the assembly:  

```bash
git clone https://github.com/krishnan-Rama/transpipeline_containerised.git
```

2. Put the raw reads in `raw_data` folder.  

3. Run the pipeline using `./deploy.sh`  

4. The prompt will ask you to (1) enter your preferred HPC partition name to submit the job (2) If you have a reference genome or not, (3) Path to the reference genome, and/or (4) the species/project identifier (e.g. Hsap or Hsap_200524 for _Homo sapiens_), simply type the name and return (it doesn't really matter).

> **Note:** 
>- You can run the pipeline multiple times simultaneously with different raw reads, simply repeat the installation process in a different directory and `./deploy` with a different species/project identifier name.
>- You can manually reconfigure Slurm parameters as per your HPC system (e.g memory, CPUs) by going through indivudal scripts in `modules` directory.  
>- All the relevent outputs will be stored in `outdir` folder, and outputs for every individual steps in the pipeline can be found in `workdir`.
---
