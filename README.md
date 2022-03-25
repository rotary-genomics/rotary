# nanopore-workflows
Snakemake workflows for Nanopore data analysis

## Contents
- `genome_longread`: workflow for genome assembly. Assembles long reads using Flye and polishes using short reads.

## Instructions
### `genome_longread`
1. Install
```bash
git clone https://github.com/jmtsuji/nanopore-workflows.git
cd nanopore-workflows
conda env create -n genome_longread --file=envs/conda_requirements.yaml
```
Also, you need to download the DFAST, EggNOG, and GTDB databases (TODO - not shown)

2. Copy and fill out the config (YAML) file (`genome_longread.yaml`). TODO - make this easier to use.

The key parameters to fill are the genome name and the paths to the long and (QC'ed) short reads.

3. Run
```bash
conda activate genome_longread

run_directory="." # Wherever you want to store the run files
conda_prefix="conda" # Wherever you want to store the conda envs
jobs=40

snakemake --snakefile snakemake/genome_longread.smk --configfile config.yaml --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee genome_longread.log
```

Enjoy!