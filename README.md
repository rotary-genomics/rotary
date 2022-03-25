# nanopore-workflows
Snakemake workflows for Nanopore data analysis

## Contents
- `genome_longread`: workflow for genome assembly. Assembles long reads using Flye and polishes using short reads.

## Instructions
### `genome_longread`
1. Install
```bash
git clone https://github.com/jmtsuji/nanopore-workflows.git
conda env create -n genome_longread --file=nanopore-workflows/envs/conda_requirements.yaml
```
Also, you need to download the DFAST, EggNOG, and GTDB databases (TODO - not shown)

2. Copy and fill out the config (YAML) file (`genome_longread.yaml`). TODO - make this easier to use.

The key parameters to fill are the genome name and the paths to the long and (QC'ed) short reads.

Save as something like `genome_longread_mycopy.yaml`

3. Run
```bash
conda activate genome_longread

run_directory="E_coli" # Wherever you want to store the run files
config="genome_longread_mycopy.yaml"
conda_prefix="/Data/databases/nanopore-workflows/conda_envs" # Wherever you want to store the conda envs
snakefile="nanopore-workflows/snakemake/genome_longread.smk"
jobs=40

mkdir -p "${run_directory}"

snakemake --dryrun --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee genome_longread_steps.log

snakemake --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee genome_longread.log
```

Enjoy!