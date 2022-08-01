# rotary
Assembly/annotation workflow for Nanopore genome data

## Quick start
```bash
git clone https://github.com/jmtsuji/rotary.git

conda env create -n rotary --file=rotary/envs/rotary.yaml

cp rotary/config.yaml myconfig.yaml # Edit this file, especially the first few lines

conda activate rotary

mkdir -p output_dir conda_envs

snakemake --snakefile rotary/rules/rotary.yaml \
  --configfile myconfig.yaml \
  --directory output_dir \
  --conda-prefix conda_envs \
  --jobs 10 \
  --use-conda \
  --conda-frontend mamba \
  --rerun-incomplete \
  --reason \
  --printshellcmds 2>&1 | \
  tee rotary.log
```

## Description
_rotary_ is a snakemake pipeline that can be used to assemble single microbial genomes using 
standlone Nanopore data<sup>[1](#Footnotes)</sup> or hybrid Nanopore + short read data. 
The pipeline performs long read QC, assembly, end repair, polishing, contig rotation, and genome annotation.

### Some advantages of using _rotary_:
- All databases auto-install, so you can start analyzing genomes reproducibly with limited effort!
- Snakemake checkpointing allows you to re-start a failed run from where you left off
- Circularization is handled fairly carefully. Unlike the defaults in most pipelines, _rotary_ fixes the 
  [short gap region](https://github.com/fenderglass/Flye/issues/315#issuecomment-720679812) that can occur at the ends 
  of circular contigs produced by Flye. It also polishes the circular contigs in two 
  different rotation states to try to correct errors near contig ends.

## Requirements
- OS: Runs on Linux (tested on Ubuntu 20.04)
- Software: requires `miniconda`
- Resources: The majority of the pipeline is not too resource intensive. The limiting factors are:
  - Flye requires moderate RAM (e.g., < 64 GB for typical bacterial genome runs)
  - GTDB-Tk v2 (with GTDB r207) needs ~55 GB RAM
  - EggNOG-mapper is a bit slow (e.g., ~40 minutes on 40 CPU threads)

## Usage
### 1. Install
```bash
git clone https://github.com/jmtsuji/rotary.git
conda env create -n rotary --file=rotary/envs/rotary.yaml
```

### 2. Copy and fill out the config (YAML) file (`config.yaml`).

The basic parameters to fill are the genome name, the paths to the long and (QC'ed) short reads, and the DB dir.

Other very important parameters:
- `flye_input_mode`: set to "nano-hq" if your reads have <= 5% error rate; otherwise use "nano-raw"
- `medaka_model`: set this parameter to match the flow cell version / basecalling model you used to generate the long reads. See details in the [Medaka Github repo's "Models" section](https://github.com/nanoporetech/medaka#models)
- "Post-polishing contig filter" section: set the min depth and coverage you would like for a contig to be kept (more info in that section of the YAML file)

Save as something like `myconfig.yaml`.

**Note** that short read QC is not performed, so you'll want to do short read QC prior to using _rotary_. 
I recomend using the qc module of the [ATLAS pipeline](https://github.com/metagenome-atlas/atlas) for quick/robust QC.

### 3. Run
```bash
conda activate rotary

run_directory="E_coli" # Wherever you want to store the run files
config="myconfig.yaml"
conda_prefix="/Data/databases/rotary/conda_envs" # Wherever you want to store the conda envs, which can be re-used between runs
snakefile="rotary/rules/rotary.smk"
jobs=40

mkdir -p "${run_directory}"

# Optional command to just see what rules will be run
snakemake --dryrun --snakefile "${snakefile}" --configfile "${config}" \
  --directory "${run_directory}" --use-conda --conda-frontend mamba \
  --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee genome_longread_steps.log

# Run the pipeline
snakemake --snakefile "${snakefile}" --configfile "${config}" \
  --directory "${run_directory}" --use-conda --conda-frontend mamba \
  --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee genome_longread.log
```

### 4. Post-run tips
To double check that everything ran correctly, I recommend to briefly check all files in the `log` and `stats` folders 
once the run is finished.

In addition, please check the following:
- QC results (`logs/qc_long.log`) -- shows how many reads were retained vs. discarded during the QC filter step (this is technically in the `logs` folder mentioned above, but I wanted to add it here again to stress that it is worth checking)
- Assembly quality (`assembly/assembly_info.txt`) -- see how many contigs you got and circular vs linear status
- Detailed end repair log (`assembly/end_repair/verbose.log`) -- this is hard to read, but it's helpful to look for any 
  signs of errors or warnings. You can also see how contig rotation went after the stitch.
- Filtering of contigs by coverage (`polish/cov_filter/filtered_contigs.list`) -- did anything get removed that you wanted?
- Identification of the start marker gene (e.g., _dnaA_) - see `circularize/identify/hmmsearch_hits.txt` 
  and `start_genes.ffn` in the same folder. Did you get the hits you expected? Were they filtered properly into the `.ffn` file?
- Contig rotation (`circularize/circlator/rotated.log`) -- was the start gene reasonably far off from the contig end 
  (so that the second round of polishing will actually improve things?)
- Re-polishing (`stats/circularize/polypolish_changes.log`) - there should be few to zero changes if everything went smoothly. 
  If you see more than about 20 changes, it means your genome might have some odd difficult-to-correct regions.

## Final comments
Enjoy! I hope to continue to update/improve this pipeline (and remove some of the below caveats) over time, but for now, 
please feel free to use this basic working version.

## Appendix

### _rotary_ workflow summary
1. Perform simple long read QC using 
   [reformat.sh from bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/reformat-guide/)
2. Assemble long reads using [Flye](https://github.com/fenderglass/Flye)
3. Perform end repair of any circular contigs (using a heavily modified version of the 
   [circlator](https://github.com/sanger-pathogens/circlator) workflow)
4. Polish contigs using long reads via [medaka](https://github.com/nanoporetech/medaka)
5. If short read are provided, polishes using [Polypolish](https://github.com/rrwick/Polypolish) 
   and [POLCA](https://github.com/alekseyzimin/masurca)
6. Filters resulting contigs by a user-provided coverage threshold
7. Rotates any circular contigs to start at a marker gene of your choice (_dnaA_ by default), with help from 
   [circlator](https://github.com/sanger-pathogens/circlator)
8. Performs one more round of polishing on circular contigs, either using 
   Polypolish (if short reads were provided) or medaka (if only long reads were provided)
9. Gene prediction via [DFAST](https://github.com/nigyta/dfast_core)
10. Functional and taxonomic annotation via [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) 
    and [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

### Known issues
- Short read QC is not performed (you have to do it in advance)
- Limited flags are exposed for some key tools in the pipeline (e.g., Flye)
- Edge case: If you get "unlucky" and genome rotation is not substantial (e.g., _dnaA_ is already at the end of the contig, re-polishing will have effect in improving assembly quality). Ideally, I should add a test of how the contig was rotated after finding _dnaA_.
- Very rare edge case: Similarly, for short contigs, the error-prone region from end repair might end up near the end of short contigs (e.g., < 100kb long). This means that it could miss the benefits of long read polishing. It will still receive short read polishing after the circularization module, but sometimes short read polishing is less efficient if long read polishing has not been performed properly in advance. Ideally, I should add a check of how short contigs have been rotated after end repair.

### Future ideas
- Add proper CLI (including auto generation of config)
- Create a conda install (e.g., in bioconda)
- Add summary reports
- Add other assembly options
- Auto classify contigs as chromosomes vs. plasmids
- Move the end repair script from Bash to Python

## Footnotes
<sup>1</sup> Currently, Nanopore data generated with R9.4.1 flow cells (or earlier) requires short reads for error correction. A [recent paper](https://doi.org/10.1038/s41592-022-01539-7) demonstrated that microbial genome assemblies from R10.4 flow cells no longer benefit from short read error correction compared to just using long reads (very cool!).
