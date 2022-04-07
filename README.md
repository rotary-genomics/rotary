# _tōyako_
Assembly/annotation workflow for Nanopore genome data

## Quick start
```bash
git clone https://github.com/jmtsuji/toyako.git

conda env create -n toyako --file=toyako/envs/toyako.yaml

cp toyako/config.yaml myconfig.yaml # Edit this file, especially the first few lines

conda activate toyako

mkdir -p output_dir conda_envs

snakemake --snakefile toyako/rules/toyako.yaml \
  --configfile myconfig.yaml \
  --directory output_dir \
  --conda-prefix conda_envs \
  --jobs 10 \
  --use-conda \
  --conda-frontend mamba \
  --rerun-incomplete \
  --reason \
  --printshellcmds 2>&1 | \
  tee toyako.log
```

## Description
_tōyako_ is a snakemake pipeline that can be used to assemble single microbial genomes using 
standlone Nanopore data<sup>[1](#Footnotes)</sup> or hybrid Nanopore + short read data. 
The pipeline performs long read QC, assembly, end repair, polishing, contig rotation, and genome annotation.

### Some advantages of using _tōyako_:
- All databases auto-install, so you can start analyzing genomes reproducibly with limited effort!
- Snakemake checkpointing allows you to re-start a failed run from where you left off
- Circularization is handled fairly carefully. Unlike the defaults in most pipelines, _tōyako_ fixes the 
  [short gap region](https://github.com/fenderglass/Flye/issues/315#issuecomment-720679812) that can occur at the ends 
  of circular contigs produced by Flye. If short reads are provided, it then polishes the circular contigs in two 
  different rotation states to try to correct errors near contig ends.

## Requirements
- OS: Runs on Linux (tested on Ubuntu 20.04)
- Software: requires `miniconda`
- Resources: The majority of the pipeline is not too resource intensive. The limiting factors are:
  - Flye requires moderate RAM (e.g., < 64 GB for typical bacterial genome runs)
  - GTDB-Tk needs a huge amount of RAM (~220 GB!)
  - EggNOG-mapper is a bit slow (e.g., ~40 minutes on 40 CPU threads)

## Usage
1. Install
```bash
git clone https://github.com/jmtsuji/toyako.git
conda env create -n toyako --file=toyako/envs/toyako.yaml
```

2. Copy and fill out the config (YAML) file (`config.yaml`).

The key parameters to fill are the genome name, the paths to the long and (QC'ed) short reads, and the DB dir.

Save as something like `genome_longread_mycopy.yaml`.

**Note** that short read QC is not performed, so you'll want to do short read QC prior to using _tōyako_. 
I recomend using the qc module of the [ATLAS pipeline](https://github.com/metagenome-atlas/atlas) for quick/robust QC.

3. Run
```bash
conda activate genome_longread

run_directory="E_coli" # Wherever you want to store the run files
config="genome_longread_mycopy.yaml"
conda_prefix="/Data/databases/toyako/conda_envs" # Wherever you want to store the conda envs, which can be re-used between runs
snakefile="toyako/rules/toyako.smk"
jobs=40

mkdir -p "${run_directory}"

# Optional command to just see what will be run
snakemake --dryrun --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee genome_longread_steps.log

# Run the pipeline
snakemake --snakefile "${snakefile}" --configfile "${config}" --directory "${run_directory}" \
  --use-conda --conda-frontend mamba --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete \
  --reason --printshellcmds 2>&1 | tee genome_longread.log
```

4. Post-run tips
This pipeline is in its early development stages, so I recommend to check all files in the `log` and `stats` folders 
once the run is complete to confirm everything worked normally.

In addition, please check the following:
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
Enjoy! I hope to continue to update this pipeline (and remove some of the above caveats) over time, but for now, 
please feel free to use this basic working version.

## Appendix

### _tōyako_ workflow summary
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
8. Performs one more round of short read polishing via [Polypolish](https://github.com/rrwick/Polypolish), 
   if short reads are provided, on circular contigs
9. Gene prediction via [DFAST](https://github.com/nigyta/dfast_core)
10. Functional and taxonomic annotation via [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) 
    and [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

### Known issues
- Re-polishing is not performed if short reads are not provided
- Short read QC is not performed (you have to do it in advance)
- If you get "unlucky" and genome rotation is not substantial (e.g., _dnaA_ is already at the end of the contig, re-polishing will have effect in improving assembly quality)
- Similarly, for short contigs, the error-prone region from end repair might end up near the end of short contigs (e.g., < 100kb long)
- Limited flags are exposed for some key tools in the pipeline (e.g., Flye)

### Future ideas
- Add proper CLI (including auto generation of config)
- Make conda install
- Add summary reports
- Add other assembly options
- Auto classify contigs as chromosomes vs. plasmids

## Footnotes
<sup>1</sup> Currently, most Nanopore data requires short reads for error correction, although there are some reports that r10.4 assemblies might not need short read error correction.
