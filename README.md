# rotary

[![GitHub release](https://img.shields.io/badge/Version-0.2.0--beta4-lightgrey.svg)](https://github.com/jmtsuji/rotary/releases)
[![DOI](https://zenodo.org/badge/473891963.svg)](https://zenodo.org/badge/latestdoi/473891963)

Assembly/annotation workflow for Nanopore-based microbial genome data containing circular DNA elements

## Quick start

### Install
```bash
git clone https://github.com/rotary-genomics/rotary.git

conda env create -n rotary --file=rotary/enviroment.yaml

conda activate rotary

cd rotary

pip install --editable .
```

### Run One Sample
```bash
mkdir output_dir
mkdir rotary_db_dir

cd output_dir

rotary run_one -l s1_long.fastq.gz -r1 s1_R1.fastq.gz -r2 s1_R2.fastq.gz -d ../rotary_db_dir
```
**Note**: Multiple samples can also be run in batch using the `rotary init` and `rotary run` commands. 

## Description

_rotary_ is a snakemake pipeline that can be used to assemble single microbial genomes using
standalone Nanopore data<sup>[1](#Footnotes)</sup> or hybrid Nanopore + short read data.
The pipeline performs short read qc, short read decontamination, long read QC, assembly, 
end repair, polishing, contig rotation, and genome annotation.

### Some advantages of using _rotary_:

- All databases auto-install, so you can start analyzing genomes reproducibly with limited effort!
- Snakemake checkpointing allows you to restart a failed run from where you left off
- Circularization is handled fairly carefully. Unlike the defaults in most pipelines, _rotary_ fixes the
  [short gap region](https://github.com/fenderglass/Flye/issues/315#issuecomment-720679812) that can occur at the ends
  of circular contigs produced by Flye. It also polishes the circular contigs in two
  different rotation states to try to correct errors near contig ends.
- A robust annotation pipeline that covers gene annotation, GTDB taxonomy prediction, and completeness and 
  contamination estimation.

## Requirements

- OS: Runs on Linux (tested on Ubuntu 20.04)
- Software: requires `miniconda`
- Resources: The majority of the pipeline is not too resource intensive. The limiting factors are:
    - bbduk requires a lot of RAM when decontaminating short reads with large host genomes as reference 
      (e.g., the human genome requires 148 GB). The decontaminating with large genomes step can be skipped 
       in most use cases.
    - Flye requires moderate RAM (e.g., < 64 GB for typical bacterial genome runs)
    - GTDB-Tk v2 (with GTDB r214) needs ~55 GB RAM
    - EggNOG-mapper is a bit slow (e.g., ~40 minutes on 40 CPU threads)

## Usage

### 1. Install

```bash
git clone https://github.com/jmtsuji/rotary.git
conda env create -n rotary --file=rotary/enviroment.yaml
conda activate rotary
cd rotary
pip install --editable .
```

### 2a. Run One Sample

```bash
mkdir output_dir
mkdir rotary_db_dir

cd output_dir

rotary run_one -l s1_long.fastq.gz -r1 s1_R1.fastq.gz -r2 s1_R2.fastq.gz -d ../rotary_db_dir
```
**Note**: If you are using older nanopore flow cells you should stop the run, modify the config file 
(see **Advanced Usage** below) and restart the run using the `rotary run` command.

### 2b. Run Multiple Samples

Rotary can target a directory containing numerous FASTQ files derived from various samples.
It automatically organizes these files into sets corresponding to each sample and constructs a project 
directory that includes the necessary configuration files for executing the `rotary run` command. The 
`rotary run` command can then be used to run the workflow on an entire batch of samples.

```bash
mkdir output_dir
mkdir rotary_db_dir

rotary init -d rotary_db_dir -i sample_fastq_dir -o output_dir

cd output_dir

rotary run
```

## Advanced Usage

### Modifying the config (YAML) file (`config.yaml`).

Both `rotary run_one` and `rotary init` generate a `config.yaml` file in the output directory. Running `rotary run` in
the directory will utilize this file for run parameters. The config file can be modified to change rotary's behavior.

Very Important Parameters:

- `flye_input_mode`: set to "nano-hq" if your reads have <= 5% error rate; otherwise use "nano-raw"
- `medaka_model`: set this parameter to match the flow cell version / basecalling model you used to generate the long
  reads. See details in the [Medaka Github repo's "Models" section](https://github.com/nanoporetech/medaka#models)
- "Post-polishing contig filter" section: set the min depth and coverage you would like for a contig to be kept (more
  info in that section of the YAML file)

There are other advanced parameters that you can also edit if you'd like.

Lastly, make sure you set the threads and memory to values that make sense for your server.

### Post-run tips

To double check that everything ran correctly, I recommend to briefly check all files in the `log` and `stats` folders
once the run is finished.

In addition, please check the following in each sample folder:

- QC results (`logs/qc/qc_long.log`) -- shows how many reads were retained vs. discarded during the QC filter step (this
  is technically in the `logs` folder mentioned above, but I wanted to add it here again to stress that it is worth
  checking)
- Assembly quality (`assembly/flye/[SAMPLE_ID]_assembly_info.txt`) -- see how many contigs you got and circular vs linear status
- End repair results (`assembly/[SAMPLE_ID]_circular_info.tsv`) -- you can see if any circular contigs could not be 
  repaired successfully at their ends
- Filtering of contigs by coverage (`polish/cov_filter/[SAMPLE_ID]_filtered_contigs.list`) -- did anything get removed that you
  wanted?
- Identification of the start marker gene (e.g., _dnaA_) - see `circularize/identify/[SAMPLE_ID]_hmmsearch_hits.txt`
  and `[SAMPLE_ID]_start_genes.ffn` in the same folder. Did you get the hits you expected? Were they appropriately filtered into 
  the `.ffn` file?
- Contig rotation (`circularize/circlator/rotated.log`) -- was the start gene reasonably far off from the contig end
  (so that the second round of polishing will actually improve things?)
- Re-polishing (`stats/circularize/polypolish_changes.log`) - there should be few to zero changes if everything went
  smoothly. If you see more than about 20 changes, it means your genome might have some odd difficult-to-correct regions.

## Citation

_rotary_ is currently described in the methods of a bioRxiv pre-print. Please cite this pre-print if you use _rotary_:
> Tsuji JM, Shaw NA, Nagashima S, Venkiteswaran JJ, Schiff SL, Watanabe T, Fukui M, Hanada S, Tank M, Neufeld JD (2023).
> Anoxygenic phototrophic _Chloroflexota_ member uses a Type I reaction center. _bioRxiv_, DOI:10.1101/2020.07.07.190934

## Final comments

Enjoy! I hope to continue to update/improve this pipeline (and remove some of the below caveats) over time, but for now,
please feel free to use this basic working version.

## Appendix

### _rotary_ workflow summary

1. Perform short read QC (e.g., score, length, and adapter trimming) using [bbduk.sh from bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
2. Remove short read contamination using [bbduk.sh from bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
3. Perform simple long read QC using
   [reformat.sh from bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/reformat-guide/)
4. Assemble long reads using [Flye](https://github.com/fenderglass/Flye)
5. Perform end repair of any circular contigs (using a heavily modified version of the
   [circlator](https://github.com/sanger-pathogens/circlator) workflow)
6. Polish contigs using long reads via [medaka](https://github.com/nanoporetech/medaka)
7. If short reads are provided, polishes using [Polypolish](https://github.com/rrwick/Polypolish)
   and [POLCA](https://github.com/alekseyzimin/masurca)
8. Filters resulting contigs by a user-provided coverage threshold
9. Rotates any circular contigs to start at a marker gene of your choice (_dnaA_ by default), with help from
   [circlator](https://github.com/sanger-pathogens/circlator)
10. Performs one more round of polishing on circular contigs, either using
    Polypolish (if short reads were provided) or medaka (if only long reads were provided)
11. Gene prediction via [DFAST](https://github.com/nigyta/dfast_core)
12. Functional and taxonomic annotation via [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)
    and [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

### Known issues

- Limited flags are exposed for some key tools in the pipeline (e.g., Flye)
- Minimum supported length for circular contigs is 100 kb (it should be a fairly easy fix in a future version to
  decrease this threshold)
- Edge case: If you get "unlucky" and genome rotation is not substantial (e.g., _dnaA_ is already at the end of the
  contig), re-polishing will have limited effect on improving assembly quality. Ideally, I should add a test of how the 
  contig was rotated after finding _dnaA_.
- Very rare edge case: Similarly, for short contigs, the error-prone region from end repair might end up near the end of
  short contigs (e.g., < 100kb long). This means that it could miss the benefits of long read polishing. It will still
  receive short read polishing after the circularization module, but sometimes short read polishing is less efficient if
  long read polishing has not been performed properly in advance. Ideally, I should add a check of how short contigs
  have been rotated after end repair.
- If rotary is run in batch mode, then all samples in one batch must have only long reads, or all samples in one batch 
  must have long and short reads. In the future, we hope to make the pipeline more flexible in this area.

### Future to-do's and ideas

- Create a conda install (e.g., in bioconda)
- Add summary reports
- Add other assembly options
- Auto classify contigs as chromosomes vs. plasmids

## Footnotes

<sup>1</sup> Currently, Nanopore data generated with R9.4.1 flow cells (or earlier) requires short reads for error
correction. A [recent paper](https://doi.org/10.1038/s41592-022-01539-7) demonstrated that microbial genome assemblies
from R10.4 flow cells no longer benefit from short read error correction compared to just using long reads (very cool!).
Thus, long-read only assembly and polishing might be sufficient if you have R10.4 (or higher) data - probably should be
called with a super-high accuracy Guppy model.
