# rotary: utilities and workflow for long-read DNA assemblies including circular elements
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os
import sys
import pandas as pd
import itertools
from snakemake.utils import min_version

from rotary.sample import parse_sample_tsv
from rotary.utils import symlink_or_compress


ZENODO_VERSION = "10087395"
VERSION_POLYPOLISH="0.5.0"
VERSION_DFAST="1.2.18"
VERSION_EGGNOG="5.0.0" # See http://eggnog5.embl.de/#/app/downloads
START_HMM_NAME = os.path.splitext(os.path.basename(config.get("hmm_url")))[0]
VERSION_GTDB_COMPLETE= "214.1" # See https://data.gtdb.ecogenomic.org/releases/
VERSION_GTDB_MAIN=VERSION_GTDB_COMPLETE.split('.')[0] # Remove subversion
DB_DIR_PATH = config.get('db_dir')

SAMPLE_TSV_PATH = 'samples.tsv'
SAMPLES = parse_sample_tsv(SAMPLE_TSV_PATH)

SAMPLE_NAMES = list(SAMPLES.keys())

CONTAMINATION_NCBI_ACCESSIONS = config.get("contamination_references_ncbi_accessions")
CUSTOM_CONTAMINATION_FILEPATHS = config.get("contamination_references_custom_filepaths")

if CUSTOM_CONTAMINATION_FILEPATHS:
    CUSTOM_CONTAMINATION_FILE_NAMES = [os.path.basename(path).split('.')[0] for path in CUSTOM_CONTAMINATION_FILEPATHS]
else:
    CUSTOM_CONTAMINATION_FILE_NAMES = []

if not CONTAMINATION_NCBI_ACCESSIONS and not CUSTOM_CONTAMINATION_FILEPATHS:
    CONTAMINANT_REFERENCE_GENOMES = False
else:
    CONTAMINANT_REFERENCE_GENOMES = True

# Specify the minimum snakemake version allowable
min_version("7.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")

if str(config.get('polish_with_short_reads')).lower() == 'true':
    POLISH_WITH_SHORT_READS = True
else:
    POLISH_WITH_SHORT_READS = False

rule all:
    input:
        "checkpoints/qc",
        "checkpoints/assembly",
        "checkpoints/polish",
        "checkpoints/circularize",
        "checkpoints/annotation"


rule download_short_read_adapters:
    """
    Downloads sequence adapters to trim from short reads (for now, uses the ATLAS version)
    """
    output:
        os.path.join(DB_DIR_PATH, "adapters.fasta")
    log:
        "logs/download/download_short_read_adapters.log"
    benchmark:
        "benchmarks/download/download_short_read_adapters.benchmark.txt"
    params:
        url = f"https://zenodo.org/records/{ZENODO_VERSION}/files/adapters.fasta"
    shell:
        """
        wget -O {output} {params.url} > {log} 2>&1
        """


rule download_ncbi_contamination_reference:
    """
    Downloads references genome for contamination screening using its NCBI genome accession
    """
    output:
        genome = os.path.join(DB_DIR_PATH, 'contamination_references', 'ncbi', '{accession}.fna.gz'),
        zip = temp(os.path.join(DB_DIR_PATH, 'contamination_references', 'ncbi', '{accession}.zip')),
        zip_dir = temp(directory(os.path.join(DB_DIR_PATH, 'contamination_references', 'ncbi', '{accession}')))
    conda:
        "../envs/download.yaml"
    log:
        "logs/download/download_contamination_reference_{accession}.log"
    benchmark:
        "benchmarks/download/download_contamination_reference_{accession}.benchmark.txt"
    params:
        accession = "{accession}"
    threads:
        min(config.get("threads", 1), 4)
    shell:
        """        
        echo "### Downloading genome: {params.accession} ###" > {log}
        datasets download genome accession {params.accession} --include genome --filename {output.zip} 2>> {log}
        unzip -d {output.zip_dir} {output.zip} > /dev/null
        
        # Confirm that there is only one genome file matching the expected pattern in the unzipped folder
        genome_file=($(find {output.zip_dir}/ncbi_dataset/data/{params.accession} -type f -name "{params.accession}_*_genomic.fna"))
        if [[ "${{#genome_file[@]}}" == 1 ]]; then
          pigz -c -p {threads} "${{genome_file[0]}}" > {output.genome}
        else
          echo "ERROR: more than 1 genome file (or no genome file) in dir {output.zip_dir}/ncbi_dataset/data/{params.accession}"
          exit 1
        fi
        """

rule set_up_custom_contamination_references:
    output:
        expand(os.path.join(DB_DIR_PATH, 'contamination_references', 'custom', '{contaminant_name}.fna.gz'),
            contaminant_name=CUSTOM_CONTAMINATION_FILE_NAMES),
    run:
        for name, path in zip(CUSTOM_CONTAMINATION_FILE_NAMES,CUSTOM_CONTAMINATION_FILEPATHS):
            symlink_or_compress(path, os.path.join(DB_DIR_PATH, 'contamination_references', 'custom', f'{name}.fna.gz'))

rule install_polypolish:
    output:
        polypolish_filter=os.path.join(DB_DIR_PATH,"polypolish_" + VERSION_POLYPOLISH,"polypolish_insert_filter.py"),
        polypolish=os.path.join(DB_DIR_PATH,"polypolish_" + VERSION_POLYPOLISH,"polypolish"),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","polypolish_" + VERSION_POLYPOLISH)
    log:
        "logs/download/install_polypolish.log"
    benchmark:
        "benchmarks/download/install_polypolish.txt"
    params:
        db_dir=os.path.join(DB_DIR_PATH,"polypolish_" + VERSION_POLYPOLISH),
        url="https://github.com/rrwick/Polypolish/releases/download/v" + VERSION_POLYPOLISH + "/polypolish-linux-x86_64-musl-v" + VERSION_POLYPOLISH + ".tar.gz"
    shell:
        """
        mkdir -p {params.db_dir}
        wget -O - {params.url} 2> {log} | tar -C {params.db_dir} -xzf - >> {log} 2>&1
        touch {output.install_finished}
        """


# TODO - does not check the HMM version, only ID. If the HMM version updates, it won't automatically re-download
rule download_hmm:
    output:
        hmm=os.path.join(DB_DIR_PATH,"hmm",START_HMM_NAME + ".hmm")
    log:
        "logs/download/hmm_download.log"
    benchmark:
        "benchmarks/download/hmm_download.txt"
    params:
        db_dir=os.path.join(DB_DIR_PATH,"hmm"),
        url=config.get("hmm_url")
    shell:
        """
        mkdir -p {params.db_dir}
        wget -O {output.hmm} {params.url} 2> {log}
        """


rule download_dfast_db:
    output:
        db=directory(os.path.join(DB_DIR_PATH,"dfast_" + VERSION_DFAST)),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","dfast_" + VERSION_DFAST)
    conda:
        "../envs/annotation_dfast.yaml"
    log:
        "logs/download/dfast_db_download.log"
    benchmark:
        "benchmarks/download/dfast_db_download.txt"
    params:
        db_dir=os.path.join(DB_DIR_PATH,"dfast_" + VERSION_DFAST)
    shell:
        """
        mkdir -p {params.db_dir}
        dfast_file_downloader.py --protein dfast --dbroot {params.db_dir} > {log} 2>&1
        dfast_file_downloader.py --cdd Cog --hmm TIGR --dbroot {params.db_dir} >> {log} 2>&1
        touch {output.install_finished}
        """


rule download_eggnog_db:
    output:
        db=directory(os.path.join(DB_DIR_PATH,"eggnog_" + VERSION_EGGNOG)),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","eggnog_" + VERSION_EGGNOG)
    conda:
        "../envs/eggnog.yaml"
    log:
        "logs/download/eggnog_db_download.log"
    benchmark:
        "benchmarks/download/eggnog_db_download.txt"
    params:
        db_dir=os.path.join(DB_DIR_PATH,"eggnog_" + VERSION_EGGNOG)
    shell:
        """
        mkdir -p {params.db_dir}
        download_eggnog_data.py -y -M --data_dir {params.db_dir} > {log} 2>&1
        touch {output.install_finished}
        """


# TODO - if there is an error during download, the initial_download_dir is not deleted during cleanup
rule download_gtdb_db:
    output:
        db=directory(os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE)),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_download")
    log:
        "logs/download/gtdb_db_download.log"
    benchmark:
        "benchmarks/download/gtdb_db_download.txt"
    params:
        db_dir_root=os.path.join(DB_DIR_PATH),
        initial_download_dir=os.path.join(DB_DIR_PATH,"release" + VERSION_GTDB_MAIN),
        db_dir=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE),
        url="https://data.gtdb.ecogenomic.org/releases/release" + VERSION_GTDB_MAIN + "/" + VERSION_GTDB_COMPLETE + "/auxillary_files/gtdbtk_r" + VERSION_GTDB_MAIN + "_data.tar.gz"
    shell:
        """
        mkdir -p {params.db_dir_root}

        if [[ -d {params.initial_download_dir} ]]; then
          echo "GTDB cannot be downloaded because of pre-existing dir: {params.initial_download_dir}"
          exit 1
        fi
        
        wget -O - {params.url} 2> {log} | tar -C {params.db_dir_root} -xzf - >> {log} 2>&1
        mv {params.initial_download_dir} {params.db_dir}

        touch {output.install_finished}
        """


rule validate_gtdb:
    input:
        os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_download")
    output:
        os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_validate")
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/download/gtdb_db_validate.log"
    benchmark:
        "benchmarks/download/gtdb_db_validate.txt"
    params:
        db_dir=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE)
    shell:
        """
        # Split the "progress bar" style output into multiple lines with sed
        # See: https://stackoverflow.com/a/60786606 (accessed 2022.08.02)
        GTDBTK_DATA_PATH={params.db_dir} gtdbtk check_install | sed 's/\\r/\\n/g' > {log} 2>&1

        touch {output}
        """


rule build_gtdb_mash_ref_database:
    input:
        os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_validate")
    output:
        ref_msh_file=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE + '_mash','gtdb_ref_sketch.msh'),
        ref_genome_path_list=temp(os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE + '_mash','ref_mash_genomes.txt'))
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/download/gtdb_db_mash.log"
    benchmark:
        "benchmarks/download/gtdb_mash.txt"
    threads:
        config.get("threads",1)
    params:
        fast_ani_genomes_dir=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE,'fastani','database')
    shell:
        """
        find {params.fast_ani_genomes_dir} -name *_genomic.fna.gz -type f > {output.ref_genome_path_list}
        mash sketch -l {output.ref_genome_path_list} -p {threads} -o {output.ref_msh_file} -k 16 -s 5000 > {log} 2>&1   
        """

rule download_checkm_db:
    output:
        db=directory(os.path.join(DB_DIR_PATH,"checkm2")),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","checkm2")
    conda:
        "../envs/checkm2.yaml"
    log:
        "logs/download/checkm2_db_download.log"
    benchmark:
        "benchmarks/download/checkm2_db_download.txt"
    shell:
        """
        mkdir -p {output.db}
        checkm2 database --download --path {output.db} > {log} 2>&1
        touch {output.install_finished}
        """

rule set_up_sample_directories:
    input:
        SAMPLE_TSV_PATH
    output:
        long_reads = expand("{sample}/raw/{sample}_long.fastq.gz", sample=SAMPLE_NAMES),
        short_R1_reads = expand("{sample}/raw/{sample}_R1.fastq.gz", sample=SAMPLE_NAMES),
        short_R2_reads = expand("{sample}/raw/{sample}_R2.fastq.gz", sample=SAMPLE_NAMES),
    run:
        for sample in SAMPLES.values():
            identifier = sample.identifier
            symlink_or_compress(sample.long_read_path,f'{identifier}/raw/{identifier}_long.fastq.gz')
            symlink_or_compress(sample.short_read_left_path,f'{identifier}/raw/{identifier}_R1.fastq.gz')
            symlink_or_compress(sample.short_read_right_path,f'{identifier}/raw/{identifier}_R2.fastq.gz')


rule nanopore_qc_filter:
    input:
        "{sample}/raw/{sample}_long.fastq.gz"
    output:
        "{sample}/qc/long/{sample}_nanopore_qc.fastq.gz"
    conda:
        "../envs/mapping.yaml"
    log:
        "{sample}/logs/qc/long/qc_long.log"
    benchmark:
        "{sample}/benchmarks/qc/long/qc_long.benchmark.txt"
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    params:
        minlength = config.get("minlength"),
        minavgquality = config.get("minavgquality")
    shell:
        """
        reformat.sh in={input} out={output} minlength={params.minlength} minavgquality={params.minavgquality} \
          interleaved=f qin=33 threads={threads} -Xmx{resources.mem}G > {log} 2>&1
        """


rule qc_long_length_hist:
    input:
        "{sample}/qc/long/{sample}_nanopore_qc.fastq.gz"
    output:
        "{sample}/qc/long/{sample}_length_hist.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "{sample}/logs/qc/long/qc_long_length_hist.log"
    benchmark:
        "{sample}/benchmarks/qc/long/qc_long_length_hist.benchmark.txt"
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        reformat.sh in={input} lhist={output} maxhistlen=10000000 \
          interleaved=f qin=33 threads={threads} -Xmx{resources.mem}G > {log} 2>&1
        """


rule qc_long_length_stats:
    input:
        "{sample}/qc/long/{sample}_length_hist.tsv"
    output:
        "{sample}/stats/{sample}_qc_long_length_stats.txt"
    log:
        "{sample}/logs/qc/long/qc_long_length_stats.log"
    benchmark:
        "{sample}/benchmarks/qc/long/qc_long_length_stats.benchmark.txt"
    run:
        length_hist = pd.read_csv(input[0], sep='\t')

        lengths = []

        for index, row in length_hist.iterrows():
            length, count = row
            lengths = lengths + list(itertools.repeat(length,count))

        lengths = pd.Series(lengths)

        length_stats = pd.DataFrame({'Total reads':   [lengths.shape[0]],
                                     'Mean length':   [round(lengths.mean(),2)],
                                     'Median length': [lengths.median()],
                                     'Min length':    [lengths.min()],
                                     'Max length':    [lengths.max()]})\
          .transpose()

        length_stats.to_csv(output[0], sep='\t', header=None, index=True)


rule finalize_qc_long:
    input:
        "{sample}/qc/long/{sample}_nanopore_qc.fastq.gz"
    output:
        "{sample}/qc/{sample}_qc_long.fastq.gz"
    run:
        source_relpath = os.path.relpath(str(input), os.path.dirname(str(output)))
        os.symlink(source_relpath, str(output))


rule qc_long:
    input:
        expand("{sample}/qc/{sample}_qc_long.fastq.gz", sample=SAMPLE_NAMES),
        expand("{sample}/stats/{sample}_qc_long_length_stats.txt", sample=SAMPLE_NAMES)
    output:
        temp(touch("checkpoints/qc_long"))


rule short_read_reformat:
    """
    Makes the input format of the short reads consistent (e.g., by trimming read descriptions and tossing invalid 
    nucleotide characters). Threads are locked at a maximum of 4 because this code is IO limited.
    """
    input:
        short_r1 = "{sample}/raw/{sample}_R1.fastq.gz",
        short_r2 = "{sample}/raw/{sample}_R2.fastq.gz"
    output:
        short_reformat_r1 = temp("{sample}/qc/short/{sample}_reformat_R1.fastq.gz"),
        short_reformat_r2 = temp("{sample}/qc/short/{sample}_reformat_R2.fastq.gz"),
        quality_histogram = "{sample}/qc/short/{sample}_reformat_qhist.tsv"
    conda:
        "../envs/qc.yaml"
    log:
        "{sample}/logs/qc/short/short_read_reformat.log"
    benchmark:
        "{sample}/benchmarks/qc/short/short_read_reformat.benchmark.txt"
    threads:
        min(config.get("threads", 1), 4)
    resources:
        mem = config.get("memory")
    shell:
        """
        reformat.sh -Xmx{resources.mem}g threads={threads} in={input.short_r1} in2={input.short_r2} \
          out={output.short_reformat_r1} out2={output.short_reformat_r2} qhist={output.quality_histogram} \
          overwrite=t interleaved=f qin=33 verifypaired=t trimreaddescription=t tossjunk=t pigz=t unpigz=t \
          2> {log}
        """


rule short_read_adapter_trimming:
    """
    Trims 3' adapters off the reads (e.g., that are caused by having a short insert).
    """
    input:
        short_r1 = "{sample}/qc/short/{sample}_reformat_R1.fastq.gz",
        short_r2 = "{sample}/qc/short/{sample}_reformat_R2.fastq.gz",
        adapters = os.path.join(DB_DIR_PATH, "adapters.fasta")
    output:
        short_r1 = temp("{sample}/qc/short/{sample}_adapter_trim_R1.fastq.gz"),
        short_r2 = temp("{sample}/qc/short/{sample}_adapter_trim_R2.fastq.gz"),
        adapter_trim_stats= "{sample}/stats/{sample}_short_read_adapter_trimming.txt"
    conda:
        "../envs/qc.yaml"
    log:
        "{sample}/logs/qc/short/short_read_adapter_trimming.log"
    benchmark:
        "{sample}/benchmarks/qc/short/short_read_adapter_trimming.benchmark.txt"
    params:
        adapter_trimming_kmer_length = config.get("adapter_trimming_kmer_length"),
        minimum_detectable_adapter_length_on_read_end = config.get("minimum_detectable_adapter_length_on_read_end"),
        trim_adapters_by_overlap = "t" if str(config.get("overlap_based_trimming")).lower() == "true" else "f",
        min_read_length = config.get("minimum_read_length_adapter_trim")
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("memory")
    shell:
        """
        bbduk.sh -Xmx{resources.mem}g threads={threads} in={input.short_r1} in2={input.short_r2} \
          out={output.short_r1} out2={output.short_r2} ref={input.adapters} \
          k={params.adapter_trimming_kmer_length} ktrim=r mink={params.minimum_detectable_adapter_length_on_read_end} \
          rcomp=t trimbyoverlap={params.trim_adapters_by_overlap} minoverlap=14 mininsert=40 \
          minlength={params.min_read_length} stats={output.adapter_trim_stats} overwrite=t interleaved=f qin=33 \
          pigz=t unpigz=t \
          2> {log}
        """


rule short_read_quality_trimming:
    """
    Performs quality trimming after trimming adapters.
    The conditional statements in the input field select whether to skip adapter trimming based on the config file.
    Note: if quality_trim_direction is set as "f", the reads are just passed through this rule and not trimmed.
    """
    input:
        short_r1 = "{sample}/qc/short/{sample}_adapter_trim_R1.fastq.gz" if str(config.get("perform_adapter_trimming")).lower() == 'true' else "{sample}/qc/short/{sample}_reformat_R1.fastq.gz",
        short_r2 = "{sample}/qc/short/{sample}_adapter_trim_R2.fastq.gz" if str(config.get("perform_adapter_trimming")).lower() == 'true' else "{sample}/qc/short/{sample}_reformat_R2.fastq.gz"
    output:
        short_r1 = temp("{sample}/qc/short/{sample}_quality_trim_R1.fastq.gz"),
        short_r2 = temp("{sample}/qc/short/{sample}_quality_trim_R2.fastq.gz"),
        quality_histogram = "{sample}/qc/short/{sample}_quality_trim_qhist.tsv"
    conda:
        "../envs/qc.yaml"
    log:
        "{sample}/logs/qc/short/short_read_quality_trimming.log"
    benchmark:
        "{sample}/benchmarks/qc/short/short_read_quality_trimming.benchmark.txt"
    params:
        quality_trim_direction = config.get("quality_trim_direction"),
        quality_trim_cutoff = config.get("quality_trim_score_cutoff"),
        min_read_length = config.get("minimum_read_length_quality_trim"),
        min_average_quality = config.get("minimum_average_quality_score_post_trim") if config.get("quality_trim_direction") != "f" else 0
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("memory")
    shell:
        """
        bbduk.sh -Xmx{resources.mem}g threads={threads} in={input.short_r1} in2={input.short_r2} \
          out={output.short_r1} out2={output.short_r2} qhist={output.quality_histogram} \
          qtrim={params.quality_trim_direction} trimq={params.quality_trim_cutoff} minlength={params.min_read_length} \
          minavgquality={params.min_average_quality} overwrite=t interleaved=f qin=33 pigz=t unpigz=t \
          2> {log}
        """


rule short_read_contamination_filter:
    """
    Filters short reads based on match to a reference. 
    Note that contamination_references in the input is only there to trigger download of the genome files; the actual 
    text to be input into the shell (with proper comma separation) is prepared in all_contaminant_references in params.
    """
    input:
        short_r1 = "{sample}/qc/short/{sample}_quality_trim_R1.fastq.gz",
        short_r2 = "{sample}/qc/short/{sample}_quality_trim_R2.fastq.gz",
        ncbi_contamination_references = expand(os.path.join(DB_DIR_PATH, 'contamination_references', 'ncbi',
            '{accession}.fna.gz'), accession=CONTAMINATION_NCBI_ACCESSIONS),
        custom_contamination_references = expand(os.path.join(DB_DIR_PATH, 'contamination_references', 'custom',
            '{contaminant_name}.fna.gz'), contaminant_name=CUSTOM_CONTAMINATION_FILE_NAMES)
    output:
        short_r1 = temp("{sample}/qc/short/{sample}_filter_R1.fastq.gz"),
        short_r2 = temp("{sample}/qc/short/{sample}_filter_R2.fastq.gz"),
        quality_histogram = "{sample}/qc/short/{sample}_filter_qhist.tsv",
        filter_stats= "{sample}/stats/{sample}_short_read_contamination_filter.txt"
    conda:
        "../envs/qc.yaml"
    log:
        "{sample}/logs/qc/short/short_read_contamination_filter.log"
    benchmark:
        "{sample}/benchmarks/qc/short/short_read_contamination_filter.benchmark.txt"
    params:
        contamination_filter_kmer_length = config.get("contamination_filter_kmer_length"),
        contamination_references = lambda wildcards, input: ','.join(input.ncbi_contamination_references +
                                                                     input.custom_contamination_references)
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("memory")
    shell:
        """
        bbduk.sh -Xmx{resources.mem}g threads={threads} in={input.short_r1} in2={input.short_r2} \
          ref={params.contamination_references} out={output.short_r1} out2={output.short_r2} \
          stats={output.filter_stats} qhist={output.quality_histogram} \
          k={params.contamination_filter_kmer_length} ktrim=f rcomp=t \
          overwrite=t interleaved=f qin=33 pigz=t unpigz=t \
          2> {log}
        """


rule finalize_qc_short:
    """
    The conditional statements in this rule control whether or not contaminant filtration is performed.
    """
    input:
        short_r1 = "{sample}/qc/short/{sample}_filter_R1.fastq.gz" if CONTAMINANT_REFERENCE_GENOMES == True else "{sample}/qc/short/{sample}_quality_trim_R1.fastq.gz",
        short_r2 = "{sample}/qc/short/{sample}_filter_R2.fastq.gz" if CONTAMINANT_REFERENCE_GENOMES == True else "{sample}/qc/short/{sample}_quality_trim_R2.fastq.gz"
    output:
        short_r1 = "{sample}/qc/short/{sample}_qc_R1.fastq.gz",
        short_r2 = "{sample}/qc/short/{sample}_qc_R2.fastq.gz"
    shell:
        """
        cp {input.short_r1} {output.short_r1}
        cp {input.short_r2} {output.short_r2}
        """


rule symlink_qc_short:
    input:
        short_r1 = "{sample}/qc/short/{sample}_qc_R1.fastq.gz",
        short_r2 = "{sample}/qc/short/{sample}_qc_R2.fastq.gz"
    output:
        short_final_r1 = "{sample}/qc/{sample}_qc_R1.fastq.gz",
        short_final_r2 = "{sample}/qc/{sample}_qc_R2.fastq.gz"
    run:
        source_relpath = os.path.relpath(str(input.short_r1), os.path.dirname(str(output.short_final_r1)))
        os.symlink(source_relpath, str(output.short_final_r1))

        source_relpath = os.path.relpath(str(input.short_r2), os.path.dirname(str(output.short_final_r2)))
        os.symlink(source_relpath, str(output.short_final_r2))


rule qc_short:
    input:
        expand("{sample}/qc/{sample}_qc_R1.fastq.gz", sample=SAMPLE_NAMES),
        expand("{sample}/qc/{sample}_qc_R2.fastq.gz", sample=SAMPLE_NAMES)
    output:
        temp(touch("checkpoints/qc_short"))


rule qc:
    input:
        ["checkpoints/qc_long"] if POLISH_WITH_SHORT_READS == False else ["checkpoints/qc_long", "checkpoints/qc_short"]
    output:
        temp(touch("checkpoints/qc"))


rule assembly_flye:
    input:
        "{sample}/qc/{sample}_qc_long.fastq.gz"
    output:
        assembly="{sample}/assembly/flye/{sample}_assembly.fasta",
        info="{sample}/assembly/flye/{sample}_assembly_info.txt",
        output_dir=directory("{sample}/assembly/flye")
    conda:
        "../envs/assembly_flye.yaml"
    log:
        "{sample}/logs/assembly/assembly_flye.log"
    benchmark:
        "{sample}/benchmarks/assembly/assembly_flye.benchmark.txt"
    params:
        input_mode=config.get("flye_input_mode"),
        read_error="" if config.get("flye_read_error") == "auto" else "--read-error " + config.get("flye_read_error"),
        meta_mode="--meta" if config.get("flye_meta_mode") == "True" else "",
        polishing_rounds=config.get("flye_polishing_rounds")
    threads:
        config.get("threads",1)
    shell:
        """
        flye --{params.input_mode} {input} {params.read_error} --out-dir {output.output_dir} {params.meta_mode} \
          --iterations {params.polishing_rounds} -t {threads} > {log} 2>&1 
        mv {output.output_dir}/assembly.fasta {output.assembly} 
        mv {output.output_dir}/assembly_info.txt {output.info} 
        """


# TODO - eventually make this rule more generic so that outputs from other assemblers can go to the same output files
# TODO - the math for memory per thread should really be done somewhere other than 'resources' - what is best practice?
rule assembly_end_repair:
    input:
        qc_long_reads="{sample}/qc/{sample}_qc_long.fastq.gz",
        assembly="{sample}/assembly/flye/{sample}_assembly.fasta",
        info="{sample}/assembly/flye/{sample}_assembly_info.txt"
    output:
        assembly="{sample}/assembly/end_repair/{sample}_repaired.fasta",
        info="{sample}/assembly/end_repair/{sample}_repaired_info.tsv",
        output_dir=directory("{sample}/assembly/end_repair"),
    log:
        "{sample}/logs/assembly/end_repair.log"
    benchmark:
        "{sample}/benchmarks/assembly/end_repair.txt"
    params:
        flye_input_mode=config.get("flye_input_mode"),
        flye_read_error="0" if config.get("flye_read_error") == "auto" else config.get("flye_read_error"),
        min_id=config.get("circlator_merge_min_id"),
        min_length=config.get("circlator_merge_min_length"),
        ref_end=config.get("circlator_merge_ref_end"),
        reassemble_end=config.get("circlator_merge_reassemble_end"),
        keep_going="--keep_going_with_failed_contigs" if config.get("keep_unrepaired_contigs") == "True" else "",
    threads:
        config.get("threads",1)
    resources:
        mem=int(config.get("memory") / config.get("threads",1))
    shell:
        """
        rotary-repair --long_read_filepath {input.qc_long_reads} \
          --assembly_fasta_filepath {input.assembly} \
          --assembly_info_filepath {input.info} \
          --output_dir {output.output_dir} \
          --flye_read_mode {params.flye_input_mode} \
          --flye_read_error {params.flye_read_error} \
          --circlator_min_id {params.min_id} \
          --circlator_min_length {params.min_length} \
          --circlator_ref_end {params.ref_end} \
          --circlator_reassemble_end {params.reassemble_end} \
          --threads {threads} \
          --threads_mem {resources.mem} \
          --overwrite \
          {params.keep_going} \
          > {log} 2>&1
        mv {output.output_dir}/repaired.fasta {output.assembly} 
        mv {output.output_dir}/repaired_info.tsv {output.info}
        """


rule finalize_assembly:
    input:
        assembly="{sample}/assembly/end_repair/{sample}_repaired.fasta",
        info="{sample}/assembly/end_repair/{sample}_repaired_info.tsv"
    output:
        assembly="{sample}/assembly/{sample}_assembly.fasta",
        info="{sample}/assembly/{sample}_circular_info.tsv"
    run:
        source_relpath = os.path.relpath(str(input.assembly),os.path.dirname(str(output.assembly)))
        os.symlink(source_relpath,str(output.assembly))

        source_relpath = os.path.relpath(str(input.info),os.path.dirname(str(output.info)))
        os.symlink(source_relpath,str(output.info))


rule assembly:
    input:
        expand("{sample}/assembly/{sample}_assembly.fasta",sample=SAMPLE_NAMES),
        expand("{sample}/assembly/{sample}_circular_info.tsv", sample=SAMPLE_NAMES)
    output:
        temp(touch("checkpoints/assembly"))


rule prepare_medaka_polish_input:
    input:
        "{sample}/assembly/{sample}_assembly.fasta"
    output:
        temp("{sample}/polish/medaka_input/{sample}_input.fasta")
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule polish_medaka:
    input:
        qc_long_reads="{sample}/qc/{sample}_qc_long.fastq.gz",
        contigs="{sample}/{step}/medaka_input/{sample}_input.fasta"
    output:
        dir=directory("{sample}/{step}/medaka"),
        contigs="{sample}/{step}/medaka/{sample}_consensus.fasta"
    conda:
        "../envs/medaka.yaml"
    log:
        "{sample}/logs/{step}/medaka.log"
    benchmark:
        "{sample}/benchmarks/{step}/medaka.txt"
    params:
        medaka_model=config.get("medaka_model"),
        batch_size=config.get("medaka_batch_size")
    threads:
        config.get("threads",1)
    shell:
        """
        medaka_consensus -i {input.qc_long_reads} -d {input.contigs} -o {output.dir} \
          -m {params.medaka_model} -t {threads} -b {params.batch_size} > {log} 2>&1
        mv {output.dir}/consensus.fasta {output.contigs}
        """


rule prepare_polypolish_polish_input:
    input:
        "{sample}/polish/medaka/{sample}_consensus.fasta"
    output:
        "{sample}/polish/polypolish/input/{sample}_input.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule polish_polypolish:
    input:
        qc_short_r1 = "{sample}/qc/{sample}_qc_R1.fastq.gz",
        qc_short_r2 = "{sample}/qc/{sample}_qc_R2.fastq.gz",
        contigs="{sample}/{step}/polypolish/input/{sample}_input.fasta",
        polypolish_filter=os.path.join(DB_DIR_PATH,"polypolish_" + VERSION_POLYPOLISH,"polypolish_insert_filter.py"),
        polypolish=os.path.join(DB_DIR_PATH,"polypolish_" + VERSION_POLYPOLISH,"polypolish"),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","polypolish_" + VERSION_POLYPOLISH)
    output:
        mapping_r1=temp("{sample}/{step}/polypolish/{sample}_R1.sam"),
        mapping_r2=temp("{sample}/{step}/polypolish/{sample}_R2.sam"),
        mapping_clean_r1=temp("{sample}/{step}/polypolish/{sample}_R1.clean.sam"),
        mapping_clean_r2=temp("{sample}/{step}/polypolish/{sample}_R2.clean.sam"),
        polished="{sample}/{step}/polypolish/{sample}_polypolish.fasta",
        debug="{sample}/{step}/polypolish/polypolish.debug.log",
        debug_stats="{sample}/stats/{step}/polypolish_changes.log"
    conda:
        "../envs/mapping.yaml"
    log:
        "{sample}/logs/{step}/polypolish.log"
    benchmark:
        "{sample}/benchmarks/{step}/polypolish.txt"
    threads:
        config.get("threads",1)
    shell:
        """
        printf "\n\n### Read mapping ###\n" > {log}
        bwa index {input.contigs} 2>> {log}
        bwa mem -t {threads} -a {input.contigs} {input.qc_short_r1} > {output.mapping_r1} 2>> {log}
        bwa mem -t {threads} -a {input.contigs} {input.qc_short_r2} > {output.mapping_r2} 2>> {log}
        
        printf "\n\n### Polypolish insert filter ###\n" >> {log}
        {input.polypolish_filter} --in1 {output.mapping_r1} --in2 {output.mapping_r2} \
          --out1 {output.mapping_clean_r1} --out2 {output.mapping_clean_r2} 2>> {log}
          
        printf "\n\n### Polypolish ###\n" >> {log}
        {input.polypolish} {input.contigs} --debug {output.debug} \
          {output.mapping_clean_r1} {output.mapping_clean_r2} 2>> {log} | 
          seqtk seq -A -l 0 | 
          awk \'{{ if ($0 ~ /^>/) {{ gsub("_polypolish", ""); print }} else {{ print }} }}\' | 
          seqtk seq -l 60 > {output.polished} 2>> {log}
          
        head -n 1 {output.debug} > {output.debug_stats}
        grep changed {output.debug} >> {output.debug_stats}
        
        printf "\n\n### Done. ###\n"
        """


# TODO - the relative path workarounds in the shell here are a bit odd because polca outputs files in the present working directory
rule polish_polca:
    input:
        qc_short_r1 = "{sample}/qc/{sample}_qc_R1.fastq.gz",
        qc_short_r2 = "{sample}/qc/{sample}_qc_R2.fastq.gz",
        polished = "{sample}/polish/polypolish/{sample}_polypolish.fasta"
    output:
        polca_output = "{sample}/polish/polca/{sample}_polca.fasta",
        polypolish_sam = temp("{sample}/polish/polca/{sample}_polypolish.fasta.unSorted.sam"),
        polypolish_bam = temp("{sample}/polish/polca/{sample}_polypolish.fasta.alignSorted.bam")
    conda:
        "../envs/masurca.yaml"
    log:
        "{sample}/logs/polish/polca.log"
    benchmark:
        "{sample}/benchmarks/polish/polca.txt"
    params:
        outdir="{sample}/polish/polca"
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        cd {params.outdir}
        polca.sh -a ../../../{input.polished} -r "../../../{input.qc_short_r1} ../../../{input.qc_short_r2}" -t {threads} -m {resources.mem}G > ../../../{log} 2>&1
        ln -s "{wildcards.sample}_polypolish.fasta.PolcaCorrected.fa" "{wildcards.sample}_polca.fasta"
        cd ../../../
        """


# Conditional based on whether short read polishing was performed
rule pre_coverage_filter:
    input:
        "{sample}/polish/medaka/{sample}_consensus.fasta" if POLISH_WITH_SHORT_READS == False else "{sample}/polish/polca/{sample}_polca.fasta"
    output:
        "{sample}/polish/cov_filter/{sample}_pre_filtered.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


# Different coverage methods can be used for the filter: short, long, a combo, or neither (bypass)
filtration_method = []

if (POLISH_WITH_SHORT_READS == True) & \
        ((config.get("meandepth_cutoff_short_read") != "None") | (config.get("evenness_cutoff_short_read") != "None")):
    filtration_method.append("short_read")

    # TODO - consider mapping to medaka polished contigs instead
    rule calculate_short_read_coverage:
        input:
            qc_short_r1 = "{sample}/qc/{sample}_qc_R1.fastq.gz",
            qc_short_r2 = "{sample}/qc/{sample}_qc_R2.fastq.gz",
            contigs = "{sample}/polish/cov_filter/{sample}_pre_filtered.fasta"
        output:
            mapping=temp("{sample}/polish/cov_filter/{sample}_short_read.bam"),
            mapping_index=temp("{sample}/polish/cov_filter/{sample}_short_read.bam.bai"),
            coverage="{sample}/polish/cov_filter/{sample}_short_read_coverage.tsv"
        conda:
            "../envs/mapping.yaml"
        log:
            "{sample}/logs/polish/calculate_short_read_coverage.log"
        benchmark:
            "{sample}/benchmarks/polish/calculate_short_read_coverage.txt"
        threads:
            config.get("threads",1)
        resources:
            mem=int(config.get("memory") / config.get("threads",1))
        shell:
            """
            # Note that -F 4 removes unmapped reads
            bwa index {input.contigs} 2> {log}
            bwa mem -t {threads} {input.contigs} {input.qc_short_r1} {input.qc_short_r2} 2>> {log} | \
              samtools view -b -F 4 -@ {threads} 2>> {log} | \
              samtools sort -@ {threads} -m {resources.mem}G 2>> {log} \
              > {output.mapping}
            samtools index -@ {threads} {output.mapping}
            samtools coverage {output.mapping} > {output.coverage}
            """


if (config.get("meandepth_cutoff_long_read") != "None") | (config.get("evenness_cutoff_long_read") != "None"):
    filtration_method.append("long_read")

    rule calculate_long_read_coverage:
        input:
            contigs="{sample}/polish/cov_filter/{sample}_pre_filtered.fasta",
            qc_long_reads="{sample}/qc/{sample}_qc_long.fastq.gz"
        output:
            mapping=temp("{sample}/polish/cov_filter/{sample}_long_read.bam"),
            mapping_index=temp("{sample}/polish/cov_filter/{sample}_long_read.bam.bai"),
            coverage="{sample}/polish/cov_filter/{sample}_long_read_coverage.tsv"
        conda:
            "../envs/mapping.yaml"
        log:
            "{sample}/logs/polish/calculate_long_read_coverage.log"
        benchmark:
            "{sample}/benchmarks/polish/calculate_long_read_coverage.txt"
        threads:
            config.get("threads",1)
        resources:
            mem=int(config.get("memory") / config.get("threads",1))
        shell:
            """
            # Note that -F 4 removes unmapped reads
            minimap2 -t {threads} -ax map-ont {input.contigs} {input.qc_long_reads} 2> {log} | \
              samtools view -b -F 4 -@ {threads} 2>> {log} | \
              samtools sort -@ {threads} -m {resources.mem}G 2>> {log} \
              > {output.mapping}
            samtools index -@ {threads} {output.mapping}
            samtools coverage {output.mapping} > {output.coverage}
            """


rule summarize_contigs_by_coverage:
    input:
        expand("{{sample}}/polish/cov_filter/{{sample}}_{type}_coverage.tsv",
          type=filtration_method)
    output:
        "{sample}/polish/cov_filter/{sample}_filtered_contigs.list"
    params:
        meandepth_short=config.get("meandepth_cutoff_short_read"),
        evenness_short=config.get("evenness_cutoff_short_read"),
        meandepth_long=config.get("meandepth_cutoff_long_read"),
        evenness_long=config.get("evenness_cutoff_long_read")
    run:
        # Filter a samtools coverage file by meandepth and evenness. Returns a pandas series of the contig names.
        def filter_coverage_data(coverage_file, meandepth, evenness):
            coverage_data = pd.read_csv(coverage_file, sep='\t')

            coverage_filtered = coverage_data[ \
                (coverage_data['meandepth'] >= meandepth) & \
                (coverage_data['coverage'] >= evenness)]

            return(coverage_filtered['#rname'])

        input_list = list(input)

        short_read_coverage_tsv_path = f"{wildcards.sample}/polish/cov_filter/{wildcards.sample}_short_read_coverage.tsv"
        long_read_coverage_tsv_path = f"{wildcards.sample}/polish/cov_filter/{wildcards.sample}_long_read_coverage.tsv"

        if len(input_list) == 1:
            if input_list[0] == short_read_coverage_tsv_path:
                contigs = filter_coverage_data(input_list[0], params.meandepth_short, params.evenness_short)

            elif input_list[0] == long_read_coverage_tsv_path:
                contigs = filter_coverage_data(input_list[0], params.meandepth_long, params.evenness_long)

            else:
                sys.exit("One unexpected coverage file detected in 'polish/cov_filter'.")

        elif len(input_list) == 2:
            input_list.sort()

            if (input_list[0] != long_read_coverage_tsv_path) |\
                    (input_list[1] != short_read_coverage_tsv_path):
                sys.exit("At least one unexpected coverage file detected in 'polish/cov_filter'.")

            set1 = set(filter_coverage_data(input_list[0], params.meandepth_long, params.evenness_long))
            set2 = set(filter_coverage_data(input_list[1], params.meandepth_short, params.evenness_short))

            contigs = pd.Series(set1.union(set2))

        else:
            sys.exit("More than 2 coverage files detected in 'polish/cov_filter'.")

        contigs.to_csv(output[0],header=None,index=False)


if (config.get("meandepth_cutoff_short_read") == "None") & (config.get("evenness_cutoff_short_read") == "None") & \
        (config.get("meandepth_cutoff_long_read") == "None") & (config.get("evenness_cutoff_long_read") == "None"):

    rule bypass_coverage_filter:
        input:
            "{sample}/polish/cov_filter/{sample}_pre_filtered.fasta"
        output:
            "{sample}/polish/cov_filter/{sample}_filtered_contigs.fasta"
        run:
            source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
            os.symlink(source_relpath,str(output))

else:
    rule filter_contigs_by_coverage:
        input:
            contigs="{sample}/polish/medaka/{sample}_consensus.fasta" if POLISH_WITH_SHORT_READS == False else "{sample}/polish/polca/{sample}_polca.fasta",
            filter_list="{sample}/polish/cov_filter/{sample}_filtered_contigs.list"
        output:
            "{sample}/polish/cov_filter/{sample}_filtered_contigs.fasta"
        conda:
            "../envs/mapping.yaml"
        shell:
            """
            seqtk subseq -l 60 {input.contigs} {input.filter_list} > {output}
            """


rule symlink_polish:
        input:
            "{sample}/polish/cov_filter/{sample}_filtered_contigs.fasta"
        output:
            "{sample}/polish/{sample}_polish.fasta"
        run:
            source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
            os.symlink(source_relpath,str(output))


rule polish:
    input:
        expand("{sample}/polish/{sample}_polish.fasta",sample=SAMPLE_NAMES)
    output:
        temp(touch("checkpoints/polish"))


# Writes circular.list with the names of circular contigs if there are any circular contigs
# Writes linear.list with the names of linear contigs if there are any linear contigs
# Then, the DAG is re-evaluated. Circularization is only run if there are circular contigs.
# Based on clustering tutorial at https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html (accessed 2022.3.31)
checkpoint split_circular_and_linear_contigs:
    input:
        assembly_stats="{sample}/assembly/{sample}_circular_info.tsv",
        filter_list="{sample}/polish/cov_filter/{sample}_filtered_contigs.list"
    output:
        directory("{sample}/circularize/filter/lists")
    run:
        coverage_filtered_contigs = pd.read_csv(input.filter_list, header=None)[0]

        assembly_info = pd.read_csv(input.assembly_stats, sep='\t')
        assembly_info_filtered = assembly_info[assembly_info['contig'].isin(coverage_filtered_contigs)]

        circular_contigs = assembly_info_filtered[assembly_info_filtered['circular'] == 'Y']
        linear_contigs = assembly_info_filtered[assembly_info_filtered['circular'] == 'N']

        os.makedirs(output[0], exist_ok=True)

        # Only output files if there is >=1 entry
        if circular_contigs.shape[0] >= 1:
            circular_contigs['contig'].to_csv(os.path.join(output[0], 'circular.list'), header=None, index=False)

        if linear_contigs.shape[0] >= 1:
            linear_contigs['contig'].to_csv(os.path.join(output[0], 'linear.list'), header=None, index=False)


# Makes separate files for circular and linear contigs as needed
rule get_polished_contigs:
    input:
        contigs="{sample}/polish/{sample}_polish.fasta",
        list="{sample}/circularize/filter/lists/{status}.list"
    output:
        "{sample}/circularize/filter/{sample}_{status}.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        seqtk subseq -l 60 {input.contigs} {input.list} > {output}
        """


rule search_contig_start:
    input:
        contigs="{sample}/circularize/filter/{sample}_circular.fasta",
        hmm=os.path.join(DB_DIR_PATH,"hmm",START_HMM_NAME + ".hmm")
    output:
        orf_predictions=temp("{sample}/circularize/identify/{sample}_circular.faa"),
        gene_predictions=temp("{sample}/circularize/identify/{sample}_circular.ffn"),
        annotation_gff=temp("{sample}/circularize/identify/{sample}_circular.gff"),
        search_hits="{sample}/circularize/identify/{sample}_hmmsearch_hits.txt",
        search_hits_no_comments=temp("{sample}/circularize/identify/{sample}_hmmsearch_hits_no_comments.txt")
    conda:
        "../envs/mapping.yaml"
    log:
        "{sample}/logs/circularize/search_contig_start.log"
    benchmark:
        "{sample}/benchmarks/circularize/search_contig_start.txt"
    params:
        hmmsearch_evalue=config.get("hmmsearch_evalue")
    threads:
        config.get("threads",1)
    shell:
        """
        printf "### Predict genes ###\n" > {log}
        prodigal -i {input.contigs} -a {output.orf_predictions} -d {output.gene_predictions} \
          -f gff -o {output.annotation_gff} 2>> {log}
        
        printf "\n\n### Find HMM hits ###\n\n" >> {log}
        hmmsearch --cpu {threads} -E {params.hmmsearch_evalue} --tblout {output.search_hits} -o /dev/stdout \
          {input.hmm} {output.orf_predictions} >> {log} 2>&1
          
        grep -v "^#" {output.search_hits} > {output.search_hits_no_comments}
        
        printf "\n\n### Done. ###\n" >> {log}
        """


rule process_start_genes:
    input:
        "{sample}/circularize/identify/{sample}_hmmsearch_hits_no_comments.txt"
    output:
        "{sample}/circularize/identify/{sample}_start_genes.list"
    log:
        "{sample}/logs/circularize/process_start_genes.log"
    run:
        with open(input[0], 'r') as hmmsearch_results_raw:
            line_count = len(hmmsearch_results_raw.readlines())

        if line_count == 0:
            # Write empty output if there were no HMM hits
            start_orf_ids = []

        else:
            # Load HMM search results
            hmmsearch_results = pd.read_csv(input[0], sep='\s+', header=None)[[0, 2, 3, 4]]

            hmmsearch_results.columns = ['orf', 'hmm_name', 'hmm_accession', 'evalue']

            hmmsearch_results['contig'] = hmmsearch_results['orf'] \
                .str.split('_',expand=False) \
                .apply(lambda x: x[:-1]) \
                .str.join('_')

            contig_counts = pd.DataFrame(hmmsearch_results['contig'].value_counts()) \
                .reset_index() \
                .rename(columns={'index': 'contig', 'contig': 'count'})

            # If one contig has multiple HMM hits, throw a warning and take the hit with lowest e-value
            start_orf_ids = []

            for index, row in contig_counts.iterrows():
                contig, count = row

                if count == 1:
                    start_id = hmmsearch_results[hmmsearch_results['contig'] == contig]['orf'].to_list()[0]
                    start_orf_ids.append(start_id)

                elif count > 1:
                    multistart_orfs = hmmsearch_results[hmmsearch_results['contig'] == contig] \
                        .sort_values(by='evalue',ascending=True) \
                        .reset_index(drop=True)

                    multistart_orf_ids = multistart_orfs['orf'].to_list()

                    print("WARNING: More than one possible start gene on contig '" + str(contig) + "': '" + \
                          ', '.join(multistart_orf_ids) + "'. Will take the hit with lowest e-value.")

                    # Warn if there are ties
                    if multistart_orfs['evalue'][0] == multistart_orfs['evalue'][1]:

                        tie_orf_ids = multistart_orfs[multistart_orfs['evalue'] == multistart_orfs['evalue'][0]][
                            'orf'].to_list()

                        print('WARNING: Two or more genes tied for the lowest e-value for this contig. ' + \
                              'This implies there could be something unusual going on with your data.')
                        print('         All genes with lowest evalue will be output for use by circlator for this contig.')
                        print('         Gene IDs: ' + ', '.join(tie_orf_ids))

                        start_orf_ids = start_orf_ids + tie_orf_ids

                    else:
                        # Relies on the dataframe having already been sorted by e-value above
                        start_orf_ids.append(multistart_orf_ids[0])

        pd.Series(start_orf_ids).to_csv(output[0], header=None, index=False)


rule get_start_genes:
    input:
        gene_predictions="{sample}/circularize/identify/{sample}_circular.ffn",
        start_gene_list="{sample}/circularize/identify/{sample}_start_genes.list"
    output:
        "{sample}/circularize/identify/{sample}_start_gene.ffn"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        seqtk subseq {input.gene_predictions} {input.start_gene_list} > {output}
        """


rule run_circlator:
    input:
        contigs="{sample}/circularize/filter/{sample}_circular.fasta",
        start_gene="{sample}/circularize/identify/{sample}_start_gene.ffn"
    output:
        rotated="{sample}/circularize/circlator/{sample}_rotated.fasta",
        circlator_dir=directory('{sample}/circularize/circlator')
    conda:
        "../envs/circlator.yaml"
    log:
        "{sample}/logs/circularize/circlator.log"
    benchmark:
        "{sample}/benchmarks/circularize/circlator.txt"
    params:
        min_id=90,
    shell:
        """
        if [[ -s {input.start_gene} ]]; then
        
          printf "## Running circlator with custom start gene:\n" > {log}
          circlator fixstart --min_id {params.min_id} --genes_fa {input.start_gene}\
            {input.contigs} {output.circlator_dir}/rotated >> {log} 2>&1
            
        else
        
          ## TODO - I think by default circlator might search for DnaA via its own builtins.
          #         It would be better to skip that step and just rotate everything around to a gene start on the other side of the contig.
          printf "## Start gene file is empty, so will use circlator defaults\n\n" > {log}
          circlator fixstart --min_id {params.min_id} \
            {input.contigs} {output.circlator_dir}/rotated >> {log} 2>&1
            
        fi
        
        mv {output.circlator_dir}/rotated.fasta {output.rotated}
        
        printf "### Circlator log output ###\n" >> {log}
        cat "{wildcards.sample}/circularize/circlator/rotated.log" >> {log}
        """


# Points to the main medaka rule (polish_medaka) above
rule prepare_medaka_circularize_input:
    input:
        "{sample}/circularize/circlator/{sample}_rotated.fasta"
    output:
        temp("{sample}/circularize/medaka_input/{sample}_input.fasta")
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


# Points to the main polypolish rule (polish_polypolish) above
rule prepare_polypolish_circularize_input:
    input:
        "{sample}/circularize/circlator/{sample}_rotated.fasta"
    output:
        "{sample}/circularize/polypolish/input/{sample}_input.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


# Determines whether a second round of long vs. short read polishing is performed
rule finalize_circular_contig_rotation:
    input:
        "{sample}/circularize/medaka/{sample}_consensus.fasta" if POLISH_WITH_SHORT_READS == False else "{sample}/circularize/polypolish/{sample}_polypolish.fasta"
    output:
        "{sample}/circularize/combine/{sample}_circular.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule bypass_circularization:
    input:
        "{sample}/circularize/filter/{sample}_linear.fasta"
    output:
        "{sample}/circularize/combine/{sample}_linear.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


# Gets the names of all lists in circularize/filter/lists (should be either circular.list or linear.list or both)
# Then outputs the expected paths of the finalized circular and linear FastA files
# This function allows the DAG to figure out whether to run the circular / linear specific processing steps
#   based on the split_circular_and_linear_contigs checkpoint made earlier.
def aggregate_contigs(wildcards):
    # TODO - I do not further use this variable, but checkpoints needs to be called to trigger the checkpoint.
    # Am I doing something wrong?
    checkpoint_output = checkpoints.split_circular_and_linear_contigs.get(**wildcards).output[0]
    circularize_lists_path = f"{wildcards.sample}/circularize/filter/lists"

    return expand("{{sample}}/circularize/combine/{{sample}}_{circular_or_linear}.fasta",
                  circular_or_linear=glob_wildcards(os.path.join(circularize_lists_path, "{i}.list")).i)


# TODO - consider sorting contigs by length (or by circular and then by length)
rule combine_circular_and_linear_contigs:
    input:
        aggregate_contigs
    output:
        "{sample}/circularize/combine/{sample}_combined.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        cat {input} | seqtk seq -l 60 > {output}
        """


rule symlink_circularization:
    input:
        "{sample}/circularize/combine/{sample}_combined.fasta"
    output:
        "{sample}/circularize/{sample}_circularize.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule circularize:
    input:
        expand("{sample}/circularize/{sample}_circularize.fasta",sample=SAMPLE_NAMES)
    output:
        temp(touch("checkpoints/circularize"))


# TODO - I might need to remove special characters from strain name to use for locus tag prefix
# TODO - can I auto-predict genome completeness, names, types, topologies?
rule run_dfast:
    input:
        contigs="{sample}/circularize/{sample}_circularize.fasta",
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","dfast_" + VERSION_DFAST)
    output:
        dfast_genome = "{sample}/annotation/dfast/{sample}_genome.fna",
        dfast_proteins = "{sample}/annotation/dfast/{sample}_protein.faa",
        outdir = directory("{sample}/annotation/dfast")
    conda:
        "../envs/annotation_dfast.yaml"
    log:
        "{sample}/logs/annotation/annotation_dfast.log"
    benchmark:
        "{sample}/benchmarks/annotation/annotation_dfast.txt"
    params:
        db=directory(os.path.join(DB_DIR_PATH,"dfast_" + VERSION_DFAST)),
        strain='{sample}'
    threads:
        config.get("threads",1)
    shell:
        """
        dfast --force \
          --dbroot {params.db} \
          -g {input.contigs} \
          -o {output.outdir} \
          --strain {params.strain} \
          --locus_tag_prefix {params.strain} \
          --cpu {threads} > {log} 2>&1
         # --complete t \
         # --seq_names "Chromosome,unnamed" \
         # --seq_types "chromosome,plasmid" \
         # --seq_topologies "circular,circular" 
         
         mv {output.outdir}/genome.fna {output.dfast_genome}
         mv {output.outdir}/protein.faa {output.dfast_proteins}
         """


# TODO - add option to control whether --dbmem flag is set (uses more RAM but does faster analysis)
rule run_eggnog:
    input:
        protein="{sample}/annotation/dfast/{sample}_protein.faa",
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","eggnog_" + VERSION_EGGNOG)
    output:
        eggnog_annotations="{sample}/annotation/eggnog/{sample}.emapper.annotations",
        outdir= directory("{sample}/annotation/eggnog")
    conda:
        "../envs/eggnog.yaml"
    log:
        "{sample}/logs/annotation/eggnog.log"
    benchmark:
        "{sample}/benchmarks/annotation/eggnog.txt"
    params:
        prefix='{sample}',
        db=directory(os.path.join(DB_DIR_PATH,"eggnog_" + VERSION_EGGNOG)),
        sensmode=config.get("eggnog_sensmode"),
        search_tool=config.get('eggnog_search_tool')
    threads:
        config.get("threads",1)
    shell:
        """
        mkdir -p {output.outdir}/tmp
        emapper.py --cpu {threads} -i {input.protein} --itype proteins -m {params.search_tool} \
          --sensmode {params.sensmode} --dbmem --output eggnog --output_dir {output.outdir} \
          --temp_dir {output.outdir}/tmp --output {params.prefix} \
          --data_dir {params.db} --override > {log} 2>&1
        rm -r {output.outdir}/tmp
        """


rule run_gtdbtk:
    input:
        genome="{sample}/annotation/dfast/{sample}_genome.fna",
        setup_finished=os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_validate"),
        ref_msh_file=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE + '_mash','gtdb_ref_sketch.msh')
    output:
        batchfile=temp("{sample}/annotation/gtdbtk/batchfile.tsv"),
        annotation="{sample}/annotation/gtdbtk/{sample}_gtdbtk.summary.tsv",
        outdir=directory("{sample}/annotation/gtdbtk/run_files")
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "{sample}/logs/annotation/gtdbtk.log"
    benchmark:
        "{sample}/benchmarks/annotation/gtdbtk.txt"
    params:
        db_dir=directory(os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE)),
        genome_id='{sample}', # Get sample name from wildcards.
        gtdbtk_mode="--full_tree" if config.get("gtdbtk_mode") == "full_tree" else ""
    threads:
        config.get("threads",1)
    shell:
        """
        printf "{input.genome}\t{params.genome_id}\n" > {output.batchfile}
        GTDBTK_DATA_PATH={params.db_dir} \
          gtdbtk classify_wf --batchfile {output.batchfile} --out_dir {output.outdir} {params.gtdbtk_mode} \
            --mash_db {input.ref_msh_file} --cpus {threads} --pplacer_cpus {threads} > {log} 2>&1
        head -n 1 {output.outdir}/gtdbtk.*.summary.tsv | sort -u > {output.annotation}
        tail -n +2 {output.outdir}/gtdbtk.*.summary.tsv >> {output.annotation}
        """

rule run_checkm2:
    input:
        genome="{sample}/annotation/dfast/{sample}_genome.fna",
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","checkm2")
    output:
        quality_report="{sample}/annotation/checkm/{sample}_checkm_quality_report.tsv",
        outdir=directory("{sample}/annotation/checkm/")
    log:
        "{sample}/logs/annotation/checkm.log"
    benchmark:
        "{sample}/benchmarks/annotation/checkm.txt"
    conda:
        "../envs/checkm2.yaml"
    params:
        db=directory(os.path.join(DB_DIR_PATH,"checkm2"))
    threads:
        config.get("threads",1)
    shell:
        """
        mkdir -p {output.outdir}
        checkm2 predict --threads {threads} --input {input.genome} --output-directory {output.outdir} > {log} 2>&1
        mv {output.outdir}/quality_report.tsv {output.quality_report}
        """


if POLISH_WITH_SHORT_READS == True:

    # TODO - clarify name compared to previous mapping step
    rule calculate_final_short_read_coverage:
        input:
            qc_short_r1 = "{sample}/qc/{sample}_qc_R1.fastq.gz",
            qc_short_r2 = "{sample}/qc/{sample}_qc_R2.fastq.gz",
            dfast_genome = "{sample}/annotation/dfast/{sample}_genome.fna"
        output:
            mapping="{sample}/annotation/coverage/{sample}_short_read.bam",
            index="{sample}/annotation/coverage/{sample}_short_read.bam.bai",
            coverage="{sample}/annotation/coverage/{sample}_short_read_coverage.tsv"
        conda:
            "../envs/mapping.yaml"
        log:
            "{sample}/logs/annotation/calculate_final_short_read_coverage.log"
        benchmark:
            "{sample}/benchmarks/annotation/calculate_final_short_read_coverage.txt"
        threads:
            config.get("threads",1)
        resources:
            mem=int(config.get("memory") / config.get("threads",1))
        shell:
            """
            bwa index {input.dfast_genome} 2> {log}
            bwa mem -t {threads} {input.dfast_genome} {input.qc_short_r1} {input.qc_short_r2} 2>> {log} | \
              samtools view -b -@ {threads} 2>> {log} | \
              samtools sort -@ {threads} -m {resources.mem}G 2>> {log} \
              > {output.mapping}
            samtools index -@ {threads} {output.mapping}
            samtools coverage {output.mapping} > {output.coverage}
            """


rule calculate_final_long_read_coverage:
    input:
        contigs="{sample}/annotation/dfast/{sample}_genome.fna",
        qc_long_reads="{sample}/qc/{sample}_qc_long.fastq.gz"
    output:
        mapping="{sample}/annotation/coverage/{sample}_long_read.bam",
        index="{sample}/annotation/coverage/{sample}_long_read.bam.bai",
        coverage="{sample}/annotation/coverage/{sample}_long_read_coverage.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "{sample}/logs/annotation/calculate_final_long_read_coverage.log"
    benchmark:
        "{sample}/benchmarks/annotation/calculate_final_long_read_coverage.txt"
    threads:
        config.get("threads",1)
    resources:
        mem=int(config.get("memory") / config.get("threads",1))
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.contigs} {input.qc_long_reads} 2> {log} | \
          samtools view -b -@ {threads} 2>> {log} | \
          samtools sort -@ {threads} -m {resources.mem}G 2>> {log} \
          > {output.mapping}
        samtools index -@ {threads} {output.mapping}
        samtools coverage {output.mapping} > {output.coverage}
        """


rule symlink_logs:
    input:
        long_read_coverage="{sample}/annotation/coverage/{sample}_long_read_coverage.tsv",
    output:
        logs=temp(directory("{sample}/annotation/logs")),
        stats=temp(directory("{sample}/annotation/stats"))
    params:
        logs="{sample}/logs",
        stats="{sample}/stats"
    run:
        source_relpath = os.path.relpath(str(params.logs),os.path.dirname(str(output.logs)))
        os.symlink(source_relpath, str(output.logs))

        source_relpath = os.path.relpath(str(params.stats),os.path.dirname(str(output.stats)))
        os.symlink(source_relpath, str(output.stats))


# TODO - can I remove the use of cd?
# TODO - I currently avoid writing to the log directory during zip because the folder is being zipped,
#        but the resulting code seems a bit unnatural
rule summarize_annotation:
    input:
        "{sample}/annotation/dfast/{sample}_genome.fna",
        "{sample}/annotation/eggnog/{sample}.emapper.annotations",
        "{sample}/annotation/gtdbtk/{sample}_gtdbtk.summary.tsv",
        "{sample}/annotation/checkm/",
        expand("{{sample}}/annotation/coverage/{{sample}}_{type}_coverage.tsv",
            type=["short_read", "long_read"] if POLISH_WITH_SHORT_READS == True else ["long_read"]),
        "{sample}/annotation/logs",
        "{sample}/annotation/stats"
    output:
        "{sample}/{sample}_annotation_summary.zip"
    log:
        "{sample}/logs/annotation/summarize_annotation.log"
    params:
        zipdir="{sample}/annotation"
    shell:
        """
        cd {params.zipdir}
        zip -r ../../{output} * -x \*.bam\* gtdbtk/run_files/\* > "../../summarize_annotation.log" 2>&1
        cd ../../
        mv "summarize_annotation.log" {log}
        """


rule annotation:
    input:
        expand("{sample}/{sample}_annotation_summary.zip",sample=SAMPLE_NAMES),
    output:
        temp(touch("checkpoints/annotation"))


# TODO - add nice summaries of reads removed during QC, polishing stats, assembly stats, circular vs. non.
#        These aren't in the final ZIP summary at the moment.
