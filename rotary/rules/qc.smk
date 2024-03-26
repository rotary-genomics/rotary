# qc: rules for long and short read quality control.
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os
import itertools
import pandas as pd
from pungi.utils import symlink_or_compress, is_config_parameter_true

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

ZENODO_VERSION = "10087395"

DB_DIR_PATH = config.get('db_dir')

KEEP_QC_READ_FILES = is_config_parameter_true(config,'keep_final_qc_read_files')

# SAMPLE_NAMES and POLISH_WITH_SHORT_READS are instantiated in rotary.smk

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
        genome=os.path.join(DB_DIR_PATH,'contamination_references','ncbi','{accession}.fna.gz'),
        zip=temp(os.path.join(DB_DIR_PATH,'contamination_references','ncbi','{accession}.zip')),
        zip_dir=temp(directory(os.path.join(DB_DIR_PATH,'contamination_references','ncbi','{accession}')))
    conda:
        "../envs/download.yaml"
    log:
        "logs/download/download_contamination_reference_{accession}.log"
    benchmark:
        "benchmarks/download/download_contamination_reference_{accession}.benchmark.txt"
    params:
        accession="{accession}"
    threads:
        min(config.get("threads",1),4)
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

rule nanopore_qc_filter:
    input:
        "{sample}/raw/{sample}_long.fastq.gz"
    output:
        temp("{sample}/qc/long/{sample}_nanopore_qc.fastq.gz")
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
          interleaved=f qin=33 threads={threads} -Xmx{resources.mem}G pigz=t unpigz=t > {log} 2>&1
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
        reformat.sh in={input} lhist={output} maxhistlen=10000000 pigz=t unpigz=t \
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
        "{sample}/qc/{sample}_qc_long.fastq.gz" if KEEP_QC_READ_FILES else temp("{sample}/qc/{sample}_qc_long.fastq.gz")
    shell:
        """
        cp {input} {output}
        """


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
        short_final_r1="{sample}/qc/{sample}_qc_R1.fastq.gz" if KEEP_QC_READ_FILES else temp("{sample}/qc/{sample}_qc_R1.fastq.gz"),
        short_final_r2="{sample}/qc/{sample}_qc_R2.fastq.gz" if KEEP_QC_READ_FILES else temp("{sample}/qc/{sample}_qc_R2.fastq.gz")
    shell:
        """
        cp {input.short_r1} {output.short_final_r1}
        cp {input.short_r2} {output.short_final_r2}
        """

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
