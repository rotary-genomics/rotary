# polish: rules for long read polishing.
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os
import sys

import pandas as pd

VERSION_POLYPOLISH="0.5.0"

DB_DIR_PATH = config.get('db_dir')

# SAMPLE_NAMES and POLISH_WITH_SHORT_READS are instantiated in rotary.smk

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
        temp("{sample}/polish/polypolish/input/{sample}_input.fasta")
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule polish_polypolish:
    input:
        qc_short_r1="{sample}/qc/{sample}_qc_R1.fastq.gz",
        qc_short_r2="{sample}/qc/{sample}_qc_R2.fastq.gz",
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
        debug=temp("{sample}/{step}/polypolish/polypolish.debug.log"),
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
        polca_output = temp("{sample}/polish/polca/{sample}_polca.fasta")
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
        temp("{sample}/polish/cov_filter/{sample}_pre_filtered.fasta")
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
            qc_short_r1="{sample}/qc/{sample}_qc_R1.fastq.gz",
            qc_short_r2="{sample}/qc/{sample}_qc_R2.fastq.gz",
            contigs="{sample}/polish/cov_filter/{sample}_pre_filtered.fasta"
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
        #
        def filter_coverage_data(coverage_file, mean_depth, evenness):
            """
            Filter a samtools coverage file by mean depth and evenness.
            Returns a pandas series of passing contig names.

            :param coverage_file: Path to the coverage file (CSV format).
            :param mean_depth: Minimum required mean depth of coverage.
            :param evenness: Minimum required evenness of coverage.
            :return: Filtered coverage data containing only passing contig names.
            """
            coverage_data = pd.read_csv(coverage_file,sep='\t')

            coverage_filtered = coverage_data[ \
                (coverage_data['meandepth'] >= mean_depth) & \
                (coverage_data['coverage'] >= evenness)]

            return (coverage_filtered['#rname'])


        input_list = list(input)

        short_read_coverage_tsv_path = f"{wildcards.sample}/polish/cov_filter/{wildcards.sample}_short_read_coverage.tsv"
        long_read_coverage_tsv_path = f"{wildcards.sample}/polish/cov_filter/{wildcards.sample}_long_read_coverage.tsv"

        if len(input_list) == 1:
            if input_list[0] == short_read_coverage_tsv_path:
                contigs = filter_coverage_data(input_list[0],params.meandepth_short,params.evenness_short)

            elif input_list[0] == long_read_coverage_tsv_path:
                contigs = filter_coverage_data(input_list[0],params.meandepth_long,params.evenness_long)

            else:
                sys.exit("One unexpected coverage file detected in 'polish/cov_filter'.")

        elif len(input_list) == 2:
            input_list.sort()

            if (input_list[0] != long_read_coverage_tsv_path) | \
                    (input_list[1] != short_read_coverage_tsv_path):
                sys.exit("At least one unexpected coverage file detected in 'polish/cov_filter'.")

            set1 = set(filter_coverage_data(input_list[0],params.meandepth_long,params.evenness_long))
            set2 = set(filter_coverage_data(input_list[1],params.meandepth_short,params.evenness_short))

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
