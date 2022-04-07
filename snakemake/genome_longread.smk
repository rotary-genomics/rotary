# Snakemake rules long read-based genome assembly workflow
# Copyright Jackson M. Tsuji, ILTS, Hokkaido University, 2022

import os
import sys
import shutil
import pandas as pd
import itertools
from snakemake.utils import logger, min_version, update_config

VERSION="0.1.0"
VERSION_POLYPOLISH="0.5.0"
VERSION_DFAST="1.2.15"
VERSION_EGGNOG="5.0.0" # See http://eggnog5.embl.de/#/app/downloads
VERSION_GTDB="202" # See https://data.gtdb.ecogenomic.org/releases/

# Specify the minimum snakemake version allowable
min_version("6.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")


rule all:
    input:
        "checkpoints/qc_long",
        "checkpoints/assembly",
        "checkpoints/polish",
        "checkpoints/circularize",
        "checkpoints/annotation"


rule install_internal_scripts:
    output:
        end_repair=os.path.join(config.get("db_dir"), "nanopore-workflows-" + VERSION, "scripts", "flye_end_repair.sh"),
        end_repair_utils=os.path.join(config.get("db_dir"), "nanopore-workflows-" + VERSION, "scripts", "flye_end_repair_utils.py"),
        install_finished=os.path.join(config.get("db_dir"), "checkpoints", "internal_scripts_" + VERSION)
    log:
        "logs/download/install_internal_scripts.log"
    benchmark:
        "benchmarks/download/install_internal_scripts.txt"
    params:
        db_dir=config.get("db_dir"),
        url="https://github.com/jmtsuji/nanopore-workflows/archive/refs/tags/" + VERSION + ".tar.gz"
    shell:
        """
        mkdir -p {params.db_dir}
        wget -O - {params.url} 2> {log} | tar -C {params.db_dir} -xzf - >> {log} 2>&1
        touch {output.install_finished}
        """


rule install_polypolish:
    output:
        polypolish_filter=os.path.join(config.get("db_dir"), "polypolish_" + VERSION_POLYPOLISH, "polypolish_insert_filter.py"),
        polypolish=os.path.join(config.get("db_dir"), "polypolish_" + VERSION_POLYPOLISH, "polypolish"),
        install_finished=os.path.join(config.get("db_dir"), "checkpoints", "polypolish_" + VERSION_POLYPOLISH)
    log:
        "logs/download/install_polypolish.log"
    benchmark:
        "benchmarks/download/install_polypolish.txt"
    params:
        db_dir=os.path.join(config.get("db_dir"), "polypolish_" + VERSION_POLYPOLISH),
        url="https://github.com/rrwick/Polypolish/releases/download/v0.5.0/polypolish-linux-x86_64-musl-v" + VERSION_POLYPOLISH + ".tar.gz"
    shell:
        """
        mkdir -p {params.db_dir}
        wget -O - {params.url} 2> {log} | tar -C {params.db_dir} -xzf - >> {log} 2>&1
        touch {output.install_finished}
        """


# TODO - does not check the HMM version, only ID. If the HMM version updates, it won't automatically re-download
rule download_hmm:
    output:
        hmm=os.path.join(config.get("db_dir"), "hmm", config.get("start_hmm_pfam_id") + ".hmm")
    log:
        "logs/download/hmm_download.log"
    benchmark:
        "benchmarks/download/hmm_download.txt"
    params:
        db_dir=os.path.join(config.get("db_dir"), "hmm"),
        url="https://pfam.xfam.org/family/" + config.get("start_hmm_pfam_id") + "/hmm"
    shell:
        """
        mkdir -p {params.db_dir}
        wget -O {output.hmm} {params.url} 2> {log}
        # --no-check-certificate
        """


rule download_dfast_db:
    output:
        db=directory(os.path.join(config.get("db_dir"), "dfast_" + VERSION_DFAST)),
        install_finished=os.path.join(config.get("db_dir"), "checkpoints", "dfast_" + VERSION_DFAST)
    conda:
        "../envs/annotation_dfast.yaml"
    log:
        "logs/download/dfast_db_download.log"
    benchmark:
        "benchmarks/download/dfast_db_download.txt"
    params:
        db_dir=os.path.join(config.get("db_dir"), "dfast_" + VERSION_DFAST)
    shell:
        """
        mkdir -p {params.db_dir}
        dfast_file_downloader.py --protein dfast --dbroot {params.db_dir} > {log} 2>&1
        dfast_file_downloader.py --cdd Cog --hmm TIGR --dbroot {params.db_dir} >> {log} 2>&1
        touch {output.install_finished}
        """


rule download_eggnog_db:
    output:
        db=directory(os.path.join(config.get("db_dir"), "eggnog_" + VERSION_EGGNOG)),
        install_finished=os.path.join(config.get("db_dir"), "checkpoints", "eggnog_" + VERSION_EGGNOG)
    conda:
        "../envs/eggnog.yaml"
    log:
        "logs/download/eggnog_db_download.log"
    benchmark:
        "benchmarks/download/eggnog_db_download.txt"
    params:
        db_dir=os.path.join(config.get("db_dir"), "eggnog_" + VERSION_EGGNOG)
    shell:
        """
        mkdir -p {params.db_dir}
        download_eggnog_data.py -y --data_dir {params.db_dir} > {log} 2>&1
        touch {output.install_finished}
        """


rule download_gtdb_db:
    output:
        db=directory(os.path.join(config.get("db_dir"), "GTDB_" + VERSION_GTDB)),
        install_finished=os.path.join(config.get("db_dir"), "checkpoints", "GTDB_" + VERSION_GTDB)
    log:
        "logs/download/gtdb_db_download.log"
    benchmark:
        "benchmarks/download/gtdb_db_download.txt"
    params:
        db_dir_root=os.path.join(config.get("db_dir")),
        initial_download_dir=os.path.join(config.get("db_dir"), "release" + VERSION_GTDB),
        db_dir=os.path.join(config.get("db_dir"), "GTDB_" + VERSION_GTDB),
        url="https://data.gtdb.ecogenomic.org/releases/release" + VERSION_GTDB + "/" + VERSION_GTDB + ".0/auxillary_files/gtdbtk_r" + VERSION_GTDB + "_data.tar.gz"
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


rule nanopore_qc_filter:
    input:
        config.get("longreads")
    output:
        "qc_long/nanopore_qc.fastq.gz"
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/qc_long.log"
    benchmark:
        "benchmarks/qc_long.benchmark.txt"
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
        "qc_long/nanopore_qc.fastq.gz"
    output:
        "qc_long/length_hist.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/qc_long_length_hist.log"
    benchmark:
        "benchmarks/qc_long_length_hist.benchmark.txt"
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
        "qc_long/length_hist.tsv"
    output:
        "stats/qc_long_length_stats.txt"
    log:
        "logs/qc_long_length_stats.log"
    benchmark:
        "benchmarks/qc_long_length_stats.benchmark.txt"
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


rule qc_long:
    input:
        "stats/qc_long_length_stats.txt"
    output:
        temp(touch("checkpoints/qc_long"))


rule assembly_flye:
    input:
        "qc_long/nanopore_qc.fastq.gz"
    output:
        "assembly/flye/assembly.fasta",
        "assembly/flye/assembly_info.txt"
    conda:
        "../envs/assembly_flye.yaml"
    log:
        "logs/assembly_flye.log"
    benchmark:
        "benchmarks/assembly_flye.benchmark.txt"
    params:
        output_dir="assembly/flye",
        input_mode=config.get("flye_input_mode"),
        meta_mode="--meta" if config.get("flye_meta_mode") == "True" else "",
        polishing_rounds=config.get("flye_polishing_rounds")
    threads:
        config.get("threads",1)
    shell:
        """
        flye --{params.input_mode} {input} --out-dir {params.output_dir} {params.meta_mode} \
          --iterations {params.polishing_rounds} -t {threads} > {log} 2>&1
        """


rule assembly_end_repair:
    input:
        qc_long_reads="qc_long/nanopore_qc.fastq.gz",
        assembly="assembly/flye/assembly.fasta",
        info="assembly/flye/assembly_info.txt",
        end_repair=os.path.join(config.get("db_dir"),"nanopore-workflows-" + VERSION,"scripts","flye_end_repair.sh"),
        end_repair_utils=os.path.join(config.get("db_dir"),"nanopore-workflows-" + VERSION,"scripts","flye_end_repair_utils.py"),
        install_finished=os.path.join(config.get("db_dir"),"checkpoints","internal_scripts_" + VERSION)
    output:
        assembly="assembly/end_repair/repaired.fasta",
        info="assembly/end_repair/assembly_info.txt"
    conda:
        "../envs/circlator.yaml"
    log:
        "logs/end_repair.log"
    benchmark:
        "benchmarks/end_repair.txt"
    params:
        output_dir="assembly/end_repair",
        flye_input_mode=config.get("flye_input_mode"),
        min_id=config.get("circlator_merge_min_id"),
        min_length=config.get("circlator_merge_min_length"),
        ref_end=config.get("circlator_merge_ref_end"),
        reassemble_end=config.get("circlator_merge_reassemble_end")
    threads:
        config.get("threads",1)
    resources:
        mem=int(config.get("memory") / config.get("threads",1))
    shell:
        """
        {input.end_repair} -f {params.flye_input_mode} -i {params.min_id} -l {params.min_length} \
          -e {params.ref_end} -E {params.reassemble_end} -t {threads} -m {resources.mem} \
          {input.qc_long_reads} {input.assembly} {input.info} {params.output_dir} > {log} 2>&1
        cp {input.info} {output.info}
        """


# TODO - eventually make this rule more generic so that outputs from other assemblers can go to the same output files
#        In particular, the format of assembly_info.txt will need to be standardized (also in assembly_end_repair).
rule finalize_assembly:
    input:
        assembly="assembly/end_repair/repaired.fasta",
        info="assembly/end_repair/assembly_info.txt"
    output:
        assembly="assembly/assembly.fasta",
        info="assembly/assembly_info.txt"
    run:
        source_relpath = os.path.relpath(str(input.assembly),os.path.dirname(str(output.assembly)))
        os.symlink(source_relpath,str(output.assembly))

        source_relpath = os.path.relpath(str(input.info),os.path.dirname(str(output.info)))
        os.symlink(source_relpath,str(output.info))


rule assembly:
    input:
        "assembly/assembly.fasta",
        "assembly/assembly_info.txt"
    output:
        temp(touch("checkpoints/assembly"))


rule polish_medaka:
    input:
        qc_long_reads="qc_long/nanopore_qc.fastq.gz",
        contigs="assembly/assembly.fasta"
    output:
        "polish/medaka/consensus.fasta"
    conda:
        "../envs/medaka.yaml"
    log:
        "logs/medaka.log"
    benchmark:
        "benchmarks/medaka.txt"
    params:
        output_dir="polish/medaka",
        medaka_model=config.get("medaka_model")
    threads:
        config.get("threads",1)
    shell:
        """
        medaka_consensus -i {input.qc_long_reads} -d {input.contigs} -o {params.output_dir} \
          -m {params.medaka_model} -t {threads} > {log} 2>&1
        """


rule prepare_polypolish_polish_input:
    input:
        "polish/medaka/consensus.fasta"
    output:
        "polish/polypolish/input.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule polish_polypolish:
    input:
        contigs="{step}/polypolish/input.fasta",
        polypolish_filter=os.path.join(config.get("db_dir"),"polypolish_" + VERSION_POLYPOLISH,"polypolish_insert_filter.py"),
        polypolish=os.path.join(config.get("db_dir"),"polypolish_" + VERSION_POLYPOLISH,"polypolish"),
        install_finished=os.path.join(config.get("db_dir"),"checkpoints","polypolish_" + VERSION_POLYPOLISH)
    output:
        mapping_r1=temp("{step}/polypolish/R1.sam"),
        mapping_r2=temp("{step}/polypolish/R2.sam"),
        mapping_clean_r1=temp("{step}/polypolish/R1.clean.sam"),
        mapping_clean_r2=temp("{step}/polypolish/R2.clean.sam"),
        polished="{step}/polypolish/polypolish.fasta",
        debug="{step}/polypolish/polypolish.debug.log",
        debug_stats="stats/{step}/polypolish_changes.log"
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/polypolish_{step}.log"
    benchmark:
        "benchmarks/polypolish_{step}.txt"
    params:
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2")
    threads:
        config.get("threads",1)
    shell:
        """
        printf "\n\n### Read mapping ###\n" > {log}
        bwa index {input} 2>> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r1} > {output.mapping_r1} 2>> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r2} > {output.mapping_r2} 2>> {log}
        
        printf "\n\n### Polypolish insert filter ###\n" >> {log}
        {input.polypolish_filter} --in1 {output.mapping_r1} --in2 {output.mapping_r2} \
          --out1 {output.mapping_clean_r1} --out2 {output.mapping_clean_r2} 2>> {log}
          
        printf "\n\n### Polypolish ###\n" >> {log}
        {input.polypolish} {input} --debug {output.debug} \
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
        "polish/polypolish/polypolish.fasta"
    output:
        "polish/polca/polca.fasta",
        temp("polish/polca/polypolish.fasta.unSorted.sam"),
        temp("polish/polca/polypolish.fasta.alignSorted.bam")
    conda:
        "../envs/masurca.yaml"
    log:
        "logs/polca.log"
    benchmark:
        "benchmarks/polca.txt"
    params:
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2"),
        outdir="polish/polca"
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        cd {params.outdir}
        polca.sh -a ../../{input} -r "{params.qc_short_r1} {params.qc_short_r2}" -t {threads} -m {resources.mem}G > ../../{log} 2>&1
        ln -s "polypolish.fasta.PolcaCorrected.fa" "polca.fasta"
        cd ../..
        """


# Conditional based on whether short read polishing was performed
rule pre_coverage_filter:
    input:
        "polish/medaka/consensus.fasta" if config.get("qc_short_r1") == "None" else "polish/polca/polca.fasta"
    output:
        "polish/cov_filter/pre_filtered.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


# Different coverage methods can be used for the filter: short, long, a combo, or neither (bypass)
filtration_method = []

if (config.get("qc_short_r1") != "None") & \
        ((config.get("meandepth_cutoff_short_read") != "None") | (config.get("evenness_cutoff_short_read") != "None")):
    filtration_method.append("short_read")

    # TODO - consider mapping to medaka polished contigs instead
    rule calculate_short_read_coverage:
        input:
            "polish/cov_filter/pre_filtered.fasta"
        output:
            mapping=temp("polish/cov_filter/short_read.bam"),
            mapping_index=temp("polish/cov_filter/short_read.bam.bai"),
            coverage="polish/cov_filter/short_read_coverage.tsv"
        conda:
            "../envs/mapping.yaml"
        log:
            "logs/calculate_short_read_coverage.log"
        benchmark:
            "benchmarks/calculate_short_read_coverage.txt"
        params:
            qc_short_r1=config.get("qc_short_r1"),
            qc_short_r2=config.get("qc_short_r2")
        threads:
            config.get("threads",1)
        resources:
            mem=int(config.get("memory") / config.get("threads",1))
        shell:
            """
            # Note that -F 4 removes unmapped reads
            bwa index {input} 2> {log}
            bwa mem -t {threads} {input} {params.qc_short_r1} {params.qc_short_r2} 2>> {log} | \
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
            contigs="polish/cov_filter/pre_filtered.fasta",
            qc_long_reads="qc_long/nanopore_qc.fastq.gz"
        output:
            mapping=temp("polish/cov_filter/long_read.bam"),
            mapping_index=temp("polish/cov_filter/long_read.bam.bai"),
            coverage="polish/cov_filter/long_read_coverage.tsv"
        conda:
            "../envs/mapping.yaml"
        log:
            "logs/calculate_long_read_coverage.log"
        benchmark:
            "benchmarks/calculate_long_read_coverage.txt"
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
        expand("polish/cov_filter/{type}_coverage.tsv",
          type=filtration_method)
    output:
        "polish/cov_filter/filtered_contigs.list"
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

        if len(input_list) == 1:
            if input_list[0] == "polish/cov_filter/short_read_coverage.tsv":
                contigs = filter_coverage_data(input_list[0], params.meandepth_short, params.evenness_short)

            elif input_list[0] == "polish/cov_filter/long_read_coverage.tsv":
                contigs = filter_coverage_data(input_list[0], params.meandepth_long, params.evenness_long)

            else:
                sys.exit("One unexpected coverage file detected in 'polish/cov_filter'.")

        elif len(input_list) == 2:
            input_list.sort()

            if (input_list[0] != "polish/cov_filter/long_read_coverage.tsv") |\
                    (input_list[1] != "polish/cov_filter/short_read_coverage.tsv"):
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
            "polish/cov_filter/pre_filtered.fasta"
        output:
            "polish/cov_filter/filtered_contigs.fasta"
        run:
            source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
            os.symlink(source_relpath,str(output))

else:

    rule filter_contigs_by_coverage:
        input:
            contigs="polish/medaka/consensus.fasta" if config.get("qc_short_r1") == "None" else "polish/polca/polca.fasta",
            filter_list="polish/cov_filter/filtered_contigs.list"
        output:
            "polish/cov_filter/filtered_contigs.fasta"
        conda:
            "../envs/mapping.yaml"
        shell:
            """
            seqtk subseq -l 60 {input.contigs} {input.filter_list} > {output}
            """


rule symlink_polish:
        input:
            "polish/cov_filter/filtered_contigs.fasta"
        output:
            "polish/polish.fasta"
        run:
            source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
            os.symlink(source_relpath,str(output))


rule polish:
    input:
        "polish/polish.fasta"
    output:
        temp(touch("checkpoints/polish"))


# Writes circular.list with the names of circular contigs if there are any circular contigs
# Writes linear.list with the names of linear contigs if there are any linear contigs
# Then, the DAG is re-evaluated. Circularization is only run if there are circular contigs.
# Based on clustering tutorial at https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html (accessed 2022.3.31)
checkpoint split_circular_and_linear_contigs:
    input:
        assembly_stats="assembly/assembly_info.txt",
        filter_list="polish/cov_filter/filtered_contigs.list"
    output:
        directory("circularize/filter/lists")
    run:
        coverage_filtered_contigs = pd.read_csv(input.filter_list, header=None)[0]

        assembly_info = pd.read_csv(input.assembly_stats, sep='\t')
        assembly_info_filtered = assembly_info[assembly_info['#seq_name'].isin(coverage_filtered_contigs)]

        circular_contigs = assembly_info_filtered[assembly_info_filtered['circ.'] == 'Y']
        linear_contigs = assembly_info_filtered[assembly_info_filtered['circ.'] == 'N']

        os.makedirs(output[0], exist_ok=True)

        # Only output files if there is >1 entry
        if circular_contigs.shape[0] > 1:
            circular_contigs['#seq_name'].to_csv(os.path.join(output[0], 'circular.list'), header=None, index=False)

        if linear_contigs.shape[0] > 1:
            linear_contigs['#seq_name'].to_csv(os.path.join(output[0], 'linear.list'), header=None, index=False)


# Makes separate files for circular and linear contigs as needed
rule get_polished_contigs:
    input:
        contigs="polish/polish.fasta",
        list="circularize/filter/lists/{status}.list"
    output:
        "circularize/filter/{status}.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        seqtk subseq -l 60 {input.contigs} {input.list} > {output}
        """


rule search_contig_start:
    input:
        contigs="circularize/filter/circular.fasta",
        hmm=os.path.join(config.get("db_dir"), "hmm", config.get("start_hmm_pfam_id") + ".hmm")
    output:
        orf_predictions=temp("circularize/identify/circular.faa"),
        gene_predictions=temp("circularize/identify/circular.ffn"),
        annotation_gff=temp("circularize/identify/circular.gff"),
        search_hits="circularize/identify/hmmsearch_hits.txt"
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/search_contig_start.log"
    benchmark:
        "benchmarks/search_contig_start.txt"
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
        
        printf "\n\n### Done. ###\n" >> {log}
        """


rule process_start_genes:
    input:
        "circularize/identify/hmmsearch_hits.txt"
    output:
        "circularize/identify/start_genes.list"
    log:
        "logs/process_start_genes.log"
    run:
        # Load HMM search results
        hmmsearch_results = pd.read_csv(input[0], sep='\s+', header=None, comment='#')[[0, 2, 3, 4]]
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
        gene_predictions="circularize/identify/circular.ffn",
        start_gene_list="circularize/identify/start_genes.list"
    output:
        "circularize/identify/start_gene.ffn"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        seqtk subseq {input.gene_predictions} {input.start_gene_list} > {output}
        """


rule run_circlator:
    input:
        contigs="circularize/filter/circular.fasta",
        start_gene="circularize/identify/start_gene.ffn"
    output:
        "circularize/circlator/rotated.fasta"
    conda:
        "../envs/circlator.yaml"
    log:
        "logs/circlator.log"
    benchmark:
        "benchmarks/circlator.txt"
    params:
        min_id=90,
        run_name="circularize/circlator/rotated"
    shell:
        """
        circlator fixstart --min_id {params.min_id} --genes_fa {input.start_gene}\
          {input.contigs} {params.run_name} > {log} 2>&1
        
        printf "### Circlator log output ###\n" >> {log}
        cat "circularize/circlator/rotated.log" >> {log}
        """


# Points to the main polypolish rule (polish_polypolish) above
rule prepare_polypolish_circularize_input:
    input:
        "circularize/circlator/rotated.fasta"
    output:
        "circularize/polypolish/input.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule finalize_circular_contig_rotation:
    input:
        "circularize/polypolish/polypolish.fasta"
    output:
        "circularize/combine/circular.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule bypass_circularization:
    input:
        "circularize/filter/linear.fasta"
    output:
        "circularize/combine/linear.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


# Gets the names of all lists in circularize/filter/lists (should be either circular.list or linear.list or both)
# Then outputs the expected paths of the finalized circular and linear FastA files
# This function allows the DAG to figure out whether to run the circular / linear specific processing steps
#   based on the split_circular_and_linear_contigs checkpoint made earlier.
def aggregate_contigs(wildcards):
    # TODO - I do not further use this variable, but checkpoints needs to be called to trigger the checkpoint. Am I doing something wrong?
    checkpoint_output = checkpoints.split_circular_and_linear_contigs.get(**wildcards).output[0]

    return expand("circularize/combine/{circular_or_linear}.fasta",
                  circular_or_linear=glob_wildcards(os.path.join("circularize/filter/lists", "{i}.list")).i)


# TODO - consider sorting contigs by length (or by circular and then by length)
rule combine_circular_and_linear_contigs:
    input:
        aggregate_contigs
    output:
        "circularize/combine/combined.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        cat {input} | seqtk seq -l 60 > {output}
        """


rule symlink_circularization:
        input:
            "circularize/combine/combined.fasta"
        output:
            "circularize/circularize.fasta"
        run:
            source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
            os.symlink(source_relpath,str(output))


rule circularize:
    input:
        "circularize/circularize.fasta"
    output:
        temp(touch("checkpoints/circularize"))


# TODO - I might need to remove special characters from strain name to use for locus tag prefix
# TODO - can I auto-predict genome completeness, names, types, topologies?
rule run_dfast:
    input:
        contigs="circularize/circularize.fasta",
        install_finished=os.path.join(config.get("db_dir"),"checkpoints","dfast_" + VERSION_DFAST)
    output:
        "annotation/dfast/genome.fna",
        "annotation/dfast/protein.faa"
    conda:
        "../envs/annotation_dfast.yaml"
    log:
        "logs/annotation_dfast.log"
    benchmark:
        "benchmarks/annotation_dfast.txt"
    params:
        outdir="annotation/dfast",
        db=directory(os.path.join(config.get("db_dir"),"dfast_" + VERSION_DFAST)),
        strain=config.get("sample_id")
    threads:
        config.get("threads",1)
    shell:
        """
        dfast --force \
          --dbroot {params.db} \
          -g {input.contigs} \
          -o {params.outdir} \
          --strain {params.strain} \
          --locus_tag_prefix {params.strain} \
          --cpu {threads} > {log} 2>&1
         # --complete t \
         # --seq_names "Chromosome,unnamed" \
         # --seq_types "chromosome,plasmid" \
         # --seq_topologies "circular,circular" \
         """


# TODO - add option to control whether --dbmem flag is set (uses more RAM but does faster analysis)
rule run_eggnog:
    input:
        protein="annotation/dfast/protein.faa",
        install_finished=os.path.join(config.get("db_dir"),"checkpoints","eggnog_" + VERSION_EGGNOG)
    output:
        "annotation/eggnog/eggnog.emapper.annotations"
    conda:
        "../envs/eggnog.yaml"
    log:
        "logs/eggnog.log"
    benchmark:
        "benchmarks/eggnog.txt"
    params:
        outdir = "annotation/eggnog",
        tmpdir="annotation/eggnog/tmp",
        db=directory(os.path.join(config.get("db_dir"), "eggnog_" + VERSION_EGGNOG)),
        sensmode=config.get("eggnog_sensmode")
    threads:
        config.get("threads",1)
    shell:
        """
        mkdir -p {params.tmpdir}
        emapper.py --cpu {threads} -i {input.protein} --itype proteins -m diamond --sensmode ultra-sensitive \
          --dbmem --output eggnog --output_dir {params.outdir} --temp_dir {params.tmpdir} \
          --data_dir {params.db} > {log} 2>&1
        rm -r {params.tmpdir}
        """


# TODO - if adding support for multiple genomes, this could be run once with all genomes together to save time
rule run_gtdbtk:
    input:
        genome="annotation/dfast/genome.fna",
        install_finished=os.path.join(config.get("db_dir"),"checkpoints","GTDB_" + VERSION_GTDB)
    output:
        batchfile=temp("annotation/gtdbtk/batchfile.tsv"),
        annotation="annotation/gtdbtk/gtdbtk.summary.tsv"
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/gtdbtk.log"
    benchmark:
        "benchmarks/gtdbtk.txt"
    params:
        outdir="annotation/gtdbtk/run_files",
        db=directory(os.path.join(config.get("db_dir"), "GTDB_" + VERSION_GTDB)),
        genome_id=config.get("sample_id")
    threads:
        config.get("threads",1)
    shell:
        """
        GTDBTK_DATA_PATH={params.db}
        printf "{input.genome}\t{params.genome_id}\n" > {output.batchfile}
        gtdbtk classify_wf --batchfile {output.batchfile} --out_dir {params.outdir} \
          --cpus {threads} --pplacer_cpus {threads} > {log} 2>&1
        head -n 1 {params.outdir}/gtdbtk.*.summary.tsv | sort -u > {output.annotation}
        tail -n +2 {params.outdir}/gtdbtk.*.summary.tsv >> {output.annotation}
        """


if config.get("qc_short_r1") != "None":

    # TODO - clarify name compared to previous mapping step
    rule calculate_final_short_read_coverage:
        input:
            "annotation/dfast/genome.fna"
        output:
            mapping="annotation/coverage/short_read.bam",
            index="annotation/coverage/short_read.bam.bai",
            coverage="annotation/coverage/short_read_coverage.tsv"
        conda:
            "../envs/mapping.yaml"
        log:
            "logs/calculate_final_short_read_coverage.log"
        benchmark:
            "benchmarks/calculate_final_short_read_coverage.txt"
        params:
            qc_short_r1=config.get("qc_short_r1"),
            qc_short_r2=config.get("qc_short_r2")
        threads:
            config.get("threads",1)
        resources:
            mem=int(config.get("memory") / config.get("threads",1))
        shell:
            """
            bwa index {input} 2> {log}
            bwa mem -t {threads} {input} {params.qc_short_r1} {params.qc_short_r2} 2>> {log} | \
              samtools view -b -@ {threads} 2>> {log} | \
              samtools sort -@ {threads} -m {resources.mem}G 2>> {log} \
              > {output.mapping}
            samtools index -@ {threads} {output.mapping}
            samtools coverage {output.mapping} > {output.coverage}
            """


rule calculate_final_long_read_coverage:
    input:
        contigs="annotation/dfast/genome.fna",
        qc_long_reads="qc_long/nanopore_qc.fastq.gz"
    output:
        mapping="annotation/coverage/long_read.bam",
        index="annotation/coverage/long_read.bam.bai",
        coverage="annotation/coverage/long_read_coverage.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/calculate_final_long_read_coverage.log"
    benchmark:
        "benchmarks/calculate_final_long_read_coverage.txt"
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
        "annotation/coverage/long_read_coverage.tsv"
    output:
        logs=temp(directory("annotation/logs")),
        stats=temp(directory("annotation/stats"))
    params:
        logs="logs",
        stats="stats"
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
        "annotation/dfast/genome.fna",
        "annotation/eggnog/eggnog.emapper.annotations",
        "annotation/gtdbtk/gtdbtk.summary.tsv",
        expand("annotation/coverage/{type}_coverage.tsv",
            type=["short_read", "long_read"] if config.get("qc_short_r1") != "None" else ["long_read"]),
        "annotation/logs",
        "annotation/stats"
    output:
        "summary.zip"
    log:
        "logs/summarize_annotation.log"
    params:
        zipdir="annotation"
    shell:
        """
        cd {params.zipdir}
        zip -r ../{output} * -x \*.bam\* gtdbtk/run_files/\* > "../summarize_annotation.log" 2>&1
        cd ..
        mv "summarize_annotation.log" {log}
        """


rule annotation:
    input:
        "summary.zip"
    output:
        temp(touch("checkpoints/annotation"))


# TODO - add nice summaries of reads removed during QC, polishing stats, assembly stats, circular vs. non.
#        These aren't in the final ZIP summary at the moment.
