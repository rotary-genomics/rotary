# Snakemake rules long read-based genome assembly workflow
# Copyright Jackson M. Tsuji, ILTS, Hokkaido University, 2022

import os
import pandas as pd
import itertools
from snakemake.utils import logger, min_version, update_config

# Specify the minimum snakemake version allowable
min_version("6.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")


rule all:
    input:
        "{genome}/checkpoints/qc_long",
        "{genome}/checkpoints/assembly_long",
        "{genome}/checkpoints/short_read_polish",
        "{genome}/checkpoints/circularize",
        "{genome}/checkpoints/annotation"


rule nanopore_qc_filter:
    input:
        lambda wildcards: config["genomes"][wildcards.genome]
    output:
        "{genome}/qc_long/{genome}.fastq.gz"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/qc_long.log"
    benchmark:
        "{genome}/benchmarks/qc_long.benchmark.txt"
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
        "{genome}/qc_long/{genome}.fastq.gz"
    output:
        "{genome}/qc_long/length_hist.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/qc_long_length_hist.log"
    benchmark:
        "{genome}/benchmarks/qc_long_length_hist.benchmark.txt"
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        reformat.sh in={input} lhist={output} maxhistlen=10000000 
          interleaved=f qin=33 threads={threads} -Xmx{resources.mem}G > {log} 2>&1
        """


rule qc_long_length_stats:
    input:
        "{genome}/qc_long/length_hist.tsv"
    output:
        "{genome}/stats/qc_long_length_stats.txt"
    log:
        "{genome}/logs/qc_long_length_stats.log"
    benchmark:
        "{genome}/benchmarks/qc_long_length_stats.benchmark.txt"
    run:
        length_hist = pd.read_csv(input, sep='\t')

        lengths = []

        for index, row in length_hist.iterrows():
            length, count = row
            lengths = lengths + list(itertools.repeat(length,count))

        lengths = pd.Series(lengths)

        length_stats = pd.DataFrame({'Total reads':   lengths.shape[0],
                                     'Mean length':   round(lengths.mean(),2),
                                     'Median length': lengths.median(),
                                     'Min length':    lengths.min(),
                                     'Max length':    lengths.max()})\
          .transpose()

        length_stats.to_csv(output, index=True)


rule qc_long:
    input:
        "{genome}/stats/qc_long_length_stats.txt"
    output:
        temp(touch("{genome}/checkpoints/qc_long"))


rule assembly_flye:
    input:
        "{genome}/qc_long/{genome}.fastq.gz"
    output:
        "{genome}/assembly/flye/assembly.fasta",
        "{genome}/assembly/flye/assembly_info.txt"
    conda:
        "../envs/assembly_flye.yaml"
    log:
        "{genome}/logs/assembly_flye.log"
    benchmark:
        "{genome}/benchmarks/assembly_flye.benchmark.txt"
    params:
        output_dir="{genome}/assembly/flye",
        flye_mode=config.get("flye_mode")
    threads:
        config.get("threads",1)
    shell:
        """
        flye --{params.flye_mode} {input} --out-dir {params.output_dir} -t {threads} > {log} 2>&1
        """


rule polish_medaka:
    input:
        qc_long_reads="{genome}/qc_long/{genome}.fastq.gz",
        contigs="{genome}/assembly/flye/assembly.fasta"
    output:
        "{genome}/polish/medaka/consensus.fasta"
    conda:
        "../envs/medaka.yaml"
    log:
        "{genome}/logs/medaka.log"
    benchmark:
        "{genome}/benchmarks/medaka.txt"
    params:
        output_dir="{genome}/polish/medaka",
        medaka_model=config.get("medaka_model")
    threads:
        config.get("threads",1)
    shell:
        """
        medaka_consensus -i {input.qc_long_reads} -d {input.contigs} -o {params.output_dir} \
          -m {params.medaka_model} -t {threads} > {log} 2>&1
        """


rule assembly_long:
    input:
        "{genome}/polish/medaka/consensus.fasta"
    output:
        temp(touch("{genome}/checkpoints/assembly_long"))


# TODO: consider splitting into multiple functions; make sam files temporary;
#   consider moving polypolish download to separate rule at start of pipeline so internet is not needed in middle
rule polish_polypolish:
    input:
        "{genome}/polish/medaka/consensus.fasta"
    output:
        polypolish_filter="{genome}/polish/polypolish/polypolish_insert_filter.py",
        polypolish="{genome}/polish/polypolish/polypolish",
        mapping_r1="{genome}/polish/polypolish/R1.sam",
        mapping_r2="{genome}/polish/polypolish/R2.sam",
        mapping_clean_r1="{genome}/polish/polypolish/R1.clean.sam",
        mapping_clean_r2="{genome}/polish/polypolish/R2.clean.sam",
        polished="{genome}/polish/polypolish/polypolish.fasta",
        debug="{genome}/polish/polypolish/polypolish.debug.log"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/polypolish.log"
    benchmark:
        "{genome}/benchmarks/polypolish.txt"
    params:
        polypolish_url="https://github.com/rrwick/Polypolish/releases/download/v0.5.0/polypolish-linux-x86_64-musl-v0.5.0.tar.gz",
        output_dir="{genome}/polish/polypolish",
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2")
    threads:
        config.get("threads",1)
    shell:
        """
        wget -O - {params.polypolish_url} | tar -xvzf -C {params.output_dir} -
        bwa index {input} 2> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r1} > {output.mapping_r1} 2>> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r2} > {output.mapping_r2} 2>> {log}
        {output.polypolish_filter} --in1 {output.mapping_r1} --in2 {output.mapping_r2} \
          --out1 {output.mapping_clean_r1} --out2 {output.mapping_clean_r2} 2>> {log} 
        {output.polypolish} {input} --debug {output.debug} \
          {output.mapping_clean_r1} {output.mapping_clean_r2} 2>> {log} | 
          seqtk seq -A -l 0 | 
          awk "{{ if ($0 ~ /^>/) {{ gsub("_polypolish", ""); print }} else {{ print }} }}" | 
          seqtk seq -l 60 > {output.polished} 2>> {log}
        """


rule polish_polca:
    input:
        "{genome}/polish/polypolish/polypolish.fasta"
    output:
        raw="{genome}/polish/polca/polypolish.fasta.PolcaCorrected.fa",
        renamed="{genome}/polish/polca/polca.fasta"
    conda:
        "../envs/masurca.yaml"
    log:
        "{genome}/logs/polca.log"
    benchmark:
        "{genome}/benchmarks/polca.txt"
    params:
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2")
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        polca.sh -a {input} -r {params.qc_short_r1} {params.qc_short_r2} -t {threads} -m {resources.mem}G > {log} 2>&1
        ln -s {output.raw} {output.renamed}
        """


# TODO - consider mapping to medaka polished contigs instead
rule calculate_short_read_coverage:
    input:
        "{genome}/polish/polca/polca.fasta"
    output:
        mapping=temp("{genome}/polish/cov_filter/short_read.bam"),
        coverage="{genome}/polish/cov_filter/coverage.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/calculate_short_read_coverage.log"
    benchmark:
        "{genome}/benchmarks/calculate_short_read_coverage.txt"
    params:
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2")
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
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


rule summarize_contigs_by_coverage:
    input:
        "{genome}/polish/cov_filter/coverage.tsv"
    output:
        "{genome}/polish/cov_filter/filtered_contigs.list"
    params:
        meandepth_cutoff=config.get("meandepth_cutoff"),
        evenness_cutoff=config.get("evenness_cutoff")
    run:
        coverage_data = pd.read_csv({input}, sep='\t')
        coverage_filtered = coverage_data[
            (coverage_data['meandepth'] >= {params.meandepth_cutoff}) &
            (coverage_data['coverage'] >= {params.evenness_cutoff})]
        coverage_filtered['#rname'].to_csv({output}, header=None, index=False)


# TODO - consider writing log and benchmark files
rule filter_contigs_by_coverage:
    input:
        contigs="{genome}/polish/polca/polca.fasta",
        filter_list="{genome}/polish/cov_filter/filtered_contigs.list"
    output:
        "{genome}/polish/cov_filter/filtered_contigs.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        seqtk subseq -l 60 {input.contigs} {input.filter_list} > {output}
        """


rule symlink_short_read_polish:
    input:
        "{genome}/polish/cov_filter/filtered_contigs.fasta"
    output:
        "{genome}/polish/polish.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule short_read_polish:
    input:
        "{genome}/polish/polish.fasta"
    output:
        temp(touch("{genome}/checkpoints/short_read_polish"))


rule find_circular_contigs:
    input:
        assembly_stats="{genome}/assembly/flye/assembly_info.txt",
        filter_list="{genome}/polish/cov_filter/filtered_contigs.list"
    output:
        circular_list="{genome}/circularize/filter/circular.list",
        linear_list="{genome}/circularize/filter/linear.list"
    run:
        coverage_filtered_contigs = pd.read_csv({input.filter_list}, header=None)[0]

        assembly_info = pd.read_csv({input.assembly_stats}, sep='\t')
        assembly_info_filtered = assembly_info[assembly_info['#seq_name'].isin(coverage_filtered_contigs)]

        circular_contigs = assembly_info_filtered[assembly_info_filtered['circ.'] == 'Y']
        linear_contigs = assembly_info_filtered[assembly_info_filtered['circ.'] == 'N']

        circular_contigs['#seq_name'].to_csv({output.circular_list}, header=None, index=False)
        linear_contigs['#seq_name'].to_csv({output.linear_list}, header=None, index=False)


# TODO - consider writing log and benchmark files
rule filter_circular_contigs:
    input:
        contigs="{genome}/polish/polish.fasta",
        circular_list="{genome}/circularize/filter/circular.list",
        linear_list="{genome}/circularize/filter/linear.list"
    output:
        circular="{genome}/circularize/filter/circular.fasta",
        linear="{genome}/circularize/filter/linear.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        # Note: command still works (and outputs an empty file) if the input list is empty
        seqtk subseq -l 60 {input.contigs} {input.circular_list} > {output.circular}
        seqtk subseq -l 60 {input.contigs} {input.linear_list} > {output.linear}
        """


# TODO - only run the following commands if circular contigs actually exist
# TODO - move HMM download to its own rule
# TODO - check the number of hmmsearch hits. Right now, all hits in the evalue threshold are taken
#        (but could take the top with head -n 1, for example). Also flag if there are no hits.
#        Indeed, if multiple circular contigs from different organisms, might be multiple start genes.
#        Not sure how fixstart would handle two closely matching genes in the same circular contig
#        **Thus, the ideal solution would be to take the top hit (and throw a warning) for each contig if there are multiple hits.
#        Could get this using some python code and the GFF file from prodigal
rule find_genome_start:
    input:
        "{genome}/circularize/filter/circular.fasta"
    output:
        hmm="{genome}/circularize/identify/start.hmm",
        orf_predictions="{genome}/circularize/identify/circular.faa",
        gene_predictions="{genome}/circularize/identify/circular.ffn",
        search_hits="{genome}/circularize/identify/hmmsearch_hits.txt",
        start_gene_list="{genome}/circularize/identify/start_gene.list",
        start_gene="{genome}/circularize/identify/start_gene.ffn"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/find_genome_start.log"
    benchmark:
        "{genome}/benchmarks/find_genome_start.txt"
    params:
        start_hmm_pfam_id=config.get("start_hmm_pfam_id"),
        hmmsearch_evalue=config.get("hmmsearch_evalue")
    threads:
        config.get("threads",1)
    shell:
        """
        wget -O {output.hmm} --no-check-certificate \
          https://pfam.xfam.org/family/{params.start_hmm_pfam_id}/hmm 2> {log}
        prodigal -i {input} -a {output.orf_predictions} -d {output.gene_predictions} >/dev/null 2>> {log}
        hmmsearch --cpu {threads} -E {params.hmmsearch_evalue} --tblout /dev/stdout -o /dev/null \
          {output.hmm} {output.orf_predictions} > {output.search_hits} 2>> {log}
        cat {output.search_hits} | grep -v "^#" | cut -d " " -f 1 > {output.start_gene_list}
        seqtk subseq {output.gene_predictions} {output.start_gene_list} > {output.start_gene}.ffn 2>> {log}
        """


rule run_circlator:
    input:
        contigs="{genome}/circularize/filter/circular.fasta",
        start_gene="{genome}/circularize/identify/start_gene.ffn"
    output:
        "{genome}/circularize/circlator/rotated/rotated.fasta"
    conda:
        "../envs/circlator.yaml"
    log:
        "{genome}/logs/circlator.log"
    benchmark:
        "{genome}/benchmarks/circlator.txt"
    params:
        min_id=config.get("circlator_min_id"),
        run_name="{genome}/circularize/circlator/rotated"
    shell:
        """
        circlator fixstart --min_id {params.min_id} --genes_fa {input.start_gene}\
          {input.contigs} {params.run_name} > {log} 2>&1
        """


# TODO: can I merge this code with the first run of polypolish? The only difference is to change the input/output names and to add the summarizer code at the end
rule repolish_polypolish:
    input:
        "{genome}/circularize/circlator/rotated/rotated.fasta"
    output:
        polypolish_filter="{genome}/circularize/polypolish/polypolish_insert_filter.py",
        polypolish="{genome}/circularize/polypolish/polypolish",
        mapping_r1="{genome}/circularize/polypolish/R1.sam",
        mapping_r2="{genome}/circularize/polypolish/R2.sam",
        mapping_clean_r1="{genome}/circularize/polypolish/R1.clean.sam",
        mapping_clean_r2="{genome}/circularize/polypolish/R2.clean.sam",
        polished="{genome}/circularize/polypolish/polypolish.fasta",
        debug="{genome}/circularize/polypolish/polypolish.debug.log"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/repolypolish.log"
    benchmark:
        "{genome}/benchmarks/repolypolish.txt"
    params:
        polypolish_url="https://github.com/rrwick/Polypolish/releases/download/v0.5.0/polypolish-linux-x86_64-musl-v0.5.0.tar.gz",
        output_dir="{genome}/circularize/polypolish",
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2")
    threads:
        config.get("threads",1)
    shell:
        """
        wget -O - {params.polypolish_url} | tar -xvzf -C {params.output_dir} -
        bwa index {input} 2> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r1} > {output.mapping_r1} 2>> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r2} > {output.mapping_r2} 2>> {log}
        {output.polypolish_filter} --in1 {output.mapping_r1} --in2 {output.mapping_r2} \
          --out1 {output.mapping_clean_r1} --out2 {output.mapping_clean_r2} 2>> {log} 
        {output.polypolish} {input} --debug {output.debug} \
          {output.mapping_clean_r1} {output.mapping_clean_r2} 2>> {log} | 
          seqtk seq -A -l 0 | 
          awk "{{ if ($0 ~ /^>/) {{ gsub("_polypolish", ""); print }} else {{ print }} }}" | 
          seqtk seq -l 60 > {output.polished} 2>> {log}
          
        printf \n\nSummary of polishing changes:\n" >> {log}
        (head -n 1 polypolish.debug.log; grep changed polypolish.debug.log; ) | column -t >> {log}
        """


# TODO - consider adding log and benchmark
# TODO - consider sorting contigs by length (or by circular and then by length)
rule combine_circular_and_linear_contigs:
    input:
        circular_rotated="{genome}/circularize/polypolish/polypolish.fasta",
        linear="{genome}/circularize/filter/linear.fasta"
    output:
        "{genome}/circularize/combine/combined.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        cat {input.circular_rotated} {input.linear} | seqtk seq -l 60 > {output}
        """


rule symlink_circularization:
    input:
        "{genome}/circularize/combine/combined.fasta"
    output:
        "{genome}/circularize/{genome}.fna"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule circularize:
    input:
        "{genome}/circularize/{genome}.fna"
    output:
        temp(touch("{genome}/checkpoints/circularize"))


# TODO - I might need to remove special characters from strain name to use for locus tag prefix
# TODO - can I auto-predict genome completeness, names, types, topologies?
rule run_dfast:
    input:
        "{genome}/circularize/{genome}.fna"
    output:
        "{genome}/annotation/dfast/genome.fna",
        "{genome}/annotation/dfast/protein.faa"
    conda:
        "../envs/annotation_dfast.yaml"
    log:
        "{genome}/logs/annotation_dfast.log"
    benchmark:
        "{genome}/benchmarks/annotation_dfast.txt"
    params:
        outdir="{genome}/annotation/dfast",
        dfast_db=config.get("dfast_db"),
        strain="{genome}"
    threads:
        config.get("threads",1)
    shell:
        """
        dfast --dbroot {params.dfast_db} \
          -g {input} \
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
        "{genome}/annotation/dfast/protein.faa"
    output:
        "{genome}/annotation/eggnog/eggnog.emapper.annotations"
    conda:
        "../envs/eggnog.yaml"
    log:
        "{genome}/logs/eggnog.log"
    benchmark:
        "{genome}/benchmarks/eggnog.txt"
    params:
        outdir = "{genome}/annotation/eggnog",
        tmpdir="{genome}/annotation/eggnog/tmp",
        db_dir=config.get("eggnog_db"),
        sensmode=config.get("eggnog_sensmode")
    threads:
        config.get("threads",1)
    shell:
        """
        mkdir -p {params.tmpdir}
        emapper.py --cpu {threads} -i {input} --itype proteins -m diamond --sensmode ultra-sensitive \
          --dbmem --output eggnog --output_dir {params.outdir} --temp_dir {params.tmpdir} \
          --data_dir {params.db_dir} > {log} 2>&1
        rm -r {params.tmpdir}
        """


# TODO - if adding support for multiple genomes, this could be run once with all genomes together to save time
rule run_gtdbtk:
    input:
        "{genome}/annotation/dfast/genome.fna"
    output:
        batchfile="{genome}/annotation/gtdbtk/batchfile.tsv",
        annotation="{genome}/annotation/gtdbtk/gtdbtk.summary.tsv"
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "{genome}/logs/gtdbtk.log"
    benchmark:
        "{genome}/benchmarks/gtdbtk.txt"
    params:
        outdir="{genome}/annotation/gtdbtk/run_files",
        db_dir=config.get("gtdbtk_db"),
        genome_id="{genome}"
    threads:
        config.get("threads",1)
    shell:
        """
        GTDBTK_DATA_PATH={params.db_dir}
        printf "{input}\t{params.genome_id}\n" > {output.batchfile}
        gtdbtk classify_wf --batchfile {output.batchfile} --out_dir {params.outdir} \
          --cpus {threads} --pplacer_cpus {threads} > {log} 2>&1
        head -n 1 {params.outdir}/gtdbtk.*.summary.tsv | sort -u > {output.annotation}
        tail -n +2 {params.outdir}/gtdbtk.*.summary.tsv >> {output.annotation}
        """


# TODO - clarify name compared to previous mapping step
rule calculate_final_short_read_coverage:
    input:
        "{genome}/annotation/dfast/genome.fna"
    output:
        mapping="{genome}/annotation/coverage/short_read.bam",
        coverage="{genome}/annotation/coverage/short_read_coverage.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/calculate_final_short_read_coverage.log"
    benchmark:
        "{genome}/benchmarks/calculate_final_short_read_coverage.txt"
    params:
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2")
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
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
        contigs="{genome}/annotation/dfast/genome.fna",
        qc_long_reads="{genome}/qc_long/{genome}.fastq.gz"
    output:
        mapping="{genome}/annotation/coverage/long_read.bam",
        coverage="{genome}/annotation/coverage/long_read_coverage.tsv"
    conda:
        "../envs/mapping.yaml"
    log:
        "{genome}/logs/calculate_final_long_read_coverage.log"
    benchmark:
        "{genome}/benchmarks/calculate_final_long_read_coverage.txt"
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.contigs} {input.qc_long_reads} 2> {log} | \
          samtools view -b -@ {threads} 2>> {log} | \
          samtools sort -@ {threads} -m {resources.mem}G 2>> {log} \
          > {output.mapping}
        samtools index -@ {threads} {output.mapping}
        samtools coverage {output.mapping} > {output.coverage}
        """


rule summarize_annotation:
    input:
        "{genome}/annotation"
    output:
        "{genome}/{genome}_summary.zip"
    params:
        bam_short="{genome}/annotation/coverage/short_read.bam",
        bam_long="{genome}/annotation/coverage/long_read.bam"
    shell:
        """
        zip -r -x {params.bam_short} -x {params.bam_long} {output} {input}
        """

rule annotation:
    input:
        "{genome}/annotation/dfast/genome.fna",
        "{genome}/annotation/eggnog/eggnog.emapper.annotations",
        "{genome}/annotation/gtdbtk/gtdbtk.summary.tsv",
        "{genome}/annotation/coverage/short_read_coverage.tsv",
        "{genome}/annotation/coverage/long_read_coverage.tsv",
        "{genome}/{genome}_summary.zip"
    output:
        temp(touch("{genome}/checkpoints/annotation"))


# TODO - add nice summaries of reads removed during QC, polishing stats, assembly stats, circular vs. non.
#        These aren't in the final ZIP summary at the moment.
