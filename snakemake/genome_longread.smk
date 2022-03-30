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
        "checkpoints/qc_long",
        "checkpoints/assembly_long",
        "checkpoints/short_read_polish",
        "checkpoints/circularize",
        "checkpoints/annotation"


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
        flye_mode=config.get("flye_mode")
    threads:
        config.get("threads",1)
    shell:
        """
        flye --{params.flye_mode} {input} --out-dir {params.output_dir} -t {threads} > {log} 2>&1
        """


rule polish_medaka:
    input:
        qc_long_reads="qc_long/nanopore_qc.fastq.gz",
        contigs="assembly/flye/assembly.fasta"
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


rule assembly_long:
    input:
        "polish/medaka/consensus.fasta"
    output:
        temp(touch("checkpoints/assembly_long"))


# TODO: consider moving polypolish download to separate rule at start of pipeline so internet is not needed in middle
rule polish_polypolish:
    input:
        "polish/medaka/consensus.fasta"
    output:
        polypolish_filter="polish/polypolish/polypolish_insert_filter.py",
        polypolish="polish/polypolish/polypolish",
        mapping_r1=temp("polish/polypolish/R1.sam"),
        mapping_r2=temp("polish/polypolish/R2.sam"),
        mapping_clean_r1=temp("polish/polypolish/R1.clean.sam"),
        mapping_clean_r2=temp("polish/polypolish/R2.clean.sam"),
        polished="polish/polypolish/polypolish.fasta",
        debug="polish/polypolish/polypolish.debug.log",
        debug_stats="stats/polypolish_changes_round1.log"
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/polypolish.log"
    benchmark:
        "benchmarks/polypolish.txt"
    params:
        polypolish_url="https://github.com/rrwick/Polypolish/releases/download/v0.5.0/polypolish-linux-x86_64-musl-v0.5.0.tar.gz",
        output_dir="polish/polypolish",
        qc_short_r1=config.get("qc_short_r1"),
        qc_short_r2=config.get("qc_short_r2")
    threads:
        config.get("threads",1)
    shell:
        """
        printf "### File download ###\n" > {log}
        wget -O - {params.polypolish_url} 2>> {log} | tar -C {params.output_dir} -xvzf - >> {log} 2>&1
        
        printf "\n\n### Read mapping ###\n" >> {log}
        bwa index {input} 2>> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r1} > {output.mapping_r1} 2>> {log}
        bwa mem -t {threads} -a {input} {params.qc_short_r2} > {output.mapping_r2} 2>> {log}
        
        printf "\n\n### Polypolish insert filter ###\n" >> {log}
        {output.polypolish_filter} --in1 {output.mapping_r1} --in2 {output.mapping_r2} \
          --out1 {output.mapping_clean_r1} --out2 {output.mapping_clean_r2} 2>> {log}
          
        printf "\n\n### Polypolish ###\n" >> {log}
        {output.polypolish} {input} --debug {output.debug} \
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
        "polish/polca/polca.fasta"
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


# TODO - consider mapping to medaka polished contigs instead
rule calculate_short_read_coverage:
    input:
        "polish/polca/polca.fasta"
    output:
        mapping=temp("polish/cov_filter/short_read.bam"),
        mapping_index=temp("polish/cov_filter/short_read.bam.bai"),
        coverage="polish/cov_filter/coverage.tsv"
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


rule summarize_contigs_by_coverage:
    input:
        "polish/cov_filter/coverage.tsv"
    output:
        "polish/cov_filter/filtered_contigs.list"
    params:
        meandepth_cutoff=config.get("meandepth_cutoff"),
        evenness_cutoff=config.get("evenness_cutoff")
    run:
        coverage_data = pd.read_csv(input[0], sep='\t')
        coverage_filtered = coverage_data[\
            (coverage_data['meandepth'] >= params.meandepth_cutoff) & \
            (coverage_data['coverage'] >= params.evenness_cutoff)]
        coverage_filtered['#rname'].to_csv(output[0], header=None, index=False)


rule filter_contigs_by_coverage:
    input:
        contigs="polish/polca/polca.fasta",
        filter_list="polish/cov_filter/filtered_contigs.list"
    output:
        "polish/cov_filter/filtered_contigs.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        seqtk subseq -l 60 {input.contigs} {input.filter_list} > {output}
        """


rule symlink_short_read_polish:
    input:
        "polish/cov_filter/filtered_contigs.fasta"
    output:
        "polish/polish.fasta"
    run:
        source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
        os.symlink(source_relpath,str(output))


rule short_read_polish:
    input:
        "polish/polish.fasta"
    output:
        temp(touch("checkpoints/short_read_polish"))


rule find_circular_contigs:
    input:
        assembly_stats="assembly/flye/assembly_info.txt",
        filter_list="polish/cov_filter/filtered_contigs.list"
    output:
        circular_list="circularize/filter/circular.list",
        linear_list="circularize/filter/linear.list"
    run:
        coverage_filtered_contigs = pd.read_csv(input.filter_list, header=None)[0]

        assembly_info = pd.read_csv(input.assembly_stats, sep='\t')
        assembly_info_filtered = assembly_info[assembly_info['#seq_name'].isin(coverage_filtered_contigs)]

        circular_contigs = assembly_info_filtered[assembly_info_filtered['circ.'] == 'Y']
        linear_contigs = assembly_info_filtered[assembly_info_filtered['circ.'] == 'N']

        circular_contigs['#seq_name'].to_csv(output.circular_list, header=None, index=False)
        linear_contigs['#seq_name'].to_csv(output.linear_list, header=None, index=False)


rule filter_circular_contigs:
    input:
        contigs="polish/polish.fasta",
        circular_list="circularize/filter/circular.list",
        linear_list="circularize/filter/linear.list"
    output:
        circular="circularize/filter/circular.fasta",
        linear="circularize/filter/linear.fasta"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        # Note: command still works (and outputs an empty file) if the input list is empty
        seqtk subseq -l 60 {input.contigs} {input.circular_list} > {output.circular}
        seqtk subseq -l 60 {input.contigs} {input.linear_list} > {output.linear}
        """


# Get number of circular contigs
circular_contig_count = 0
for line in open("circularize/filter/circular.list"):
    circular_contig_count = circular_contig_count + 1

if circular_contig_count > 0:

    rule hmm_download:
        output:
            "circularize/identify/start.hmm"
        log:
            "logs/hmm_download.log"
        benchmark:
            "benchmarks/hmm_download.txt"
        params:
            pfam_id=config.get("start_hmm_pfam_id")
        shell:
            """
            wget -O {output} --no-check-certificate \
              https://pfam.xfam.org/family/{params.start_hmm_pfam_id}/hmm 2> {log}
            """


    rule search_genome_start:
        input:
            contigs="circularize/filter/circular.fasta",
            hmm="circularize/identify/start.hmm"
        output:
            orf_predictions=temp("circularize/identify/circular.faa"),
            gene_predictions=temp("circularize/identify/circular.ffn"),
            annotation_gff=temp("circularize/identify/circular.gff"),
            search_hits="circularize/identify/hmmsearch_hits.txt"
        conda:
            "../envs/mapping.yaml"
        log:
            "logs/search_genome_start.log"
        benchmark:
            "benchmarks/search_genome_start.txt"
        params:
            hmmsearch_evalue=config.get("hmmsearch_evalue")
        threads:
            config.get("threads",1)
        shell:
            """
            printf "\n\n### Predict genes ###\n" > {log}
            prodigal -i {input.contigs} -a {output.orf_predictions} -d {output.gene_predictions} \
              -f gff -o {output.annotation_gff} 2>> {log}
            
            printf "\n\n### Find HMM hits ###\n" >> {log}
            hmmsearch --cpu {threads} -E {params.hmmsearch_evalue} --tblout /dev/stdout -o /dev/stderr \
              {input.hmm} {output.orf_predictions} > {output.search_hits} 2>> {log}
            
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


    rule filter_start_genes:
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
            min_id=config.get("circlator_min_id"),
            run_name="circularize/circlator/rotated"
        shell:
            """
            circlator fixstart --min_id {params.min_id} --genes_fa {input.start_gene}\
              {input.contigs} {params.run_name} > {log} 2>&1
            
            printf "\n\n### Circlator log output ###\n" >> {log}
            cat "circularize/circlator/rotated.log" >> {log}
            """


    # TODO: can I merge this code with the first run of polypolish? The only difference is to change the input/output names and to add the summarizer code at the end
    rule repolish_polypolish:
        input:
            "circularize/circlator/rotated.fasta"
        output:
            polypolish_filter="circularize/polypolish/polypolish_insert_filter.py",
            polypolish="circularize/polypolish/polypolish",
            mapping_r1=temp("circularize/polypolish/R1.sam"),
            mapping_r2=temp("circularize/polypolish/R2.sam"),
            mapping_clean_r1=temp("circularize/polypolish/R1.clean.sam"),
            mapping_clean_r2=temp("circularize/polypolish/R2.clean.sam"),
            polished="circularize/polypolish/polypolish.fasta",
            debug="circularize/polypolish/polypolish.debug.log",
            debug_stats="stats/polypolish_changes_round2.log"
        conda:
            "../envs/mapping.yaml"
        log:
            "logs/repolypolish.log"
        benchmark:
            "benchmarks/repolypolish.txt"
        params:
            polypolish_url="https://github.com/rrwick/Polypolish/releases/download/v0.5.0/polypolish-linux-x86_64-musl-v0.5.0.tar.gz",
            output_dir="circularize/polypolish",
            qc_short_r1=config.get("qc_short_r1"),
            qc_short_r2=config.get("qc_short_r2")
        threads:
            config.get("threads",1)
        shell:
            """
            printf "### File download ###\n" > {log}
            wget -O - {params.polypolish_url} 2>> {log} | tar -C {params.output_dir} -xvzf - >> {log} 2>&1
            
            printf "\n\n### Read mapping ###\n" >> {log}
            bwa index {input} 2>> {log}
            bwa mem -t {threads} -a {input} {params.qc_short_r1} > {output.mapping_r1} 2>> {log}
            bwa mem -t {threads} -a {input} {params.qc_short_r2} > {output.mapping_r2} 2>> {log}
            
            printf "\n\n### Polypolish insert filter ###\n" >> {log}
            {output.polypolish_filter} --in1 {output.mapping_r1} --in2 {output.mapping_r2} \
              --out1 {output.mapping_clean_r1} --out2 {output.mapping_clean_r2} 2>> {log}
              
            printf "\n\n### Polypolish ###\n" >> {log}
            {output.polypolish} {input} --debug {output.debug} \
              {output.mapping_clean_r1} {output.mapping_clean_r2} 2>> {log} | 
              seqtk seq -A -l 0 | 
              awk \'{{ if ($0 ~ /^>/) {{ gsub("_polypolish", ""); print }} else {{ print }} }}\' | 
              seqtk seq -l 60 > {output.polished} 2>> {log}
            
            head -n 1 {output.debug} > {output.debug_stats}
            grep changed {output.debug} >> {output.debug_stats}
            
            printf "\n\n### Done. ###\n"
            """


    # TODO - consider sorting contigs by length (or by circular and then by length)
    rule combine_circular_and_linear_contigs:
        input:
            circular_rotated="circularize/polypolish/polypolish.fasta",
            linear="circularize/filter/linear.fasta"
        output:
            "circularize/combine/combined.fasta"
        conda:
            "../envs/mapping.yaml"
        shell:
            """
            cat {input.circular_rotated} {input.linear} | seqtk seq -l 60 > {output}
            """


    rule symlink_circularization:
        input:
            "circularize/combine/combined.fasta"
        output:
            "circularize/circularize.fasta"
        run:
            source_relpath = os.path.relpath(str(input),os.path.dirname(str(output)))
            os.symlink(source_relpath,str(output))

else:

    rule bypass_circularization:
        input:
            "circularize/filter/linear.fasta"
        output:
            "circularize/circularize.fasta"
        run:
            print('No circular contigs, so will bypass circularization.')
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
        "circularize/circularize.fasta"
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
        dfast_db=config.get("dfast_db"),
        strain=config.get("sample_id")
    threads:
        config.get("threads",1)
    shell:
        """
        dfast --force \
          --dbroot {params.dfast_db} \
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
        "annotation/dfast/protein.faa"
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
        "annotation/dfast/genome.fna"
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
        db_dir=config.get("gtdbtk_db"),
        genome_id=config.get("sample_id")
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
        "annotation/dfast/genome.fna"
    output:
        mapping="annotation/coverage/short_read.bam",
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


rule summarize_annotation:
    input:
        "annotation/dfast/genome.fna",
        "annotation/eggnog/eggnog.emapper.annotations",
        "annotation/gtdbtk/gtdbtk.summary.tsv",
        "annotation/coverage/short_read_coverage.tsv",
        "annotation/coverage/long_read_coverage.tsv"
    output:
        "summary.zip"
    params:
        zipdir="annotation"
    shell:
        """
        zip -r {output} {params.zipdir} -x \*.bam\* \*/gtdbtk/run_files/\*
        """

rule annotation:
    input:
        "summary.zip"
    output:
        temp(touch("checkpoints/annotation"))


# TODO - add nice summaries of reads removed during QC, polishing stats, assembly stats, circular vs. non.
#        These aren't in the final ZIP summary at the moment.
