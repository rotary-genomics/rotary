# circularize: rules for circularizing genomes.
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os
import pandas as pd

START_HMM_NAME = os.path.splitext(os.path.basename(config.get("hmm_url")))[0]

DB_DIR_PATH = config.get('db_dir')

# SAMPLE_NAMES and POLISH_WITH_SHORT_READS are instantiated in rotary.smk

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

# Writes circular.list with the names of circular contigs if there are any circular contigs
# Writes linear.list with the names of linear contigs if there are any linear contigs
# Then, the DAG is re-evaluated. Circularization is only run if there are circular contigs.
# Based on clustering tutorial at https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html (accessed 2022.3.31)
checkpoint split_circular_and_linear_contigs:
    input:
        "{sample}/polish/{sample}_contig_info.tsv"
    output:
        directory("{sample}/circularize/filter/lists")
    run:
        contig_info = pd.read_csv(input[0], sep='\t')
        contig_info_filtered = contig_info[contig_info['pass_coverage_filter'] == 'Y']

        circular_contigs = contig_info_filtered[contig_info_filtered['circular'] == 'Y']
        linear_contigs = contig_info_filtered[contig_info_filtered['circular'] == 'N']

        os.makedirs(output[0],exist_ok=True)

        # Only output files if there is >=1 entry
        if circular_contigs.shape[0] >= 1:
            circular_contigs['contig'].to_csv(os.path.join(output[0],'circular.list'),header=None,index=False)

        if linear_contigs.shape[0] >= 1:
            linear_contigs['contig'].to_csv(os.path.join(output[0],'linear.list'),header=None,index=False)


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
        with open(input[0],'r') as hmmsearch_results_raw:
            line_count = len(hmmsearch_results_raw.readlines())

        if line_count == 0:
            # Write empty output if there were no HMM hits
            start_orf_ids = []

        else:
            # Load HMM search results
            hmmsearch_results = pd.read_csv(input[0],sep='\s+',header=None)[[0, 2, 3, 4]]

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

        pd.Series(start_orf_ids).to_csv(output[0],header=None,index=False)


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
        circlator_dir=directory('{sample}/circularize/circlator'),
        contigs_with_ends=temp('{sample}/circularize/circlator/rotated.promer.contigs_with_ends.fa')
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
        temp("{sample}/circularize/polypolish/input/{sample}_input.fasta")
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
        circular_or_linear=glob_wildcards(os.path.join(circularize_lists_path,"{i}.list")).i)


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
