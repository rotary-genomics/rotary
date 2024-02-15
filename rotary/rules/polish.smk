# polish: rules for long read polishing.
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os
import sys


DB_DIR_PATH = config.get('db_dir')

READ_MAPPING_FILE_EXTENSIONS = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']

# SAMPLE_NAMES and POLISH_WITH_SHORT_READS are instantiated in rotary.smk

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
        contigs="{sample}/{step}/medaka/{sample}_consensus.fasta",
        calls_to_draft=temp(multiext("{sample}/{step}/medaka/calls_to_draft", '.bam', '.bam.bai')),
        consensus_probs=temp("{sample}/{step}/medaka/consensus_probs.hdf")
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

rule map_short_reads_for_polishing:
    input:
        qc_short_r1 = "{sample}/qc/{sample}_qc_R1.fastq.gz",
        qc_short_r2 = "{sample}/qc/{sample}_qc_R2.fastq.gz",
        contigs = "{sample}/{step}/polypolish/input/{sample}_input.fasta",
    output:
        mapping_r1 = temp("{sample}/{step}/polypolish/{sample}_R1.sam"),
        mapping_r2 = temp("{sample}/{step}/polypolish/{sample}_R2.sam"),
        read_mapping_files= temp(multiext("{sample}/{step}/polypolish/input/{sample}_input.fasta",
            *READ_MAPPING_FILE_EXTENSIONS))
    conda:
        "../envs/mapping.yaml"
    log:
        "{sample}/logs/{step}/bwa_mem.log"
    benchmark:
        "{sample}/benchmarks/{step}/bwa_mem.txt"
    threads:
        config.get("threads", 1)
    shell:
        """
        bwa-mem2 index {input.contigs} 2>> {log}
        bwa-mem2 mem -t {threads} -a {input.contigs} {input.qc_short_r1} > {output.mapping_r1} 2>> {log}
        bwa-mem2 mem -t {threads} -a {input.contigs} {input.qc_short_r2} > {output.mapping_r2} 2>> {log}
        """

rule polish_polypolish:
    input:
        contigs="{sample}/{step}/polypolish/input/{sample}_input.fasta",
        mapping_r1 = "{sample}/{step}/polypolish/{sample}_R1.sam",
        mapping_r2 = "{sample}/{step}/polypolish/{sample}_R2.sam"
    output:
        mapping_clean_r1 = temp("{sample}/{step}/polypolish/{sample}_R1.clean.sam"),
        mapping_clean_r2 = temp("{sample}/{step}/polypolish/{sample}_R2.clean.sam"),
        polished = temp("{sample}/{step}/polypolish/{sample}_polypolish.fasta"),
        debug = temp("{sample}/{step}/polypolish/polypolish.debug.log"),
        debug_stats = "{sample}/stats/{step}/polypolish_changes.log"
    conda:
        "../envs/polypolish.yaml"
    log:
        "{sample}/logs/{step}/polypolish.log"
    benchmark:
        "{sample}/benchmarks/{step}/polypolish.txt"
    params:
        careful = "--careful" if str(config.get("careful_short_read_polishing")).lower() == "true" else ""
    shell:
        """
        printf "### Polypolish insert filter ###\n" >> {log}
        polypolish filter --in1 {input.mapping_r1} --in2 {input.mapping_r2} \
          --out1 {output.mapping_clean_r1} --out2 {output.mapping_clean_r2} 2>> {log}
          
        printf "\n\n### Polypolish ###\n" >> {log}
        polypolish polish --debug {output.debug} {params.careful} {input.contigs}  \
          {output.mapping_clean_r1} {output.mapping_clean_r2} 2>> {log} |
          seqtk seq -A -C -l 60 > {output.polished} 2>> {log}
          
        head -n 1 {output.debug} > {output.debug_stats}
        grep changed {output.debug} >> {output.debug_stats}
        
        printf "\n\n### Done. ###\n"
        """


rule polish_pypolca:
    input:
        qc_short_r1 = "{sample}/qc/{sample}_qc_R1.fastq.gz",
        qc_short_r2 = "{sample}/qc/{sample}_qc_R2.fastq.gz",
        polished = "{sample}/polish/polypolish/{sample}_polypolish.fasta"
    output:
        polished = temp("{sample}/polish/pypolca/{sample}_corrected.fasta"),
        report = "{sample}/stats/polish/pypolca_changes.report",
        variant_calling = temp("{sample}/polish/pypolca/{sample}.vcf"),
        intermediate_logs = temp(directory("{sample}/polish/pypolca/logs"))
    conda:
        "../envs/pypolca.yaml"
    log:
        "{sample}/logs/polish/pypolca.log"
    benchmark:
        "{sample}/benchmarks/polish/pypolca.txt"
    params:
        outdir = "{sample}/polish/pypolca",
        report= "{sample}/polish/pypolca/{sample}.report",
        careful = "--careful" if str(config.get("careful_short_read_polishing")).lower() == "true" else "",
        mem_per_thread = int(config.get("memory") / config.get("threads",1))
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        pypolca run -a {input.polished} -1 {input.qc_short_r1} -2 {input.qc_short_r2} -o {params.outdir} \
            -p {wildcards.sample} {params.careful} -f -t {threads} -m {params.mem_per_thread}G > {log} 2>&1
        
        mv {params.report} {output.report}
        
        # Delete the logfile generated by pypolca that has a timestamp on the end - we already save another log file
        default_logfiles=($(find {params.outdir} -name "pypolca_*.log"))
        if [[ "${{#default_logfiles[@]}}" == 1 ]]; then
            rm "${{default_logfiles[0]}}"
        else
            echo "ERROR: more than one pypolca log found: ${{default_logfiles[@]}}"
            exit 1
        fi
        """


# Conditional based on whether short read polishing was performed
rule pre_coverage_filter:
    input:
        "{sample}/polish/medaka/{sample}_consensus.fasta" if POLISH_WITH_SHORT_READS == False else "{sample}/polish/pypolca/{sample}_corrected.fasta"
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
            coverage="{sample}/polish/cov_filter/{sample}_short_read_coverage.tsv",
            read_mapping_files= temp(multiext("{sample}/polish/cov_filter/{sample}_pre_filtered.fasta",
                *READ_MAPPING_FILE_EXTENSIONS))
        conda:
            "../envs/mapping.yaml"
        log:
            "{sample}/logs/polish/calculate_short_read_coverage.log"
        benchmark:
            "{sample}/benchmarks/polish/calculate_short_read_coverage.txt"
        params:
            mem_per_thread = int(config.get("memory") / config.get("threads",1))
        threads:
            config.get("threads",1)
        resources:
            mem=config.get("memory")
        shell:
            """
            # Note that -F 4 removes unmapped reads
            bwa-mem2 index {input.contigs} 2> {log}
            bwa-mem2 mem -t {threads} {input.contigs} {input.qc_short_r1} {input.qc_short_r2} 2>> {log} | \
              samtools view -b -F 4 -@ {threads} 2>> {log} | \
              samtools sort -@ {threads} -m {params.mem_per_thread}G 2>> {log} \
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
            coverage="{sample}/polish/cov_filter/{sample}_long_read_coverage.tsv",
            long_read_map=temp("{sample}/assembly/end_repair/{sample}_repaired.fasta.map-ont.mmi")
        conda:
            "../envs/mapping.yaml"
        log:
            "{sample}/logs/polish/calculate_long_read_coverage.log"
        benchmark:
            "{sample}/benchmarks/polish/calculate_long_read_coverage.txt"
        params:
            mem_per_thread = int(config.get("memory") / config.get("threads",1))
        threads:
            config.get("threads",1)
        resources:
            mem=config.get("memory")
        shell:
            """
            # Note that -F 4 removes unmapped reads
            minimap2 -t {threads} -ax map-ont {input.contigs} {input.qc_long_reads} 2> {log} | \
              samtools view -b -F 4 -@ {threads} 2>> {log} | \
              samtools sort -@ {threads} -m {params.mem_per_thread}G 2>> {log} \
              > {output.mapping}
            samtools index -@ {threads} {output.mapping}
            samtools coverage {output.mapping} > {output.coverage}
            """


if ((POLISH_WITH_SHORT_READS == False) |
        (config.get("meandepth_cutoff_short_read") == "None") & (config.get("evenness_cutoff_short_read") == "None")) & \
        (config.get("meandepth_cutoff_long_read") == "None") & (config.get("evenness_cutoff_long_read") == "None"):

    rule bypass_coverage_filter:
        input:
            contigs="{sample}/polish/cov_filter/{sample}_pre_filtered.fasta",
            contig_info="{sample}/assembly/{sample}_circular_info.tsv"
        output:
            contigs="{sample}/polish/cov_filter/{sample}_filtered_contigs.fasta",
            contig_info="{sample}/polish/cov_filter/{sample}_filtration_summary.tsv"
        run:
            import shutil
            import pandas as pd

            # Copy the pre-filtered contig file rather than perform filtration
            shutil.copyfile(str(input.contigs), str(output.contigs), follow_symlinks=True)

            # Write output summary file based on the assembly guide file
            contig_info = pd.read_csv(input.contig_info, sep='\t')
            contig_info['pass_coverage_filter'] = 'Y'

            contig_info.to_csv(output.contig_info, sep='\t', index=False)


else:
    rule summarize_contigs_by_coverage:
        input:
            coverage_info=expand("{{sample}}/polish/cov_filter/{{sample}}_{type}_coverage.tsv",
                type=filtration_method),
            contig_info="{sample}/assembly/{sample}_circular_info.tsv"
        output:
            contig_info_with_coverage = "{sample}/polish/cov_filter/{sample}_filtration_summary.tsv",
            filtered_contig_list="{sample}/polish/cov_filter/{sample}_filtered_contigs.list"
        params:
            meandepth_short=config.get("meandepth_cutoff_short_read"),
            evenness_short=config.get("evenness_cutoff_short_read"),
            meandepth_long=config.get("meandepth_cutoff_long_read"),
            evenness_long=config.get("evenness_cutoff_long_read")
        run:
            import pandas as pd

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

                # Set to zero if the user did not specify a filtration threshold
                mean_depth = 0 if mean_depth == "None" else mean_depth
                evenness = 0 if evenness == "None" else evenness

                coverage_filtered = coverage_data[ \
                    (coverage_data['meandepth'] >= mean_depth) & \
                    (coverage_data['coverage'] >= evenness)]

                return coverage_filtered['#rname']


            coverage_info_filepaths = list(input.coverage_info)

            # These are the expected filepaths of the coverage files that could be generated if short vs. long read
            #  coverage filtration is requested
            short_read_coverage_tsv_path = f"{wildcards.sample}/polish/cov_filter/{wildcards.sample}_short_read_coverage.tsv"
            long_read_coverage_tsv_path = f"{wildcards.sample}/polish/cov_filter/{wildcards.sample}_long_read_coverage.tsv"

            if len(coverage_info_filepaths) == 1:
                # Perform coverage filtration just based on short reads or just based on long reads if only one
                #  coverage file is provided
                if coverage_info_filepaths[0] == short_read_coverage_tsv_path:
                    contigs = filter_coverage_data(coverage_info_filepaths[0],params.meandepth_short,params.evenness_short)

                elif coverage_info_filepaths[0] == long_read_coverage_tsv_path:
                    contigs = filter_coverage_data(coverage_info_filepaths[0],params.meandepth_long,params.evenness_long)

                else:
                    sys.exit(f"Unexpected coverage file detected in 'polish/cov_filter': {coverage_info_filepaths[0]}")

            elif len(coverage_info_filepaths) == 2:
                # Sort using information from both short and long read coverage if both coverage files are provided
                coverage_info_filepaths.sort()

                if (coverage_info_filepaths[0] != long_read_coverage_tsv_path) | \
                        (coverage_info_filepaths[1] != short_read_coverage_tsv_path):
                    sys.exit("At least one unexpected coverage file detected in 'polish/cov_filter'.")

                set1 = set(filter_coverage_data(coverage_info_filepaths[0],params.meandepth_long,params.evenness_long))
                set2 = set(filter_coverage_data(coverage_info_filepaths[1],params.meandepth_short,params.evenness_short))

                contigs = pd.Series(list(set1.intersection(set2)))

            elif len(coverage_info_filepaths) == 0:
                sys.exit("No coverage files detected in 'polish/cov_filter'")

            else:
                sys.exit(f"More than 2 coverage files detected in 'polish/cov_filter': {coverage_info_filepaths}.")

            # Add the filtered contig information onto a contig info file
            contig_info = pd.read_csv(input.contig_info, sep='\t')

            passed_coverage_filter = []
            for contig_id in contig_info['contig']:
                if contig_id in contigs.values:
                    passed_coverage_filter.append('Y')
                else:
                    passed_coverage_filter.append('N')

            contig_info['pass_coverage_filter'] = passed_coverage_filter

            contig_info.to_csv(output.contig_info_with_coverage, sep='\t', index=False)
            contigs.to_csv(output.filtered_contig_list, header=None, index=False)


    rule filter_contigs_by_coverage:
        input:
            contigs="{sample}/polish/medaka/{sample}_consensus.fasta" if POLISH_WITH_SHORT_READS == False else "{sample}/polish/pypolca/{sample}_corrected.fasta",
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
        contigs="{sample}/polish/cov_filter/{sample}_filtered_contigs.fasta",
        contig_info="{sample}/polish/cov_filter/{sample}_filtration_summary.tsv"
    output:
        contigs="{sample}/polish/{sample}_polish.fasta",
        contig_info="{sample}/polish/{sample}_contig_info.tsv"
    run:
        source_relpath_contigs = os.path.relpath(str(input.contigs),os.path.dirname(str(output.contigs)))
        os.symlink(source_relpath_contigs,str(output.contigs))

        source_relpath_contig_info = os.path.relpath(str(input.contig_info),os.path.dirname(str(output.contig_info)))
        os.symlink(source_relpath_contig_info,str(output.contig_info))


rule polish:
    input:
        expand("{sample}/polish/{sample}_polish.fasta",sample=SAMPLE_NAMES),
        expand("{sample}/polish/{sample}_contig_info.tsv",sample=SAMPLE_NAMES)
    output:
        temp(touch("checkpoints/polish"))
