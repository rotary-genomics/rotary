#!/usr/bin/env python
# flye_end_repair.py
# Fixes ends of circular contigs produced by Flye
# Jackson M. Tsuji, Hokkaido University & JAMSTEC, 2023

import os
import sys
import time
import argparse
import logging
import shutil
import shlex
import subprocess
from Bio import SeqIO
import pandas as pd

# GLOBAL VARIABLES
SCRIPT_VERSION = '0.2.0'
DEPENDENCY_NAMES = ['flye', 'minimap2', 'samtools', 'circlator', 'seqtk']

# Set up the logger
logging.basicConfig(format='[ %(asctime)s UTC ]: %(levelname)s: %(message)s')
logging.Formatter.converter = time.gmtime
logger = logging.getLogger(__name__)


def check_dependency(dependency_name, dependency_type='shell'):
    """
    Checks if a required dependency is present

    :param dependency_name: name of the dependency
    :param dependency_type: specify the dependency type - currently only 'shell' command is available
    :return: string of the dependency type and name
    """

    if dependency_type == 'shell':

        if shutil.which(dependency_name) is None:

            logger.error(f'Dependency not found: {dependency_name}')
            raise RuntimeError

    else:

        logger.error(f'Available dependency types are "shell" and "python" - you provided "{dependency_type}"')
        raise RuntimeError

    dependency_message = f'{dependency_type} : {dependency_name}'

    # TODO: Consider reporting the path to the dependency for downstream analysis
    return dependency_message


def parse_assembly_info_file(assembly_info_filepath, return_type='circular'):
    """
    List circular and linear contigs from the Flye assembly info file

    :param assembly_info_filepath: path to assembly_info.txt output by Flye
    :param return_type: whether to return a list of 'circular' or 'linear' contigs
    :return: list of contig names, either circular or linear depending on return_type
    """

    logger.debug('Loading assembly info file')
    assembly_info = pd.read_csv(assembly_info_filepath, sep='\t')

    circular_contigs = assembly_info[assembly_info['circ.'] == 'Y'][['#seq_name', 'length', 'circ.']]
    linear_contigs = assembly_info[assembly_info['circ.'] == 'N'][['#seq_name', 'length', 'circ.']]

    if return_type == 'circular':
        output_list = list(circular_contigs['#seq_name'])
    elif return_type == 'linear':
        output_list = list(linear_contigs['#seq_name'])
    else:
        logger.error(f'return_type must be "circular" or "linear"; you provided "{return_type}"')
        raise RuntimeError

    return output_list


def subset_sequences(input_fasta_filepath, subset_sequence_ids):
    """
    Given an input FastA file, subsets the file to the provided sequence IDs. The output can be saved as a FastA file
    or saved in RAM as a SeqRecord. The order of records can optionally be sorted to be in the same order as the input
    IDs.

    :param input_fasta_filepath: Path(s) to the input FastA file(s); provide a list if multiple FastAs are desired
    :param subset_sequence_ids: list of names of the sequences to keep. If any names are in the list that aren't in the
                                input file, the function will not return them. In the case of conflicting IDs between
                                the multiple input FastA's, the script will take the first record with the correct name
                                and WILL NOT throw a warning.
    :return: generator of a SeqRecord object for the subset sequences
    """

    # Get the original (pre-stitched) contig sequence as a SeqRecord
    with open(input_fasta_filepath) as fasta_handle:

        for record in SeqIO.parse(fasta_handle, 'fasta'):

            if len(subset_sequence_ids) == 0:
                break

            if record.name in subset_sequence_ids:

                subset_sequence_ids.remove(record.name)

                yield record


def generate_bed_file(contig_seqrecord, bed_filepath, length_threshold=100000):
    """
    Generates a BED file for the desired contig for the desired bp threshold around that contig's ends

    :param contig_seqrecord: SeqRecord for the contig of interest
    :param bed_filepath: desired output filepath for the BED file
    :param length_threshold: length (bp) around the contig end to target in the BED file
    :return: writes BED file to the bed_filepath
    """

    contig_name = contig_seqrecord.name
    contig_length = len(contig_seqrecord.seq)

    logger.debug('Making BED file with end proximity threshold of ' + str(length_threshold))

    if contig_length < length_threshold:
        logger.warning(f'Contig length ({contig_length}) is less than the supplied length threshold '
                       f'({length_threshold}), so the BED file will be for the whole contig.')

        contigs = [contig_name]
        starts = [0]
        stops = [contig_length]

    else:
        half_threshold = int(length_threshold / 2)

        contigs = [contig_name, contig_name]
        starts = [0, contig_length - half_threshold]
        stops = [half_threshold, contig_length]

    end_regions = pd.DataFrame({'contig': contigs, 'start': starts, 'stop': stops})

    end_regions.to_csv(bed_filepath, sep='\t', header=None, index=False)


def map_long_reads(contig_filepath, long_read_filepath, output_bam_filepath, log_filepath, append_log: bool = True,
                   threads=1, thread_mem=1):
    """
    Maps long reads (via minimap2) to contigs and sorts/indexes the resulting BAM file

    :param contig_filepath: file containing contigs to map the reads to (FastA)
    :param long_read_filepath: file containing long reads to map to the contigs (FastQ; compression is OK)
    :param output_bam_filepath: path to the BAM file to be saved
    :param log_filepath: path to the log file to be saved
    :param append_log: whether to log should append onto an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :param thread_mem: memory in GB per thread (to use for samtools)
    :return: output_bam_filepath and log_filepath are saved to disk; nothing is returned.
    """

    if append_log is True:
        write_mode = 'a'
    elif append_log is False:
        write_mode = 'w'
    else:
        raise ValueError

    with open(log_filepath, write_mode) as logfile_handle:
        with open(output_bam_filepath, 'w') as bam_handle:

            # TODO - add support for different flags like -ax for pacbio
            minimap_args = ['minimap2', '-t', threads, '-ax', 'map-ont', contig_filepath, long_read_filepath]
            logger.debug(f'{shlex.join(minimap_args)} | \\')
            minimap = subprocess.run(minimap_args, check=True, stdout=subprocess.PIPE, stderr=logfile_handle)

            samtools_view_args = ['samtools', 'view', '-b', '-@', threads]
            logger.debug(f'{shlex.join(minimap_args)} | \\')
            samtools_view = subprocess.run(samtools_view_args, check=True, input=minimap.stdout,
                                           stdout=subprocess.PIPE, stderr=logfile_handle)

            samtools_sort_args = ['samtools', 'sort', '-@', threads, '-m', f'{thread_mem}G']
            logger.debug(shlex.join(samtools_sort_args))
            samtools_sort = subprocess.run(samtools_sort_args, check=True, input=samtools_view.stdout,
                                           stdout=bam_handle, stderr=logfile_handle)

        samtools_index_args = ['samtools', 'index', '-@', threads, output_bam_filepath]
        logger.debug(shlex.join(samtools_index_args))
        samtools_index = subprocess.run(samtools_index_args, check=True, stderr=logfile_handle)

    logger.debug('Read mapping finished')


def get_selected_reads(bam_filepath, bed_filepath, output_ends_fastq_filepath, log_filepath,
                       append_log: bool = True, threads=1):
    """
    Pulls mapped reads from a BAM file that were mapped to regions defined in a BED file and saves to FastQ

    :param bam_filepath: BAM file containing mapped reads to a reference; needs to be sorted and indexed I think
    :param bed_filepath: BED file containing the regions of contigs to recover reads for
    :param output_ends_fastq_filepath: path to the FastQ file to be saved (.fastq.gz extension saves as Gzipped FastQ)
    :param log_filepath: path to the log file to be saved
    :param append_log: whether to log should append onto an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :return: output_ends_fastq_filepath is saved to disk; nothing is returned.
    """

    if append_log is True:
        write_mode = 'a'
    elif append_log is False:
        write_mode = 'w'
    else:
        raise ValueError

    with open(log_filepath, write_mode) as logfile_handle:

        # TODO - consider adding option to split long reads in half if they go around a short circular contig, like in circlator
        samtools_view_args = ['samtools', 'view', '-@', threads, '-L', bed_filepath, '-b', bam_filepath]
        logger.debug(f'{shlex.join(samtools_view_args)} | \\')
        samtools_view = subprocess.run(samtools_view_args, check=True, stdout=subprocess.PIPE, stderr=logfile_handle)

        samtools_fastq_args = ['samtools', 'fastq', '-0', output_ends_fastq_filepath, '-n', '-@', threads]
        logger.debug(shlex.join(samtools_fastq_args))
        samtools_fastq = subprocess.run(samtools_fastq_args, check=True, input=samtools_view.stdout,
                                        stderr=logfile_handle)

        # TODO - delete this CLI code reference once I confirm the Python code works
        # samtools view -@ "${threads}" -L "${regions}" -b "${bam_file}" 2>> "${verbose}" | \
        #   samtools fastq -0 "${fastq}" -n -@ "${threads}" 2>> "${verbose}"


def run_flye(fastq_filepath, flye_outdir, flye_read_mode, flye_read_error, log_filepath,
             append_log: bool = True, threads=1):
    """
    Pulls mapped reads from a BAM file that were mapped to regions defined in a BED file and saves to FastQ

    :param fastq_filepath: path to the input read FastQ file (gzipped is OK)
    :param flye_outdir: directory to save Flye output to
    :param flye_read_mode: type of raw reads, either 'nano-raw' or 'nano-hq'
    :param flye_read_error: expected error rate of reads as a proportion; specify 0 to use default Flye settings
    :param log_filepath: path to the log file to be saved
    :param append_log: whether to log should append onto an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :return: flye output is saved to disk at flye_outdir; nothing is returned.
    """

    if append_log is True:
        write_mode = 'a'
    elif append_log is False:
        write_mode = 'w'
    else:
        raise ValueError

    flye_args = ['flye', f'--{flye_read_mode}', fastq_filepath, '-o', flye_outdir, '-t', threads]

    if flye_read_error != 0:

        flye_args.append('--read_error')
        flye_args.append(flye_read_error)

    with open(log_filepath, write_mode) as logfile_handle:
        # TODO - add error handling if flye fails
        logger.debug(shlex.join(flye_args))
        subprocess.run(flye_args, check=True, stderr=logfile_handle)

        # TODO - delete CLI code example once I confirm Python is working
        # flye "--${flye_read_mode}" "${fastq}" -o "${flye_length_outdir}" --read_error "${flye_read_error}" \
        #           -t "${threads}" >> "${verbose}" 2>&1


def run_circlator_merge(original_contig_filepath, patch_contig_filepath, merge_outdir, circlator_min_id,
                        circlator_min_length, circlator_ref_end, circlator_reassemble_end,
                        log_filepath, append_log: bool = True):
    """
    Pulls mapped reads from a BAM file that were mapped to regions defined in a BED file and saves to FastQ

    :param original_contig_filepath: path to the original (non-stitched) contig file (FastA)
    :param patch_contig_filepath: path to the contig(s) that should span the circularization point in the original
                                  contig (FastA)
    :param merge_outdir: directory to save merge output to
    :param circlator_min_id: Percent identity threshold for circlator merge
    :param circlator_min_length: Minimum required overlap (bp) between original and merge contigs
    :param circlator_ref_end: Minimum required overlap (bp) between original and merge contigs
    :param circlator_reassemble_end: Minimum distance (bp) between end of merge contig and nucmer hit
    :param log_filepath: path to the log file to be saved
    :param append_log: whether to log should append onto an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :return: circlator merge output is saved to disk at merge_outdir; nothing is returned.
    """

    if append_log is True:
        write_mode = 'a'
    elif append_log is False:
        write_mode = 'w'
    else:
        raise ValueError

    os.makedirs(merge_outdir, exist_ok=True)

    with open(log_filepath, write_mode) as logfile_handle:

        circlator_merge_args = ['circlator', 'merge', '--verbose', '--mid_id', circlator_min_id, '--min_length',
                                circlator_min_length, '--ref_end', circlator_ref_end, '--reassemble_end',
                                circlator_reassemble_end, original_contig_filepath, patch_contig_filepath,
                                os.path.join(merge_outdir, 'merge')]
        logger.debug(shlex.join(circlator_merge_args))
        subprocess.run(circlator_merge_args, check=True, stderr=logfile_handle)

        # TODO - delete CLI example once I confirm Python is working
        # circlator merge --verbose --min_id "${circlator_min_id}" --min_length "${circlator_min_length}" \
        #         --ref_end "${circlator_ref_end}" --reassemble_end "${circlator_reassemble_end}" \
        #         "${assembly}" "${reassembly_dir}/assembly.fasta" "${merge_dir}/merge" >> "${verbose}" 2>&1


def check_circlator_success(circlator_logfile):
    """
    Checks the circlator log file to see if contig stitching was successful

    :param circlator_logfile: path to the circlator log file
    :return: Boolean of whether the contigs were successfully stitched or not
    """

    logger.debug('Loading circlator logfile')
    circlator_info = pd.read_csv(circlator_logfile, sep='\t')[['#Contig', 'stitched']]

    # TODO - check I understood this properly when making edits (look at a real file as an example)
    # If a row is '1', it means it was stitched properly, but if '0', it means it was not stitched.
    # So if all rows are 1 (i.e., the sum of rows / # of rows is 1), it means everything was stiched properly.
    if circlator_info['stitched'].sum() / circlator_info.shape[0] == 1:
        logger.debug('Everything is stitched.')
        result = True

    elif circlator_info['stitched'].sum() >= 0:
        logger.debug('Not everything is stitched.')
        result = False

    else:
        logger.error('File processing error. # of non-stitched contigs is not >=0. Exiting...')
        raise RuntimeError

    return result


def rotate_contig_to_midpoint(contig_fasta_filepath, output_filepath):
    """
    Rotates an input (circular) contig to its approximate midpoint

    :param contig_fasta_filepath: FastA file containing a single circular contig
    :param output_filepath: Filepath for the output rotated FastA file
    :return: writes FastA file to output_filepath, 60 nt per line
    """

    # Read contig FastA
    contig_count = 0

    with open(contig_fasta_filepath) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):

            if contig_count == 0:
                contig_record = record

            elif contig_count > 0:
                logger.error('More than one contig in input FastA file')
                raise RuntimeError

            contig_count = contig_count + 1

    contig_length = len(contig_record.seq)
    contig_midpoint = int(contig_length / 2)

    logger.debug(f'Rotating contig to midpoint at {contig_midpoint} bp')
    contig_sequence_rotated_front = contig_record.seq[contig_midpoint:contig_length]
    contig_sequence_rotated_back = contig_record.seq[0:contig_midpoint]
    contig_sequence_rotated = contig_sequence_rotated_front + contig_sequence_rotated_back

    # Update SeqRecord
    contig_record.seq = contig_sequence_rotated
    contig_record.description = contig_record.name  # trim off description

    # Write
    with open(output_filepath, 'w') as output_handle:
        SeqIO.write(contig_record, output_handle, 'fasta')


def run_end_repair(long_read_filepath, assembly_fasta_filepath, assembly_info_filepath, output_dir,
                   flye_read_mode, flye_read_error, length_cutoffs, circlator_min_id, circlator_min_length,
                   circlator_ref_end, circlator_reassemble_end, threads, thread_mem):

    # Define core file names and directory structures
    # These will be the files and folders in the main output directory:
    end_repaired_contigs_filepath = os.path.join(output_dir, 'repaired.fasta')
    verbose_logfile = os.path.join(output_dir, 'verbose.log')
    bam_filepath = os.path.join(output_dir, 'long_read.bam')
    circlator_logdir = os.path.join(output_dir, 'circlator_logs')
    reassembly_outdir_base = os.path.join(output_dir, 'contigs')
    circular_contig_tmp_fasta = os.path.join(reassembly_outdir_base, 'circular_input.fasta')

    # Check output dir
    if os.path.isdir(output_dir):
        logger.warning(f'Output directory already exists: "{output_dir}"; files may be overwritten.')

    os.makedirs(output_dir, exist_ok=True)

    # Get lists of circular contigs
    circular_contig_names = parse_assembly_info_file(assembly_info_filepath, return_type='circular')

    # No need to run the pipeline if there are no circular contigs
    if len(circular_contig_names) == 0:

        logger.info('No circular contigs. Will copy the input file and finish early.')
        shutil.copyfile(assembly_fasta_filepath, end_repaired_contigs_filepath)
        logger.info('Pipeline finished.')
        sys.exit(0)

    # Subset circular contigs from the full file and save to disk
    with open(circular_contig_tmp_fasta, 'w') as output_handle:
        for record in subset_sequences(assembly_fasta_filepath, circular_contig_names):
            SeqIO.write(record, output_handle, 'fasta')

    # Initialize the repaired contigs FastA file (so it will overwrite an old file rather than just append later)
    with open(end_repaired_contigs_filepath, 'w') as output_handle:
        output_handle.write('')

    os.makedirs(circlator_logdir, exist_ok=True)
    os.makedirs(reassembly_outdir_base, exists_ok=True)

    # Start the main workflow from here:
    logger.info('Mapping reads to all contigs')
    map_long_reads(contig_filepath=assembly_fasta_filepath, long_read_filepath=long_read_filepath,
                   output_bam_filepath=bam_filepath, log_filepath=verbose_logfile, append_log=False,
                   threads=threads, thread_mem=thread_mem)

    failed_contigs = 0
    for contig_name in circular_contig_names:

        logger.info(f'End repair: {contig_name}')

        # Define temp files and folders that will be generated for this contig during stitching
        reassembly_outdir = os.path.join(reassembly_outdir_base, contig_name)
        log_dir_base = os.path.join(reassembly_outdir, 'logs')
        original_contig_fasta_filepath = os.path.join(reassembly_outdir, f'{contig_name}.fasta')
        reassembly_outfile = os.path.join(reassembly_outdir, f'{contig_name}_stitched.fasta')

        os.makedirs(reassembly_outdir, exist_ok=True)
        os.makedirs(log_dir_base, exist_ok=True)

        # Get the original (pre-stitched) contig sequence as a SeqRecord and save as FastA
        contig_records = []
        with open(original_contig_fasta_filepath, 'w') as output_handle:
            for record in subset_sequences(circular_contig_tmp_fasta, [contig_name]):
                contig_records.append(record)
                SeqIO.write(record, output_handle, 'fasta')

        if len(contig_records) > 1:
            logger.error(f'Got more than one contig for {contig_name}')
            raise RuntimeError

        contig_record = contig_records[0]

        # Keep trying to stitch the contig, iterating over different sizes of the region around the contig ends to
        #   target, until stitching works
        assembly_attempts = 0
        linked_ends = False

        while linked_ends is False:

            # End the script if all length cutoffs have been tried
            if assembly_attempts > len(length_cutoffs)-1:

                logger.warning(f'End repair: {contig_name}: FAILED to linked contig ends')
                failed_contigs = failed_contigs + 1
                os.makedirs(os.path.join(output_dir, 'troubleshooting'), exist_ok=True)
                shutil.move(log_dir_base, os.path.join(output_dir, 'troubleshooting', contig_name))
                shutil.rmtree(reassembly_outdir)

                break

            # Get the length cutoff for this attempt
            length_cutoff = length_cutoffs[assembly_attempts]

            if len(contig_record.seq) <= length_cutoff:
                logger.debug(f'Skipping length cutoff of {length_cutoff} because '
                             f'contig is shorter than this ({len(contig_record.seq)} bp)')
                continue

            logger.debug(f'Starting reassembly with a length cutoff of {length_cutoff} bp')

            # Define folders and filenames that will be used for this length cutoff iteration
            length_outdir = os.path.join(reassembly_outdir, f'L{length_cutoff}')
            bed_filepath = os.path.join(length_outdir, 'ends.bed')
            ends_fastq_filepath = os.path.join(length_outdir, 'ends.fastq.gz')
            flye_length_outdir = os.path.join(length_outdir, 'assembly')
            merge_dir = os.path.join(length_outdir, 'merge')
            log_dir = os.path.join(log_dir_base, f'L{length_cutoff}')

            os.makedirs(length_outdir, exist_ok=True)
            os.makedirs(log_dir, exist_ok=True)

            # Get reads for the selected region around the contig ends
            generate_bed_file(contig_record, bed_filepath, length_threshold=length_cutoff)
            get_selected_reads(bam_filepath=bam_filepath, bed_filepath=bed_filepath,
                               output_ends_fastq_filepath=ends_fastq_filepath, log_filepath=verbose_logfile,
                               append_log=True, threads=threads)

            # Assemble the reads to get (hopefully) a joined contig end
            run_flye(fastq_filepath=ends_fastq_filepath, flye_outdir=flye_length_outdir,
                     flye_read_mode=flye_read_mode, flye_read_error=flye_read_error, log_filepath=verbose_logfile,
                     append_log=True, threads=threads)
            shutil.copy(os.path.join(flye_length_outdir, 'assembly_info.txt'), log_dir)

            # Stitch the joined contig end onto the original assembly
            run_circlator_merge(original_contig_filepath=original_contig_fasta_filepath,
                                patch_contig_filepath=os.path.join(flye_length_outdir, 'assembly.fasta'),
                                merge_outdir=merge_dir, circlator_min_id=circlator_min_id,
                                circlator_min_length=circlator_min_length, circlator_ref_end=circlator_ref_end,
                                circlator_reassemble_end=circlator_reassemble_end, log_filepath=verbose_logfile,
                                append_log=True)

            shutil.copy(os.path.join(merge_dir, 'merge.circularise.log'), log_dir)
            shutil.copy(os.path.join(merge_dir, 'merge.circularise_details.log'), log_dir)

            if check_circlator_success(os.path.join(merge_dir, 'merge.circularise.log')) is True:

                logger.info('End repair: successfully linked contig ends')

                # Rotate to midpoint so that the stitched points around the midpoint can be polished more effectively
                # TODO - sometimes small contigs are already rotated far from original origin. Stitch point is
                #  hard to find. Does circlator report stitch point?
                rotate_contig_to_midpoint(os.path.join(merge_dir, 'merge.fasta'), reassembly_outfile)

                # Append the contig onto the main repaired contig output file
                with open(reassembly_outfile) as input_handle:
                    with open(end_repaired_contigs_filepath, 'a') as append_handle:
                        append_handle.write(input_handle.read())

                # Cleanup
                shutil.copy(os.path.join(merge_dir, 'merge.circularise_details.log'),
                            os.path.join(circlator_logdir, f'{contig_name}.log'))
                shutil.rmtree(reassembly_outdir)

                linked_ends = True

            assembly_attempts = assembly_attempts + 1

    if failed_contigs > 0:

        logger.error(f'{failed_contigs} contigs could not be circularized. A partial output file including '
                     f'successfully circularized contigs (and no linear contigs) is available at '
                     f'{end_repaired_contigs_filepath} for debugging. Exiting with error status. See temporary files '
                     f'and verbose logs for more details.')
        sys.exit(1)

    linear_contig_names = parse_assembly_info_file(assembly_info_filepath, return_type='linear')

    with open(end_repaired_contigs_filepath, 'a') as append_handle:

        for record in subset_sequences(assembly_fasta_filepath, linear_contig_names):
            SeqIO.write(record, append_handle, 'fasta')

    # Clean up temp files
    os.remove(bam_filepath)
    os.remove(f'{bam_filepath}.bai')
    shutil.rmtree(reassembly_outdir_base)

    logger.info(f'End repair finished. Output contigs saved at {end_repaired_contigs_filepath}.')


def main(args):
    # Set user variables
    long_read_filepath = args.long_reads
    assembly_fasta_filepath = args.assembly_fasta
    assembly_info_filepath = args.assembly_info
    output_dir = args.output_dir
    flye_read_mode = args.flye_read_mode
    flye_read_error = args.flye_read_error
    length_cutoffs = [int(x) for x in args.length_cutoffs.split(',')]
    circlator_min_id = args.circlator_min_id
    circlator_min_length = args.circlator_min_length
    circlator_ref_end = args.circlator_ref_end
    circlator_reassemble_end = args.circlator_reassemble_end
    threads = args.threads
    thread_mem = args.thread_mem
    verbose = args.verbose

    # Startup checks
    if verbose is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Check dependencies
    for dependency_name in DEPENDENCY_NAMES:
        check_dependency(dependency_name)

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.info('Version: ' + SCRIPT_VERSION)
    logger.info('### SETTINGS ###')
    logger.info(f'Long read filepath: {long_read_filepath}')
    logger.info(f'Assembly FastA filepath: {assembly_fasta_filepath}')
    logger.info(f'Assembly info filepath: {assembly_info_filepath}')
    logger.info(f'Output directory: {output_dir}')
    logger.info(f'Flye read mode: {flye_read_mode}')
    logger.info(f'Flye read error: {flye_read_error}')
    logger.info(f'Length cutoffs to test (bp): {length_cutoffs}')
    logger.info(f'Circlator min. ID: {circlator_min_id}')
    logger.info(f'Circlator min. length: {circlator_min_length}')
    logger.info(f'Circlator ref. end: {circlator_ref_end}')
    logger.info(f'Circlator reassembly end: {circlator_reassemble_end}')
    logger.info(f'Threads: {threads}')
    logger.info(f'Memory per thread (GB): {thread_mem}')
    logger.info(f'Verbose logging: {verbose}')
    logger.info('################')

    run_end_repair(long_read_filepath, assembly_fasta_filepath, assembly_info_filepath, output_dir,
                   flye_read_mode, flye_read_error, length_cutoffs, circlator_min_id, circlator_min_length,
                   circlator_ref_end, circlator_reassemble_end, threads, thread_mem)

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'{os.path.basename(sys.argv[0])}: pipeline to repair ends of circular contigs from Flye. \n'
                    f'Copyright Jackson M. Tsuji, Hokkaido University / JAMSTEC, 2023. \n'
                    f'Version: {SCRIPT_VERSION}')
    parser.add_argument('longreads.fastq.gz', dest='long_read_filepath', help='QC-passing Nanopore reads')
    parser.add_argument('assembly.fasta', dest='assembly_fasta_filepath', help='Contigs output from Flye')
    parser.add_argument('assembly_info.txt', dest='assembly_info_filepath', help='Assembly info file from Flye')
    parser.add_argument('output_dir',
                        help='Output directory path (might overwrite contents if the dir already exists!)')
    parser.add_argument('-f', '--flye_read_mode', required=False, default='nano-hq', choices=['nano-hq','nano-raw'],
                        help='Type of input reads for Flye')
    parser.add_argument('-F', '--flye_read_error', required=False, default=0, type=float,
                        help='Expected error rate of input reads, expressed as proportion (e.g., 0.03) '
                             '[0 = use flye defaults]')
    parser.add_argument('-c', '--length_cutoffs', required=False,
                        default='100000,90000,80000,70000,60000,50000,20000,10000,5000,2000', type=str,
                        help='Comma-separated list of length thresholds for reassembly around the contig ends (bp)')
    parser.add_argument('-i', '--circlator_min_id', required=False, default=99, type=float,
                        help='Percent identity threshold for circlator merge')
    parser.add_argument('-l', '--circlator_min_length', required=False, default=10000, type=int,
                        help='Minimum required overlap (bp) between original and merge contigs')
    parser.add_argument('-e', '--circlator_ref_end', required=False, default=100, type=int,
                        help='Minimum distance (bp) between end of original contig and nucmer hit')
    parser.add_argument('-E', '--circlator_reassemble_end', required=False, default=100, type=int,
                        help='Minimum distance (bp) between end of merge contig and nucmer hit')
    parser.add_argument('-t', '--threads', required=False, default=1, type=int,
                        help='Number of processors threads to use')
    parser.add_argument('-m', '--threads_mem', required=False, default=1, type=float,
                        help='Memory (GB) to use **per thread** for samtools sort')
    parser.add_argument('-v', '--verbose', required=False, action='store_false',
                        help='Enable for verbose logging.')
    command_line_args = parser.parse_args()
    main(command_line_args)

