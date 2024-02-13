#!/usr/bin/env python
# repair.py
# Fixes ends of circular contigs produced by Flye
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2023

import argparse
import logging
import os
import shlex
import shutil
import subprocess
import sys

import pandas as pd
from Bio import SeqIO

# GLOBAL VARIABLES
DEPENDENCY_NAMES = ['flye', 'minimap2', 'samtools', 'circlator']

# Set up the logger
logger = logging.getLogger(__name__)
formatter = logging.Formatter('[ %(asctime)s ]: %(levelname)s: %(funcName)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


def main():
    """
    Collects input arguments and runs the end repair workflow
    """

    """
    When installing through pip with pyprojects.toml a new python script is generated
    that calls main() out of rotary.py. Run parser inside main() so it can be called
    externally as a function.
    """
    parser = parse_cli()
    args = parser.parse_args()

    # Startup checks
    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Check output dir
    output_dir_exists = os.path.isdir(args.output_dir)

    if (output_dir_exists is True) & (args.overwrite is False):

        logger.error(f'Output directory already exists: "{args.output_dir}". Will not continue. Set the '
                     f'--overwrite flag at your own risk if you want to use an existing directory.')
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Start log file in the output dir
    file_handler = logging.FileHandler(filename=os.path.join(args.output_dir, 'repaired.log'), mode='w')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

    # Check dependencies
    dependency_paths = []
    for dependency_name in DEPENDENCY_NAMES:
        dependency_paths.append(check_dependency(dependency_name))

    dependency_dict = dict(zip(DEPENDENCY_NAMES, dependency_paths))

    # Parse some command line inputs further
    assembly_info_type = 'custom' if args.custom_assembly_info_file is True else 'flye'
    # TODO - improve error handling if a string is provided instead of a real length
    length_thresholds = [int(x) for x in args.length_thresholds.split(',')]
    cli_tool_settings_dict = {'flye_read_mode': args.flye_read_mode,
                              'flye_read_error': args.flye_read_error,
                              'circlator_min_id': args.circlator_min_id,
                              'circlator_min_length': args.circlator_min_length,
                              'circlator_ref_end': args.circlator_ref_end,
                              'circlator_reassemble_end': args.circlator_reassemble_end}
    threads_mem_mb = int(args.threads_mem * 1024)

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.debug('### SETTINGS ###')
    logger.debug(f'Long read filepath: {args.long_read_filepath}')
    logger.debug(f'Assembly FastA filepath: {args.assembly_fasta_filepath}')
    logger.debug(f'Assembly info filepath: {args.assembly_info_filepath}')
    logger.debug(f'Use a custom assembly info file: {args.custom_assembly_info_file}')
    logger.debug(f'Output directory: {args.output_dir}')
    logger.debug(f'Overwrite output dir if needed?: {args.overwrite}')
    logger.debug(f'Length thresholds to test (bp): {length_thresholds}')
    logger.debug(f'Keep going if some contigs cannot be re-circularized?: {args.keep_going_with_failed_contigs}')
    logger.debug(f'Flye read mode: {cli_tool_settings_dict["flye_read_mode"]}')
    logger.debug(f'Flye read error: {cli_tool_settings_dict["flye_read_error"]}')
    logger.debug(f'Circlator min. ID: {cli_tool_settings_dict["circlator_min_id"]}')
    logger.debug(f'Circlator min. length: {cli_tool_settings_dict["circlator_min_length"]}')
    logger.debug(f'Circlator ref. end: {cli_tool_settings_dict["circlator_ref_end"]}')
    logger.debug(f'Circlator reassembly end: {cli_tool_settings_dict["circlator_reassemble_end"]}')
    logger.debug(f'Threads: {args.threads}')
    logger.debug(f'Memory per thread (GB = MB): {args.threads_mem} = {threads_mem_mb}')
    logger.debug(f'Verbose logging: {args.verbose}')
    logger.debug('### DEPENDENCIES ###')
    for key, value in dependency_dict.items():
        logger.debug(f'{key}: {value}')
    logger.debug('################')

    # Record a warning if the output dir was already there before the script started
    if output_dir_exists is True:
        logger.warning(f'Output directory already exists: "{args.output_dir}". Files may be overwritten.')

    run_end_repair(args.long_read_filepath, args.assembly_fasta_filepath, args.assembly_info_filepath,
                   assembly_info_type, args.output_dir, length_thresholds, args.keep_going_with_failed_contigs,
                   cli_tool_settings_dict, args.threads, threads_mem_mb)

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


def check_dependency(dependency_name: str):
    """
    Checks if a required shell dependency is present
    :param dependency_name: name of the dependency
    :return: path to the dependency
    """

    dependency_path = shutil.which(dependency_name)

    if dependency_path is None:

        logger.error(f'Dependency not found: {dependency_name}')
        raise RuntimeError

    return dependency_path


def parse_assembly_info_file(assembly_info_filepath: str, info_type: str, return_type: str = 'circular'):
    """
    List circular and linear contigs from a Flye (or custom format) assembly info file
    :param assembly_info_filepath: path to assembly_info.txt output by Flye
    :param info_type: whether the info file is in 'flye' format or is a 'custom' format.
                      'flye' format refers to the 'assembly_info.txt' format output by Flye after a successful assembly.
                      'custom' info files are tab-separated, have no headers, and have two columns: contig name and
                      contig type, either 'circular' or 'linear'.
    :param return_type: whether to return a list of 'circular' or 'linear' contigs
    :return: list of contig names (either circular or linear depending on return_type)
    """

    if info_type == 'flye':

        logger.debug('Loading Flye-format assembly info file')

        assembly_info = pd.read_csv(assembly_info_filepath, sep='\t')

        circular_contigs = assembly_info[assembly_info['circ.'] == 'Y'][['#seq_name', 'length', 'circ.']]
        linear_contigs = assembly_info[assembly_info['circ.'] == 'N'][['#seq_name', 'length', 'circ.']]

    elif info_type == 'custom':

        logger.debug('Loading custom format assembly info file')

        assembly_info = pd.read_csv(assembly_info_filepath, sep='\t', header=None)
        assembly_info.columns = ['#seq_name', 'status']

        expected_statuses = {'circular', 'linear'}
        actual_statuses = set(assembly_info['status'])
        unexpected_statuses = actual_statuses.difference(expected_statuses)

        if len(unexpected_statuses) != 0:
            logger.warning(f'Some entries in the assembly info file had unexpected contig status names, i.e.: '
                           f'{", ".join(unexpected_statuses)}')
            logger.warning('These entries will be treated as linear contigs... they will not be rotated and will be '
                           'returned as-is at the end of the script. Please make sure you did not make a typo or '
                           'include a header for your custom assembly info file.')

        circular_contigs = assembly_info[assembly_info['status'] == 'circular']
        linear_contigs = assembly_info[assembly_info['status'] != 'circular']

    else:

        logger.error(f'Assembly info file format should be "flye" or "custom"; you provided {info_type}')
        raise ValueError

    # Check for duplicate sequence IDs
    if assembly_info['#seq_name'].drop_duplicates().shape[0] < assembly_info['#seq_name'].shape[0]:
        duplicated_ids = set(assembly_info['#seq_name'][assembly_info['#seq_name'].duplicated() == True])

        logger.error(f'Some sequence IDs are duplicated in the input assembly info file: {", ".join(duplicated_ids)}')
        raise RuntimeError

    if return_type == 'circular':

        output_list = list(circular_contigs['#seq_name'])
        logger.debug(f'Found {len(output_list)} circular sequences')

    elif return_type == 'linear':

        output_list = list(linear_contigs['#seq_name'])
        logger.debug(f'Found {len(output_list)} linear sequences')

    else:

        logger.error(f'Return_type must be "circular" or "linear"; you provided "{return_type}"')
        raise RuntimeError

    return output_list


def subset_sequences(input_fasta_filepath: str, subset_sequence_ids: list):
    """
    Given an input FastA file, subsets the file to the provided sequence IDs
    :param input_fasta_filepath: Path to the input FastA file
    :param subset_sequence_ids: list of the names of sequences to keep. If any names in the list are not in the
                                input file, the function will not return anything for those names. The function will
                                raise an error if any duplicate sequence names are detected in the input FastA file.
    :return: generator of a SeqRecord object for the subset sequences
    """

    sequence_names = []

    # Get the original (pre-stitched) contig sequence as a SeqRecord
    with open(input_fasta_filepath) as fasta_handle:

        for record in SeqIO.parse(fasta_handle, 'fasta'):

            sequence_names.append(record.name)

            if record.name in subset_sequence_ids:

                yield record

    # Raise an error if there are duplicate sequence names
    if len(set(sequence_names)) < len(sequence_names):
        sequence_names_series = pd.Series(sequence_names)
        duplicates_names = set(sequence_names_series[sequence_names_series.duplicated() == True])

        logger.error(f'Duplicate sequence IDs were detected in the input FastA file "{input_fasta_filepath}": '
                     f'{", ".join(duplicates_names)}')
        raise RuntimeError


def set_write_mode(append_log: bool):
    """
    Converts the boolean append_log to 'w' or 'a' write modes
    :param append_log: boolean of whether to append to an existing log file (True) or to overwrite an existing log
                       file (False)
    :return: string of either 'a' (append mode) or 'w' (write mode)
    """

    if append_log is True:
        write_mode = 'a'

    elif append_log is False:
        write_mode = 'w'

    else:

        logger.error(f'append_log should be a boolean True or False; you provided {append_log}')
        raise ValueError

    return write_mode


def generate_bed_file(contig_seqrecord: SeqIO.SeqRecord, bed_filepath: str, length_threshold: int = 100000):
    """
    Generates a BED file for a desired region around the ends of a contig
    :param contig_seqrecord: SeqRecord of the contig
    :param bed_filepath: desired output filepath for the BED file
    :param length_threshold: length (bp) around the contig end to target in the BED file.
                             Half of this length will be returned around each end of the contig.
    :return: writes the BED file to the bed_filepath
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


def map_long_reads(contig_filepath: str, long_read_filepath: str, output_bam_filepath: str, log_filepath: str,
                   append_log: bool = True, threads: int = 1, threads_mem_mb: float = 1):
    """
    Maps long reads (via minimap2) to contigs and sorts/indexes the resulting BAM file
    :param contig_filepath: path to the FastA file containing the reference contigs
    :param long_read_filepath: path to the FastQ file containing long reads to map (compressed is OK)
    :param output_bam_filepath: path to the BAM file to be saved
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :param threads_mem_mb: memory in MB per thread (to use for samtools); must be an integer
    :return: output_bam_filepath and log_filepath are saved to disk
    """

    write_mode = set_write_mode(append_log)

    with open(log_filepath, write_mode) as logfile_handle:
        with open(output_bam_filepath, 'w') as bam_handle:
            # TODO - add support for different flags like -ax for pacbio
            minimap_args = ['minimap2', '-t', str(threads), '-ax', 'map-ont', contig_filepath, long_read_filepath]
            minimap = run_pipeline_subcommand(command_args=minimap_args, stdout=subprocess.PIPE, stderr=logfile_handle)

            samtools_view_args = ['samtools', 'view', '-b', '-@', str(threads)]
            samtools_view = run_pipeline_subcommand(command_args=samtools_view_args, stdin=minimap,
                                                    stdout=subprocess.PIPE, stderr=logfile_handle)

            samtools_sort_args = ['samtools', 'sort', '-@', str(threads), '-m', f'{threads_mem_mb}M']
            run_pipeline_subcommand(command_args=samtools_sort_args, stdin=samtools_view, stdout=bam_handle,
                                    stderr=logfile_handle)

        samtools_index_args = ['samtools', 'index', '-@', str(threads), output_bam_filepath]
        run_pipeline_subcommand(command_args=samtools_index_args, stderr=logfile_handle)

    logger.debug('Read mapping finished')


def run_pipeline_subcommand(command_args, stdin=None, stdout=None, stderr=None, check=True):
    """
    Wrapper function for running subcommands.

    :param command_args: The command line arguments of the subcommand (e.g., ['samtools', '-h'])
    :param stdin: A subprocess.PIPE or None if stdin is not to be used.
    :param stdout: Where to send stdout or None if stdout is not to be used.
    :param stderr: Where to send stderr or None if stderr is not to be used.
    :param check: Cause a runtime error if the subcommand fails.
    :return: The output of the subcommand.
    """
    logger.debug(shlex.join(command_args))

    if stdin:
        stdin = stdin.stdout

    return subprocess.run(command_args, check=check, input=stdin, stdout=stdout, stderr=stderr)


def subset_reads_from_bam(bam_filepath: str, bed_filepath: str, subset_fastq_filepath: str, log_filepath: str,
                          append_log: bool = True, threads: int = 1):
    """
    Subsets reads from a BAM file that were mapped to regions defined in a BED file; saves reads to a FastQ file
    :param bam_filepath: path to a BAM file containing reads mapped to a reference; BAM needs to be sorted and indexed
    :param bed_filepath: path to a BED file containing the regions of reference contigs to subset reads for
    :param subset_fastq_filepath: path to the FastQ file to be saved (.fastq.gz extension saves as Gzipped FastQ)
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :return: subset_fastq_filepath is saved to disk
    """

    write_mode = set_write_mode(append_log)

    with open(log_filepath, write_mode) as logfile_handle:

        # TODO - consider adding option to split long reads in half if they go around a short circular contig,
        #  like in circlator
        samtools_view_args = ['samtools', 'view', '-@', str(threads), '-L', bed_filepath, '-b', bam_filepath]
        samtools_view = run_pipeline_subcommand(command_args=samtools_view_args, stdout=subprocess.PIPE,
                                                stderr=logfile_handle)

        samtools_fastq_args = ['samtools', 'fastq', '-0', subset_fastq_filepath, '-n', '-@', str(threads)]
        run_pipeline_subcommand(command_args=samtools_fastq_args, stdin=samtools_view, stderr=logfile_handle)


def run_flye(fastq_filepath: str, flye_outdir: str, flye_read_mode: str, flye_read_error: float, log_filepath: str,
             append_log: bool = True, threads: int = 1):
    """
    Runs Flye to assemble the reads in the input FastQ file. This function allows Flye to fail without raising an error
    :param fastq_filepath: path to a FastQ file containing the input reads (gzipped is OK)
    :param flye_outdir: directory to save Flye output to
    :param flye_read_mode: type of long reads, either 'nano-raw' or 'nano-hq'
    :param flye_read_error: expected error rate of reads as a proportion; specify 0 to use default Flye settings
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :return: return code of Flye (0 if it finished successfully).
    """

    write_mode = set_write_mode(append_log)

    # TODO - add support for PacBio reads
    if (flye_read_mode != 'nano-raw') & (flye_read_mode != 'nano-hq'):

        logger.error(f'flye_read_mode must be "nano-raw" or "nano-hq"; you provided {flye_read_mode}')
        raise ValueError

    flye_args = ['flye', f'--{flye_read_mode}', fastq_filepath, '-o', flye_outdir, '-t', str(threads)]

    if flye_read_error != 0:

        flye_args.append('--read_error')
        flye_args.append(flye_read_error)

    with open(log_filepath, write_mode) as logfile_handle:

        flye_run = run_pipeline_subcommand(command_args=flye_args, check=False, stderr=logfile_handle)

    if flye_run.returncode != 0:

        logger.warning(f'Flye did not finish successfully; see log for details at "{log_filepath}"')

    return flye_run.returncode


def run_circlator_merge(circular_contig_filepath: str, patch_contig_filepath: str, merge_outdir: str,
                        circlator_min_id: float, circlator_min_length: int, circlator_ref_end: int,
                        circlator_reassemble_end: int, log_filepath: str, append_log: bool = True):
    """
    Runs the 'circlator merge' module to stitch a gap-spanning contig onto the ends of a circular contig to confirm and
    repair the circularization of the contig

    :param circular_contig_filepath: path to a FastA file containing the original (non-stitched) circular contig
    :param patch_contig_filepath: path to a FastA file containing a linear contig that should span the 'ends' of the
                                  circular contig
    :param merge_outdir: directory to save circlator merge output to
    :param circlator_min_id: Percent identity threshold for circlator merge to stitch the contigs
    :param circlator_min_length: Minimum required overlap (bp) between the circular contig and the patch contig
    :param circlator_ref_end: Minimum distance (bp) between end of circular contig and the nucmer hit
    :param circlator_reassemble_end: Minimum distance (bp) between end of patch contig and the nucmer hit
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :return: circlator merge output is saved to disk at merge_outdir
    """

    write_mode = set_write_mode(append_log)

    os.makedirs(merge_outdir, exist_ok=True)

    with open(log_filepath, write_mode) as logfile_handle:

        circlator_merge_args = ['circlator', 'merge', '--verbose', '--min_id', str(circlator_min_id),
                                '--min_length', str(circlator_min_length), '--ref_end', str(circlator_ref_end),
                                '--reassemble_end', str(circlator_reassemble_end), circular_contig_filepath,
                                patch_contig_filepath, os.path.join(merge_outdir, 'merge')]
        run_pipeline_subcommand(command_args=circlator_merge_args, stdout=logfile_handle, stderr=subprocess.STDOUT)


def check_circlator_success(circlator_logfile: str):
    """
    Checks the circlator log file to see if contig stitching was successful
    :param circlator_logfile: path to the circlator log file
    :return: Boolean of whether the contigs were successfully stitched or not
    """

    logger.debug('Loading circlator logfile')
    circlator_info = pd.read_csv(circlator_logfile, sep='\t')[['#Contig', 'circularised']]

    # TODO: I think I might be able to safely just look for if the sum is 1 (>1 might mean something odd happened)
    # If a row is '1', it means it was stitched properly, but if '0', it means it was not stitched.
    # So if all rows are 1 (i.e., the sum of rows / # of rows is 1), it means everything was stitched properly.
    if circlator_info['circularised'].sum() / circlator_info.shape[0] == 1:

        logger.debug('Everything is stitched.')
        result = True

    elif circlator_info['circularised'].sum() >= 0:

        logger.debug('Not everything is stitched.')
        result = False

    else:

        logger.error('File processing error. # of non-stitched contigs is not >=0. Exiting...')
        raise RuntimeError

    return result


def rotate_contig_to_midpoint(contig_fasta_filepath: str, output_filepath: str, append: bool = False):
    """
    Rotates an input (circular) contig to its approximate midpoint
    :param contig_fasta_filepath: path to a FastA file containing a single circular contig
    :param output_filepath: path where the FastA file containing the output rotated contig should be saved
    :param append: whether to append the output FastA onto an existing file (True) or overwrite (False)
    :return: writes FastA file to output_filepath, 60 bp per line
    """

    write_mode = set_write_mode(append)

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
    with open(output_filepath, write_mode) as output_handle:
        SeqIO.write(contig_record, output_handle, 'fasta')


def link_contig_ends(contig_record: SeqIO.SeqRecord, bam_filepath: str, length_outdir: str, length_threshold: int,
                     cli_tool_settings_dict: dict, verbose_logfile: str, threads: int = 1):
    """
    Attempt to stitch the ends of an input circular contig via assembling reads mapped within x bp of the contig ends

    :param contig_record: SeqRecord of the circular contig
    :param bam_filepath: path to a BAM file with mapping information of long reads to the contig (it is OK if this file
                         also contains mappings to other contigs outside the contig of interest)
    :param length_outdir: output directory for the analysis; will be created by the function, although it can exist
                          before the function is run
    :param length_threshold: bp region around the contig ends to subset for the assembly (half of this distance will be
                             targeted from each contig end)
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param verbose_logfile: path to a log file where shell script logs will be saved
    :param threads: parallel processor threads to use for the analysis
    :return: exit status code for Flye
    """

    # Define folder and file names
    bed_filepath = os.path.join(length_outdir, 'ends.bed')
    ends_fastq_filepath = os.path.join(length_outdir, 'ends.fastq.gz')
    flye_length_outdir = os.path.join(length_outdir, 'assembly')
    merge_dir = os.path.join(length_outdir, 'merge')

    os.makedirs(length_outdir, exist_ok=True)

    # Get reads for the selected region around the contig ends
    generate_bed_file(contig_record, bed_filepath, length_threshold=length_threshold)
    subset_reads_from_bam(bam_filepath=bam_filepath, bed_filepath=bed_filepath,
                          subset_fastq_filepath=ends_fastq_filepath, log_filepath=verbose_logfile,
                          append_log=True, threads=threads)

    # Assemble the reads to get (hopefully) a joined contig end
    flye_exit_status = run_flye(fastq_filepath=ends_fastq_filepath, flye_outdir=flye_length_outdir,
                                flye_read_mode=cli_tool_settings_dict['flye_read_mode'],
                                flye_read_error=cli_tool_settings_dict['flye_read_error'],
                                log_filepath=verbose_logfile, append_log=True, threads=threads)

    if flye_exit_status == 0:

        # Write the original contig sequence to a file
        circular_contig_filepath = os.path.join(length_outdir, 'original.fasta')

        with open(circular_contig_filepath, 'w') as output_handle:
            SeqIO.write(contig_record, output_handle, 'fasta')

        # Stitch the joined contig end onto the original assembly
        # TODO - sometimes small contigs are already rotated far from original origin. Stitch point is
        #        hard to find. Does circlator report stitch point?
        run_circlator_merge(circular_contig_filepath=circular_contig_filepath,
                            patch_contig_filepath=os.path.join(flye_length_outdir, 'assembly.fasta'),
                            merge_outdir=merge_dir, circlator_min_id=cli_tool_settings_dict['circlator_min_id'],
                            circlator_min_length=cli_tool_settings_dict['circlator_min_length'],
                            circlator_ref_end=cli_tool_settings_dict['circlator_ref_end'],
                            circlator_reassemble_end=cli_tool_settings_dict['circlator_reassemble_end'],
                            log_filepath=verbose_logfile, append_log=True)

    else:
        logger.warning('Flye assembly FAILED')

    return flye_exit_status


def iterate_linking_contig_ends(contig_record: SeqIO.SeqRecord, bam_filepath: str, linking_outdir: str,
                                length_thresholds: list, cli_tool_settings_dict: dict,  verbose_logfile: str,
                                threads: int = 1):
    """
    Iterate link_contig_ends to try to stitch the ends of a circular contig using multiple length thresholds

    :param contig_record: SeqRecord of the circular contig
    :param bam_filepath: path to a BAM file with mapping information of long reads to the contig (it is OK if this file
                         also contains mappings to other contigs outside the contig of interest)
    :param linking_outdir: output directory for the analysis; will be created by the function, although it can exist
                           before the function is run
    :param length_thresholds: list of bp regions around the contig ends to attempt to subset for the assembly
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param verbose_logfile: path to a logfile where shell script logs will be saved
    :param threads: parallel processor threads to use for the analysis
    :return: boolean of whether end linkage was successful (True) or not (False)
    """

    os.makedirs(linking_outdir, exist_ok=True)

    # Where important log files will be copied
    log_dir_base = os.path.join(linking_outdir, 'logs')

    # Keep trying to link the contig ends until successful or until all length thresholds have been attempted
    assembly_attempts = 0
    linked_ends = False

    while linked_ends is False:

        # Exit the function if all length thresholds have been tried
        if assembly_attempts >= len(length_thresholds):
            break

        # Get the length threshold for this attempt
        length_threshold = length_thresholds[assembly_attempts]

        # Don't attempt to link the ends if the contig is too short
        if len(contig_record.seq) <= length_threshold:
            logger.info(f'Skipping length threshold of {length_threshold} because '
                        f'contig is shorter than this length ({len(contig_record.seq)} bp)')

            assembly_attempts = assembly_attempts + 1
            continue

        # Make sure that the circlator_min_length is shorter than the length threshold (otherwise merging is impossible)
        if cli_tool_settings_dict['circlator_min_length'] > length_threshold:

            revised_min_length = int(length_threshold * 0.9)
            cli_tool_settings_dict_revised = cli_tool_settings_dict.copy()
            cli_tool_settings_dict_revised['circlator_min_length'] = revised_min_length

            logger.warning(f'The minimum required length of alignment between the original and reassembled contig '
                           f'(specified by circlator_min_length; {cli_tool_settings_dict["circlator_min_length"]} bp) '
                           f'is longer than the length_threshold the original contig will be subset to '
                           f'({length_threshold} bp). This means that finding a matching merge will be impossible. To '
                           f'overcome this, the script will shorter the circlator_min_length for this iteration to 90% '
                           f'of the length threshold, i.e., {revised_min_length} bp.')

            cli_tool_settings_dict_used = cli_tool_settings_dict_revised

        else:
            cli_tool_settings_dict_used = cli_tool_settings_dict

        # Try to stitch the contig ends
        logger.info(f'Attempting reassembly with a length threshold of {length_threshold} bp')
        length_outdir = os.path.join(linking_outdir, 'tmp', f'L{length_threshold}')

        flye_exit_status = link_contig_ends(contig_record=contig_record, bam_filepath=bam_filepath,
                                            length_outdir=length_outdir, length_threshold=length_threshold,
                                            cli_tool_settings_dict=cli_tool_settings_dict_used,
                                            verbose_logfile=verbose_logfile, threads=threads)

        # Make the output logging directory (this is needed even if Flye fails so the pipeline can keep going)
        log_dir = os.path.join(log_dir_base, f'L{length_threshold}')
        os.makedirs(log_dir, exist_ok=True)

        # Stop early if Flye fails
        if flye_exit_status != 0:
            assembly_attempts = assembly_attempts + 1
            continue

        # If Flye was successful, copy important log files for that specific length threshold
        shutil.copy(os.path.join(length_outdir, 'assembly', 'assembly_info.txt'), log_dir)
        shutil.copy(os.path.join(length_outdir, 'merge', 'merge.circularise.log'), log_dir)
        shutil.copy(os.path.join(length_outdir, 'merge', 'merge.circularise_details.log'), log_dir)

        # See if stitching the contigs ends worked
        if check_circlator_success(os.path.join(length_outdir, 'merge', 'merge.circularise.log')) is True:
            logger.info('Successfully linked contig ends')

            # Rotate to midpoint so that the stitched points can be polished more effectively downstream
            rotate_contig_to_midpoint(os.path.join(length_outdir, 'merge', 'merge.fasta'),
                                      os.path.join(linking_outdir, 'stitched.fasta'), append=False)

            # Save a copy of the final circlator merge logfile in the main log directory
            shutil.copy(os.path.join(length_outdir, 'merge', 'merge.circularise_details.log'),
                        os.path.join(log_dir_base, f'{contig_record.name}_circlator_final.log'))

            linked_ends = True

        assembly_attempts = assembly_attempts + 1

    shutil.rmtree(os.path.join(linking_outdir, 'tmp'))

    return linked_ends


def stitch_all_contigs(circular_contig_tmp_fasta: str, bam_filepath: str, linking_outdir_base: str,
                       end_repaired_contigs_filepath: str, length_thresholds: list, cli_tool_settings_dict: dict,
                       verbose_logfile: str, threads: int):
    """
    Run the iterate_linking_contig_ends function on all contigs in an input FastA file, i.e., attempt to stitch the ends
    of all the contigs (assumed circular) in the file.

    :param circular_contig_tmp_fasta: Input FastA file of circular contigs to be stitched
    :param bam_filepath: Path to a BAM file with mapping information of long reads to the contigs (it is OK if this file
                         also contains mappings to other contigs outside those in the input file)
    :param end_repaired_contigs_filepath: Path (to be overwritten by this funtion) to save the end repaired contigs
    :param linking_outdir_base: Temp directory to save re-assembly and stitching files to
    :param length_thresholds: list of bp regions around the contig ends to attempt to subset for the assembly
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param verbose_logfile: path to a logfile where shell script logs will be added
    :param threads: parallel processor threads to use for the analysis
    :return: Writes stitched contigs to end_repaired_contigs_filepath. Returns a list of the names of any contigs that
             could not be stitched successfully (list length will be zero if all contigs stitched successfully)
    """

    # Initialize the repaired contigs FastA file (so it will overwrite an old file rather than just append later)
    with open(end_repaired_contigs_filepath, 'w') as output_handle:
        output_handle.write('')

    # Make the tmp directory for output files and the tmp directory where final log files (if any) will be moved
    os.makedirs(os.path.join(linking_outdir_base, 'log_summary'), exist_ok=True)

    failed_contig_names = []

    with open(circular_contig_tmp_fasta) as fasta_handle:

        for contig_record in SeqIO.parse(fasta_handle, 'fasta'):

            logger.info(f'Starting end repair on contig: {contig_record.name}')

            # This is the temp folder that will be used for this contig during stitching
            linking_outdir = os.path.join(linking_outdir_base, contig_record.name)

            end_linkage_complete = iterate_linking_contig_ends(contig_record, bam_filepath, linking_outdir,
                                                               length_thresholds, cli_tool_settings_dict,
                                                               verbose_logfile, threads)

            if end_linkage_complete is False:

                logger.warning(f'Contig {contig_record.name}: FAILED to linked contig ends')
                os.makedirs(os.path.join(linking_outdir_base, 'troubleshooting'), exist_ok=True)
                shutil.move(os.path.join(linking_outdir, 'logs'),
                            os.path.join(linking_outdir_base, 'troubleshooting', contig_record.name))

                failed_contig_names.append(contig_record.name)

            elif end_linkage_complete is True:

                # Append the successful contig onto the main file
                with open(os.path.join(linking_outdir, 'stitched.fasta')) as input_handle:
                    with open(end_repaired_contigs_filepath, 'a') as append_handle:
                        append_handle.write(input_handle.read())

                shutil.move(os.path.join(linking_outdir, 'logs', f'{contig_record.name}_circlator_final.log'),
                            os.path.join(linking_outdir_base, 'log_summary', f'{contig_record.name}.log'))

                shutil.rmtree(linking_outdir)

            else:

                logger.error(f'end_linkage_complete should be boolean True or False; is actually '
                             f'{end_linkage_complete}')
                raise ValueError

    return failed_contig_names


def run_end_repair(long_read_filepath: str, assembly_fasta_filepath: str, assembly_info_filepath: str,
                   assembly_info_type: str, output_dir: str, length_thresholds: list, keep_failed_contigs: bool,
                   cli_tool_settings_dict: dict, threads: int, threads_mem_mb: int):
    """
    Runs the end repair workflow

    :param long_read_filepath: path to the QC-passing Nanopore reads, FastQ format; gzipped is OK
    :param assembly_fasta_filepath: path to the assembled contigs output by Flye, FastA format
    :param assembly_info_filepath: path to assembly_info.txt output by Flye
    :param assembly_info_type: 'flye' type for assembly_info.txt from flye; 'custom' for custom format (see parse_cli
                               for custom format details)
    :param output_dir: path to the directory to save output files
    :param length_thresholds: list of bp regions around the contig ends to attempt to subset for the assembly
    :param keep_failed_contigs: boolean that defines whether to 1) continue the code even if some contigs cannot be end
                                 repaired (True) vs. to 2) exit with an error code if some contigs cannot be end
                                 repaired (False)
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param threads: parallel processor threads to use for the analysis
    :param threads_mem_mb: memory in MB per thread (to use for samtools); must be an integer
    :return: None
    """

    # Define core file names and directory structures
    # These will be the files and folders in the main output directory:
    end_repaired_contigs_filepath = os.path.join(output_dir, 'repaired.fasta')
    end_repair_status_filepath = os.path.join(output_dir, 'repaired_info.tsv')
    verbose_logfile = os.path.join(output_dir, 'verbose.log')
    bam_filepath = os.path.join(output_dir, 'long_read.bam')
    circlator_logdir = os.path.join(output_dir, 'circlator_logs')
    linking_outdir_base = os.path.join(output_dir, 'contigs')
    circular_contig_tmp_fasta = os.path.join(linking_outdir_base, 'circular_input.fasta')

    # Get lists of circular contigs
    circular_contig_names = parse_assembly_info_file(assembly_info_filepath, assembly_info_type, return_type='circular')

    # No need to run the pipeline if there are no circular contigs
    if len(circular_contig_names) == 0:
        logger.info('No circular contigs. Will copy the input file, repaired_info.tsv, and finish early.')
        shutil.copyfile(assembly_fasta_filepath, end_repaired_contigs_filepath)
        linear_contig_names = parse_assembly_info_file(assembly_info_filepath, assembly_info_type, return_type='linear')
        repair_info = pd.DataFrame({'contig': linear_contig_names, 'circular': 'N', 'repaired': 'N'})
        repair_info.to_csv(end_repair_status_filepath, sep='\t', index=False)
        logger.info('Pipeline finished.')
        sys.exit(0)

    # Subset circular contigs from the full file and save to disk
    os.makedirs(linking_outdir_base, exist_ok=True)
    with open(circular_contig_tmp_fasta, 'w') as output_handle:
        for record in subset_sequences(assembly_fasta_filepath, circular_contig_names):
            SeqIO.write(record, output_handle, 'fasta')

    # Start the main workflow from here:
    logger.info('Mapping reads to all contigs')
    map_long_reads(contig_filepath=assembly_fasta_filepath, long_read_filepath=long_read_filepath,
                   output_bam_filepath=bam_filepath, log_filepath=verbose_logfile,  append_log=False, threads=threads,
                   threads_mem_mb=threads_mem_mb)

    failed_contig_names = stitch_all_contigs(circular_contig_tmp_fasta, bam_filepath, linking_outdir_base,
                                             end_repaired_contigs_filepath, length_thresholds, cli_tool_settings_dict,
                                             verbose_logfile, threads)

    os.makedirs(circlator_logdir, exist_ok=True)
    shutil.move(os.path.join(linking_outdir_base, 'log_summary'), circlator_logdir)

    # Check for contigs that could not be circularized
    if len(failed_contig_names) != 0:

        if os.path.isdir(os.path.join(linking_outdir_base, 'troubleshooting')):
            shutil.move(os.path.join(linking_outdir_base, 'troubleshooting'),
                        os.path.join(output_dir, 'troubleshooting'))

        if keep_failed_contigs is False:

            logger.error(f'{len(failed_contig_names)} contigs could not be circularized. A partial output file '
                         f'including successfully circularized contigs (and no linear contigs) is available at '
                         f'{end_repaired_contigs_filepath} for debugging. Exiting with error status. See temporary '
                         f'files and verbose logs for more details.')
            sys.exit(1)

        elif keep_failed_contigs is True:

            logger.warning(f'{len(failed_contig_names)} contigs could not be circularized. The original (non-repaired) '
                           f'versions of these contigs will be included in the final output file')
            logger.warning(f'Names of contigs that could not be circularized: {", ".join(failed_contig_names)}')

            # Get the non-repaired circular contigs and append them to the repaired contigs file
            with open(end_repaired_contigs_filepath, 'a') as append_handle:

                for record in subset_sequences(assembly_fasta_filepath, failed_contig_names):
                    SeqIO.write(record, append_handle, 'fasta')

        else:

            logger.error(f'keep_failed_contigs should be a boolean True or False; you provided {keep_failed_contigs}')
            raise ValueError

    # Get the linear contigs and append them to the repaired contigs file
    # TODO - consider just getting all non-circular contigs regardless of whether they are in the assembly info file.
    #   This edit might require modification of the subset_sequences function and refactoring of code in several places,
    #   but I imagine it is closer to the behaviour that the user would expect.
    linear_contig_names = parse_assembly_info_file(assembly_info_filepath, assembly_info_type, return_type='linear')

    with open(end_repaired_contigs_filepath, 'a') as append_handle:

        for record in subset_sequences(assembly_fasta_filepath, linear_contig_names):
            SeqIO.write(record, append_handle, 'fasta')

    # Write info file of how contigs were repaired
    repaired_contig_names = list(set(circular_contig_names).difference(set(failed_contig_names)))
    repair_info = pd.concat([pd.DataFrame({'contig': repaired_contig_names, 'circular': 'Y', 'repaired': 'Y'}),
                             pd.DataFrame({'contig': failed_contig_names, 'circular': 'Y', 'repaired': 'N'}),
                             pd.DataFrame({'contig': linear_contig_names, 'circular': 'N', 'repaired': 'N'})], axis=0)
    repair_info.to_csv(end_repair_status_filepath, sep='\t', index=False)

    # Clean up temp files
    os.remove(bam_filepath)
    os.remove(f'{bam_filepath}.bai')
    shutil.rmtree(linking_outdir_base)

    logger.info('End repair finished!')
    logger.info('#### Final stats: ####')
    logger.info(f'Circular, repaired: {len(repaired_contig_names)} contig(s)')
    logger.info(f'Circular, not repairable: {len(failed_contig_names)} contig(s)')
    logger.info(f'Linear: {len(linear_contig_names)} contig(s)')
    logger.info('######################')
    logger.info(f'Output contigs are saved at {end_repaired_contigs_filepath}. '
                f'A summary of repair work is saved at {end_repair_status_filepath}.')


def parse_cli():
    """
    Parses the CLI arguments
    :return: An argparse parser object
    """

    parser = argparse.ArgumentParser(
        description=f'{os.path.basename(sys.argv[0])}: pipeline to repair ends of circular contigs from Flye. \n'
                    '  Copyright Jackson M. Tsuji and Lee Bergstrand, 2023',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    required_settings = parser.add_argument_group('Required')
    read_settings = parser.add_argument_group('Input read options')
    merge_settings = parser.add_argument_group('Merge options')
    workflow_settings = parser.add_argument_group('Workflow options')

    required_settings.add_argument('-l', '--long_read_filepath', required=True, type=str,
                                   help='QC-passing Nanopore reads')
    required_settings.add_argument('-a', '--assembly_fasta_filepath', required=True, type=str,
                                   help='Contigs to be end-repaired')
    required_settings.add_argument('-i', '--assembly_info_filepath', required=True, type=str,
                                   help='assembly_info.txt file from Flye showing which assembled contigs are circular '
                                        'vs. linear (or a custom guide file; see -c below)')
    required_settings.add_argument('-o', '--output_dir', required=True, type=str, help='Output directory path')

    read_settings.add_argument('-f', '--flye_read_mode', required=False, default='nano-hq',
                               choices=['nano-hq', 'nano-raw'],
                               help='Type of long read reads provided by -l, to be used for reassembly by Flye. See '
                                    'details on these settings in the Flye documentation. (default: nano-hq)')
    read_settings.add_argument('-F', '--flye_read_error', required=False, default=0, type=float,
                               help='Expected error rate of input reads, expressed as proportion (e.g., 0.03). '
                                    'If "0", then have flye set the read error automatically (default: 0)')

    merge_settings.add_argument('-I', '--circlator_min_id', required=False, default=99, type=float,
                                help='Percent identity threshold for circlator merge (default: 99)')
    merge_settings.add_argument('-L', '--circlator_min_length', required=False, default=10000, type=int,
                                help='Minimum required overlap (bp) between original and merge contigs '
                                     '(default: 10000)')
    merge_settings.add_argument('-e', '--circlator_ref_end', required=False, default=100, type=int,
                                help='Minimum distance (bp) between end of original contig and nucmer hit '
                                     '(default: 100)')
    merge_settings.add_argument('-E', '--circlator_reassemble_end', required=False, default=100, type=int,
                                help='Minimum distance (bp) between end of merge contig and nucmer hit '
                                     '(default: 100)')

    workflow_settings.add_argument('-k', '--keep_going_with_failed_contigs', required=False, action='store_true',
                                   help='Set this flag to continue running this script even if some contigs '
                                        'cannot be circularized and end-repaired. For any non-repaired contigs, an '
                                        'exact copy of the original contig will be output.')
    workflow_settings.add_argument('-c', '--custom_assembly_info_file', required=False, action='store_true',
                                   help='Set this flag if you provide a custom tab-separated file to specify the '
                                        'circular vs. linear status of each contig in the assembly, rather than using '
                                        'the assembly_info.txt file output by Flye (default), as your argument to -i. '
                                        'The custom file must have the following format: no headers; tab-separated; '
                                        'first column is the contig names; second column is the status of the contigs, '
                                        'either "circular" or "linear". Any contig names not in this file will be '
                                        'dropped by the script!')
    workflow_settings.add_argument('-O', '--overwrite', required=False, action='store_true',
                                   help='Set this flag to overwrite an existing output directory and its contents; by '
                                        'default this code will throw an error if the output directory already exists. '
                                        'By setting this flag, you risk erasing old data and/or running into errors.')
    workflow_settings.add_argument('-T', '--length_thresholds', required=False,
                                   default='100000,75000,50000,25000,5000,2500,1000', type=str,
                                   help='Comma-separated list of length thresholds for reassembly around the contig '
                                        'ends (bp) (default: 100000,75000,50000,25000,5000,2500,1000)')
    workflow_settings.add_argument('-t', '--threads', required=False, default=1, type=int,
                                   help='Number of processors threads to use (default: 1)')
    workflow_settings.add_argument('-m', '--threads_mem', required=False, default=1, type=float,
                                   help='Memory (GB) to use **per thread** for samtools sort (default: 1)')
    workflow_settings.add_argument('-v', '--verbose', required=False, action='store_true',
                                   help='Enable verbose logging')

    return parser


if __name__ == '__main__':
    main()
