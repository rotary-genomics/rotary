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
import subprocess
from Bio import SeqIO
import pandas as pd

# GLOBAL VARIABLES
SCRIPT_VERSION = '0.2.0'

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

    return dependency_message


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


def sort_fasta(input_multi_fasta, output_fasta, sort_order, low_memory=True):
    """
    Sorts entries in a multi-FastA file in the specified sequence order

    :param input_multi_fasta: filepath to a FastA file (1 or more entries)
    :param output_fasta: path to write the output sorted FastA file
    :param sort_order: list of contig names for sorting.
                       If sort_order is missing entries in input_multi_fasta, then those entries are not returned.
                       If sort_order has entries not in input_multi_fasta, then those entries are not returned.
    :param low_memory: whether to run in low memory mode (uses more disk).
                       In high memory mode, will load the whole FastA file into RAM
    :return: writes the output FastA file to disk
    """

    if low_memory is False:

        contig_names = []
        contig_records = []

        with open(input_multi_fasta) as fasta_handle:
            for record in SeqIO.parse(fasta_handle, 'fasta'):

                contig_names.append(record.name)
                contig_records.append(record)

        # Load the names and positions of the sequences as a DataFrame
        sequence_positions = pd.DataFrame({'contig_name': contig_names}).reset_index()

        # Sort/subset the DataFrame based on the provided sort order list
        sequence_positions_sorted = sequence_positions\
            .set_index('contig_name')\
            .reindex(sort_order)\
            .reset_index()
        sequence_positions_sorted = sequence_positions_sorted[
            sequence_positions_sorted['index'].isna() == False]

        # Write output in the desired order
        with open(output_fasta, 'w') as output_handle:
            for record_index in sequence_positions_sorted['index'].apply(int):

                logger.debug(f'Writing record index {record_index}')
                record = contig_records[record_index]

                SeqIO.write(record, output_handle, 'fasta')

    elif low_memory is True:

        # First read only the contig names (not sequences)
        contig_names = []

        with open(input_multi_fasta) as fasta_handle:
            for record in SeqIO.parse(fasta_handle, 'fasta'):

                contig_names.append(record.name)

        # Load the names and positions of the sequences as a DataFrame
        sequence_positions = pd.DataFrame({'contig_name': contig_names}).reset_index()

        # Sort/subset the DataFrame based on the provided sort order list
        sequence_positions_sorted = sequence_positions \
            .set_index('contig_name') \
            .reindex(sort_order) \
            .reset_index()
        sequence_positions_sorted = sequence_positions_sorted[
            sequence_positions_sorted['index'].isna() == False]

        # Write output in the desired order, reading each sequence on demand
        with open(output_fasta, 'w') as output_handle:
            for record_index in sequence_positions_sorted['index'].apply(int):

                with open(input_multi_fasta) as fasta_handle:

                    record_count = 0
                    for record in SeqIO.parse(fasta_handle, 'fasta'):

                        if record_count == record_index:
                            logger.debug(f'Writing record index {record_index}')
                            SeqIO.write(record, output_handle, 'fasta')

                        record_count = record_count + 1

    else:

        logger.error(f'low_memory mode should be True or False; you provided {low_memory}')
        raise RuntimeError


def parse_assembly_info_file(assembly_info_filepath, return_type='circular', linear_contig_list_filepath=None,
                             circular_contig_list_filepath=None):
    """
    List circular and linear contigs from the Flye assembly info file

    :param assembly_info_filepath: path to assembly_info.txt output by Flye
    :param return_type: whether to return a list of 'circular' or 'linear' contigs
    :param linear_contig_list_filepath: output filepath to list the names of linear contigs
    :param circular_contig_list_filepath: output filepath to list the names of circular contigs
    :return: list of contig names, either circular or linear depending on return_type
    """

    logger.debug('Loading assembly info file')
    assembly_info = pd.read_csv(assembly_info_filepath, sep='\t')

    circular_contigs = assembly_info[assembly_info['circ.'] == 'Y'][['#seq_name', 'length', 'circ.']]
    linear_contigs = assembly_info[assembly_info['circ.'] == 'N'][['#seq_name', 'length', 'circ.']]

    if linear_contig_list_filepath is not None:
        logger.debug('Exporting linear contig list')
        linear_contigs['#seq_name'].to_csv(linear_contig_list_filepath, index=False, header=None)

    if circular_contig_list_filepath is not None:
        logger.debug('Exporting circular contig list')
        circular_contigs['#seq_name'].to_csv(circular_contig_list_filepath, index=False, header=None)

    if return_type == 'circular':
        output_list = list(circular_contigs['#seq_name'])
    elif return_type == 'linear':
        output_list = list(linear_contigs['#seq_name'])
    else:
        logger.error(f'return_type must be "circular" or "linear"; you provided "{return_type}"')
        raise RuntimeError

    return output_list


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


def check_circlator_logfile(circlator_logfile):
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


def run_end_repair(long_read_filepath, assembly_fasta_filepath, assembly_info_filepath, output_dir,
                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                   circlator_ref_end, circlator_reassemble_end, threads, thread_mem):

    # Variable names
    end_repaired_contigs_filepath = os.path.join(output_dir, 'repaired.fasta')
    end_repaired_contigs_filepath_unsorted = f'{end_repaired_contigs_filepath}.tmp'  # TODO - initialize this file
    verbose_logfile = os.path.join(output_dir, 'verbose.log')
    bam_filepath = os.path.join(output_dir, 'long_read.bam')
    circlator_logdir = os.path.join(output_dir, 'circlator_logs')
    linear_contig_list_filepath = os.path.join(output_dir, 'linear.list')
    reassembly_outdir_base = os.path.join(output_dir, 'contigs')

    # Check output dir
    if os.path.isdir(output_dir):
        logger.warning(f'Output directory already exists: "{output_dir}"; files may be overwritten.')

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(circlator_logdir, exist_ok=True)

    # Get lists of circular contigs in RAM while saving linear contig list to file
    circular_contig_names = parse_assembly_info_file(assembly_info_filepath, return_type='circular',
                                                     linear_contig_list_filepath=linear_contig_list_filepath)

    # No need to run the pipeline if there are no circular contigs
    if len(circular_contig_names) == 0:

        logger.info('No circular contigs. Will copy the input file and finish early.')
        shutil.copyfile(assembly_fasta_filepath, end_repaired_contigs_filepath)
        os.remove(linear_contig_list_filepath)
        logger.info('Pipeline finished.')
        sys.exit(0)

    logger.info('Mapping reads to all contigs')
    # TODO - fix this command; some ideas at https://stackoverflow.com/a/13332300 (accessed 2023.8.29)
    subprocess.run(['minimap2', '-t', threads, '-ax', 'map-ont', assembly_fasta_filepath, long_read_filepath,
                    '2>>', verbose_logfile, '|',
                    'samtools', 'view', '-b', '-@', threads, '2>>', verbose_logfile, '|',
                    'samtools', 'sort', '-@', threads, '-m', str(thread_mem) + 'G', '2>>', verbose_logfile,
                    '>', bam_filepath])
    # # TODO - add support for different flags like -ax for pacbio
    #   minimap2 -t "${threads}" -ax map-ont "${all_contigs}" "${qc_long}" 2>> "${verbose}" | \
    #     samtools view -b -@ "${threads}" 2>> "${verbose}" | \
    #     samtools sort -@ "${threads}" -m "${thread_mem}G" 2>> "${verbose}" \
    #     > "${bam_file}"
    # # TODO - add this command
    # samtools index -@ "${threads}" "${bam_file}" 2>> "${verbose}"

    failed_contigs = 0

    for contig_name in circular_contig_names:

        logger.info(f'End repair: {contig_name}')

        # Output dir setup
        reassembly_outdir = os.path.join(reassembly_outdir_base, contig_name)
        reassembly_outfile = os.path.join(reassembly_outdir, f'{contig_name}.fasta')
        os.makedirs(reassembly_outdir, exist_ok=True)

        # Get contig sequence as a SeqRecord
        with open(assembly_fasta_filepath) as fasta_handle:
            for record in SeqIO.parse(fasta_handle):
                if record.name == contig_name:
                    contig_record = record

        # TODO - allow user to customize these lengths
        length_cutoffs = [100000, 90000, 80000, 70000, 60000, 50000, 20000, 10000, 5000, 2000]
        assembly_attempts = 0
        linked_ends = False

        while linked_ends is False:

            # End the script if all length cutoffs have been tried
            if assembly_attempts > len(length_cutoffs)-1:

                logger.warning(f'End repair: {contig_name}: FAILED to linked contig ends')
                failed_contigs = failed_contigs + 1
                os.makedirs(os.path.join(output_dir, 'troubleshooting'))
                shutil.move(log_dir, os.path.join(output_dir, 'troubleshooting', contig_name))
                shutil.rmtree(reassembly_outdir)

                break

            # Get the length cutoff for this attempt
            length_cutoff = length_cutoffs[assembly_attempts]

            if len(contig_record.seq) <= length_cutoff:
                logger.debug(f'Skipping length cutoff of {length_cutoff} because '
                             f'contig is shorter than this ({len(contig_record.seq)} bp)')
                continue

            logger.debug(f'Starting reassembly with a length cutoff of {length_cutoff}')

            # Define some folder and filenames for this attempt
            length_outdir = os.path.join(reassembly_outdir, f'L{length_cutoff}')
            bed_filepath = os.path.join(length_outdir, 'ends.bed')
            flye_length_outdir = os.path.join(length_outdir, 'assembly')
            merge_dir = os.path.join(length_outdir, 'merge')
            log_dir = os.path.join(reassembly_outdir, 'logs', f'L{length_cutoff}')
            os.makedirs(length_outdir, exist_ok=True)
            os.makedirs(log_dir, exist_ok=True)

            # Make a BED file
            generate_bed_file(contig_record, bed_filepath, length_threshold=length_cutoff)

            # TODO - finish code for getting reads for that region
            # # TODO - consider adding option to split long reads in half if they go around a short circular contig, like in circlator
            #       samtools view -@ "${threads}" -L "${regions}" -b "${bam_file}" 2>> "${verbose}" | \
            #         samtools fastq -0 "${fastq}" -n -@ "${threads}" 2>> "${verbose}"

            if flye_read_error == 0:

                # TODO - finish command
                # flye "--${flye_read_mode}" "${fastq}" -o "${flye_length_outdir}" \
                #           -t "${threads}" >> "${verbose}" 2>&1

            else:

                # TODO - finish command
                # flye "--${flye_read_mode}" "${fastq}" -o "${flye_length_outdir}" --read_error "${flye_read_error}" \
                #           -t "${threads}" >> "${verbose}" 2>&1

            shutil.copy(os.path.join(flye_length_outdir, 'assembly_info.txt'), log_dir)

            os.makedirs(merge_dir, exist_ok=True)

            # TODO - finish command
            # circlator merge --verbose --min_id "${circlator_min_id}" --min_length "${circlator_min_length}" \
            #         --ref_end "${circlator_ref_end}" --reassemble_end "${circlator_reassemble_end}" \
            #         "${assembly}" "${reassembly_dir}/assembly.fasta" "${merge_dir}/merge" >> "${verbose}" 2>&1
            shutil.copy(os.path.join(merge_dir, 'merge.circularise.log'), log_dir)
            shutil.copy(os.path.join(merge_dir, 'merge.circularise_details.log'), log_dir)

            if check_circlator_logfile(os.path.join(merge_dir, 'merge.circularise.log')) is True:

                logger.info('End repair: successfully linked contig ends')

                # Rotate to midpoint so that the stitched ends can be polished more effectively (especially stich points)
                # TODO - sometimes small contigs are already rotated far from original origin. Stitch point hards to find. Does circlator report stitch point?
                rotate_contig_to_midpoint(os.path.join(merge_dir, 'merge.fasta'), reassembly_outfile)

                # TODO - fix this code and also make sure to initialize the unsorted outfile first - consider if this can all be done within Python
                subprocess.run(['cat', reassembly_outfile, '>>', end_repaired_contigs_filepath_unsorted])

                # Cleanup
                # TODO - clean up path "${outdir}/circlator_logs/${contig}.log"
                shutil.copy(os.path.join(merge_dir, 'merge.circularise_details.log'),
                            os.path.join(circlator_logdir, f'{contig_name}.log'))
                shutil.rmtree(reassembly_outdir)

                linked_ends = True

            assembly_attempts = assembly_attempts + 1

    if failed_contigs > 0:

        # TODO - allow user to optionally keep going
        shutil.move(end_repaired_contigs_filepath_unsorted, end_repaired_contigs_filepath)

        logger.error(f'{failed_contigs} contigs could not be circularized. A partial output file including '
                     f'successfully circularized contigs (and no linear contigs) is available at '
                     f'{end_repaired_contigs_filepath} for debugging. Exiting with error status. See temporary files '
                     f'and verbose logs for more details.')
        sys.exit(1)

    # TODO - the final contigs will not be sorted
    # with open(end_repaired_contigs_filepath_unsorted) as circular_handle:
    #     with open(assembly_fasta_filepath) as linear_handle:
    #         with open(end_repaired_contigs_filepath, 'w') as output_handle:
    #             shutil.copyfileobj(circular_handle, output_handle)

    shutil.move(end_repaired_contigs_filepath_unsorted, end_repaired_contigs_filepath)
    # TODO - fix this code
    subprocess.run(['seqtk','subseq',assembly_fasta_filepath,linear_contig_list_filepath,'>>',end_repaired_contigs_filepath])

    # Clean up temp files
    os.remove(bam_filepath)
    os.remove(f'{bam_filepath}.bai')
    os.remove(linear_contig_list_filepath)
    os.remove(end_repaired_contigs_filepath_unsorted)
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
    logger.info(f'Circlator min. ID: {circlator_min_id}')
    logger.info(f'Circlator min. length: {circlator_min_length}')
    logger.info(f'Circlator ref. end: {circlator_ref_end}')
    logger.info(f'Circlator reassembly end: {circlator_reassemble_end}')
    logger.info(f'Threads: {threads}')
    logger.info(f'Memory per thread (GB): {thread_mem}')
    logger.info(f'Verbose logging: {verbose}')
    logger.info('################')

    run_end_repair(long_read_filepath, assembly_fasta_filepath, assembly_info_filepath, output_dir,
                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
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

