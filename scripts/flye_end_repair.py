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
SCRIPT_NAME = sys.argv[0]

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


def parse_assembly_info_file(assembly_info_filepath, linear_contig_list_filepath=None,
                             circular_contig_list_filepath=None):
    """
    List circular and linear contigs from the Flye assembly info file

    :param assembly_info_filepath: path to assembly_info.txt output by Flye
    :param linear_contig_list_filepath: output filepath to list the names of linear contigs
    :param circular_contig_list_filepath: output filepath to list the names of circular contigs
    :return: list of circular contig names
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

    return list(circular_contigs['#seq_name'])


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

    # Check output dir
    if os.path.isdir(output_dir):
        logger.warning(f'Output directory already exists: "{output_dir}"; files may be overwritten.')

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(circlator_logdir, exist_ok=True)

    # Get lists of circular and linear contigs
    circular_contig_names = parse_assembly_info_file(assembly_info_filepath)

    # No need to run the pipeline if there are no circular contigs
    if len(circular_contig_names) == 0:

        logger.info('No circular contigs. Will copy the input file and finish early.')
        shutil.copyfile(assembly_fasta_filepath, end_repaired_contigs_filepath)
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
        reassembly_outdir = os.path.join(output_dir, 'contigs', contig_name)
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
                rotate_contig_to_midpoint(os.path.join(merge_dir, 'merge.fasta'),
                                          os.path.join(length_outdir, 'rotate.fasta'))

                # TODO - fix this code and also make sure to initialize the unsorted outfile first - consider if this can all be done within Python
                subprocess.run(['cat', os.path.join(length_outdir, 'rotate.fasta'), '>>', end_repaired_contigs_filepath_unsorted])

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

    # TODO - stopped here - add final FastA sort code


def main(args):
    # Set user variables
    assembly_info_filepath = args.assembly_info
    linear_list_outfile = args.output_linear_contig_list
    circular_list_outfile = args.output_circular_contig_list
    end_regions_outfile = args.output_bed_file
    length_threshold = args.length_threshold
    contig_id = args.contig_id
    circlator_logfile = args.circlator_merge_log
    verbose = args.verbose

    # Startup checks
    if verbose is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if (assembly_info_filepath is not None) & (circlator_logfile is not None):
        logger.error('Specified --assembly_info and --circlator_merge_log, ' +
                     'but both options cannot be set at the same time. Exiting...')
        sys.exit(1)

    if (assembly_info_filepath is not None) & (linear_list_outfile is None) & (circular_list_outfile is None) & \
            (end_regions_outfile is None):
        logger.error('You specified --assembly_info, but you did not specify an output file ' +
                     '(output_circular_contig_list and/or output_linear_contig_list and/or output_bed_file). ' +
                     'Exiting...')
        sys.exit(1)

    # Startup messages
    logger.debug('Running ' + os.path.basename(sys.argv[0]))
    logger.debug('Version: ' + SCRIPT_VERSION)
    logger.debug('### SETTINGS ###')
    logger.debug('Assembly info filepath: ' + str(assembly_info_filepath))
    logger.debug('Output linear contig list filepath: ' + str(linear_list_outfile))
    logger.debug('Output circular contig list filepath: ' + str(circular_list_outfile))
    logger.debug('Output BED file of circular contig end regions: ' + str(end_regions_outfile))
    logger.debug('End region threshold (bp): ' + str(length_threshold))
    logger.debug('Contig ID for subsetting: ' + str(contig_id))
    logger.debug('Logfile from circlator merge: ' + str(circlator_logfile))
    logger.debug('Verbose logging: ' + str(verbose))
    logger.debug('################')

    if assembly_info_filepath is not None:
        summarize_assembly_info(assembly_info_filepath, linear_contig_list_filepath=linear_list_outfile,
                                circular_contig_list_filepath=circular_list_outfile,
                                bed_filepath=end_regions_outfile, length_threshold=length_threshold,
                                contig_id=contig_id)

    elif circlator_logfile is not None:
        check_circlator_logfile(circlator_logfile)

    else:
        logger.error('Neither --assembly_info nor --circlator_merge_log was set, ' +
                     'but choosing one of these options is required. Exiting...')
        logger.error('Were you looking for the help statement? If so, use the -h flag.')
        sys.exit(1)

    logger.debug(os.path.basename(sys.argv[0]) + ': done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Utilities for the flye_end_repair.sh pipeline. '
                    'Copyright Jackson M. Tsuji, ILTS, Hokkaido University, 2022. '
                    'Version: ' + SCRIPT_VERSION)
    parser.add_argument('-i', '--assembly_info', required=False, default=None,
                        help='The path to the assembly_info.txt file output by Flye.')
    parser.add_argument('-j', '--output_linear_contig_list', required=False, default=None,
                        help='Outputs a list of linear contigs to the designated filepath based on -i.')
    parser.add_argument('-k', '--output_circular_contig_list', required=False, default=None,
                        help='Outputs a list of circular contigs to the designated filepath based on -i.')
    parser.add_argument('-b', '--output_bed_file', required=False, default=None,
                        help='Outputs a BED file to the designated filepath specifying the ends of circular contigs. '
                             'Contig information is determined from -i. Proximity to ends is determined from -l.')
    parser.add_argument('-l', '--length_threshold', required=False, default=100000, type=int,
                        help='Proximity (bp) to the contig ends to be used in the output BED file (-b) [100000]')
    parser.add_argument('-n', '--contig_id', required=False, default=None,
                        help='Optional contig_id (matching the #seq_name header in -i). '
                             'Only output information for that contig.')
    parser.add_argument('-c', '--circlator_merge_log', required=False, default=None,
                        help='Path to [prefix].circularize.log file from the circlator merge command. '
                             'Will output "true" if all contigs were circularized and "false" if not (to STDOUT). '
                             'This function is the only function independent of -i.')
    parser.add_argument('-v', '--verbose', required=False, action='store_true',
                        help='Enable for verbose logging.')
    command_line_args = parser.parse_args()
    main(command_line_args)

