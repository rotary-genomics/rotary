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


def summarize_assembly_info(assembly_info_filepath, linear_contig_list_filepath=None,
                            circular_contig_list_filepath=None, bed_filepath=None,
                            length_threshold=100000, contig_id=None):

    logger.debug('Loading assembly info file')
    assembly_info = pd.read_csv(assembly_info_filepath, sep='\t')

    circular_contigs = assembly_info[assembly_info['circ.'] == 'Y'][['#seq_name', 'length', 'circ.']]
    linear_contigs = assembly_info[assembly_info['circ.'] == 'N'][['#seq_name', 'length', 'circ.']]

    # TODO - add multi-contig support
    if contig_id is not None:
        logger.debug('Filtering by contig ID')
        circular_contigs = circular_contigs[circular_contigs['#seq_name'] == contig_id]

    if linear_contig_list_filepath is not None:
        logger.debug('Exporting linear contig list')
        linear_contigs['#seq_name'].to_csv(linear_contig_list_filepath, index=False, header=None)

    if circular_contig_list_filepath is not None:
        logger.debug('Exporting circular contig list')
        circular_contigs['#seq_name'].to_csv(circular_contig_list_filepath, index=False, header=None)

    if bed_filepath is not None:
        logger.debug('Making BED file with end proximity threshold of ' + str(length_threshold))

        contigs = []
        starts = []
        stops = []

        for index, row in circular_contigs.iterrows():
            contig, length, circ = row

            if length < length_threshold:
                contigs.append(contig)
                starts.append(0)
                stops.append(length)

            else:
                half_threshold = int(round(length_threshold / 2, 0))

                contigs.append(contig)
                starts.append(0)
                stops.append(half_threshold)

                contigs.append(contig)
                starts.append(length - half_threshold)
                stops.append(length)
    
        end_regions = pd.DataFrame({'contig': contigs, 'start': starts, 'stop': stops})
        end_regions.to_csv(bed_filepath, sep='\t', header=None, index=False)


def check_circlator_logfile(circlator_logfile):

    logger.debug('Loading circlator logfile')
    circlator_info = pd.read_csv(circlator_logfile, sep='\t')[['#Contig', 'circularised']]

    # If a row is '1', it means it was circularized, but if '0', it means it was not circularized.
    # So if all rows are 1 (i.e., the sum of rows / # of rows is 1), it means everything was circularized.
    if circlator_info['circularised'].sum() / circlator_info.shape[0] == 1:
        logger.debug('Everything is circularized.')
        print('true')

    elif circlator_info['circularised'].sum() >= 0:
        logger.debug('Not everything is circularized.')
        print('false')

    else:
        logger.error('File processing error. # of non-circularized contigs is not >=0. Exiting...')
        sys.exit(1)


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

