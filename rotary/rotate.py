#!/usr/bin/env python
# rotate.py
# Utility for rotating circular DNA elements
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024

import argparse
import logging
import os
import sys

import pandas as pd
from Bio import SeqIO


# Set up the logger
logger = logging.getLogger(__name__)
formatter = logging.Formatter('[ %(asctime)s ]: %(levelname)s: %(funcName)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


def main():
    """
    Collects input arguments and runs the rotate utility
    """

    """
    When installing through pip with pyprojects.toml a new python script is generated
    that calls main() out of rotate.py. Run parser inside main() so it can be called
    externally as a function.
    """
    parser = parse_cli()
    args = parser.parse_args()

    # Startup checks
    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Check output file
    # TODO - add overwrite flag
    # TODO - make this info a function and check all other output files
    # (TODO - add prefix flag)
    # TODO - consider also checking if the directory where the output file will be saved exists and is writeable
    # TODO - include a warning like logger.warning(f'Output file already exists: "{args.output_fasta}". Files may be overwritten.') to the log
    output_file_exists = os.path.isfile(args.output_fasta)

    if (output_file_exists is True) & (args.overwrite is False):

        logger.error(f'Output file already exists: "{args.output_fasta}". Will not continue. Set the '
                     f'--overwrite flag at your own risk if you want to overwrite files.')
        sys.exit(1)

    # Parse some command line inputs further
    # TODO - improve error handling if a string is provided instead of a real length
    # TODO - confirm what the default value of args.contig_names will be - might need to specify
    if args.contig_names is not None:
        contig_names = [int(x) for x in args.contig_names.split(',')]

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.debug('### SETTINGS ###')
    logger.debug(f'Input FastA: {args.input_fasta}')
    logger.debug(f'Output FastA: {args.output_fasta}')
    # TODO - finish the rest of the startup messages

    logger.debug(f'Threads: {args.threads}')
    logger.debug(f'Verbose logging: {args.verbose}')
    logger.debug('################')

    # TODO - add rotate code

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


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


def parse_cli():
    """
    Parses the CLI arguments
    :return: An argparse parser object
    """

    parser = argparse.ArgumentParser(
        description=f'{os.path.basename(sys.argv[0])}: utility for rotating circular contigs. \n'
                    '  Copyright Jackson M. Tsuji and Lee Bergstrand, 2024',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    required_settings = parser.add_argument_group('Required')
    rotate_settings = parser.add_argument_group('Rotate options')
    workflow_settings = parser.add_argument_group('Workflow options')

    required_settings.add_argument('-i', '--input_fasta', required=True, type=str,
                                   help='Input contig fasta file')
    required_settings.add_argument('-o', '--output_fasta', required=True, type=str,
                                   help='Output contig fasta file. (Note that all input contigs will always be '
                                        'written to output. If you select no rotate options, the input contigs will '
                                        'be output with no sequence changes.)')

    rotate_settings.add_argument('-m', '--midpoint', required=False, action='store_true',
                                 help='Rotate all contigs to their midpoint (overrides -p, -P, -t, and -T)')
    rotate_settings.add_argument('-p', '--rotate_position', required=False, type=int,
                                 help='Base pair position to rotate contigs to (to specify a unique rotate '
                                      'position for each contig, see -t). Incompatible with -P.')
    rotate_settings.add_argument('-t', '--rotate_position_table', required=False, type=str,
                                 help='Path to a tab-separated file that contains the exact rotate position of '
                                      'each contig. Columns (with headers) should be "contig_id" and "rotate_position". '
                                      'Overrides -p and -P. Only contigs specified in the table will be rotated, although '
                                      'all contigs in the input file will be written to the output file. Incompatible with '
                                      '-T.')
    rotate_settings.add_argument('-P', '--rotate_proportion', required=False, type=float,
                                 help='Fractional position to rotate contigs to (e.g., 0.3 to rotate 30% of total length). '
                                      'To specify a unique fractional position to rotate for each contig, see -T). '
                                      'Incompatible with -p.')
    rotate_settings.add_argument('-T', '--rotate_proportion_table', required=False, type=str,
                                 help='Path to a tab-separated file that contains the exact fractional positions to rotate '
                                      'each contig to. Columns (with headers) should be "contig_id" and "rotate_proportion". '
                                      'Overrides -p and -P. Only contigs specified in the table will be rotated, although '
                                      'all contigs in the input file will be written to the output file. Incompatible '
                                      'with -t.')
    rotate_settings.add_argument('-c', '--contig_names', required=False, type=str,
                                 help='Contigs to be rotated (comma-separated list of IDs). '
                                      'Rotate operations will only be applied to these contigs, '
                                      'although all contigs in the input file will be written to output. '
                                      'If used with -t or -T, the contigs listed in the table will be further subset '
                                      'to those that are also in the provided list of contig names.')
    rotate_settings.add_argument('-r', '--output_report', required=False, type=str,
                                 help='Path to write an optional report file that shows how contigs were rotated')
    rotate_settings.add_argument('-b', '--base_start_position', required=False, default=1, type=int,
                                 choices=[0, 1],
                                 help='Whether the first base position on a contig is considered as 1 (default) or 0.')

    workflow_settings.add_argument('-t', '--threads', required=False, default=1, type=int,
                                   help='Number of processors threads to use (default: 1)')
    workflow_settings.add_argument('-v', '--verbose', required=False, action='store_true',
                                   help='Enable verbose logging to screen')

    return parser


if __name__ == '__main__':
    main()
