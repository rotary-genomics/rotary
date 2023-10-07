#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: A command-line tool for the Rotary hybrid assembly workflow.
"""

import argparse
import os
import csv
import re

short_read_r_regex = re.compile("[Rr][12]")

def main(args):
    """
    Collects input arguments and selects a command to perform.
    :param args: The command line arguments.
    """
    output_dir_path = get_cli_arg_path(args, 'output_dir')
    database_dir_path = get_cli_arg_path(args, 'database_dir')
    config_path = get_cli_arg_path(args, 'config')


    if hasattr(args, 'run'):
        jobs = args.jobs
        snakemake_args = args.snakemake_args

        if hasattr(args, 'long'):
            print('has input files')
            long_path =  get_cli_arg_path(args, 'long')
            left_path = get_cli_arg_path(args, 'left')
            right_path = get_cli_arg_path(args, 'right')
            sample = Sample(long_path=long_path,short_path_one=left_path,short_path_two=right_path)
            create_sample_tsv(output_dir_path, [sample])
        else:
            print('has directory')
    elif hasattr(args, 'init'):
        print('entered init')
    else:
        parser.print_help()

def create_sample_tsv(output_dir_path, samples):
    """
    Generates a TSV file in the output directory with a series of CLI paths for files belonging to each sample.

    :param output_dir_path: The path to the output Rotary directory.
    :param samples: A list of Sample objects.
    """
    with open(os.path.join(output_dir_path, 'samples.tsv'), 'w') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        header = ['sample_id', 'long-read','short-read_R1', 'short-read_R2']
        tsv_writer.writerow(header)
        for current_sample in samples:
            row = [current_sample.identifier, current_sample.long_read_file_path,
                   current_sample.short_read_left_file_path, current_sample.short_read_right_file_path]
            tsv_writer.writerow(row)

class Sample(object):
    """
    An object representing a series of FASTQ files that all belong to the same sample.
    """
    def __init__(self, long_path, short_path_one, short_path_two):
        file_paths = [long_path, short_path_one, short_path_two]
        file_names = [os.path.basename(path) for path in file_paths]

        for name in file_names:
            if 'fastq' not in name:
                raise ValueError('{} is not a fastq file.'.format(name))

        # Split file names into two lists that either contain the sample identifier of
        # each file name or the remaining part of the file name after the identifier is removed.
        current_id = None
        identifiers = []
        file_name_no_identifier = []
        for file in file_names:
            if '_' in file:
                current_id, remaining = file.split('_', 1) # Split on first _ if found.
            else:
                current_id, remaining = file.split('.', 1) # Split on file extension . if _ not found.
            identifiers.append(current_id)
            file_name_no_identifier.append(remaining)

        # The identifiers should all be the same. If not raise and exception.
        if len(set(identifiers)) == 1:
            self.identifier = current_id
        else:
            raise ValueError('Sample identifiers of the input fastq files do not match: {}'.format(identifiers))

        # Check for R1 or R2 in the file names with identifier removed.
        r_matches = [short_read_r_regex.findall(name) for name in file_name_no_identifier]
        long_r_value, short_one_r_value, short_two_r_value = r_matches

        # The long read file should not have a R1 or R2.
        if not long_r_value:
            self.long_read_file_path = long_path
        else:
            raise ValueError("The long-read file ({}) should not have an R designator.".format(long_path))

        # The short read files should both have R1s or R2s.
        if short_one_r_value and short_two_r_value:
            short_one_r_value = short_one_r_value[0].upper()
            short_two_r_value = short_two_r_value[0].upper()
        else:
            raise ValueError("Short-read file '{}' or '{}' is missing its R designator.".format(short_path_one,short_path_two))

        # The short read files should not have the same R1 or R2.
        if short_one_r_value == short_two_r_value:
            matching_r_values_message = 'Left ({}) and right ({}) short-read files must have different R designators.'
            raise ValueError(matching_r_values_message.format(short_path_one, short_path_two))

        # Assign the file with R1 to the left file path and the file with R2 to the right filepath.
        if short_one_r_value == 'R1':
            self.short_read_left_file_path = short_path_one
            self.short_read_right_file_path = short_path_two
        else:
            self.short_read_left_file_path = short_path_two
            self.short_read_right_file_path = short_path_one

def get_cli_arg_path(args, argument):
    """
    Processes command line arguments with CLI paths.

    :param args: All command line arguments
    :param argument: The CLI arg to be found.
    :return: Either None or a sanitized absolute path to the file.
    """
    if hasattr(args, argument):
        cli_path = getattr(args, argument)
        if cli_path:
            cli_path = sanitize_cli_path(cli_path)
    else:
        cli_path = None

    return cli_path


def sanitize_cli_path(cli_path):
    """
    Performs expansion of '~' and shell variables such as "$HOME" into absolute paths.
    Expands the path to the absolut path.

    :param cli_path: The path to expand
    :return: An expanded path.
    """
    sanitized_path = os.path.abspath(os.path.expanduser(os.path.expandvars(cli_path)))
    return sanitized_path

if __name__ == '__main__':
    cli_title = """A command-line interface for the Rotary hybrid assembly workflow."""
    parser = argparse.ArgumentParser(description=cli_title)
    subparsers = parser.add_subparsers(help='Available Sub-commands')

    # =========================
    # Declare Run Sub-command
    run_help = """Runs the Rotary workflow either on specified sequencing files or on a Rotary project directory 
    generated by the rotary init command."""
    parser_run = subparsers.add_parser('run', help=run_help)

    parser_run.add_argument('-l', '--long', metavar='FASTQ',
                            help='path to the Oxford Nanopore long-read .fastq(.gz) file')
    parser_run.add_argument('-r1', '--left', metavar='FASTQ',
                            help='path to the left Illumina short-read .fastq(.gz) file')
    parser_run.add_argument('-r2', '--right', metavar='FASTQ',
                            help='path to the right Illumina short-read .fastq(.gz) file')
    parser_run.add_argument('-c', '--config', metavar='YAML',
                            help='path to the rotary yaml config file')
    parser_run.add_argument('-o', '--output_dir', metavar='PATH', default=os.getcwd(),
                            help='path the output/rotary project directory')
    parser_run.add_argument('-d', '--database_dir', metavar='PATH',
                            help='path the rotary database directory')
    parser_run.add_argument('-j', '--jobs', metavar='JOBS',
                            help='number of jobs that rotary should use (overrides config)')
    parser_run.add_argument('-s', '--snakemake_args', metavar='kwargs', nargs='*',
                            help='quoted string with arguments to be passed to snakemake')

    # Add attribute to tell main() what sub-command was called.
    parser_run.set_defaults(run=True)

    # =========================
    # Declare Init Sub-command
    run_help = """Creates a rotary project directory containing a default """
    parser_init = subparsers.add_parser('init', help=run_help)

    parser_init.add_argument('-o', '--output_dir', metavar='PATH',
                            help='path the output/rotary project directory', default=os.getcwd())
    parser_init.add_argument('-d', '--database_dir', metavar='PATH', required=True,
                            help='path the rotary database directory')
    parser_init.add_argument('-i', '--input_dir', metavar='PATH', required=True,
                             help='path to a directory containing Oxford Nanopore long-read and Illumina short-read .fastq(.gz) files')

    parser_run.set_defaults(init=True)

    cli_args = parser.parse_args()
    main(cli_args)