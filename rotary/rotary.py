#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: A command-line tool for the Rotary hybrid assembly workflow.
"""
import argparse
import csv
import os
import subprocess
import sys

import psutil
from ruamel.yaml import YAML

from rotary.sample import SequencingFile, Sample
from rotary.utils import check_for_files, get_cli_arg_path

yaml = YAML()
rotary_config_name = 'config.yaml'


def main():
    """
    Collects input arguments and selects a command to perform.
    """

    """
    When installing through pip with pyprojects.toml a new python script is generated
    that calls main() out of rotary.py. Run parser inside main() so it can be called 
    externally as a function.
    """
    parser = parse_cli()
    args = parser.parse_args()

    output_dir_path = get_cli_arg_path(args, 'output_dir')
    config_path = get_cli_arg_path(args, 'config')
    database_dir_path = get_cli_arg_path(args, 'database_dir')

    if hasattr(args, 'jobs'):
        jobs = args.jobs
    else:
        jobs = None

    if hasattr(args, 'snakemake_args'):
        snakemake_args = args.snakemake_args
        if snakemake_args:
            snakemake_args = snakemake_args.split()
    else:
        snakemake_args = None

    if hasattr(args, 'force'):
        override_existing_files = args.force
    else:
        override_existing_files = None

    try:
        os.makedirs(output_dir_path, exist_ok=True)
    except TypeError:
        # Catches error that occurs if you run just 'rotary' rather than 'rotary -h' in the CLI.
        parser.print_help()
        sys.exit()

    if not config_path: # If no config is specified via CLI.
        if hasattr(args, 'run'):
            # If the run sub-command was used, use the config that is in the output dir.
            config_path = os.path.join(output_dir_path, 'config.yaml')
        else:
            # Else get the default config that is saved with this python package.
            config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.yaml')

    with open(config_path) as config_file:
        config = yaml.load(config_file)

    if database_dir_path:
        config['db_dir'] = database_dir_path
    else:
        config['db_dir'] = os.path.join(output_dir_path, 'databases')

    config['threads'] = os.cpu_count()
    # 95% of the machines memory rounded to the nearest gigabyte.
    config['memory'] = round(psutil.virtual_memory().total / (1024 ** 3) * 0.95)

    # Check for the presence of the run configuration files in the output directory.
    run_files = ['config.yaml', 'samples.tsv']
    run_file_presences = check_for_files(directory_path=output_dir_path,
                                         file_names=run_files)

    if run_file_presences == [False, False]:
        has_run_files = False
    elif run_file_presences == [True, True]:
        has_run_files = True
    else:
        raise ValueError('The output directory is missing a run configuration file: {}'.format(zip(run_files,
                                                                                                   run_file_presences)))
    existing_message = f'Existing run configuration files {run_files} in the output directory.'

    conda_env_directory = os.path.join(output_dir_path, 'conda_env')
    os.makedirs(conda_env_directory, exist_ok=True)

    # Select the sub-command to run.
    if hasattr(args, 'run'):
        if has_run_files:
            run_rotary_workflow(config_path=config_path, output_dir_path=output_dir_path, jobs=jobs,
                                conda_env_directory=conda_env_directory, snakemake_custom_args=snakemake_args)
        else:
            raise FileNotFoundError(
                f'Missing run configuration files {run_files}, run either the run_one or init subcommands.')
    elif hasattr(args, 'run_one'):
        if has_run_files:
            if override_existing_files:
                run_one(args, config, output_dir_path, conda_env_directory, jobs, snakemake_args)
            else:
                raise FileExistsError(existing_message)
        else:
            run_one(args, config, output_dir_path, conda_env_directory, jobs, snakemake_args)
    elif hasattr(args, 'init'):
        if has_run_files:
            raise FileExistsError(existing_message)
        else:
            init(args, config, output_dir_path)
    else:
        parser.print_help()


def run_one(args, config, output_dir_path, conda_env_directory, jobs, snakemake_args):
    """
    Run the Rotary workflow in a single sample.

    :param args: The command-line arguments.
    :param config: The config information as produced by ruamel.yaml parsing of a config.yaml file.
    :param output_dir_path: The path to the output Rotary directory.
    :param conda_env_directory: The path to where snakemake should store its conda environments.
    :param jobs: The number of system CPU cores to use.
    :param snakemake_args: Custom CLI args to be passed to snakemake.
    """
    long_path = get_cli_arg_path(args, 'long')
    left_path = get_cli_arg_path(args, 'left')
    right_path = get_cli_arg_path(args, 'right')
    sequencing_files = [SequencingFile(path) for path in [long_path, left_path, right_path]]
    sample = Sample(long_file=sequencing_files[0],
                    short_file_one=sequencing_files[1],
                    short_file_two=sequencing_files[2])
    create_sample_tsv(output_dir_path, [sample])
    config['sample_id'] = sample.identifier
    config['longreads'] = sample.long_read_path
    config['qc_short_r1'] = sample.short_read_left_path
    config['qc_short_r2'] = sample.short_read_right_path
    config_path = write_config_file(output_dir_path=output_dir_path, config=config)
    run_rotary_workflow(config_path=config_path, output_dir_path=output_dir_path, jobs=jobs,
                        conda_env_directory=conda_env_directory, snakemake_custom_args=snakemake_args)


def run_rotary_workflow(config_path, output_dir_path, conda_env_directory, jobs=None, snakemake_custom_args=None):
    """
    Run the rotary snakemake workflow.

    :param conda_env_directory: Path to the conda environment directory.
    :param output_dir_path: Path to the output directory.
    :param config_path: The path to the config file.
    :param jobs: The number of system CPU cores to use.
    :param snakemake_custom_args: Custom CLI args to be passed to snakemake.
    """
    snake_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rules', 'rotary.smk')

    snakemake_args = ['snakemake', '--snakefile', snake_file_path, '--configfile', config_path, '--directory',
                      output_dir_path, '--conda-prefix', conda_env_directory, '--conda-frontend',
                      'mamba', '--use-conda', '--reason', '--rerun-incomplete', '--printshellcmds']

    if jobs:
        snakemake_args = snakemake_args + ['--jobs', jobs]
    else:
        snakemake_args = snakemake_args + ['--jobs', str(os.cpu_count())]

    if snakemake_custom_args:
        snakemake_args = snakemake_args + snakemake_custom_args

    cmd = ' '.join(snakemake_args)

    print(f'Executing: {cmd}')
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        sys.exit(1)


def init(args, config, output_dir_path):
    """
    Runs code for the init command-line mode. Sets up a rotary project directory.

    :param args: The command-line arguments.
    :param config: The config information as produced by ruamel.yaml parsing of a config.yaml file.
    :param output_dir_path: The path to the output Rotary directory.
    """
    input_path = get_cli_arg_path(args, 'input_dir')
    fastq_files = []

    for file in os.listdir(input_path):
        if 'fastq' in file:
            file = SequencingFile(file_path=(os.path.join(input_path, file)))
            fastq_files.append(file)

    samples_files = {}
    for file in fastq_files:
        identifier = file.identifier
        if identifier in samples_files.keys():
            samples_files[identifier].append(file)
        else:
            samples_files[identifier] = [file]

    samples = []
    for identifier, sequencing_files in samples_files.items():
        if len(sequencing_files) != 3:
            raise ValueError(f'Sample {identifier} should have three sequencing files')

        long_file = None
        left_short_file = None
        right_short_file = None
        for file in sequencing_files:
            if file.r_value == 'R1':
                left_short_file = file
            elif file.r_value == 'R2':
                right_short_file = file
            else:
                long_file = file

        sample = Sample(long_file, left_short_file, right_short_file)
        samples.append(sample)

    create_sample_tsv(output_dir_path, samples)
    write_config_file(config, output_dir_path)


def create_sample_tsv(output_dir_path, samples):
    """
    Generates a TSV file in the output directory with a series of CLI paths for files belonging to each sample.

    :param output_dir_path: The path to the output Rotary directory.
    :param samples: A list of Sample objects.
    """
    with open(os.path.join(output_dir_path, 'samples.tsv'), 'w') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        header = ['sample_id', 'long-read', 'short-read_R1', 'short-read_R2']
        tsv_writer.writerow(header)
        for current_sample in samples:
            tsv_writer.writerow(current_sample.sample_file_row)


def write_config_file(config, output_dir_path):
    """
    Write the updated config to the output directory.

    :param config: The config information as produced by ruamel.yaml parsing of a config.yaml file.
    :param output_dir_path: The path to the output directory.
    :return: The absolute path the output config file.
    """
    output_directory_config_path = os.path.join(output_dir_path, 'config.yaml')
    with open(output_directory_config_path, 'w') as config_file:
        yaml.dump(config, config_file)

    return output_directory_config_path


def parse_cli():
    """
    Parses the CLI arguments.
    :return: An argparse parser object.
    """
    cli_title = """A command-line interface for the Rotary hybrid assembly workflow."""
    parser = argparse.ArgumentParser(description=cli_title)
    subparsers = parser.add_subparsers(help='Available Sub-commands')

    # =========================
    # Declare Run Sub-command
    run_help = """Runs the Rotary workflow on a Rotary project directory generated by the rotary init command."""
    parser_run = subparsers.add_parser('run', help=run_help)
    parser_run.add_argument('-c', '--config', metavar='YAML',
                            help='path to the rotary yaml config file')
    parser_run.add_argument('-o', '--output_dir', metavar='PATH', default=os.getcwd(),
                            help='path the output/rotary project directory')
    parser_run.add_argument('-j', '--jobs', metavar='JOBS',
                            help='number of threads that rotary should use (overrides config)')
    parser_run.add_argument('-s', '--snakemake_args', metavar='',
                            help="quoted string with arguments to be passed to snakemake. i.e., -s'--dag' (no space after -s)")
    # Add attribute to tell main() what sub-command was called.
    parser_run.set_defaults(run=True)

    # =========================
    # Declare Run_One Sub-command
    run_help = """Runs the Rotary workflow on specified sequencing files"""
    parser_run_one = subparsers.add_parser('run_one', help=run_help)
    parser_run_one.add_argument('-l', '--long', metavar='FASTQ',
                                help='path to the Oxford Nanopore long-read .fastq(.gz) file')
    parser_run_one.add_argument('-r1', '--left', metavar='FASTQ',
                                help='path to the left Illumina short-read .fastq(.gz) file')
    parser_run_one.add_argument('-r2', '--right', metavar='FASTQ',
                                help='path to the right Illumina short-read .fastq(.gz) file')
    parser_run_one.add_argument('-c', '--config', metavar='YAML',
                                help='path to the rotary yaml config file')
    parser_run_one.add_argument('-o', '--output_dir', metavar='PATH', default=os.getcwd(),
                                help='path the output/rotary project directory')
    parser_run_one.add_argument('-d', '--database_dir', metavar='PATH',
                                help='path the rotary database directory')
    parser_run_one.add_argument('-j', '--jobs', metavar='JOBS',
                                help='number of threads that rotary should use (overrides config)')
    parser_run_one.add_argument('-s', '--snakemake_args', metavar='',
                                help="quoted string with arguments to be passed to snakemake. i.e., -s'--dag' (no space after -s)")
    parser_run_one.add_argument('-f', '--force', action='store_true',
                                help="override existing run configuration files.")
    # Add attribute to tell main() what sub-command was called.
    parser_run_one.set_defaults(run_one=True)

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
    parser_init.set_defaults(init=True)
    return parser


if __name__ == '__main__':
    main()
