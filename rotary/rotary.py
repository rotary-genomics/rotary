#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: A command-line interface for the Rotary hybrid assembly workflow.
"""
import argparse
import os
import subprocess
import sys

import psutil
from ruamel.yaml import YAML

from rotary.sample import SequencingFile, Sample, is_fastq_file, create_sample_tsv
from rotary.utils import get_cli_arg_path, get_cli_arg, check_for_files

yaml = YAML()
rotary_config_name = 'config.yaml'
run_files = [rotary_config_name, 'samples.tsv']


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

    # Select the sub-command to run.
    if hasattr(args, 'run'):
        run(args)
    elif hasattr(args, 'run_one'):
        run_one(args)
    elif hasattr(args, 'init'):
        init(args)
    else:
        parser.print_help()


def run(args):
    """
    Run the Rotary workflow a Rotary project directory.

    :param args: The command-line arguments.
    """
    output_dir_path = get_cli_arg_path(args, 'output_dir')
    # Check for the presence of the run configuration files in the output directory.

    if not check_for_files(output_dir_path, run_files):
        raise FileNotFoundError(
            f'Missing run configuration files {run_files}, run either the run_one or init subcommands.')

    jobs = get_cli_arg(args, 'jobs')
    snakemake_args = get_snakemake_args(args)

    config_path = get_cli_arg_path(args, 'config')
    if not config_path:
        config_path = os.path.join(output_dir_path, 'config.yaml')

    with open(config_path) as config_file:
        config = yaml.load(config_file)

    conda_env_directory = os.path.join(config['db_dir'], 'rotary_conda_env')
    os.makedirs(conda_env_directory, exist_ok=True)

    run_rotary_workflow(config_path=config_path, output_dir_path=output_dir_path, jobs=jobs,
                        conda_env_directory=conda_env_directory, snakemake_custom_args=snakemake_args)


def run_one(args):
    """
    Run the Rotary workflow in a single sample.

    :param args: The command-line arguments.
    """
    output_dir_path = get_cli_arg_path(args, 'output_dir')

    # Checks for existing run files and sets up the run directory.
    conda_env_directory, config_path = setup_run_directory(args, output_dir_path)

    jobs = get_cli_arg(args, 'jobs')
    snakemake_args = get_snakemake_args(args)

    sequencing_files = [SequencingFile(path) for path in [(get_cli_arg_path(args, 'long')),
                                                          (get_cli_arg_path(args, 'left')),
                                                          (get_cli_arg_path(args, 'right'))]]

    sample = Sample(long_file=sequencing_files[0],
                    short_file_one=sequencing_files[1],
                    short_file_two=sequencing_files[2],
                    integrity_check=False)  # Don't do integrity check on user specified files.

    create_sample_tsv(output_dir_path, [sample])

    run_rotary_workflow(config_path=config_path, output_dir_path=output_dir_path, jobs=jobs,
                        conda_env_directory=conda_env_directory, snakemake_custom_args=snakemake_args)


def init(args):
    """
    Runs code for the init command-line mode. Sets up a rotary project directory.

    :param args: The command-line arguments.
    """
    output_dir_path = get_cli_arg_path(args, 'output_dir')

    # Checks for existing run files and sets up the run directory.
    setup_run_directory(args, output_dir_path)

    input_path = get_cli_arg_path(args, 'input_dir')
    fastq_files = []

    for file_path in os.listdir(input_path):
        filename = os.path.basename(file_path)
        if is_fastq_file(filename):
            fastq_files.append(SequencingFile(file_path=os.path.join(input_path, file_path)))

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


def setup_run_directory(args, output_dir_path):
    """
    Sets up the basic configuration files inside output directory.

    :param args: The command-line arguments.
    :param output_dir_path: The path to the output Rotary directory.
    :return: Path to the conda environment directory and the path to the configuration file.
    """
    # Check if run files exist in the output.
    if not get_cli_arg(args, 'force') and check_for_files(output_dir_path, run_files):
        raise FileExistsError(
            f'The output directory ({output_dir_path}) already has configuration files {run_files} from a prior run.')

    os.makedirs(output_dir_path, exist_ok=True)

    input_config_path = get_cli_arg_path(args, 'config')
    if not input_config_path:  # If no config is specified via CLI grab the default config.yaml.
        input_config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.yaml')

    with open(input_config_path) as config_file:
        config = yaml.load(config_file)

    # modify config to use available CPU and RAM.
    config = modify_config_with_available_computational_resources(config)

    database_dir_path = get_cli_arg_path(args, 'database_dir')
    if database_dir_path:
        config['db_dir'] = database_dir_path
    else:
        parent_dir = os.path.split(output_dir_path)[0]
        config['db_dir'] = os.path.join(parent_dir, 'rotary_database')

    output_config_path = os.path.join(output_dir_path, 'config.yaml')
    conda_env_path = os.path.join(config['db_dir'], 'rotary_conda_env')
    os.makedirs(conda_env_path, exist_ok=True)

    with open(output_config_path, 'w') as file:
        yaml.dump(config, file)

    return conda_env_path, output_config_path


def get_snakemake_args(args):
    """
    Gets the snakemake arguments from the command line arguments.

    :param args: The command line arguments.
    :return: The snakemake arguments.
    """
    snakemake_args = get_cli_arg(args, 'snakemake_args')

    if snakemake_args:
        snakemake_args = snakemake_args.split()

    return snakemake_args


def modify_config_with_available_computational_resources(config):
    """
    Modify the config to use available CPU and RAM.

    :param config: The config object.
    """
    config['threads'] = os.cpu_count()
    # 95% of the machines memory rounded to the nearest gigabyte.
    config['memory'] = round(psutil.virtual_memory().total / (1024 ** 3) * 0.95)

    return config


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
    run_one_help = """Runs the Rotary workflow on specified sequencing files"""
    parser_run_one = subparsers.add_parser('run_one', help=run_one_help)
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
    init_help = """Creates a rotary project directory containing a default """
    parser_init = subparsers.add_parser('init', help=init_help)
    parser_init.add_argument('-o', '--output_dir', metavar='PATH',
                             help='path the output/rotary project directory', default=os.getcwd())
    parser_init.add_argument('-d', '--database_dir', metavar='PATH', required=True,
                             help='path the rotary database directory')
    parser_init.add_argument('-i', '--input_dir', metavar='PATH', required=True,
                             help='path to a directory containing Oxford Nanopore long-read and Illumina short-read .fastq(.gz) files')
    parser_init.add_argument('-f', '--force', action='store_true',
                             help="override existing run configuration files.")
    parser_init.set_defaults(init=True)
    return parser


if __name__ == '__main__':
    main()
