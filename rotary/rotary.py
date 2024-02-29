#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: A command-line interface for the Rotary hybrid assembly workflow.
"""
import argparse
import os

from pungi.dataset import generate_dataset_from_fastq_directory, Dataset
from pungi.run import setup_run_directory, run_snakemake_workflow, load_yaml_config, get_snakemake_args
from pungi.sample import SequencingFile, LongReadSampleWithPairedPolishingShortReads
from pungi.utils import get_cli_arg_path, get_cli_arg, check_for_files, get_config_path

rotary_config_name = 'config.yaml'
run_files = [rotary_config_name, 'samples.tsv']

snake_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rules', 'rotary.smk')
provided_config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.yaml')
sample_tsv_header_fields = ['sample_id', 'long-read', 'short-read_R1', 'short-read_R2']


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

    config_path = get_config_path(args, default_config_path=os.path.join(output_dir_path, 'config.yaml'))

    config = load_yaml_config(config_path)
    conda_env_dir = os.path.join(config['db_dir'], 'rotary_conda_envs')

    run_snakemake_workflow(config_path=config_path, snake_file_path=snake_file_path, output_dir_path=output_dir_path,
                           jobs=jobs, conda_env_dir_path=conda_env_dir, snakemake_custom_args=snakemake_args)


def run_one(args):
    """
    Run the Rotary workflow in a single sample.

    :param args: The command-line arguments.
    """
    config_path = get_config_path(args, default_config_path=provided_config_path)
    config = load_yaml_config(config_path)
    output_dir_path = get_cli_arg_path(args, 'output_dir')

    # Checks for existing run files and sets up the run directory.
    setup_run_directory(args, output_dir_path, run_files, config)

    jobs = get_cli_arg(args, 'jobs')
    snakemake_args = get_snakemake_args(args)

    sequencing_files = [SequencingFile(path) for path in [(get_cli_arg_path(args, 'long')),
                                                          (get_cli_arg_path(args, 'left')),
                                                          (get_cli_arg_path(args, 'right'))]]

    sample = LongReadSampleWithPairedPolishingShortReads(long_file=sequencing_files[0],
                                                         short_file_left=sequencing_files[1],
                                                         short_file_right=sequencing_files[2],
                                                         # Don't do an identifier check on user-specified files.
                                                         identifier_check=False,
                                                         # Don't do integrity check on user-specified files.
                                                         integrity_check=False)

    dataset = Dataset(sample)

    dataset.create_sample_tsv(output_dir_path, header=sample_tsv_header_fields)

    conda_env_dir = os.path.join(config['db_dir'], 'rotary_conda_envs')

    run_snakemake_workflow(config_path=config_path, snake_file_path=snake_file_path, output_dir_path=output_dir_path,
                           jobs=jobs, conda_env_dir_path=conda_env_dir, snakemake_custom_args=snakemake_args)


def init(args):
    """
    Runs code for the init command-line mode. Sets up a rotary project directory.

    :param args: The command-line arguments.
    """
    config = load_yaml_config(get_config_path(args, default_config_path=provided_config_path))
    output_dir_path = get_cli_arg_path(args, 'output_dir')

    # Checks for existing run files and sets up the run directory.
    setup_run_directory(args, output_dir_path, run_files, config)

    input_path = get_cli_arg_path(args, 'input_dir')

    dataset = generate_dataset_from_fastq_directory(input_path, expected_files_per_sample=3)

    dataset.create_sample_tsv(output_dir_path, header=sample_tsv_header_fields)


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
