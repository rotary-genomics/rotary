#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: Code for building CLI programs that call Snakemake.
"""
import os
import subprocess
import sys

import psutil
from ruamel.yaml import YAML

from rotary.utils import get_cli_arg, check_for_files, get_cli_arg_path

yaml = YAML()


def load_yaml_config(config_path):
    """
    Load YAML configuration file.

    :param config_path: The path to the YAML configuration file.
    :return: A ruamel.yaml object with the configuration data.
    """
    with open(config_path) as config_file:
        config = yaml.load(config_file)
    return config


def dump_yaml_config(config, output_config_path):
    """
    Dumps the given configuration dictionary into a YAML file.

    :param config: A ruamel.yaml object containing the configuration data to be dumped into YAML format.
    :param output_config_path: The path of the output YAML file.
    """
    with open(output_config_path, 'w') as file:
        yaml.dump(config, file)


def setup_run_directory(args, output_dir_path, run_files):
    """
    Sets up the basic configuration files inside the output directory.

    :param args: The command-line arguments.
    :param output_dir_path: The path to the output Rotary directory.
    :param run_files: A list of output files to check for once the snakemake run directory is created.
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

    config = load_yaml_config(input_config_path)

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

    dump_yaml_config(config, output_config_path)

    return conda_env_path, output_config_path


def modify_config_with_available_computational_resources(config):
    """
    Modify the config to use available CPU and RAM.

    :param config: The config object.
    """
    config['threads'] = os.cpu_count()
    # 95% of the machine's memory rounded to the nearest gigabyte.
    config['memory'] = round(psutil.virtual_memory().total / (1024 ** 3) * 0.95)

    return config


def run_snakemake_workflow(config_path, snake_file_path, output_dir_path, conda_env_directory, jobs=None,
                           snakemake_custom_args=None):
    """
    Run a workflow using Snakemake.

    :param config_path: Path to the configuration file.
    :param snake_file_path: Path to the Snakefile.
    :param output_dir_path: Path to the output directory.
    :param conda_env_directory: Path to the Conda environment directory.
    :param jobs: Number of parallel jobs to run (default is the number of CPU cores).
    :param snakemake_custom_args: Additional Snakemake command-line arguments.
    """

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
