#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: Utilities for the Rotary hybrid assembly workflow.
"""
import gzip
import os

from rotary.sample import file_is_gzipped


def check_for_files(output_dir_path, files_to_check):
    """
    Checks for existing configuration files from and existing run.

    :param output_dir_path: The path to the output directory.
    :param files_to_check: The files to check for.
    :return: True if the files are all present and False if they are not.
    """
    paths = [os.path.join(output_dir_path, file_name) for file_name in files_to_check]
    presences = [True for path in paths if os.path.exists(path)]

    if len(presences) == len(files_to_check):
        all_files_present = True
    elif len(presences) == 0:
        all_files_present = False
    else:
        raise FileNotFoundError(
            f'The directory ({output_dir_path}) has a mixture of files {files_to_check} present and absent.')

    return all_files_present


def get_cli_arg(args, argument):
    """
    Gets the command line argument if present or not.

    :param args: The command line arguments.
    :param argument: The name of the CLI arg.
    :return:
    """
    if hasattr(args, argument):
        result = getattr(args, argument)
    else:
        result = None
    return result


def get_cli_arg_path(args, argument):
    """
    Processes command line arguments with CLI paths.

    :param args: All command line arguments
    :param argument: The CLI arg to be found.
    :return: Either None or a sanitized absolute path to the file.
    """
    cli_path = get_cli_arg(args, argument)

    if cli_path:
        cli_path = sanitize_cli_path(cli_path)

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


def gzip_file(in_file_path, out_file_path):
    """
    Gzips a file.

    :param in_file_path: The path to the input file to be compressed.
    :param out_file_path: The path to the gzipped output file.
    """
    with open(in_file_path, 'rb') as file_in:
        with gzip.open(out_file_path, 'wb') as file_out:
            file_out.writelines(file_in)


def symlink_or_compress(in_file_path, out_file_path):
    """
    Symlinks the input file to the new path if it is already compressed or
    compresses and saves it at the new path otherwise.

    :param in_file_path: The path to the input file.
    :param out_file_path: The path to the output file.
    """
    if file_is_gzipped(in_file_path):
        os.symlink(in_file_path, out_file_path)
    else:
        gzip_file(in_file_path, out_file_path)


def get_contamination_reference_files(config_entry: list, db_path: str, phix_accession: str, human_accession: str):
    """
    Returns a list of file paths for contamination references based on short IDs

    :param config_entry: value of config.get('contamination_references'); should be a list
    :param db_path: path to the rotary database folder
    :param phix_accession: NCBI assembly accession code for the PhiX genome assembly
    :param human_accession: NCBI assembly accession code for the human genome assembly
    :return: list of file paths for the contamination reference sequences
    """

    contamination_reference_paths = []

    if isinstance(config_entry, list) is True:

        for contamination_entry in config_entry:

            if contamination_entry.lower() == 'phix':
                contamination_reference_paths.append(os.path.join(db_path, 'contamination_references',
                    f'PhiX_{phix_accession}.fna.gz'))

            elif contamination_entry.lower() == 'human':
                contamination_reference_paths.append(os.path.join(db_path,'contamination_references',
                    f'human_{human_accession}.fna.gz'))

            else:
                raise ValueError(f'Contaminant genome must be PhiX or human, but you provided {contamination_entry}')

    else:
        raise TypeError(f'Input type must be a list, but you provided type {type(config_entry)}')

    return contamination_reference_paths
