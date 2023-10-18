#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: Utilities for the Rotary hybrid assembly workflow.
"""
import os


def check_for_files(directory_path, file_names: list = None):
    """
    Checks if there are specific files in a given directory.

    :param directory_path: The path to the directory to check.
    :param file_names: A list of file names to check for in the output directory.
    :return: A list of true/false values for if each file is present.
    """

    run_paths = [os.path.join(directory_path, file_name) for file_name in file_names]

    return [os.path.exists(path) for path in run_paths]


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
