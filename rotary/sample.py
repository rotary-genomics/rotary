#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: Classes and methods for representing Rotary sampling files.
"""
import os
import re

short_read_r_regex = re.compile("[Rr][12]")


class SequencingFile(object):
    """
    An object representing a FASTQ file.
    """

    def __init__(self, file_path):
        self.path = file_path
        self.name = os.path.basename(self.path)

        if not is_fastq_file(self.name):
            raise ValueError(f'{self.name} is not a fastq file.')

        if '_' in self.name:
            identifier, remaining = self.name.split('_', 1)  # Split on first _ if found.
        else:
            identifier, remaining = self.name.split('.', 1)  # Split on file extension . if _ not found.

        self.identifier = identifier

        # Check the remainder of the file name (e.g. _blam_blam_blam_R1.fastq.gz) for R1 or R2
        r_value = short_read_r_regex.findall(remaining)
        if r_value:
            self.r_value = r_value[0].upper()
        else:
            self.r_value = None


class Sample(object):
    """
    An object representing a series of FASTQ files that all belong to the same sample.
    """

    def __init__(self, long_file: SequencingFile, short_file_one: SequencingFile, short_file_two: SequencingFile,
                 integrity_check: bool = True):

        if integrity_check:
            self.assign_variables_with_integrity_check(long_file, short_file_one, short_file_two)
        else:
            self.identifier = long_file.identifier
            self.long_read_path = long_file.path
            self.short_read_left_path = short_file_one.path
            self.short_read_right_path = short_file_two.path

    def assign_variables_with_integrity_check(self, long_file, short_file_one, short_file_two):
        """
        Builds the sample object with integrity checks to ensure the proper files are mapping to a singel sample.

        :param long_file: A SequencingFile object representing the long read file.
        :param short_file_one: A SequencingFile object representing the first short read file.
        :param short_file_two: A SequencingFile object representing the second short read file.
        """
        sequencing_files = [long_file, short_file_one, short_file_two]
        identifiers = [file.identifier for file in sequencing_files]

        # Check for R1 or R2 in the file names with identifier removed.
        r_matches = [file.r_value for file in sequencing_files]
        long_r_value, short_one_r_value, short_two_r_value = r_matches

        long_path = long_file.path
        short_path_one = short_file_one.path
        short_path_two = short_file_two.path

        # The identifiers should all be the same. If not raise and exception.
        if len(set(identifiers)) == 1:
            self.identifier = identifiers[0]
        else:
            raise ValueError(f'Sample identifiers of the input fastq files do not match: {identifiers}')

        # The long read file should not have a R1 or R2.
        if not long_r_value:
            self.long_read_path = long_path
        else:
            raise ValueError(f"The long-read file ({long_path}) should not have an R designator.")

        # The short read files should both have R1s or R2s.
        if not short_one_r_value or not short_two_r_value:
            raise ValueError(f"Short-read file '{short_path_one}' or '{short_path_two}' is missing its R designator.")

        # The short read files should not have the same R1 or R2.
        if short_one_r_value == short_two_r_value:
            raise ValueError(
                f'Left ({short_path_one}) and right ({short_path_two}) short-read files must have different R designators.')

        # Assign the file with R1 to the left file path and the file with R2 to the right filepath.
        if short_one_r_value == 'R1':
            self.short_read_left_path = short_path_one
            self.short_read_right_path = short_path_two
        else:
            self.short_read_left_path = short_path_two
            self.short_read_right_path = short_path_one

    @property
    def sample_file_row(self):
        """
        :return: Returns a list containing the sample identifier and the paths to sample's the input files.
        """
        return [self.identifier, self.long_read_path, self.short_read_left_path, self.short_read_right_path]


def is_fastq_file(file_name):
    """
    Determines if file is a fastq file based in its extension.
    :param file_name: The name of the file.
    :return: True if the file contains 'fastq' or 'fq'.
    """
    file_extension = file_name.split(os.path.extsep, 1)[1]

    if 'fastq' in file_extension or 'fq' in file_extension:
        is_fastq = True
    else:
        is_fastq = False

    return is_fastq

def make_sample_from_sample_tsv_row(row):
    """
    Parses a row from a sample TSV file and returns a Sample object.

    :param row: A row from a sample tsv file as parsed by the CSV module.
    :return: A Sample object representing the row.
    """
    sample_identifier = row[0]
    long = SequencingFile(file_path=(row[1]))
    short_left = SequencingFile(file_path=(row[2]))
    short_right = SequencingFile(file_path=(row[3]))

    sample = Sample(long, short_left, short_right, integrity_check=False)
    sample.identifier = sample_identifier

    return sample
