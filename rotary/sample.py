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

        if 'fastq' not in self.name:
            raise ValueError('{} is not a fastq file.'.format(self.name))

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

    def __init__(self, long_file: SequencingFile, short_file_one: SequencingFile, short_file_two: SequencingFile):
        sequencing_files = [long_file, short_file_one, short_file_two]
        identifiers = [file.identifier for file in sequencing_files]

        # The identifiers should all be the same. If not raise and exception.
        if len(set(identifiers)) == 1:
            self.identifier = identifiers[0]
        else:
            raise ValueError('Sample identifiers of the input fastq files do not match: {}'.format(identifiers))

        # Check for R1 or R2 in the file names with identifier removed.
        r_matches = [file.r_value for file in sequencing_files]
        long_r_value, short_one_r_value, short_two_r_value = r_matches

        long_path = long_file.path
        short_path_one = short_file_one.path
        short_path_two = short_file_two.path

        # The long read file should not have a R1 or R2.
        if not long_r_value:
            self.long_read_path = long_path
        else:
            raise ValueError("The long-read file ({}) should not have an R designator.".format(long_path))

        # The short read files should both have R1s or R2s.
        if not short_one_r_value or not short_two_r_value:
            raise ValueError("Short-read file '{}' or '{}' is missing its R designator.".format(short_path_one,
                                                                                                short_path_two))

        # The short read files should not have the same R1 or R2.
        if short_one_r_value == short_two_r_value:
            matching_r_values_message = 'Left ({}) and right ({}) short-read files must have different R designators.'
            raise ValueError(matching_r_values_message.format(short_path_one, short_path_two))

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
