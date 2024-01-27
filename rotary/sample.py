#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: Classes and methods for representing Rotary sampling files.
"""
import csv
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
            identifier, remaining = self.name.split('.', 1)  # Split on the file extension . if _ not found.

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


def make_sample_from_sample_tsv_row(row, integrity_check=False):
    """
    Parses a row from a sample TSV file and returns a Sample object.

    :param row: A row from a sample tsv file as parsed by the CSV module.
    :param integrity_check: Perform checking that ensures that each sample object maps to
                            files that are all from the same physical sample.
    :return: A Sample object representing the row.
    """
    sample_identifier = row[0]
    long = SequencingFile(file_path=(row[1]))
    short_left = SequencingFile(file_path=(row[2]))
    short_right = SequencingFile(file_path=(row[3]))

    sample = Sample(long, short_left, short_right, integrity_check=integrity_check)
    sample.identifier = sample_identifier

    return sample


def parse_sample_tsv(sample_tsv_path, integrity_check=False):
    """
    Parses a sample tsv file and returns a dictionary of sample objects representing each sample.

    :param sample_tsv_path: The path to the sample TSV file.
    :param integrity_check: Perform checking that ensures that the sample object maps
                            files that are all from the same physical sample.
    :return: A dictionary of all the sample identifiers mapped to sample objects.
    """
    sample_dict = {}
    with open(sample_tsv_path) as sample_file:
        tsv_reader = csv.reader(sample_file, delimiter="\t")
        next(tsv_reader)  # Skip header row.
        for row in tsv_reader:
            sample = make_sample_from_sample_tsv_row(row, integrity_check=integrity_check)
            sample_dict[sample.identifier] = sample

    return sample_dict


def create_sample_tsv(output_dir_path, samples):
    """
    Generates a TSV file in the output directory with a series of CLI paths for files belonging to each sample.

    :param output_dir_path: The path to the output Rotary directory.
    :param samples: A list of Sample objects.
    """
    sample_tsv_path = os.path.join(output_dir_path, 'samples.tsv')
    with open(sample_tsv_path, 'w') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        header = ['sample_id', 'long-read', 'short-read_R1', 'short-read_R2']
        tsv_writer.writerow(header)
        for current_sample in samples:
            tsv_writer.writerow(current_sample.sample_file_row)

    return sample_tsv_path

def file_is_gzipped(file_path):
    """
    Determine if a file is gzipped based in the file extension.

    :param file_path: The path to the file to be checked.
    :return: True if the file is gzipped, False otherwise.
    """
    extension = os.path.splitext(file_path)[1]

    if extension == '.gz':
        return True
    else:
        return False


def find_samples_in_fastq_directory(input_path):
    """
    Find samples in a directory containing fastq files.

    :param input_path: Path to the directory containing fastq files.
    :return: A list of sample objects representing the samples found in the directory.
    :raises ValueError: If a sample does not have exactly three sequencing files.
    """
    fastq_files = get_fastq_files_in_directory(input_path)

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
    return samples


def get_fastq_files_in_directory(input_path):
    """
    Get a list of fastq files in a given directory.

    :param input_path: The path to the directory containing the fastq files.
    :return: A list of SequencingFile objects representing the fastq files.
    """
    fastq_files = []
    for file_path in os.listdir(input_path):
        filename = os.path.basename(file_path)
        if is_fastq_file(filename):
            fastq_files.append(SequencingFile(file_path=os.path.join(input_path, file_path)))
    return fastq_files
