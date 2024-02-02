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
    An object representing a single FASTQ sequencing file.
    """

    def __init__(self, file_path):
        self.path = file_path
        self.name = os.path.basename(self.path)

        if not is_fastq_file(self.name):
            raise ValueError(f'{self.name} is not a fastq file.')

        if '_' in self.name:
            identifier, remaining = self.name.split('_', 1)  # Split on first _ if found.
        else:
            identifier, remaining = self.name.split('.', 1)  # Split on the file extension. if _ not found.

        self.identifier = identifier

        # Check the remainder of the file name (e.g. _blam_blam_blam_R1.fastq.gz) for R1 or R2
        r_value = short_read_r_regex.findall(remaining)
        if r_value:
            self.r_value = r_value[0].upper()
            self.paired = True
            if self.r_value == 'R1':
                self.contains_left_reads = True
                self.contains_right_reads = False
            else:
                self.contains_left_reads = False
                self.contains_right_reads = True
        else:
            self.contains_left_reads = False
            self.contains_right_reads = False
            self.r_value = None
            self.paired = False


class SequencingFilePair(object):
    """
    An object representing a set of paired end FASTQ sequencing files.
    """

    def __init__(self, file_one: SequencingFile, file_two: SequencingFile):
        if file_one.contains_left_reads:
            self.short_read_file_left = file_one
            self.short_read_file_right = file_two
        else:
            self.short_read_file_left = file_two
            self.short_read_file_right = file_one

    def check_integrity(self):
        """
        Check the integrity of the short read files.

        :raises ValueError: if the short read files have different R designators.
        """

        short_file_left = self.short_read_file_left
        short_file_right = self.short_read_file_right
        if not short_file_left.paired or not short_file_right.paired:
            raise ValueError(
                f"Short-read file '{short_file_left.path}' or '{short_file_right.path}' is missing its R designator.")

        # The short read files should not have the same R1 or R2.
        if short_file_left.r_value == short_file_right.r_value:
            raise ValueError(
                f'Left ({short_file_left.path}) and right ({short_file_right.path}) short-read files must have different R designators.')

        return True


class Sample(object):
    """
    An object representing a series of sequencing files representing a single sample.
    """

    def __init__(self, *files):
        self.files = files

    @property
    def sample_file_paths(self):
        """
        Returns a list containing the paths of the files.

        :return: A list of strings representing the paths of the files.
        """
        return [file.path for file in self.files]

    @property
    def identifier(self):
        """
        Checks the first file to get the sample's identifier.

        :return: The sample's identifier.
        """
        return self.files[0].identifier

    def check_file_identifiers_match(self):
        """
        Check if the sample identifiers of the input FASTQ files match.

        :raises ValueError: If the sample identifiers do not match.
        """
        identifiers = [file.identifier for file in self.files]
        if not len(set(identifiers)) == 1:
            raise ValueError(f'Sample identifiers of the input fastq files do not match: {identifiers}')
        else:
            return True


class ShortReadSample(Sample):
    """
    An object representing a pair of short read sequencing files representing a single sample.
    """

    def __init__(self, short_file_one: SequencingFile, short_file_two: SequencingFile, identifier_check: bool = True,
                 integrity_check: bool = True):

        self.short_read_pair = SequencingFilePair(short_file_one, short_file_two)

        self.short_read_file_left = self.short_read_pair.short_read_file_left
        self.short_read_file_right = self.short_read_pair.short_read_file_right

        super().__init__(self.short_read_file_left, self.short_read_file_right)

        if identifier_check:
            self.check_file_identifiers_match()

        if integrity_check:
            self.check_short_integrity()

    def check_short_integrity(self):
        """
        Check the integrity of the short read pair.
        """
        self.short_read_pair.check_integrity()

    @property
    def short_read_left_path(self):
        """
        Returns the path of the left short read file.

        :return: The path of the left short read file.
        """
        return self.short_read_file_left.path

    @property
    def short_read_right_path(self):
        """
        Returns the path of the short read file on the right side.

        :return: The path of the short read file on the right side.
        """
        return self.short_read_file_right.path


class LongReadSample(Sample):
    def __init__(self, long_read_file: SequencingFile, integrity_check: bool = True):
        self.long_read_file = long_read_file

        super().__init__(self.long_read_file)

        if integrity_check:
            self.check_long_integrity()

    @property
    def long_read_path(self):
        """
        Return the path of the long read file.

        :return: The path of the long read file.
        """
        return self.long_read_file.path

    def check_long_integrity(self):
        """
        Check the integrity of the long-read file.

        :raises ValueError: If the long-read file has an R designator (R1 or R2) in its name.
        """
        if self.long_read_file.r_value:
            raise ValueError(f"The long-read file ({self.long_read_file.path}) should not have an R designator.")
        else:
            return True


class LongReadSampleWithPairedPolishingShortReads(LongReadSample, ShortReadSample):
    """
    An object representing a long read sequence with paired short reads for polishing.

    We use multiple inheritance here to combine the attributes and methods of both LongReadSample and ShortReadSample.
    """

    def __init__(self, long_file: SequencingFile, short_file_left: SequencingFile, short_file_right: SequencingFile,
                 identifier_check: bool = True, integrity_check: bool = True):

        self.long_read_file = long_file
        self.short_read_pair = SequencingFilePair(short_file_left, short_file_right)
        self.short_read_file_left = self.short_read_pair.short_read_file_left
        self.short_read_file_right = self.short_read_pair.short_read_file_right

        self.files = [self.long_read_file, self.short_read_file_left, self.short_read_file_right]

        if identifier_check:
            self.check_file_identifiers_match() # Checks that all the file identifiers match.

        if integrity_check:
            self.check_long_integrity()  # Checks the long read file is not paired.
            self.check_short_integrity()  # Checks the integrity of the short read files (they are paired and opposites)


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


def make_sample_from_sample_tsv_row(row, identifier_check=False, integrity_check=False,
                                    first_meta_column_position=None):
    """
    Parses a row from a sample TSV file and returns a Sample object.

    :param identifier_check: Check if the sample identifiers of the input FASTQ files match.
    :param first_meta_column_position: The column number where metadata starts.
    :param row: A row from a sample tsv file as parsed by the CSV module.
    :param integrity_check: Perform checking that ensures that each sample object maps to
                            files that are all from the same physical sample.
    :return: A Sample object representing the row.
    """
    sample_identifier = row[0]
    if first_meta_column_position:
        sample_paths = row[1:first_meta_column_position]
    else:
        sample_paths = row[1:]

    sequencing_files = [SequencingFile(file_path=path) for path in sample_paths if path]

    if len(sequencing_files) == 1:
        sample = LongReadSample(sequencing_files[0], integrity_check=integrity_check)
    elif len(sample_paths) == 2:
        sample = ShortReadSample(sequencing_files[0], sequencing_files[1], identifier_check=identifier_check,
                                 integrity_check=integrity_check)
    elif len(sample_paths) == 3:
        sample = LongReadSampleWithPairedPolishingShortReads(sequencing_files[0], sequencing_files[1],
                                                             sequencing_files[2], identifier_check=identifier_check,
                                                             integrity_check=integrity_check)
    else:
        raise ValueError("Invalid number of sequencing files for the sample.")

    sample.identifier = sample_identifier

    return sample


def parse_sample_tsv(sample_tsv_path, identifier_check=False, integrity_check=False):
    """
    Parses a sample tsv file and returns a dictionary of sample objects representing each sample.

    :param identifier_check: Check if the sample identifiers of the input FASTQ files match.
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
            sample = make_sample_from_sample_tsv_row(row, identifier_check=identifier_check,
                                                     integrity_check=integrity_check)
            sample_dict[sample.identifier] = sample

    return sample_dict


def create_sample_tsv(output_dir_path, samples, header):
    """
    Generates a TSV file in the output directory with a series of CLI paths for files belonging to each sample.
    :param header: The header to use in the sample TSV
    :param output_dir_path: The path to the output Rotary directory.
    :param samples: A list of Sample objects.
    """
    sample_tsv_path = os.path.join(output_dir_path, 'samples.tsv')
    with open(sample_tsv_path, 'w') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerow(header)
        for current_sample in samples:
            # Ensure that each row has the same number of fields as the header
            row = current_sample.sample_file_paths
            if len(row) < len(header):
                # Add empty fields to match the length of the header
                row += [''] * (len(header) - len(row))
            tsv_writer.writerow(row)
    return sample_tsv_path


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
