#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2023)

Description: Classes and methods for representing Rotary sampling files.
"""
import os
import re

from rotary.utils import is_fastq_file

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

    def __repr__(self):
        """
        Give a useful string representation of this object
        """
        return (f'SequencingFile(path={self.path!r}, name={self.name!r}, '
                f'identifier={self.identifier!r}, paired={self.paired!r}, '
                f'contains_left_reads={self.contains_left_reads!r}, '
                f'contains_right_reads={self.contains_right_reads!r}, '
                f'r_value={self.r_value!r})')


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

    def __repr__(self):
        """
        Give a useful string representation of this object
        """
        return (f'SequencingFilePair(short_read_file_left={self.short_read_file_left!r}, '
                f'short_read_file_right={self.short_read_file_right!r})')

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
                f'Left ({short_file_left.path}) and right ({short_file_right.path}) '
                f'short-read files must have different R designators.')

        return True


class Sample(object):
    """
    An object representing a series of sequencing files representing a single sample.
    """

    def __init__(self, *files):
        self.files = files

    def __repr__(self):
        return f"Sample(identifier={self.identifier}, file_paths={self.sample_file_paths})"

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

    @identifier.setter
    def identifier(self, new_identifier):
        """
        Sets a new identifier for all the files.
        This should be used with care as it modifies the original SequencingFile objects.
        """
        for file in self.files:
            file.identifier = new_identifier

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

    def __repr__(self):
        return (f"ShortReadSample(identifier={self.identifier}, short_read_left_path={self.short_read_left_path}, "
                f"short_read_right_path={self.short_read_right_path})")

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

    def __repr__(self):
        return f"LongReadSample(identifier={self.identifier}, long_read_path={self.long_read_path})"

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
            self.check_file_identifiers_match()  # Checks that all the file identifiers match.

        if integrity_check:
            self.check_long_integrity()  # Checks the long read file is not paired.
            self.check_short_integrity()  # Checks the integrity of the short read files (they are paired and opposites)

    def __repr__(self):
        return (f"LongReadSampleWithPairedPolishingShortReads(identifier={self.identifier}, "
                f"long_read_path={self.long_read_path}, short_read_left_path={self.short_read_left_path}, "
                f"short_read_right_path={self.short_read_right_path})")


def auto_create_sample_from_files(*sequencing_files, identifier_check=False, integrity_check=False):
    """
    This method automatically creates a sample object based on the provided sequencing files.
    The number of sequencing files determines the type of sample object created.

    :param sequencing_files: List of sequencing files.
    :param identifier_check: Boolean value indicating whether to perform an identifier check.
    :param integrity_check: Boolean value indicating whether to perform an integrity check.
    :return: A sample object created from the sequencing files.

    - If there is only one sequencing file, a LongReadSample object is created.
    - If there are two sequencing files, a ShortReadSample object is created.
    - If there are three sequencing files, a LongReadSampleWithPairedPolishingShortReads object is created.

    If the number of sequencing files is not 1, 2, or 3, a ValueError is raised.

    Example usage:

    auto_create_sample_from_files(file1, file2, identifier_check=True, integrity_check=False)
    """
    if len(sequencing_files) == 1:
        sample = LongReadSample(sequencing_files[0], integrity_check=integrity_check)
    elif len(sequencing_files) == 2:
        sample = ShortReadSample(sequencing_files[0], sequencing_files[1], identifier_check=identifier_check,
                                 integrity_check=integrity_check)
    elif len(sequencing_files) == 3:
        sample = LongReadSampleWithPairedPolishingShortReads(sequencing_files[0], sequencing_files[1],
                                                             sequencing_files[2], identifier_check=identifier_check,
                                                             integrity_check=integrity_check)
    else:
        raise ValueError("Invalid number of sequencing files for the sample.")
    return sample
