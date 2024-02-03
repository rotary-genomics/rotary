"""
Created by: Lee Bergstrand (2023)

Description: Classes and methods for representing a set of sampling files.
"""
import csv
import os

from rotary.sample import SequencingFile, Sample, auto_create_sample_from_files
from rotary.utils import is_fastq_file


class Dataset:
    """
    An object representing a series of samples.
    """

    def __init__(self, *samples: Sample):
        sample_dict = {}
        for sample in samples:
            sample_dict[sample.identifier] = sample

        self.sample_dict = sample_dict

    @property
    def samples(self):
        """
        Returns a list of all samples.

        :return: A list containing all samples.
        """
        return list(self.sample_dict.values())

    @property
    def identifiers(self):
        """
        Get the list of identifiers from the sample dictionary.

        :return: A list of identifiers from the sample dictionary.
        """
        return list(self.sample_dict.keys())

    def __getitem__(self, sample_id):
        """
        Retrieve an item from the sample dictionary using the given sample_id.

        :param sample_id: The ID of the sample to retrieve.
        :return: The item associated with the sample_id.
        """
        return self.sample_dict[sample_id]

    def __setitem__(self, sample_id, sample):
        """
        Method to add a sample to the sample_dict.

        :param sample_id: The unique identifier for the sample.
        :param sample: A Sample object.
        :return: None

        """
        self.sample_dict[sample_id] = sample

    def __delitem__(self, sample_id):
        """
        Deletes a sample from the sample dict.

        :param sample_id: The ID of the sample to be deleted from the sample_dict
        """
        del self.sample_dict[sample_id]

    def __iter__(self):
        """
        Returns an iterator of all samples.
        """
        return iter(self.sample_dict.values())

    def create_sample_tsv(self, output_dir_path, header):
        """
        Generates a TSV file in the output directory with a series of CLI paths for files belonging to each sample
        in the dataset.

        :param header: The header to use in the sample TSV
        :param output_dir_path: The path to the output Rotary directory.
        """
        sample_tsv_path = os.path.join(output_dir_path, 'samples.tsv')
        with open(sample_tsv_path, 'w') as tsv_file:
            tsv_writer = csv.writer(tsv_file, delimiter='\t')
            tsv_writer.writerow(header)

            for current_sample in self.samples:
                # Ensure that each row has the same number of fields as the header
                row = [current_sample.identifier] + current_sample.sample_file_paths
                if len(row) < len(header):
                    # Add empty fields to match the length of the header
                    row += [''] * (len(header) - len(row))
                tsv_writer.writerow(row)

        return sample_tsv_path


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

    sample = auto_create_sample_from_files(*sequencing_files, identifier_check=identifier_check,
                                           integrity_check=integrity_check)

    sample.identifier = sample_identifier

    return sample


def generate_dataset_from_sample_tsv(sample_tsv_path, identifier_check=False, integrity_check=False):
    """
    Generates a Dataset from a sample TSV file.

    :param sample_tsv_path: Path to the sample TSV file.
    :param identifier_check: Optional parameter to perform identifier check. Default is False.
    :param integrity_check: Optional parameter to perform integrity check. Default is False.
    :return: A Dataset object containing a series of Sample objects.

    :Example:
    generate_dataset_from_sample_tsv("path/to/sample.tsv", identifier_check=True, integrity_check=True)
    """
    dataset = Dataset()
    with open(sample_tsv_path) as sample_file:
        tsv_reader = csv.reader(sample_file, delimiter="\t")
        next(tsv_reader)  # Skip header row.
        for row in tsv_reader:
            sample = make_sample_from_sample_tsv_row(row, identifier_check=identifier_check,
                                                     integrity_check=integrity_check)
            dataset[sample.identifier] = sample

    return dataset


def generate_dataset_from_fastq_directory(input_path):
    """
    :param input_path: The path to the directory containing the FASTQ files.
    :return: A dataset generated from the FASTQ files.

    This method generates a Dataset object from the FASTQ files located in the specified directory. The method first
    generates a series of SequencingFile objects representing the FASTQ files using the get_fastq_files_in_directory
    method.

    Then, it organizes the files into groups based on their identifiers using a dictionary called samples_files.
    Each key in the dictionary represents a unique identifier, and the corresponding value is a list of files
    associated with that identifier.

    Next, the method iterates through the samples_files dictionary to process each group of sequencing files.
    It ensures that each group has exactly three files. If not, it raises a ValueError indicating that a sample should
    have three sequencing files.

    For each group of sequencing files, the method identifies the long_file (a file with no R value), the
    left_short_file (typically an R1 file), and the right_short_file (typically an R2 file). It creates a
    Sample object with these files and adds it to the samples list.

    Finally, the method creates a Dataset object using the samples list and returns it.
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

        sample = auto_create_sample_from_files(long_file, left_short_file, right_short_file,
                                               identifier_check=True, integrity_check=True)
        samples.append(sample)
    return Dataset(*samples)


def get_fastq_files_in_directory(input_path):
    """
    Get a list of FASTQ files in a given directory.

    :param input_path: The path to the directory containing the fastq files.
    :return: A list of SequencingFile objects representing the fastq files.
    """
    fastq_files = []
    for file_path in os.listdir(input_path):
        filename = os.path.basename(file_path)
        if is_fastq_file(filename):
            fastq_files.append(SequencingFile(file_path=os.path.join(input_path, file_path)))
    return fastq_files
