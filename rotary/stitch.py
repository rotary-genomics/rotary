#!/usr/bin/env python
# stitch.py
# Aligns and stitches ends of circular contigs using an end-spanning guide contig
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2023

import argparse
import logging
import os
import subprocess
import sys
import itertools

import pandas as pd
from Bio import SeqIO

from rotary.repair import check_dependency, run_pipeline_subcommand, rotate_contig_to_midpoint, set_write_mode

# GLOBAL VARIABLES
DEPENDENCY_NAMES = ['nucmer', 'show-coords']

# Set up the logger
logger = logging.getLogger(__name__)
formatter = logging.Formatter('[ %(asctime)s ]: %(levelname)s: %(funcName)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


def run_nucmer(query_filepath: str, ref_filepath: str, output_delta_filepath: str,
               output_coords_filepath: str, log_filepath: str, min_match_length: int = 1, append_log: bool = True,
               threads: int = 1):
    """
    Maps long reads (via minimap2) to contigs and sorts/indexes the resulting BAM file
    :param query_filepath: path to the FastA file containing the query contig
    :param ref_filepath: path to the FastA file containing the reference contig
    :param output_delta_filepath: path where the output .delta file should be saved
    :param output_coords_filepath: path where the output .coords file should be saved
    :param log_filepath: path to the log file to be saved
    :param min_match_length: minimum match length (in bp) between the query and the reference for nucmer to report
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :return: None
    """

    write_mode = set_write_mode(append_log)

    with open(log_filepath, write_mode) as logfile_handle:
        nucmer_args = ['nucmer', '--threads', str(threads), '--minmatch', min_match_length, '--delta',
                       output_delta_filepath, ref_filepath, query_filepath]
        run_pipeline_subcommand(command_args=nucmer_args, stderr=logfile_handle)

        with open(output_coords_filepath, 'w') as coords_handle:
            show_coords_args = ['show-coords', '-r', '-T', '-l', output_delta_filepath]
            run_pipeline_subcommand(command_args=show_coords_args, stdout=coords_handle, stderr=logfile_handle)

    logger.debug('nucmer finished')

    # TODO - delete this command line reference once script is confirmed to work
    """
    nucmer -l 1 -p test input/processed/genome_subset.fna output1/repaired.fasta 
    show-coords -r test.delta > test.coords
    """


def run_show_snps(delta_filepath: str, coords_subset_filepath: str, output_snps_filepath: str,
                  log_filepath: str, append_log: bool = True):
    """
    Maps long reads (via minimap2) to contigs and sorts/indexes the resulting BAM file
    :param delta_filepath: path to the FastA file containing the query contig
    :param coords_subset_filepath: path to the FastA file containing the reference contig
    :param output_snps_filepath: path where the output .snps file should be saved
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :return: None
    """

    write_mode = set_write_mode(append_log)

    with open(log_filepath, write_mode) as logfile_handle:
        cat_args = ['cat', coords_subset_filepath]
        cat = run_pipeline_subcommand(command_args=cat_args, stdout=subprocess.PIPE, stderr=logfile_handle)

        with open(output_snps_filepath, 'w') as snps_handle:
            show_snps_args = ['show-snps', '-r', '-T', '-S', delta_filepath]
            run_pipeline_subcommand(command_args=show_snps_args, stdin=cat.stdout, stdout=snps_handle,
                                    stderr=logfile_handle)

    logger.debug('nucmer finished')

    # TODO - delete this command line reference once script is confirmed to work
    """
    cat subset.coords | show-snps -r -T -S test.delta > subset.snps
    """


def rename_identical_entries(entry_list: list, spacer: str = '-') -> list:
    """
    Looks for identical entries in a list and renames them to have a number counter on the end
    :param entry_list: list of entries to try renaming
    :param spacer: the spacer to put inbetween the end of the entry and the counter, for duplicate entries
    :return: list of all input entries, with duplicates renamed, in identical order to the input list
    """

    entry_names_revised = []

    for entry_name in entry_list:

        # If a duplicated entry is encountered
        if entry_name in entry_names_revised:

            entry_counter = 1
            unique_entry = False

            # Try different ways to rename the entry to make the name unique
            while unique_entry is False:

                entry_name_revised = f'{entry_name}{spacer}{entry_counter}'

                # Test if the revised name is unique
                if entry_name_revised not in entry_names_revised:
                    # TODO - remove this debug entry in future (it is too verbose)
                    logger.debug(f'Renamed {entry_name} to {entry_name_revised}')
                    entry_names_revised.append(entry_name_revised)

                    unique_entry = True

                entry_counter = entry_counter + 1

        else:
            entry_names_revised.append(entry_name)

    return entry_names_revised


def parse_mummer_stats_file_header(mummer_stats_file_path, number_identical_headers: bool = True,
                                   header_row_num: int = 3) -> list:
    """
    Parses the column names from a mummer stats file. This is tricky because TAGS includes 2 columns (reference name and
    query name), causing an error in pandas. This parser function assumes that the mummer stats file was still output
    with columns but that they were set to be tab-separated
    :param mummer_stats_file_path: path to the mummer statistics file to parse
    :param number_identical_headers: boolean of whether to add 1,2,3 counters to the end of headers with duplicate
                                     names (True) or not (False). For example, if True, then if two headers are loaded
                                     named '[SUB]' and '[SUB]', they will be renamed '[SUB1]' and '[SUB2]' (in that
                                     order)
    :param header_row_num: the row number (0-ordered) of the header in the file. In theory, you should only need to
                           change this param if the mummer file format changes. Assumes that the data immediately
                           follows the header row.
    :return: list of header names suitable for pandas
    """

    # Example stats file (.coords)
    """
    /Git/rotary/tests/repair/input/processed/genome_subset.fna /Git/rotary/tests/repair/output1/repaired.fasta
    NUCMER

    [S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[TAGS]
    1	153979	87413	241410	153979	153998	99.94	241389	241410	CP128402	CP128402
    153980	241389	1	87412	87410	87412	99.85	241389	241410	CP128402	CP128402
    """

    line_num = 0
    with open(mummer_stats_file_path) as input_handle:
        for line in input_handle:

            if line_num == header_row_num:
                header_entry = line.rstrip()
                break

            line_num = line_num + 1

    header_names = header_entry.split('\t')

    if len(header_names) == 1:
        logger.warning(f'Only 1 header was found in the mummer stats file: {header_names[0]}')
        logger.warning('This is concerning because it suggests that this function was unable to parse the header '
                       'names. Perhaps they are tab-separated, for example?')

    if len(set(header_names)) < len(header_names):
        if number_identical_headers is True:
            header_names = rename_identical_entries(header_names)

        elif number_identical_headers is False:
            logger.warning('At least two identical header names were found. Set number_identical_headers to true if '
                           'you want to add counters to these to make them distinct')

        else:
            logger.error(f'number_identical_headers must be True or False, but you specified '
                         f'{number_identical_headers}')
            raise ValueError

    # Add a custom header for the second TAG column
    header_names.append('[TAGS2]')

    return header_names


def parse_snps_file(snps_filepath, header_row_num: int = 3) -> pd.DataFrame:
    """
    Loads a .snps file from show-snps as a pandas DataFrame
    :param snps_filepath: path to the .snps file to be loaded
    :param header_row_num: the row number (0-ordered) of the header in the file. In theory, you should only need to
                           change this param if the snps file format changes. Assumes that the data immediately
                           follows the header row.
    :return:
    """

    # TODO - delete this example snippet once the code is working
    """
    /Git/rotary/tests/repair/input/processed/genome_subset.fna /Git/rotary/tests/repair/output1/repaired.fasta
    NUCMER
    
    [P1]	[SUB]	[SUB]	[P2]	[BUFF]	[DIST]	[R]	[Q]	[FRM]	[TAGS]
    627	G	A	88039	627	627	0	0	1	1	CP128402	CP128402
    1682	G	.	89093	329	1682	0	0	1	1	CP128402	CP128402
    2011	G	.	89421	329	2011	0	0	1	1	CP128402	CP128402
    2874	.	A	90285	550	2874	0	0	1	1	CP128402	CP128402
    3424	.	A	90836	550	3424	0	0	1	1	CP128402	CP128402
    5937	.	G	93350	2513	5937	0	0	1	1	CP128402	CP128402
    """

    # Get the column names from the file (this is tricky because TAGS includes 2 columns, causing an error in pandas)
    header_names = parse_mummer_stats_file_header(snps_filepath, header_row_num=header_row_num)

    # Now import the data, add on the columns, and transform them to more meaningful names
    snps_data = pd.read_csv(snps_filepath, sep='\t', skiprows=range(0, header_row_num + 1), header=None)
    snps_data.columns = header_names
    snps_data = snps_data \
        .rename(columns={'[P1]': 'ref-coord', '[SUB]': 'ref-base', '[SUB]-1': 'query-base', '[P2]': 'query-coord',
                         '[BUFF]': 'bp-to-nearest-mismatch', '[DIST]': 'bp-to-nearest-seq-end',
                         '[R]': 'n-repeat-alns-ref', '[Q]': 'n-repeat-alns-query', '[FRM]': 'reading-frame',
                         '[TAGS]': 'rseqid', '[TAGS2]': 'qseqid'}) \
        .drop(columns=['reading-frame'])

    return snps_data


def parse_coords_file(coords_filepath, add_hit_id: bool = True, header_row_num: int = 3) -> pd.DataFrame:
    """
    Loads a .coords file from show-coords as a pandas DataFrame
    :param coords_filepath: path to the .coords file to be loaded
    :param add_hit_id: boolean of whether to add a numeric hit ID to the dataframe (True) or not (False)
    :param header_row_num: the row number (0-ordered) of the header in the file. In theory, you should only need to
                           change this param if the coords file format changes. Assumes that the data immediately
                           follows the header row
    :return: pandas DataFrame of the .coords file data
    """

    # Get the column names from the file (this is tricky because TAGS includes 2 columns, causing an error in pandas)
    header_names = parse_mummer_stats_file_header(coords_filepath, header_row_num=header_row_num)

    # Now import the data, add on the columns, and transform them to more meaningful names
    coords_data = pd.read_csv(coords_filepath, sep='\t', skiprows=range(0, header_row_num + 1), header=None)
    coords_data.columns = header_names
    coords_data = coords_data \
        .rename(columns={'[S1]': 'ref-start', '[E1]': 'ref-end', '[S2]': 'query-start', '[E2]': 'query-end',
                         '[LEN 1]': 'ref-aln-len', '[LEN 2]': 'query-aln-len', '[% IDY]': 'pident',
                         '[LEN R]': 'ref-total-len', '[LEN Q]': 'query-total-len', '[TAGS]': 'rseqid',
                         '[TAGS2]': 'qseqid'})

    if add_hit_id is True:
        coords_data = coords_data.reset_index() \
            .rename(columns={'index': 'hit-id'})

    return coords_data


def identify_possible_stitch_hits(coords_data: pd.DataFrame, hit_side: str, query_end_proximity: int,
                                  ref_end_proximity: int) -> list:
    """
    Identifies possible nucmer hits between a query and reference that might be suitable for stitching. Works one
    'side' at a time (for 'sides' of the reference contig).
    :param coords_data: pandas Dataframe of the .coords file from mummer's show-coords, loaded by parse_coords_file.
                        MUST have the 'hit-id' column generated by add_hit_id
    :param hit_side: 'positive' (i.e., near the 1 position on the circular reference) or 'negative' (i.e., near the end
                     of the circular reference, which could be counted as negative from the zero point)
    :param query_end_proximity: minimum distance (bp) between the real end of the query sequence and the end of the
                                nucmer hit for the hit to be considered reliable
    :param ref_end_proximity: minimum distance (bp) between the real end of the reference sequence and the end of the
                              nucmer hit for the hit to be considered reliable
    :return: list of possible hits that could work for the specified reference contig side for stitching
    """

    possible_stitch_hit_ids = []

    if hit_side == 'negative':

        for index, rows in coords_data[['hit-id', 'query-start', 'ref-end', 'ref-total-len']].iterrows():
            hit_id, query_start, ref_end, ref_total_len = rows

            # Here, the hit should span from close to the start of the query to close to the end of reference
            if (query_start <= query_end_proximity) & (ref_total_len - ref_end <= ref_end_proximity):
                possible_stitch_hit_ids.append(hit_id)

    elif hit_side == 'positive':

        for index, rows in coords_data[['hit-id', 'ref-start', 'query-end', 'query-total-len']].iterrows():
            hit_id, ref_start, query_end, query_total_len = rows

            # Here, the hit should span from close to the start of the reference to close to the end of the query
            if (ref_start <= ref_end_proximity) & (query_total_len - query_end <= query_end_proximity):
                possible_stitch_hit_ids.append(hit_id)

    else:
        logger.error(f'hit_side should be "positive" or "negative"; you provided {hit_side}')
        raise ValueError

    logger.debug(f'Identified {len(possible_stitch_hit_ids)} possible {hit_side} side hits')

    return possible_stitch_hit_ids


def find_possible_hit_pairs(coords_data: pd.DataFrame, negative_hit_candidate_ids: list,
                            positive_hit_candidate_ids: list) -> list:
    """
    Using results from identify_possible_stitch_hits, summarize possible negative and positive-side hit pairs for the
    same query and reference contigs that could work for stitching.
    :param coords_data: pandas Dataframe of the .coords file from mummer's show-coords, loaded by parse_coords_file.
                        MUST have the 'hit-id' column generated by add_hit_id
    :param negative_hit_candidate_ids: list of hit IDs for negative-side hits to consider
    :param positive_hit_candidate_ids: list of hit IDs for positive-side hits to consider
    :return: list of tuples (negative hit id, positive hit id) of possible hit pairs to test for stitching
    """

    # Summarize the reference and query sequence IDs associated with each hit
    hit_candidate_metadata = pd.concat([pd.DataFrame({'hit-id': negative_hit_candidate_ids, 'hit-side': 'negative'}),
                                        pd.DataFrame({'hit-id': positive_hit_candidate_ids, 'hit-side': 'positive'})])\
        .merge(coords_data[['hit-id', 'qseqid', 'rseqid']], on='hit-id', how='left', validate='1:1')\
        .sort_values(by=['rseqid', 'qseqid', 'hit-id'], ascending=True)

    # Group the hits based on which reference / query sequence pair they belong to
    ref_query_pair_previous = []
    hit_candidate_group = []
    hit_candidate_group_counter = 0

    for index, row in hit_candidate_metadata.iterrows():
        hit_id, hit_side, qseqid, rseqid = row
        ref_query_pair_current = [rseqid, qseqid]

        # If the current ref ID and query ID don't match the previous one, then assign a new group ID
        # Otherwise, the same group ID will be used as in the last entry
        if ref_query_pair_previous != ref_query_pair_current:
            hit_candidate_group_counter = hit_candidate_group_counter + 1
            ref_query_pair_previous = ref_query_pair_current

        hit_candidate_group.append(hit_candidate_group_counter)

    hit_candidate_metadata['group-id'] = hit_candidate_group

    # Find possible hit pair combinations (negative and positive) for each reference / query pair
    possible_pairs = []
    for group_id in hit_candidate_metadata['group-id'].drop_duplicates():

        group_members = hit_candidate_metadata[hit_candidate_metadata['group-id'] == group_id]
        group_members_negative = list(group_members[group_members['hit-side'] == 'negative']['hit-id'])
        group_members_positive = list(group_members[group_members['hit-side'] == 'positive']['hit-id'])

        # No pairs if there are no entries on the negative or positive side
        if (len(group_members_negative) == 0) | (len(group_members_positive) == 0):
            continue

        possible_pairs_in_group = list(itertools.product(group_members_negative, group_members_positive))
        possible_pairs = possible_pairs.append(possible_pairs_in_group)

    logger.debug(f'{len(possible_pairs)} possible hit pairs identified')

    return possible_pairs


def stitch_contig_ends(circular_contig_filepath: str, guide_contig_filepath: str, output_dir: str,
                       cli_tool_settings_dict: dict, threads: int = 1):
    """

    :param circular_contig_filepath:
    :param guide_contig_filepath:
    :param output_dir:
    :param cli_tool_settings_dict:
    :param threads:
    :return:
    """

    """
    # LOGIC
    1. Run nucmer on the reference vs query (and generate the coords file). Load the coords file
    2. ID possible LHS sequences and RHS sequences (based on proximity to Q-0 (LHS-start), R-L (LHS-end), 
       R-0 (RHS-start), Q-L (RHS-end)
    3. Summarize all possible combinations of LHS:RHS pairs
    
    4. Test each pair for required params
    5. If only one possible pair, then write a summary coords file. Linking was successful (Otherwise error)
    6. Integration:
       a. Cut out the stitch part from the reference (by keeping the remainder): keep from the RHS nucmer end base to
          the LHS nucmer start
       b. Trim ends off the query (by keeping the remainder): keep from the LHS nucmer start to the RHS nucmer end
       c. Paste the trimmed q at the end of the trimmed reference
       d. Calculate the coordinate shift based on the change from reference base 1 to the RHS nucmer end base used 
          above for trimming
       e. Calculate the position (start) of the old crossover point as the negative of the above coordinate shift. 
          Note that this is an approximate value, because the new stitch part might contain indels 
       f. Calculate the position of the two stitch points. One will be at 0; the other will be at the total contig 
          length minus the length of the trimmed query contig.
       g. Determine the optimal point on the contig that is furthest away from the old crossover point and the two 
          stitch points
       h. Rotate to that point. (Warn if any of the ends up being within 10kb of the new crossover point)
       i. Recalculate all positions by shifting them by the rotation amount. Note that they need to reset to 1 after 
          passing the contig length.
    7. Reporting:
       a. Use the summary coords file to generate a SNP summary file (or I could just run nucmer again from scratch)
       b. Load the SNP summary file
       c. Find the calculated stitch region approximately (above)
       d. Search for a gap or overlap (need an algorithm for this). Report.
          Note that the region is approximate, so I should really scan a window OR do a re-alignment with the original
          Or why not just look at the original nucmer hit data? I can look at the bases between the two hit regions or 
          overlap. BUT I might need to load the original contig data
       e. Report the total other number of base changes
       f. Optionally output a cleaned up SNP file
    """

    delta_filepath_initial = os.path.join(output_dir, 'nucmer_initial.delta')
    coords_filepath_initial = os.path.join(output_dir, 'nucmer_initial.coords')

    run_nucmer(query_filepath=guide_contig_filepath, ref_filepath=circular_contig_filepath,
               output_delta_filepath=delta_filepath_initial, output_coords_filepath=coords_filepath_initial,
               log_filepath=os.path.join(output_dir, 'verbose.log'), append_log=False, threads=threads)

    coords_data = parse_coords_file(coords_filepath_initial)

    negative_hit_candidates = identify_possible_stitch_hits(coords_data, hit_side='negative',
                                                            query_end_proximity=cli_tool_settings_dict[
                                                                'circlator_reassemble_end'],
                                                            ref_end_proximity=cli_tool_settings_dict[
                                                                'circlator_ref_end'])

    positive_hit_candidates = identify_possible_stitch_hits(coords_data, hit_side='positive',
                                                            query_end_proximity=cli_tool_settings_dict[
                                                                'circlator_reassemble_end'],
                                                            ref_end_proximity=cli_tool_settings_dict[
                                                                'circlator_ref_end'])

    # List of tuples (negative hit id, positive hit id) of possible hit pairs to test for stitching
    possible_pairs = find_possible_hit_pairs(coords_data, negative_hit_candidate_ids=negative_hit_candidates,
                                             positive_hit_candidate_ids=positive_hit_candidates)

    # TODO - stopped at step 4


def main():
    """
    Collects input arguments and runs the stitch workflow
    """

    """
    When installing through pip with pyprojects.toml a new python script is generated
    that calls main() out of stitch.py. Run parser inside main() so it can be called
    externally as a function.
    """
    parser = parse_cli()
    args = parser.parse_args()

    # Startup checks
    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Check output dir
    output_dir_exists = os.path.isdir(args.output_dir)

    if (output_dir_exists is True) & (args.overwrite is False):
        logger.error(f'Output directory already exists: "{args.output_dir}". Will not continue. Set the '
                     f'--overwrite flag at your own risk if you want to use an existing directory.')
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Check dependencies
    dependency_paths = []
    for dependency_name in DEPENDENCY_NAMES:
        dependency_paths.append(check_dependency(dependency_name))

    dependency_dict = dict(zip(DEPENDENCY_NAMES, dependency_paths))

    # Parse some command line inputs further
    cli_tool_settings_dict = {'circlator_min_id': args.circlator_min_id,
                              'circlator_min_length': args.circlator_min_length,
                              'circlator_ref_end': args.circlator_ref_end,
                              'circlator_reassemble_end': args.circlator_reassemble_end}

    stitch_contig_ends(args.circular_contig_filepath, args.guide_contig_filepath, args.output_dir,
                       cli_tool_settings_dict, args.threads)

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


def parse_cli():
    """
    Parses the CLI arguments
    :return: An argparse parser object
    """

    parser = argparse.ArgumentParser(
        description=f'{os.path.basename(sys.argv[0])}: stitches ends of a circular contig using a guide contig. \n'
                    '  Copyright Jackson M. Tsuji and Lee Bergstrand, 2023',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    required_settings = parser.add_argument_group('Required')
    merge_settings = parser.add_argument_group('Merge options')
    workflow_settings = parser.add_argument_group('Workflow options')

    required_settings.add_argument('-c', '--circular_contig_filepath', required=True, type=str,
                                   help='Contig to be end repaired')
    required_settings.add_argument('-g', '--guide_contig_filepath', required=True, type=str,
                                   help='Guide contig spanning the ends (in the FastA file) of the circular contig')
    required_settings.add_argument('-o', '--output_dir', required=True, type=str, help='Output directory path')

    merge_settings.add_argument('-I', '--circlator_min_id', required=False, default=99, type=float,
                                help='Percent identity threshold for circlator merge (default: 99)')
    merge_settings.add_argument('-L', '--circlator_min_length', required=False, default=10000, type=int,
                                help='Minimum required overlap (bp) between original and merge contigs '
                                     '(default: 10000)')
    merge_settings.add_argument('-e', '--circlator_ref_end', required=False, default=100, type=int,
                                help='Minimum distance (bp) between end of original contig and nucmer hit '
                                     '(default: 100)')
    merge_settings.add_argument('-E', '--circlator_reassemble_end', required=False, default=100, type=int,
                                help='Minimum distance (bp) between end of merge contig and nucmer hit '
                                     '(default: 100)')

    workflow_settings.add_argument('-t', '--threads', required=False, default=1, type=int,
                                   help='Number of processors threads to use (default: 1)')
    workflow_settings.add_argument('-v', '--verbose', required=False, action='store_true',
                                   help='Enable verbose logging')

    return parser


if __name__ == '__main__':
    main()
