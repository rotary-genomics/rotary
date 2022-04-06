#!/usr/bin/env bash
set -euo pipefail
# flye_end_repair.sh
# Copyright Jackson M. Tsuji, ILTS, Hokkaido University, 2022

# Internal variables
readonly VERSION="0.1.0"
readonly SCRIPT_NAME="${0##*/}"
readonly SCRIPT_DIR=$(realpath "${0%/*}")

#######################################
# Sorts a FastA file in the specified order
# Arguments:
#   fasta: FastA file
#   sort_order: list of contig names for sorting
#   tmp_dir: directory for temp files. Cannot already exist.
# Returns:
#   STDOUT of the sorted FastA file, 60 nt per line
# Notes:
#   If sort_order is missing entries in fasta, then those entries are not returned
#   If sort_order has entries not in the fasta, then those entries are not returned
#######################################
function sort_fasta {

  local fasta
  fasta=$1
  local sort_order
  sort_order=$2
  local tmp_dir
  tmp_dir=$3

  entry_names=($(cat "${sort_order}"))

  if [[ -d "${tmp_dir}" ]]; then

    echo "[ $(date -u) ]: Tmp dir cannot be created because it already exists: '${tmp_dir}'." | tee -a "${verbose}" >&2
    exit 1

  fi

  mkdir -p "${tmp_dir}"

  for entry in "${entry_names[@]}"; do
    echo "${entry}" > "${tmp_dir}/${entry}.list"
    seqtk subseq -l 60 "${fasta}" "${tmp_dir}/${entry}.list" # to STDOUT
    rm "${tmp_dir}/${entry}.list"
  done

  rmdir "${tmp_dir}"

}

#######################################
# Runs the end repair pipeline
# Globals:
#   SCRIPT_DIR: the directory where this script is stored
# Arguments:
#   all command line inputs for 'main' - see help statement below
# Returns:
#   runs the end repair pipeline end-to-end
#######################################
function run_pipeline() {

  # Parse input
  local qc_long
  qc_long=$1
  local all_contigs
  all_contigs=$2
  local circular_info
  circular_info=$3
  local outdir
  outdir=$4
  local flye_read_mode
  flye_read_mode=$5
  local circlator_min_id
  circlator_min_id=$6
  local circlator_min_length
  circlator_min_length=$7
  local threads
  threads=$8
  local thread_mem
  thread_mem=$9

  # Settings
  local bam_file
  bam_file="${outdir}/long_read.bam"
  local end_repaired_contigs
  end_repaired_contigs="${outdir}/repaired.fasta"
  local verbose
  verbose="${outdir}/verbose.log"

  # Initialize log file
  mkdir "${outdir}"
  printf "" > "${verbose}"

  echo "[ $(date -u) ]: Starting Flye-nano end repair pipeline" | tee -a "${verbose}" >&2

  # Get lists of circular vs linear contigs
  "${SCRIPT_DIR}/flye_end_repair_utils.py" -v -i "${circular_info}" -j "${outdir}/linear_contigs.list" \
    -k "${outdir}/circular_contigs.list" 2>> "${verbose}"

  if [[ $(cat "${outdir}/circular_contigs.list" | wc -l) -eq 0 ]]; then

    echo "[ $(date -u) ]: No circular contigs. Will copy the input file and finish early." | tee -a "${verbose}" >&2

    # Clean up temp files
    rm "${outdir}/circular_contigs.list" "${outdir}/linear_contigs.list"

    cp "${all_contigs}" "${end_repaired_contigs}"
    exit 0

  fi

  echo "[ $(date -u) ]: Mapping reads to all contigs" | tee -a "${verbose}" >&2
  # TODO - add support for different flags like -ax for pacbio
  minimap2 -t "${threads}" -ax map-ont "${all_contigs}" "${qc_long}" 2>> "${verbose}" | \
    samtools view -b -@ "${threads}" 2>> "${verbose}" | \
    samtools sort -@ "${threads}" -m "${thread_mem}G" 2>> "${verbose}" \
    > "${bam_file}"
  samtools index -@ "${threads}" "${bam_file}" 2>> "${verbose}"

  local circular_contigs
  circular_contigs=($(cat "${outdir}/circular_contigs.list"))
  local failed_contigs
  failed_contigs=0

  mkdir -p "${outdir}/circlator_logs"
  printf "" > "${end_repaired_contigs}"

  for contig in "${circular_contigs[@]}"; do

    mkdir -p "${outdir}/contigs/${contig}"
    echo "[ $(date -u) ]: End repair: '${contig}'" | tee -a "${verbose}" >&2

    local assembly
    assembly="${outdir}/contigs/${contig}/${contig}.fasta"
    local linked_ends
    linked_ends="false"

    echo "${contig}" > "${outdir}/contigs/${contig}/${contig}.list" >&2
    seqtk subseq -l 60 "${all_contigs}" "${outdir}/contigs/${contig}/${contig}.list" > "${assembly}" 2>> "${verbose}"

    # TODO - allow user to set the length cutoff range
    for length_cutoff in 100000 90000 80000 70000 60000 50000; do

      local lenout
      lenout="${outdir}/contigs/${contig}/L${length_cutoff}"
      local regions
      regions="${lenout}/ends.bed"
      local fastq
      fastq="${lenout}/ends.fastq.gz"
      local reassembly_dir
      reassembly_dir="${lenout}/assembly"
      local merge_dir
      merge_dir="${lenout}/merge"

      mkdir -p "${lenout}"

      echo "[ $(date -u) ]: End repair: '${contig}': attempting at ${length_cutoff} bp from ends" | \
        tee -a "${verbose}" >&2
      "${SCRIPT_DIR}/flye_end_repair_utils.py" -v -i "${circular_info}" -n "${contig}" -l "${length_cutoff}" \
        -b "${regions}" 2>> "${verbose}"

      samtools view -@ "${threads}" -L "${regions}" -b "${bam_file}" 2>> "${verbose}" | \
        samtools fastq -0 "${fastq}" -n -@ "${threads}" 2>> "${verbose}"

      flye "--${flye_read_mode}" "${fastq}" -o "${reassembly_dir}" -t "${threads}" >> "${verbose}" 2>&1
      # TODO: add support for optional flags like --meta

      mkdir -p "${merge_dir}"
      circlator merge --verbose --min_id "${circlator_min_id}" --min_length "${circlator_min_length}" \
        "${assembly}" "${reassembly_dir}/assembly.fasta" "${merge_dir}/merge" >> "${verbose}" 2>&1
      # TODO - consider exposing more flags like the threshold for aligning the ends

      local pass_circularization
      pass_circularization=$("${SCRIPT_DIR}/flye_end_repair_utils.py" -c "${merge_dir}/merge.circularise.log")

      mkdir -p "${outdir}/contigs/${contig}/logs/L${length_cutoff}"
      cp "${merge_dir}/merge.circularise.log" "${merge_dir}/merge.circularise_details.log" \
        "${reassembly_dir}/assembly_info.txt" "${outdir}/contigs/${contig}/logs/L${length_cutoff}"

      # If the contig was circularized correctly, exit early
      if [[ "${pass_circularization}" == "true" ]]; then

        echo "[ $(date -u) ]: End repair: '${contig}': successfully linked contig ends" | tee -a "${verbose}" >&2

        seqtk seq -l 60 "${merge_dir}/merge.fasta" >> "${end_repaired_contigs}"
        cp "${merge_dir}/merge.circularise_details.log" "${outdir}/circlator_logs/${contig}.log"
        rm -r "${outdir:?}/contigs/${contig}"

        linked_ends="true"
        break

      fi

    done

    # If loop finishes without successful circularization
    if [[  "${linked_ends}" == "false" ]]; then

      echo "[ $(date -u) ]: End repair: '${contig}': FAILED to linked contig ends" | tee -a "${verbose}" >&2

      failed_contigs=$((failed_contigs + 1))

      mkdir -p "${outdir}/troubleshooting"
      mv "${outdir}/${contig}/logs" "${outdir}/troubleshooting/${contig}"
      rm -r "${outdir:?}/contigs/${contig}"

    fi

  done

  if [[ ${failed_contigs} -gt 0 ]]; then
    printf "[ $(date -u) ]: End repair: %s contigs could not be circularized." "${failed_contigs}" | tee -a "${verbose}" >&2
    printf " A partial output file including successfully circularized contigs (and no linear " | tee -a "${verbose}" >&2
    pritnf "contigs) is available at '%s' for debugging.\n" "${end_repaired_contigs}" | tee -a "${verbose}" >&2
    printf "Exiting with error status. See temporary files and verbose log for more details.\n" | tee -a "${verbose}" >&2
    exit 1
  fi

  echo "[ $(date -u) ]: Sorting final FastA file and adding linear contigs" | tee -a "${verbose}" >&2
  seqtk subseq -l 60 "${all_contigs}" "${outdir}/linear_contigs.list" >> "${end_repaired_contigs}"

  grep "^>" "${all_contigs}" | cut -d ">" -f 2- | cut -d " " -f 1 > "${outdir}/sort_order.list"
  sort_fasta "${end_repaired_contigs}.tmp" "${outdir}/sort_order.list" "${outdir}/sort" > "${end_repaired_contigs}"

  # Clean up temp files
  rm "${bam_file}" "${bam_file}.bai" "${outdir}/circular_contigs.list" "${outdir}/linear_contigs.list" \
    "${end_repaired_contigs}.tmp"
  rmdir "${outdir}/contigs"

  echo "[ $(date -u) ]: End repair finished. Output contigs saved at '${end_repaired_contigs}'." | tee -a "${verbose}" >&2

}

function main() {

  # If no input is provided, provide help and exit
  if [[ $# -eq 0 ]]; then
    echo "No arguments provided. Please run '-h' or '--help' to see help. Exiting..." >&2
    exit 1
  elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then

    # Help statement
    printf "%s: pipeline to repair ends of circular contigs from Flye.\n" "${SCRIPT_NAME}"
    printf "Copyright Jackson M. Tsuji, ILTS, Hokkaido University, 2022\n"
    printf "Version: %s\n\n" "${VERSION}"
    printf "Usage: %s [OPTIONS] longreads.fastq.gz assembly.fasta assembly_info.txt outdir\n\n" "${SCRIPT_NAME}"
    printf "Positional arguments:\n"
    printf "   longreads.fastq.gz:      QC-passing Nanopore reads\n"
    printf "   assembly.fasta:          Contigs output from Flye\n"
    printf "   assembly_info.txt:       Assembly info file from Flye\n"
    printf "   outdir:                  Output directory path (directory must not yet exist)\n\n"
    printf "Optional arguments:\n"
    printf "   -f flye_read_mode:       Type of input reads for Flye [nano-hq; can also set nano-raw]\n"
    printf "   -i circlator_min_id:     Percent identity threshold for circlator merge [99]\n"
    printf "   -l circlator_min_length: Required overlap (bp) between original and merge contigs [10000]\n"
    printf "   -t threads:              Number of processors to use [1]\n"
    printf "   -m memory:               Memory (GB) to use per thread for samtools sort [1]\n\n"

    # Exit
    exit 0
  fi

  # Set defaults for options
  local flye_read_mode
  flye_read_mode="nano-hq"
  local circlator_min_id
  circlator_min_id=99
  local circlator_min_length
  circlator_min_length=10000
  local threads
  threads=1
  local thread_mem
  thread_mem=1

  # Set options (help from https://wiki.bash-hackers.org/howto/getopts_tutorial; accessed March 8th, 2019)
  OPTIND=1 # reset the OPTIND counter just in case
  while getopts ":f:i:l:t:m:" opt; do
    case ${opt} in
      f)
        flye_read_mode=${OPTARG}
        ;;
      i)
        circlator_min_id=${OPTARG}
        ;;
      l)
        circlator_min_length=${OPTARG}
        ;;
      t)
        threads=${OPTARG}
        ;;
      m)
        thread_mem=${OPTARG}
        ;;
      \?)
        echo "[ $(date -u) ]: ERROR: Invalid option: '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
      :)
        echo "[ $(date -u) ]: ERROR: argument needed following '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
    esac
  done

  # Set positional arguments
  local original_arguments
  original_arguments="$*" # save for reporting later
  shift $((OPTIND - 1)) # shift to avoid flags when assigning positional arguments
  local qc_long
  qc_long=$1
  local all_contigs
  all_contigs=$2
  local circular_info
  circular_info=$3
  local outdir
  outdir=$4

  # TODO - consider printing to logfile (so move outdir setup to here)
  echo "[ $(date -u) ]: Running ${SCRIPT_NAME}" >&2
  echo "[ $(date -u) ]: Command run: ${SCRIPT_NAME} ${original_arguments}" >&2

  run_pipeline "${qc_long}" "${all_contigs}" "${circular_info}" "${outdir}" \
    "${flye_read_mode}" "${circlator_min_id}" "${circlator_min_length}" \
    "${threads}" "${thread_mem}"

}

# Only run the script if it is called from the command line
if [[ ${BASH_SOURCE[0]} = "${0}" ]]; then
  main "$@"
fi
