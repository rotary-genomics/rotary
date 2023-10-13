#!/usr/bin/env bash
set -euo pipefail
# flye_end_repair.sh
# Copyright Jackson M. Tsuji, ILTS, Hokkaido University, 2022

# Internal variables
readonly VERSION="0.1.0"
readonly SCRIPT_NAME="${0##*/}"

# Get script dir (realpath is needed)
if [[ $(which realpath > /dev/null 2>&1 ; echo $?) -ne 0 ]]; then
  echo "[ $(date -u) ]: ERROR: missing dependency: 'realpath'. Exiting..." >&2
  exit 1
fi
readonly SCRIPT_DIR=$(realpath "${0%/*}")

#######################################
# Checks if a dependency exists
# Arguments:
#   dependency: Name of dependency
#   type: Type of dependency: commandline (for testing via 'which'), or python_package (for testing via 'pip'),
#                             or internal (for testing via SCRIPT_DIR)
# Returns:
#   Exit status 0 if dependency exists; otherwise exit status 1 and error message.
#######################################
function check_dependency {

  local dependency
  dependency=$1
  local type
  type=$2

  if [[ "${type}" == "commandline" ]]; then
    set +e
    success_status=$(which "${dependency}" > /dev/null 2>&1 ; echo $?)
    set -e
  elif [[ "${type}" == "python_package" ]]; then
    set +e
    success_status=$(pip list | grep -w "${dependency}" > /dev/null 2>&1 ; echo $?)
    set -e
  elif [[ "${type}" == "internal" ]]; then
    set +e
    success_status=$(test -e "${SCRIPT_DIR}/${dependency}"; echo $?)
    set -e
  else
    echo "[ $(date -u) ]: ERROR (internal): success_status must be 'commandline' or 'python_package' or 'internal', but '${type}' was provided. Exiting..." >&2
    exit 1
  fi

  if [[ "${success_status}" -ne 0 ]]; then
    echo "[ $(date -u) ]: ERROR: missing dependency: '${dependency}'. Exiting..." >&2
    exit 1
  fi

}

#######################################
# Rotates an input (circular) contig to approximately its midpoint
# Arguments:
#   fasta: FastA file containing a single circular contig
#   tmp_dir: directory for temp files. Cannot already exist.
# Returns:
#   STDOUT of the rotated FastA file, 60 nt per line
#######################################
function rotate_contig_to_midpoint {

  local fasta
  fasta=$1
  local tmp_dir
  tmp_dir=$2

  if [[ -d "${tmp_dir}" ]]; then

    echo "[ $(date -u) ]: Tmp dir cannot be created because it already exists: '${tmp_dir}'." >&2
    exit 1

  fi

  mkdir -p "${tmp_dir}"

  # Check only one contig in FastA file
  if [[ $(grep -c "^>" "${fasta}") -ne 1 ]]; then

    echo "[ $(date -u) ]: ERROR: rotate: input FastA file does not contain a single contig" >&2
    exit 1

  fi

  # Get contig name
  contig_name=$(grep "^>" "${fasta}" | head -n 1 | cut -d ">" -f 2)
  contig_name_short=$(echo "${contig_name}" | cut -d " " -f 1) # no comments

  # Get length of contig
  contig_length=$(seqtk seq -A -l 0 "${fasta}" 2>> "${verbose}" | head -n 2 | tail -n 1 | wc -m)
  contig_length=$((contig_length - 1))

  # Determine approximate midpoint - auto rounds to integer
  midpoint=$((contig_length / 2))

  # Make BED files
  printf "%s\t0\t%s\n" "${contig_name_short}" "${midpoint}" > "${tmp_dir}/front.bed"
  printf "%s\t%s\t%s\n" "${contig_name_short}" "${midpoint}" "${contig_length}" > "${tmp_dir}/back.bed"

  # Split
  seqtk subseq -l 60 "${fasta}" "${tmp_dir}/front.bed" > "${tmp_dir}/front.fasta"
  seqtk subseq -l 60 "${fasta}" "${tmp_dir}/back.bed" > "${tmp_dir}/back.fasta"
  rm "${tmp_dir}/front.bed" "${tmp_dir}/back.bed"

  # Report to user
  front_segment=$(head -n 1 "${tmp_dir}/front.fasta" | cut -d " " -f 1 | sed "s/>${contig_name_short}://g")
  back_segment=$(head -n 1 "${tmp_dir}/back.fasta" | cut -d " " -f 1 | sed "s/>${contig_name_short}://g")
  echo "[ $(date -u) ]: Rotating to approximate midpoint at ${midpoint}: ${back_segment}, ${front_segment}" >&2

  # Combine
  printf ">%s\n" "${contig_name}" > "${tmp_dir}/rotated.fasta"
  tail -n +2 "${tmp_dir}/back.fasta" >> "${tmp_dir}/rotated.fasta"
  tail -n +2 "${tmp_dir}/front.fasta" >> "${tmp_dir}/rotated.fasta"
  rm "${tmp_dir}/front.fasta" "${tmp_dir}/back.fasta"

  # Fix FastA lines and print to STDOUT
  seqtk seq -A -l 60 "${tmp_dir}/rotated.fasta"
  rm "${tmp_dir}/rotated.fasta"

  rmdir "${tmp_dir}"

}

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

    echo "[ $(date -u) ]: Tmp dir cannot be created because it already exists: '${tmp_dir}'." >&2
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
  local flye_read_error
  flye_read_error=$6
  local circlator_min_id
  circlator_min_id=$7
  local circlator_min_length
  circlator_min_length=$8
  local circlator_ref_end
  circlator_ref_end=$9
  local circlator_reassemble_end
  circlator_reassemble_end=${10}
  local threads
  threads=${11}
  local thread_mem
  thread_mem=${12}
  local verbose
  verbose=${13}

  # Internal settings
  local bam_file
  bam_file="${outdir}/long_read.bam"
  local end_repaired_contigs
  end_repaired_contigs="${outdir}/repaired.fasta"

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
  printf "" > "${end_repaired_contigs}.tmp"

  for contig in "${circular_contigs[@]}"; do

    mkdir -p "${outdir}/contigs/${contig}"
    echo "[ $(date -u) ]: End repair: '${contig}'" | tee -a "${verbose}" >&2

    local assembly
    assembly="${outdir}/contigs/${contig}/${contig}.fasta"
    local linked_ends
    linked_ends="false"

    echo "${contig}" > "${outdir}/contigs/${contig}/${contig}.list"
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

      # TODO - consider adding option to split long reads in half if they go around a short circular contig, like in circlator
      samtools view -@ "${threads}" -L "${regions}" -b "${bam_file}" 2>> "${verbose}" | \
        samtools fastq -0 "${fastq}" -n -@ "${threads}" 2>> "${verbose}"

      if [[ "${flye_read_error}" == 0 ]]; then

        flye "--${flye_read_mode}" "${fastq}" -o "${reassembly_dir}" \
          -t "${threads}" >> "${verbose}" 2>&1

      else

        flye "--${flye_read_mode}" "${fastq}" -o "${reassembly_dir}" --read_error "${flye_read_error}" \
          -t "${threads}" >> "${verbose}" 2>&1

      fi

      mkdir -p "${merge_dir}"
      circlator merge --verbose --min_id "${circlator_min_id}" --min_length "${circlator_min_length}" \
        --ref_end "${circlator_ref_end}" --reassemble_end "${circlator_reassemble_end}" \
        "${assembly}" "${reassembly_dir}/assembly.fasta" "${merge_dir}/merge" >> "${verbose}" 2>&1

      local pass_circularization
      pass_circularization=$("${SCRIPT_DIR}/flye_end_repair_utils.py" -c "${merge_dir}/merge.circularise.log")

      mkdir -p "${outdir}/contigs/${contig}/logs/L${length_cutoff}"
      cp "${merge_dir}/merge.circularise.log" "${merge_dir}/merge.circularise_details.log" \
        "${reassembly_dir}/assembly_info.txt" "${outdir}/contigs/${contig}/logs/L${length_cutoff}"

      # If the contig was circularized correctly, exit early
      if [[ "${pass_circularization}" == "true" ]]; then

        echo "[ $(date -u) ]: End repair: '${contig}': successfully linked contig ends" | tee -a "${verbose}" >&2

        # Rotate to midpoint so that the stitched ends can be polished more effectively (especially stich points)
        # TODO - sometimes small contigs are already rotated far from original origin. Stitch point hards to find. Does circlator report stitch point?
        rotate_contig_to_midpoint "${merge_dir}/merge.fasta" "${lenout}/rotate_tmp" \
          > "${lenout}/rotate.fasta" 2>> "${verbose}"

        seqtk seq -l 60 "${lenout}/rotate.fasta" >> "${end_repaired_contigs}.tmp"
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
      mv "${outdir}/contigs/${contig}/logs" "${outdir}/troubleshooting/${contig}"
      rm -r "${outdir:?}/contigs/${contig}"

    fi

  done

  if [[ ${failed_contigs} -gt 0 ]]; then
    mv "${end_repaired_contigs}.tmp" "${end_repaired_contigs}"

    printf "[ $(date -u) ]: End repair: %s contigs could not be circularized." "${failed_contigs}" | tee -a "${verbose}" >&2
    printf " A partial output file including successfully circularized contigs (and no linear " | tee -a "${verbose}" >&2
    pritnf "contigs) is available at '%s' for debugging.\n" "${end_repaired_contigs}" | tee -a "${verbose}" >&2
    printf "Exiting with error status. See temporary files and verbose log for more details.\n" | tee -a "${verbose}" >&2
    exit 1
  fi

  echo "[ $(date -u) ]: Sorting final FastA file and adding linear contigs" | tee -a "${verbose}" >&2
  seqtk subseq -l 60 "${all_contigs}" "${outdir}/linear_contigs.list" >> "${end_repaired_contigs}.tmp"

  grep "^>" "${all_contigs}" | cut -d ">" -f 2- | cut -d " " -f 1 > "${outdir}/sort_order.list"
  sort_fasta "${end_repaired_contigs}.tmp" "${outdir}/sort_order.list" "${outdir}/sort" \
    > "${end_repaired_contigs}" 2>> "${verbose}"

  # Clean up temp files
  rm "${bam_file}" "${bam_file}.bai" "${outdir}/circular_contigs.list" "${outdir}/linear_contigs.list" \
    "${end_repaired_contigs}.tmp" "${outdir}/sort_order.list"
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
    printf "   longreads.fastq.gz:          QC-passing Nanopore reads\n"
    printf "   assembly.fasta:              Contigs output from Flye\n"
    printf "   assembly_info.txt:           Assembly info file from Flye\n"
    printf "   outdir:                      Output directory path (might overwrite contents if the dir already exists!)\n\n"
    printf "Optional arguments (general):\n"
    printf "   -f flye_read_mode:           Type of input reads for Flye [nano-hq; can also set nano-raw]\n"
    printf "   -F flye_read_error:          Expected error rate of input reads, expressed as proportion (e.g., 0.03) [0 = use flye defaults]\n"
    printf "   -t threads:                  Number of processors to use [1]\n"
    printf "   -m memory:                   Memory (GB) to use per thread for samtools sort [1]\n\n"
    printf "Optional arguments (circlator merge module, used to stitch the reassembled ends to the original contigs):\n"
    printf "   -i circlator_min_id:         Percent identity threshold for circlator merge [99]\n"
    printf "   -l circlator_min_length:     Minimum required overlap (bp) between original and merge contigs [10000]\n"
    printf "   -e circlator_ref_end:        Minimum distance (bp) between end of original contig and nucmer hit [100]\n"
    printf "   -E circlator_reassemble_end: Minimum distance (bp) between end of merge contig and nucmer hit [100]\n\n"

    # Exit
    exit 0
  fi

  # Set defaults for options
  local flye_read_mode
  flye_read_mode="nano-hq"
  local flye_read_error
  flye_read_error=0
  local circlator_min_id
  circlator_min_id=99
  local circlator_min_length
  circlator_min_length=10000
  local circlator_ref_end
  circlator_ref_end=100
  local circlator_reassemble_end
  circlator_reassemble_end=100
  local threads
  threads=1
  local thread_mem
  thread_mem=1

  # Set options (help from https://wiki.bash-hackers.org/howto/getopts_tutorial; accessed March 8th, 2019)
  OPTIND=1 # reset the OPTIND counter just in case
  while getopts ":f:F:i:l:e:E:t:m:" opt; do
    case ${opt} in
      f)
        flye_read_mode=${OPTARG}
        ;;
      F)
        flye_read_error=${OPTARG}
        ;;
      i)
        circlator_min_id=${OPTARG}
        ;;
      l)
        circlator_min_length=${OPTARG}
        ;;
      e)
        circlator_ref_end=${OPTARG}
        ;;
      E)
        circlator_reassemble_end=${OPTARG}
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

  # Test positional args
  if [[ ! -f "${qc_long}" ]]; then
    echo "[ $(date -u) ]: ERROR: provided longreads.fastq.gz is not a file: '${qc_long}'. Exiting..." >&2
    exit 1
  elif [[ ! -f "${all_contigs}" ]]; then
    echo "[ $(date -u) ]: ERROR: provided assembly.fasta is not a file: '${all_contigs}'. Exiting..." >&2
    exit 1
  elif [[ ! -f "${circular_info}" ]]; then
    echo "[ $(date -u) ]: ERROR: provided assembly_info.txt is not a file: '${circular_info}'. Exiting..." >&2
    exit 1
  fi

  # Test flags
  if [[ "${flye_read_mode}" != "nano-raw" ]] && [[ "${flye_read_mode}" != "nano-hq" ]]; then
    echo "[ $(date -u) ]: ERROR: flye_read_mode must be nano-raw or nano-hq. You provided: '${flye_read_mode}'. Exiting..." >&2
    exit 1
  elif [[ "${flye_read_error}" -gt 1 ]] || [[ "${flye_read_error}" -lt 0 ]]; then
    echo "[ $(date -u) ]: ERROR: flye_read_error must be a proportion between 0-1. You provided: '${flye_read_error}'. Exiting..." >&2
    exit 1
  elif [[ ! "${circlator_min_id}" =~ ^[0-9]+$ ]] || [[ "${circlator_min_id}" -eq 0 ]]; then
    echo "[ $(date -u) ]: ERROR: circlator_min_id must be a positive integer. You provided: '${circlator_min_id}'. Exiting..." >&2
    exit 1
  elif [[ ! "${circlator_min_length}" =~ ^[0-9]+$ ]] || [[ "${circlator_min_length}" -eq 0 ]]; then
    echo "[ $(date -u) ]: ERROR: circlator_min_length must be a positive integer. You provided: '${circlator_min_length}'. Exiting..." >&2
    exit 1
  elif [[ ! "${circlator_ref_end}" =~ ^[0-9]+$ ]] || [[ "${circlator_ref_end}" -eq 0 ]]; then
    echo "[ $(date -u) ]: ERROR: circlator_ref_end must be a positive integer. You provided: '${circlator_ref_end}'. Exiting..." >&2
    exit 1
  elif [[ ! "${circlator_reassemble_end}" =~ ^[0-9]+$ ]] || [[ "${circlator_reassemble_end}" -eq 0 ]]; then
    echo "[ $(date -u) ]: ERROR: circlator_reassemble_end must be a positive integer. You provided: '${circlator_reassemble_end}'. Exiting..." >&2
    exit 1
  elif [[ ! "${threads}" =~ ^[0-9]+$ ]] || [[ "${threads}" -eq 0 ]]; then
    echo "[ $(date -u) ]: ERROR: threads must be a positive integer. You provided: '${threads}'. Exiting..." >&2
    exit 1
  elif [[ ! "${thread_mem}" =~ ^[0-9]+$ ]] || [[ "${thread_mem}" -eq 0 ]]; then
    echo "[ $(date -u) ]: ERROR: thread_mem must be a positive integer. You provided: '${thread_mem}'. Exiting..." >&2
    exit 1
  fi

  # Check dependencies
  check_dependency python commandline
  check_dependency pip commandline
  check_dependency seqtk commandline
  check_dependency circlator commandline
  check_dependency minimap2 commandline
  check_dependency samtools commandline
  check_dependency pandas python_package
  check_dependency flye_end_repair_utils.py internal

  # Make output dir
  outdir_test=$(test -d "${outdir}"; echo $?)
  mkdir -p "${outdir}"

  # Initialize log file
  local verbose
  verbose="${outdir}/verbose.log"
  printf "" > "${verbose}"

  echo "[ $(date -u) ]: Running ${SCRIPT_NAME}" | tee -a "${verbose}" >&2
  echo "[ $(date -u) ]: Command run: ${SCRIPT_NAME} ${original_arguments}" | tee -a "${verbose}" >&2

  if [[ "${outdir_test}" -eq 0 ]]; then
    echo "[ $(date -u) ]: Warning: output dir already exists. Files could collide." | tee -a "${verbose}" >&2
  fi

  run_pipeline "${qc_long}" "${all_contigs}" "${circular_info}" "${outdir}" \
    "${flye_read_mode}" "${flye_read_error}" "${circlator_min_id}" "${circlator_min_length}" \
    "${circlator_ref_end}" "${circlator_reassemble_end}" "${threads}" "${thread_mem}" "${verbose}"

}

# Only run the script if it is called from the command line
if [[ ${BASH_SOURCE[0]} = "${0}" ]]; then
  main "$@"
fi
