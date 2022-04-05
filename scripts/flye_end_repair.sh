#!/usr/bin/env bash
set -euxo pipefail
# flye_end_repair.sh
# Copyright Jackson M. Tsuji, ILTS, Hokkaido University, 2022


# Internal variables
readonly VERSION="0.1.0"
readonly SCRIPT_NAME="${0##*/}"
readonly SCRIPT_DIR=$(realpath "${0%/*}")


readonly bam_file="${outdir}/long_read.bam"
readonly end_repaired_contigs="${outdir}/repaired.fasta"
readonly log="${outdir}/end_repair.log"
readonly verbose="${outdir}/end_repair_verbose.log"

# Initialize log files
mkdir "${outdir}"
printf "" > "${log}"
printf "" > "${verbose}"

echo "[ $(date -u) ]: Starting Flye-nano end repair pipeline" | tee -a "${log}" | tee -a "${verbose}"
echo "[ $(date -u) ]: Mapping reads to all contigs" | tee -a "${log}" | tee -a "${verbose}"

# Map reads
# TODO - add support for different flags like -ax for pacbio
minimap2 -t "${threads}" -ax map-ont "${all_contigs}" "${qc_long}" 2>> "${verbose}" | \
  samtools view -b -@ "${threads}" 2>> "${verbose}" | \
  samtools sort -@ "${threads}" -m "${thread_mem}G" 2>> "${verbose}" \
  > "${bam_file}"
samtools index -@ "${threads}" "${bam_file}" 2>> "${log}"

# Get a list of circular contigs
"${SCRIPT_DIR}/circular_contig_toolkit.py" -i "${circular_info}" -a "${outdir}/circular_contigs.list" 2>> "${verbose}"

if [[ $(cat "${outdir}/circular_contigs.list" | wc -l) == 0 ]]; then

  echo "[ $(date -u) ]: No circular contigs, so making a copy of the input file and finishing early." | \
    tee -a "${log}" | tee -a "${verbose}"

  cp "${all_contigs}" "${end_repaired_contigs}"
  exit 0

fi

circular_contigs=($(cat "${outdir}/circular_contigs.list"))
failed_contigs=0

mkdir -p "${outdir}/circlator_logs"
printf "" > "${end_repaired_contigs}"

for contig in "${circular_contigs[@]}"; do

  mkdir -p "${outdir}/${contig}"
  echo "[ $(date -u) ]: End repair: '${contig}'" | tee -a "${log}" | tee -a "${verbose}"
  
  assembly="${outdir}/${contig}/${contig}.fasta"
  
  echo "${contig}" > "${outdir}/${contig}/${contig}.list"
  seqtk subseq -l 60 "${all_contigs}" "${outdir}/${contig}/${contig}.list" > "${assembly}" 2>> "${verbose}"

  # TODO - allow user to set the length cutoff range
  for length_cutoff in 100000 90000 80000 70000 60000 50000; do
  
    lenout="${outdir}/${contig}/L${length_cutoff}"
    regions="${lenout}/ends.bed"
    fastq="${lenout}/ends.fastq.gz"
    reassembly_dir="${lenout}/assembly"
    merge_dir="${lenout}/merge"
    
    mdkir -p "${lenout}"

    echo "[ $(date -u) ]: End repair: '${contig}': attempting at ${length_cutoff} bp from ends"\
      | tee -a "${log}" | tee -a "${verbose}"
    "${SCRIPT_DIR}/circular_contig_toolkit.py" -i "${circular_info}" -n "${contig}" -l "${length_cutoff}" \
      -b "${regions}" 2>> "${verbose}"

    samtools view -@ "${threads}" -L "${regions}" -b "${bam_file}" 2>> "${verbose}" | \
      samtools fastq -0 "${fastq}" -n -@ "${threads}" 2>> "${verbose}"

    flye "--${flye_assembly_mode}" "${fastq}" -o "${reassembly_dir}" -t "${threads}" >> "${verbose}" 2>&1
    # TODO: add support for optional flags like --meta

    mkdir -p "${merge_dir}"
    circlator merge --verbose --min_id "${circlator_min_id}" --min_length "${circlator_min_length}" \
      "${assembly}" "${reassembly_dir}/assembly.fasta" "${merge_dir}/merge" >> "${verbose}" 2>&1
    
    pass_circularization=$("${SCRIPT_DIR}/check_circlator_status.py" -i "${merge_dir}/merge.circularise.log")
    
    mkdir -p "${outdir}/${contig}/logs/L${length_cutoff}"
    cp "${merge_dir}/merge.circularise_details.log" "${reassembly_dir}/assembly_info.txt" \
      "${outdir}/${contig}/logs/L${length_cutoff}"
    
    if [[ "${pass_circularization}" == "true" ]]; then

      echo "[ $(date -u) ]: End repair: '${contig}': successfully linked contig ends" | \
        tee -a "${log}" | tee -a "${verbose}"

      seqtk seq -l 60 "${merge_dir}/merge.fasta" >> "${end_repaired_contigs}"

      cp "${merge_dir}/merge.circularise_details.log" "${outdir}/circlator_logs/${contig}.log"

      rm -r "${outdir:?}/${contig}"

      break
    fi

  done

  # If loop finishes without successful circularization
  echo "[ $(date -u) ]: End repair: '${contig}': FAILED to linked contig ends" | tee -a "${log}" | tee -a "${verbose}"

  failed_contigs=$((failed_contigs + 1))

  mkdir -p "${outdir}/troubleshooting"
  mv "${outdir}/${contig}/logs" "${outdir}/troubleshooting/${contig}"
  rm -r "${outdir:?}/${contig}"

done

if [[ ${failed_contigs} -gt 0 ]]; then
  echo "[ $(date -u) ]: End repair: ${failed_contigs} could not be circularized and were not included in the final output file. Exiting with error status. See temporary files for more details." | tee -a ${log} | tee -a "${verbose}"
  exit 1

else
  echo "[ $(date -u) ]: End repair finished. Output contigs saved at '${end_repaired_contigs}'." | tee -a ${log} | tee -a "${verbose}"
  
  # Clean up temp files
  rm "${bam_file}" "${bam_file}.bai"
fi

#######################################
# Perform the 'run' command
# Globals:
#   SCRIPT_NAME: the name of this script
#   VERSION: the script version
#   SNAKEFILE: the path to the Snakefile used to run BackBLAST
# Arguments:
#   all command line inputs for the 'run' module - see help statement below
# Returns:
#   runs the 'run' command end-to-end
#######################################
function run_pipeline() {

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
    printf "Usage: %s longreads.fastq.gz assembly.fasta assembly_info.txt outdir\n\n" "${SCRIPT_NAME}"
    printf "Positional arguments:\n"
    printf "   longreads.fastq.gz:      QC-passing Nanopore reads\n"
    printf "   assembly.fasta:          Contigs output from Flye\n"
    printf "   assembly_info.txt:       Assembly info file from Flye\n"
    printf "   outdir:                  Output directory path (directory must not yet exist)\n"
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
  flye_read_mode="nano-qc"
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
