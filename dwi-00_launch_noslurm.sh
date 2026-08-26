#!/bin/bash

# ============================================================
# DWI Pipeline Runner (No-SLURM version)
#
# Processes all subjects found in a BIDS directory. The user MUST
# specify which pipeline stage to run:
#
#   preproc  — dwi preprocessing + eddy
#   tracto   — tractography + connectivity matrices
#   qc       — quality control reports
#   all      — run preproc → tracto → qc in sequence
#
# This separation enforces the workflow: preprocessing should be
# QC-ed before tractography is started.
#
# Subjects can be processed in series (default) or in parallel
# using bash background jobs with a concurrency limit.
# ============================================================

set -uo pipefail

########################
# DEFAULT CONFIGURATION
########################

MAX_PARALLEL=1
CONTAINERPATH=""
NICE_LEVEL=10
DRY_RUN=0
SINGLE_SUBJECT=""
STAGE=""

# Per-stage CPU/thread settings
PREPROC_CPUS=8
TRACTO_CPUS=16
QC_CPUS=1

# Colors
NC='\033[0m'
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'

########################
# FUNCTIONS
########################

Usage() {
  cat <<EOF
Usage: $0 <stage> [OPTIONS] <path_to_spec.json>

Run tractoprep pipelines for all subjects in a BIDS directory.
Works without SLURM — subjects can be processed in series or
in parallel.

Required argument:
  <stage>            Which pipeline stage to run:
                       preproc   — dwi preprocessing + eddy
                       tracto    — tractography + connectivity
                       qc        — quality control reports
                       all       — preproc → tracto → qc in sequence

Options:
  --parallel N       Run N subjects concurrently (default: 1, serial)
  --serial           Run subjects one at a time (same as --parallel 1)
  --container PATH   Path to apptainer .sif container
                     (default: ${CONTAINERPATH})
  --subject SUB      Process only this subject (e.g., sub-01)
  --nice N           Set nice level 0-19 (default: ${NICE_LEVEL})
  --dry-run          Print commands without executing them
  --help, -h         Show this help message

Environment variables:
  MAX_PARALLEL       Same as --parallel
  CONTAINERPATH      Same as --container

Examples:
  # Preprocess all subjects (serial)
  $0 preproc spec.json

  # Preprocess all subjects (4 in parallel)
  $0 preproc --parallel 4 spec.json

  # Run tractography for a single subject
  $0 tracto --subject sub-01 spec.json

  # Run QC for all subjects
  $0 qc spec.json

  # Run all stages in sequence (preproc → tracto → qc)
  $0 all --parallel 4 spec.json

  # Dry run (see what would be executed)
  $0 preproc --dry-run --parallel 4 spec.json

Requirements:
  - bash 4.3+ (for parallel mode; serial mode works with older bash)
  - apptainer (or singularity) installed and in PATH
  - jq installed
  - GNU find (for -printf; standard on Linux)
EOF
  exit 1
}

# ----------------------------------------------------------------
# run_stage: execute a single pipeline stage for one subject
#   $1 = stage name (for logging)
#   $2 = command string to execute
#   $3 = number of CPUs/threads for this stage
#   $4 = path to log file
#   $5 = subject ID (for logging)
# ----------------------------------------------------------------
run_stage() {
    local stage_name="$1"
    local cmd="$2"
    local cpus="$3"
    local logfile="$4"
    local subj="$5"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting ${stage_name} for ${subj}..." | tee -a "$logfile"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY RUN] Would execute:" | tee -a "$logfile"
        echo "  $cmd" | tee -a "$logfile"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${stage_name} (dry run) completed." | tee -a "$logfile"
        return 0
    fi

    OMP_NUM_THREADS="$cpus" \
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="$cpus" \
    nice -n "$NICE_LEVEL" \
        bash -c "$cmd" >> "$logfile" 2>&1

    local status=$?
    if [[ $status -ne 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: ${stage_name} failed with exit code ${status}" | tee -a "$logfile"
        return $status
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${stage_name} completed successfully." | tee -a "$logfile"
    return 0
}

# ----------------------------------------------------------------
# run_preproc: Stage 1 — dwi-preproc + eddy
#   $1 = subject ID
# ----------------------------------------------------------------
run_preproc() {
    local subj="$1"
    local session="${session:-}"

    local sessionpath="/${session:+${session}/}"
    local sessionfile="_${session:+${session}_}"

    local logdir="${outputdir}/logs/${subj}${sessionpath}"
    mkdir -p "${logdir}"

    local subjspecjson="${bidsdir}/spec_${subj}.json"
    jq --arg subj "$subj" '.subj = $subj' "${templatejson}" > "${subjspecjson}"

    # Expected output files (skip-logic)
    local preproc_nifti="${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz"
    local preproc_bvec="${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec"
    local preproc_qc="${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz"
    local preproc_anat="${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_template.nii.gz"
    local preproc_atlas="${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-400P17N_dseg.nii.gz"

    local logfile="${logdir}/dwi-preproc.log"

    if [[ -f "${preproc_nifti}" ]] && [[ -f "${preproc_bvec}" ]] && [[ -f "${preproc_qc}" ]] && [[ -f "${preproc_anat}" ]] && [[ -f "${preproc_atlas}" ]]; then
        echo "Preprocessing already done for ${subj} ${session:-}, skipping." | tee -a "$logfile"
        return 0
    fi

    local cmd="tmpdir_job=\"${workdir}/tmp/${subj}${sessionfile}\${BASHPID}\"; \
mkdir -p \"\${tmpdir_job}\"; \
trap 'rm -rf \"\${tmpdir_job}\"' EXIT; \
${RUNTIME} run --cleanenv \
  --bind ${bidsdir}:${bidsdir},${workdir}:${workdir},${outputdir}:${outputdir},\${tmpdir_job}:/scratch \
  --env TMPDIR=/scratch --env TMP=/scratch --env TEMP=/scratch \
  --env OMP_NUM_THREADS=${PREPROC_CPUS} \
  ${CONTAINERPATH} dwi-preproc ${subjspecjson}"

    if ! run_stage "dwi-preproc" "$cmd" "$PREPROC_CPUS" "$logfile" "$subj"; then
        echo -e "${RED}ERROR: dwi-preproc failed for ${subj}.${NC}" | tee -a "$logfile"
        return 1
    fi
    return 0
}

# ----------------------------------------------------------------
# run_tracto: Stage 2 — dwi-tracto + connectivity
#   $1 = subject ID
# ----------------------------------------------------------------
run_tracto() {
    local subj="$1"
    local session="${session:-}"

    local sessionpath="/${session:+${session}/}"
    local sessionfile="_${session:+${session}_}"

    local logdir="${outputdir}/logs/${subj}${sessionpath}"
    mkdir -p "${logdir}"

    local subjspecjson="${bidsdir}/spec_${subj}.json"
    jq --arg subj "$subj" '.subj = $subj' "${templatejson}" > "${subjspecjson}"

    # Expected output files (skip-logic)
    local tracto_file="${outputdir}/dwi-tracto/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck"
    local conn_file="${outputdir}/dwi-tracto/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-400P17N_desc-streams_connmatrix.csv"

    local logfile="${logdir}/dwi-tracto.log"

    if [[ -f "${tracto_file}" ]] && [[ -f "${conn_file}" ]]; then
        echo "Tractography already done for ${subj} ${session:-}, skipping." | tee -a "$logfile"
        return 0
    fi

    local cmd="tmpdir_job=\"${workdir}/tmp/${subj}${sessionfile}\${BASHPID}\"; \
mkdir -p \"\${tmpdir_job}\"; \
trap 'rm -rf \"\${tmpdir_job}\"' EXIT; \
${RUNTIME} run --cleanenv \
  --bind ${bidsdir}:${bidsdir},${workdir}:${workdir},${outputdir}:${outputdir},\${tmpdir_job}:/scratch \
  --env TMPDIR=/scratch --env TMP=/scratch --env TEMP=/scratch \
  --env OMP_NUM_THREADS=${TRACTO_CPUS} \
  ${CONTAINERPATH} dwi-tracto ${subjspecjson}"

    if ! run_stage "dwi-tracto" "$cmd" "$TRACTO_CPUS" "$logfile" "$subj"; then
        echo -e "${RED}ERROR: dwi-tracto failed for ${subj}.${NC}" | tee -a "$logfile"
        return 1
    fi
    return 0
}

# ----------------------------------------------------------------
# run_qc: Stage 3 — dwi-qc
#   $1 = subject ID
# ----------------------------------------------------------------
run_qc() {
    local subj="$1"
    local session="${session:-}"

    local sessionpath="/${session:+${session}/}"
    local sessionfile="_${session:+${session}_}"

    local logdir="${outputdir}/logs/${subj}${sessionpath}"
    mkdir -p "${logdir}"

    local subjspecjson="${bidsdir}/spec_${subj}.json"
    jq --arg subj "$subj" '.subj = $subj' "${templatejson}" > "${subjspecjson}"

    local bindcmd="${bidsdir}:${bidsdir},${workdir}:${workdir},${outputdir}:${outputdir}"
    local logfile="${logdir}/dwi-qc.log"

    local cmd="${RUNTIME} run --cleanenv --bind ${bindcmd} \
  --env OMP_NUM_THREADS=${QC_CPUS} \
  ${CONTAINERPATH} dwi-qc ${subjspecjson}"

    if ! run_stage "dwi-qc" "$cmd" "$QC_CPUS" "$logfile" "$subj"; then
        echo -e "${RED}ERROR: dwi-qc failed for ${subj}.${NC}" | tee -a "$logfile"
        return 1
    fi
    return 0
}

# ----------------------------------------------------------------
# process_subject: dispatch to the selected stage(s)
#   $1 = subject ID
#   $2 = stage name (preproc, tracto, qc, or all)
# ----------------------------------------------------------------
process_subject() {
    local subj="$1"
    local stage="$2"

    echo -e "${BLUE}Processing ${subj} (stage: ${stage})...${NC}"

    case "$stage" in
        preproc)
            run_preproc "$subj" || return 1
            ;;
        tracto)
            run_tracto "$subj" || return 1
            ;;
        qc)
            run_qc "$subj" || return 1
            ;;
        all)
            run_preproc "$subj" || return 1
            run_tracto "$subj" || return 1
            run_qc     "$subj" || return 1
            ;;
        *)
            echo "ERROR: Unknown stage '${stage}'"
            return 1
            ;;
    esac

    echo -e "${GREEN}Stage '${stage}' completed for ${subj}.${NC}"
    return 0
}

########################
# ARGUMENT PARSING
########################

# First positional argument must be the stage
if [[ $# -lt 1 ]]; then
    echo "Error: stage argument is required."
    echo ""
    Usage
fi

case "$1" in
    preproc|tracto|qc|all)
        STAGE="$1"
        shift
        ;;
    --help|-h)
        Usage
        ;;
    --*)
        echo "Error: first argument must be the stage (preproc, tracto, qc, or all)."
        echo ""
        Usage
        ;;
    *)
        echo "Error: '$1' is not a valid stage."
        echo "Choose one of: preproc, tracto, qc, all"
        echo ""
        Usage
        ;;
esac

# Parse remaining options
while [[ $# -gt 0 ]]; do
  case "$1" in
    --parallel)
      MAX_PARALLEL="$2"
      shift 2
      ;;
    --serial)
      MAX_PARALLEL=1
      shift
      ;;
    --container)
      CONTAINERPATH="$2"
      shift 2
      ;;
    --subject)
      SINGLE_SUBJECT="$2"
      shift 2
      ;;
    --nice)
      NICE_LEVEL="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --help|-h)
      Usage
      ;;
    --*)
      echo "Unknown option: $1"
      Usage
      ;;
    *)
      break
      ;;
  esac
done

# The remaining argument should be the spec.json path
if [[ $# -ne 1 ]]; then
    echo "Error: spec.json path is required."
    echo ""
    Usage
fi

templatejson="$1"

########################
# VALIDATION
########################

if ! [[ "$MAX_PARALLEL" =~ ^[0-9]+$ ]] || [[ "$MAX_PARALLEL" -lt 1 ]]; then
  echo "Error: --parallel must be a positive integer."
  exit 1
fi

if [[ "$MAX_PARALLEL" -gt 1 ]]; then
  bash_major="${BASH_VERSINFO[0]}"
  bash_minor="${BASH_VERSINFO[1]}"
  if [[ "$bash_major" -lt 4 ]] || { [[ "$bash_major" -eq 4 ]] && [[ "$bash_minor" -lt 3 ]]; }; then
    echo "Error: Parallel mode requires bash 4.3+ (for 'wait -n')."
    echo "       You have bash ${BASH_VERSION}."
    echo "       Use --serial mode or upgrade bash."
    exit 1
  fi
fi

########################
# ENVIRONMENT SETUP
########################

if command -v module &>/dev/null; then
  module load apptainer 2>/dev/null || true
fi

if command -v apptainer &>/dev/null; then
  RUNTIME="apptainer"
elif command -v singularity &>/dev/null; then
  RUNTIME="singularity"
  echo "Note: apptainer not found, using singularity instead."
else
  echo "Error: Neither apptainer nor singularity found in PATH."
  exit 1
fi

if ! command -v jq &>/dev/null; then
  echo "Error: jq is required but not found in PATH."
  exit 1
fi

if [[ ! -f "$CONTAINERPATH" ]]; then
  echo "Error: Container not found: $CONTAINERPATH"
  exit 1
fi

# Convert to absolute path
templatejson="$(cd "$(dirname "$templatejson")" && pwd)/$(basename "$templatejson")"

########################
# READ SPEC.JSON
########################

for key in $(jq -r 'keys[]' "${templatejson}"); do
  value=$(jq -r --arg k "$key" '.[$k]' "${templatejson}")
  declare "$key"="$value"
done

unset subj

for var in bidsdir outputdir workdir nstreamlines; do
  if [[ -z "${!var:-}" ]]; then
    echo "Error: Variable '$var' is not set or is empty in spec.json."
    exit 1
  fi
done

mkdir -p "${workdir}" "${outputdir}"

########################
# DISCOVER SUBJECTS
########################

cd "${bidsdir}"

mapfile -t subjects < <(find . -maxdepth 1 -mindepth 1 -type d -name 'sub-*' -printf '%f\n' | sort)

if [[ ${#subjects[@]} -eq 0 ]]; then
  echo "Error: No subjects found in ${bidsdir}"
  exit 1
fi

if [[ -n "$SINGLE_SUBJECT" ]]; then
  if [[ ! -d "${bidsdir}/${SINGLE_SUBJECT}" ]]; then
    echo "Error: Subject directory not found: ${bidsdir}/${SINGLE_SUBJECT}"
    exit 1
  fi
  subjects=("$SINGLE_SUBJECT")
fi

########################
# PRINT RUN INFO
########################

echo "============================================"
echo "DWI Pipeline Runner (local mode)"
echo "============================================"
echo "Stage:        ${STAGE}"
echo "Subjects:     ${#subjects[@]} (${subjects[*]})"
echo "Parallelism:  ${MAX_PARALLEL}"
echo "Container:    ${CONTAINERPATH}"
echo "Runtime:      ${RUNTIME}"
echo "Nice level:   ${NICE_LEVEL}"
echo "Dry run:      ${DRY_RUN}"
echo "Spec JSON:    ${templatejson}"
echo "BIDS dir:     ${bidsdir}"
echo "Output dir:   ${outputdir}"
echo "Work dir:     ${workdir}"
echo "============================================"
echo ""

########################
# MAIN EXECUTION
########################

STATUS_DIR=$(mktemp -d)
trap 'rm -rf "$STATUS_DIR"' EXIT

if [[ "$MAX_PARALLEL" -eq 1 ]]; then
    # ===== SERIAL MODE =====
    failed_subjects=()
    for subj in "${subjects[@]}"; do
        if process_subject "$subj" "$STAGE"; then
            echo -e "${GREEN}✓ ${subj} completed successfully.${NC}"
        else
            echo -e "${RED}✗ ${subj} FAILED.${NC}"
            failed_subjects+=("$subj")
        fi
        echo ""
    done

else
    # ===== PARALLEL MODE =====
    echo "Launching ${#subjects[@]} subjects with max ${MAX_PARALLEL} concurrent..."
    echo ""

    running=0
    for subj in "${subjects[@]}"; do
        (
            if process_subject "$subj" "$STAGE"; then
                touch "$STATUS_DIR/success_${subj}"
            else
                touch "$STATUS_DIR/failed_${subj}"
            fi
        ) &
        ((running++))
        if [[ $running -ge $MAX_PARALLEL ]]; then
            wait -n 2>/dev/null || true
            ((running--))
        fi
    done

    wait

    failed_subjects=()
    for subj in "${subjects[@]}"; do
        if [[ -f "$STATUS_DIR/failed_${subj}" ]]; then
            failed_subjects+=("$subj")
        fi
    done
fi

########################
# SUMMARY
########################

echo ""
echo "============================================"
echo "SUMMARY"
echo "============================================"
echo "Stage:          ${STAGE}"
echo "Total subjects: ${#subjects[@]}"
echo "Succeeded:       $((${#subjects[@]} - ${#failed_subjects[@]}))"
echo "Failed:          ${#failed_subjects[@]}"

if [[ ${#failed_subjects[@]} -gt 0 ]]; then
    echo -e "${RED}Failed subjects:${NC}"
    for s in "${failed_subjects[@]}"; do
        echo "  - $s"
    done
    echo ""
    echo "Check per-subject logs in: ${outputdir}/logs/"
    exit 1
else
    echo -e "${GREEN}All subjects completed successfully.${NC}"
    exit 0
fi
