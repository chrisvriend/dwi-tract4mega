#!/bin/bash

jsonfile=$1
scriptdir=/tracto

export FSLOUTPUTTYPE=NIFTI_GZ

# Colors (define if not present in your environment)
NC='\033[0m'
BLUE='\033[0;34m'
GREEN='\033[0;32m'
RED='\033[0;31m'

Usage() {
  echo "Usage: $0 <path_to_spec.json> [--all]"
  echo "  --all   Run QC for every subject folder found in \${outputdir}/dwi-preproc,"
  echo "          overriding any 'subj' set in the spec.json file."
  exit 1
}

# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

# Summary tracking
declare -a COMPLETED_LIST=()
declare -a FAILED_LIST=()
declare -a SKIPPED_LIST=()

all_flag=0

if [[ $# -eq 1 ]]; then
  all_flag=0
elif [[ $# -eq 2 && "$2" == "--all" ]]; then
  all_flag=1
else
  Usage
fi

# read spec.json file into variables
for key in $(jq -r 'keys[]' "${jsonfile}"); do
  value=$(jq -r --arg k "$key" '.[$k]' "${jsonfile}")
  declare "$key"="$value"
done

# check if required variables are non-empty
for var in bidsdir outputdir scriptdir; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  

# Helper: build session path and file suffix
set_session_vars() {
    local sess="$1"
    if [[ -z "${sess}" ]]; then
        sessionpath="/"
        sessionfile="_"
    else
        sessionpath="/${sess}/"
        sessionfile="_${sess}_"
    fi
}

# Helper: build argument list for generate_qc_report.py based on existing files
build_qc_args() {
    local subj="$1"
    local session="$2"
    local atlas="$3"

    set_session_vars "${session}"

    local qc_args=()

    # Define all candidate files and their corresponding flags
    # Format: "flag|full_path"
    local entries=(
        # bvals check
        "bval-check-bvals|${bidsdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.bval"
        "noise|${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_desc-noise_dwi.nii.gz"
         # topup
        "topup-before|${outputdir}/dwi-preproc/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-*_space-dwi_desc-4topup_epi.nii.gz"
        "topup-after|${outputdir}/dwi-preproc/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-unwarped_epi.nii.gz"
        "topup-acqparams|${outputdir}/dwi-preproc/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-*_desc-refparams.tsv"
        "topup-fieldmap|${outputdir}/dwi-preproc/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-topup_fieldmap.nii.gz"
        "topup-dwi-acqparams|${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}acq-dwi_desc-acqparams.tsv"
        # eddy outcomes
        "eddy-json|${outputdir}/dwi-preproc/${subj}${sessionpath}qc/qc.json"
        "eddy-rms|${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_desc-preproc.eddy_movement_rms"
        "eddy-outliers|${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_desc-preproc.eddy_outlier_report"
        "eddy-raw-dwi|${bidsdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.nii.gz"
        "eddy-preproc-dwi|${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz"
        "eddy-cnr-maps|${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz"
        # brain masks
        "brainmask-nodif|${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz"
        "brainmask-mask|${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz"
        # registraion outputs 
        "reg-t1w-dwi|${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-FS_desc-brain_T1w.nii.gz"
        "reg-nodif|${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_template.nii.gz"
        "reg-5ttvis|${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_vis.nii.gz"
        # response voxels
        "response-voxels|${outputdir}/dwi-tracto/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_desc-response_voxels.nii.gz"
        "response-underlay|${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodif-brain_dwi.nii.gz"
        # tractography outputs
        "tract-tck|${outputdir}/dwi-tracto/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_tracto-100k.tck"
        "connectivity-matrix|${outputdir}/dwi-tracto/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-${atlas}_desc-streams_connmatrix.csv"
    )

    for entry in "${entries[@]}"; do
        IFS='|' read -r flag path <<< "${entry}"

        case "$flag" in
            topup-before|topup-acqparams)
                # wildcard entries: expand glob and use first match
                shopt -s nullglob
                matches=( $path )   # unquoted for globbing
                shopt -u nullglob

                if (( ${#matches[@]} == 0 )); then
                    log "$BLUE" "Skipping --${flag}: no files match pattern (${path})"
                    continue
                fi

                qc_args+=( "--${flag}" "${matches[0]}" )
                ;;
            *)
                # normal files
                if [[ -f "$path" ]]; then
                    qc_args+=( "--${flag}" "$path" )
                else
                    log "$BLUE" "Skipping --${flag}: file not found (${path})"
                fi
                ;;
        esac
    done

    # Always append output and subject options (no file checks)
    qc_args+=( "--output" "${outputdir}/dwi-preproc/${subj}${sessionfile}qc.html" "--subject" "${subj}-${session:-}" )

    # Return via global array
    QC_ARGS=( "${qc_args[@]}" )
}

run_qc_for_subject_session() {
    local subj="$1"
    local session="$2"

    set_session_vars "${session}"

    echo
    log "$BLUE" "-- ------------------ --"
    log "$BLUE" "Starting DWI QC for subject: ${subj} ${session}"
    log "$BLUE" "-- ------------------ --"
    echo

    # Build QC_ARGS array based on existing files
    build_qc_args "${subj}" "${session}" "${atlas}"

    # If only output/subject are present and no inputs, skip
    if (( ${#QC_ARGS[@]} <= 4 )); then
        log "$RED" "No QC inputs found for ${subj} ${session}. Skipping."
        SKIPPED_LIST+=( "${subj}${sessionfile%_}" )
        return
    fi

    # Run QC script
    "${scriptdir}/generate_qc_report.py" "${QC_ARGS[@]}"
    local ret=$?

    if [[ ${ret} -ne 0 ]]; then
        log "$RED" "QC script failed for ${subj} ${session} (exit code ${ret})"
        error=1
        FAILED_LIST+=( "${subj}${sessionfile%_} (exit ${ret})" )
    else
        log "$GREEN" "-- ------------------ --"
        log "$GREEN" "DWI QC completed for subject: ${subj} ${session}"
        log "$GREEN" "-- ------------------ --"
        echo
        COMPLETED_LIST+=( "${subj}${sessionfile%_}" )
    fi
}

# Discover subjects and sessions if not specified in JSON
run_all() {
    # Subjects: look for sub-* directories in bidsdir
    for subdir in ${outputdir}/dwi-preproc/sub-*; do
        [[ -d "${subdir}" ]] || continue
        local subj
        subj=$(basename "${subdir}")

        # Discover sessions; if none, use empty session
        shopt -s nullglob
        local sess_dirs=( "${subdir}"/ses-* )
        shopt -u nullglob

        if (( ${#sess_dirs[@]} == 0 )); then
            # no sessions -> single run with session=""
            run_qc_for_subject_session "${subj}" ""
        else
            for sdir in "${sess_dirs[@]}"; do
                [[ -d "${sdir}" ]] || continue
                local session
                session=$(basename "${sdir}")
                run_qc_for_subject_session "${subj}" "${session}"
            done
        fi
    done
}

# Main logic:
# - If --all was passed on the command line, run for every subject under
#   outputdir/dwi-preproc, regardless of any 'subj' set in the JSON.
# - Else if subj is set in JSON, run only for that subject (and provided/auto-discovered sessions)
# - Else (subj not set, no --all), run for all subjects under bidsdir
if [[ ${all_flag} -eq 1 ]]; then
    log "$BLUE" "--all specified. Running QC for all subjects in ${outputdir}/dwi-preproc."
    run_all
elif [[ -z "${subj}" ]]; then
    log "$BLUE" "No subject specified in JSON. Running QC for all subjects in ${bidsdir}."
    run_all
else
    if [[ -z "${session}" ]]; then
        # discover sessions for this subject
        subdir="${outputdir}/dwi-preproc/${subj}"
        if [[ ! -d "${subdir}" ]]; then
            log "$RED" "Subject directory not found: ${subdir}"
            exit 1
        fi

        shopt -s nullglob
        sess_dirs=( "${subdir}"/ses-* )
        shopt -u nullglob

        if (( ${#sess_dirs[@]} == 0 )); then
            run_qc_for_subject_session "${subj}" ""
        else
            for sdir in "${sess_dirs[@]}"; do
                [[ -d "${sdir}" ]] || continue
                session=$(basename "${sdir}")
                run_qc_for_subject_session "${subj}" "${session}"
            done
        fi
    else
        # subj and session both supplied
        run_qc_for_subject_session "${subj}" "${session}"
    fi
fi

# Summary report
echo
log "$BLUE" "======================================"
log "$BLUE" " QC RUN SUMMARY"
log "$BLUE" "======================================"
log "$GREEN" "Completed: ${#COMPLETED_LIST[@]}"
log "$RED"   "Failed:    ${#FAILED_LIST[@]}"
log "$BLUE"  "Skipped:   ${#SKIPPED_LIST[@]}"
echo

if (( ${#FAILED_LIST[@]} > 0 )); then
    log "$RED" "-- Failed --"
    for item in "${FAILED_LIST[@]}"; do
        log "$RED" "  ${item}"
    done
    echo
fi

if (( ${#SKIPPED_LIST[@]} > 0 )); then
    log "$BLUE" "-- Skipped (no QC inputs found) --"
    for item in "${SKIPPED_LIST[@]}"; do
        log "$BLUE" "  ${item}"
    done
    echo
fi

log "$BLUE" "======================================"
echo

# Optional final message using `error` flag if you keep it
if [[ ${error:-0} -ne 1 ]]; then
  log "$GREEN" "-- ------------------ --"
  log "$GREEN" "DWI QC finished without fatal errors."
  log "$GREEN" "-- ------------------ --"
  echo
fi