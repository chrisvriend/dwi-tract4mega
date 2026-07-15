#!/bin/bash

jsonfile=$1
scriptdir=/tracto

export FSLOUTPUTTYPE=NIFTI_GZ

# Color variables
BLUE='\033[0;34m'
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

Usage() {
  echo "Usage: $0 <path_to_spec.json>"
  exit 1
}

# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

if [[ $# -ne 1 ]]; then
  Usage
fi

# read spec.json file and export all variables
for key in $(jq -r 'keys[]' ${jsonfile}); do
  value=$(jq -r --arg k "$key" '.[$k]' ${jsonfile})
  declare -x "$key"="$value"
done

# check if all variables are non-empty
for var in subj bidsdir outputdir workdir freesurferdir eddy_method nthreads; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done

echo
log "$BLUE" "-- ------------------ --"
log "$BLUE" "Starting DWI preprocessing for subject: ${subj} ${session}"
log "$BLUE" "-- ------------------ --"
echo


run_03_and_02b() {
  local session_flag="${1:-}"

  # determine nthreads
  if (( nthreads >= 8 )); then
    eddy_threads=4
    anat2dwi_threads=4
  elif (( nthreads > 4 )); then
    eddy_threads=4
    anat2dwi_threads=$((nthreads - 4))
  elif (( nthreads < 2 )); then
    log "$RED" "Error: At least 2 threads are required to run the processes."
    exit 1
  else
    eddy_threads=$nthreads
    anat2dwi_threads=1
  fi

  ${scriptdir}/dwi-03-anat2dwi_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" -f "${freesurferdir}" -c "${scriptdir}" -t "${anat2dwi_threads}" ${session_flag} > "${outputdir}/dwi-preproc/${subj}/log/${subj}_anat2dwi_$(date +"%Y-%m-%d_%H-%M").log" 2>&1 &
  pid1=$!
  ${scriptdir}/dwi-02b-eddyCPU_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" -m "${eddy_method}" -t "${eddy_threads}" ${session_flag} > "${outputdir}/dwi-preproc/${subj}/log/${subj}_eddy_$(date +"%Y-%m-%d_%H-%M").log" 2>&1 &
  pid2=$!

  wait $pid1
  wait $pid2
}


mkdir -p "${outputdir}/dwi-preproc/${subj}/log"
error=0

if [[ -z "${session}" ]]; then
  # preprocess
  log "$BLUE" "Step 02a (see log file in ${outputdir}/dwi-preproc/${subj}/log)"
  echo
  ${scriptdir}/dwi-02a-preproc_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -c "${scriptdir}" -s "${subj}" -t "${nthreads}" > "${outputdir}/dwi-preproc/${subj}/log/${subj}_preproc_$(date +"%Y-%m-%d_%H-%M").log" 2>&1
  status=$?
  if [[ $status -eq 0 ]]; then
    log "$BLUE" "Steps 02b / 03 (see log file in ${outputdir}/dwi-preproc/${subj}/log)"
    echo
    run_03_and_02b
  else
    echo "dwi-preproc failed, skipping eddy and anat2dwi"
    exit 1
  fi
else
  # preprocess
  log "$BLUE" "Step 02a (see log file in ${outputdir}/dwi-preproc/${subj}/log)"
  echo
  ${scriptdir}/dwi-02a-preproc_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -c "${scriptdir}" -s "${subj}" -z "${session}" -t "${nthreads}" > "${outputdir}/dwi-preproc/${subj}/log/${subj}_preproc_$(date +"%Y-%m-%d_%H-%M").log" 2>&1
  status=$?
  if [[ $status -eq 0 ]]; then
    log "$BLUE" "Steps 02b / 03 (see log file in ${outputdir}/dwi-preproc/${subj}/log)"
    echo
    run_03_and_02b "-z ${session}"
  else
    echo "dwi-preproc failed, skipping eddy and anat2dwi"
    exit 1
  fi
fi

## run checks before clean-up

# Set session path/file
if [[ -z "${session}" ]]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
fi

files=$(echo "
    ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz
    ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec
    ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval
    ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz
    ${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-300P7N_dseg.nii.gz
    ${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz
    ${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-gmwm_probseg.nii.gz
    ")

for file in ${files}; do
  if [ ! -f "${file}" ]; then
    echo
    log "${RED}" "!!!ERROR!!!"
    log "${RED}" "a scan was not found in the output folder"
    echo "${file}"
    error=1
  fi
done

if [[ ${error} -ne 1 ]]; then
  echo
  log "$GREEN" "-- ------------------ --"
  log "$GREEN" "DWI preprocessing completed for subject: ${subj} ${session}"
  log "$GREEN" "-- ------------------ --"
  echo
  rm -rf "${workdir}/${subj}"
fi