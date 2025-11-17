#!/bin/bash

jsonfile=$1


export FSLOUTPUTTYPE=NIFTI_GZ

usage() {
  echo "Usage: $0 <path_to_spec.json>"
  exit 1
}
if [[ $# -ne 1 ]]; then
  usage
fi

# read spec.json file
for key in $(jq -r 'keys[]' ${jsonfile}); do
  value=$(jq -r --arg k "$key" '.[$k]' ${jsonfile})
  declare "$key"="$value"
done

# check if all variables are non-empty
for var in subj bidsdir outputdir workdir scriptdir freesurferdir eddy_method; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  



# Function to allow dwi-anat and dwi-eddy to run in parallel
run_03_and_02b() {
  # 03 alleen als Freesurfer data niet bestaat
  if [[ ! -d ${freesurferdir}/${subj} || ! -f ${freesurferdir}/${subj}/mri/aseg.mgz ]]; then
    echo "Freesurfer directory does not contain data for ${subj}. Will run FastSurfer"
    ${scriptdir}/dwi-03-anat2dwi_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" $1 -f "${freesurferdir}" &
    pid1=$!
  fi

  # 02b altijd uitvoeren
  ${scriptdir}/dwi-02b-eddyCPU_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" -m "${eddy_method}" $1 &
  pid2=$!

  # Wacht op beide processen (indien gestart)
  if [[ ! -z "$pid1" ]]; then
    wait $pid1
  fi
  wait $pid2
}
###

if [[ -z "${session}" ]]; then
  # preprocess
  ${scriptdir}/dwi-02a-preproc_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -c "${scriptdir}" -s "${subj}"
  status=$?
  if [[ $status -eq 0 ]]; then
    run_03_and_02b ""
  else
    echo "dwi-preproc failed, skipping 03 and 02b"
    exit 1
  fi
else
  # preprocess
  ${scriptdir}/dwi-02a-preproc_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -c "${scriptdir}" -s "${subj}" -z "${session}"
  status=$?
  if [[ $status -eq 0 ]]; then
    run_03_and_02b "-z \"${session}\""
  else
    echo "02a-preproc failed, skipping 03 and 02b"
    exit 1
  fi
fi
