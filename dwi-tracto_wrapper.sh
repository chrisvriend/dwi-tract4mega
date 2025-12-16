#!/bin/bash

jsonfile=$1


export FSLOUTPUTTYPE=NIFTI_GZ

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


# read spec.json file
for key in $(jq -r 'keys[]' ${jsonfile}); do
  value=$(jq -r --arg k "$key" '.[$k]' ${jsonfile})
  declare "$key"="$value"
done

# check if all variables are non-empty
for var in subj bidsdir outputdir workdir scriptdir nstreamlines; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  

echo $subj
echo $session


echo
log "$BLUE" "-- ------------------ --"
log "$BLUE" "Starting DWI Tractography for subject: ${subj} ${session}"
log "$BLUE" "-- ------------------ --"
echo

 
mkdir -p "${outputdir}/dwi-tracto/"
mkdir -p "${outputdir}/dwi-tracto/${subj}/log"

if [[ -z "${session}" ]]; then
  # tractography without session
  ${scriptdir}/dwi-04a-connectome_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" -t "${nthreads}" -x "${nstreamlines}" > ${outputdir}/dwi-tracto/${subj}/log/${subj}_tracto_$(date +"%Y-%m-%d_%H-%M").log
  status=$?
  if [[ $status -eq 0 ]]; then

 ${scriptdir}/dwi-04b-tracts2conn_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" -t "${nthreads}" -x "${nstreamlines}" > ${outputdir}/dwi-tracto/${subj}/log/${subj}_conn_$(date +"%Y-%m-%d_%H-%M").log

  else
    echo "dwi-tracto failed, skipping tck2conn"
    exit 1
  fi

else
  # tractography with session
  ${scriptdir}/dwi-04a-connectome_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" -t "${nthreads}" -x "${nstreamlines}" -z ${session} > ${outputdir}/dwi-tracto/${subj}/log/${subj}_conn_$(date +"%Y-%m-%d_%H-%M").log
  
  status=$?
  if [[ $status -eq 0 ]]; then
 ${scriptdir}/dwi-04b-tracts2conn_container.sh -i "${bidsdir}" -o "${outputdir}" -w "${workdir}" -s "${subj}" -t "${nthreads}" -x "${nstreamlines}" -z ${session} > ${outputdir}/dwi-tracto/${subj}/log/${subj}_conn_$(date +"%Y-%m-%d_%H-%M").log
  else
    echo "dwi-tracto failed, skipping tck2conn"
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


    # files=$(echo "
    # ${outputdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck
    # ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz
    # ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec
    # ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval
    # ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz
    # ${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-300P7N_dseg.nii.gz
    # ${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz
    # ${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-gmwm_probseg.nii.gz
    # ")

    # for file in ${files}; do

    #   if [ ! -f ${file} ]; then
    #     log "${RED}" "!!!ERROR!!!"
    #     log "${RED}" "a scan was not found in the output folder"
    #     echo "${file}"
    #     error=1

    #   fi

    # done



 if [[ ${error} -ne 1 ]]; then
  log "$GREEN" "-- ------------------ --"
  log "$GREEN" "DWI preprocessing completed for subject: ${subj} ${session}"
  log "$GREEN" "-- ------------------ --"
  echo
  #rm -rf ${workdir}/${subj}
 fi


