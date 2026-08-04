#!/bin/bash

########################
# SLURM DIRECTIVES
########################

#SBATCH --job-name=dwipipeline
#SBATCH --partition=luna-cpu-short
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=20M
#SBATCH --time=04:00:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-1%1

set -euo pipefail

jsonfile=$1
scriptdir=/tracto

containerpath=
templatejson=template.json


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
for var in bidsdir outputdir workdir scriptdir nstreamlines; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  


cd ${bidsdir}
subjects=$(ls -d sub-*)

subj=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${subjects}")


jq --arg subj "$subj" '.subj = $subj' ${templatejson} > "spec_${subj}.json"


echo "Submitting preprocessing and tracto jobs for ${subj}"

preproc_jobs=()
tck_jobs=()
tck_submitted=()



 # Check if preprocessing is already done
    if [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" ] &&
       [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" ] &&
       [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz" ]; then

        echo "Preprocessing already done for ${subj} ${session:-}, skipping preproc and eddy."
        preproc_jobs+=("")

 # Submit preproc job
        job_id_preproc=$(sbatch --parsable \
            "${scriptdir}/dwi-02a-preproc.sh" \
            ${subjspecjson}
        )
        preproc_jobs+=("${job_id_preproc}")