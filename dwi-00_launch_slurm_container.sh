#!/bin/bash

########################
# SLURM DIRECTIVES
########################

#SBATCH --job-name=dwipipeline
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --mem=20M
#SBATCH --time=00:10:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=2-3%3
# NOTE: --array above is a fallback only. In practice this script is
# submitted via submit_pipeline.sh, which passes --array=1-N%throttle
# on the sbatch command line, and command-line flags override this.

set -euo pipefail

ml apptainer
templatejson=$1
containerpath=/net/beegfs/users/P042819/TractoFriend.sif
#hostworkdir=/scratch/users/P042819/work

#export FSLOUTPUTTYPE=NIFTI_GZ
#export APPTAINER_BINDPATH="/net/beegfs/users/P042819,/home/P042819,/scratch/users/P042819"


Usage() {
  echo "Usage: $0 <path_to_spec.json>"
  exit 1
}

log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

if [[ $# -ne 1 ]]; then
  Usage
fi

# read spec.json file
for key in $(jq -r 'keys[]' ${templatejson}); do
  value=$(jq -r --arg k "$key" '.[$k]' ${templatejson})
  declare "$key"="$value"
done

# get rid of subj variable from spec.json, since we will override it with the SLURM_ARRAY_TASK_ID
unset subj

# check if all variables are non-empty
for var in bidsdir outputdir workdir nstreamlines; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  

mkdir -p ${workdir} ${outputdir}



cd ${bidsdir}

# FIX: build a proper bash array of subject dirs instead of piping
# `ls` output through sed on a variable (fragile, breaks on odd names).
mapfile -t subjects < <(find . -maxdepth 1 -mindepth 1 -type d -name 'sub-*' -printf '%f\n' | sort)

if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${#subjects[@]}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) exceeds number of subjects (${#subjects[@]})."
    exit 1
fi

# arrays are 0-indexed, SLURM_ARRAY_TASK_ID starts at 1
subj="${subjects[$((SLURM_ARRAY_TASK_ID - 1))]}"

jq --arg subj "$subj" '.subj = $subj' ${templatejson} > "spec_${subj}.json"
subjspecjson="spec_${subj}.json"

echo "Submitting preprocessing and tracto jobs for ${subj} ${session:-}..."

sessionpath="${session:+/${session}/}"
sessionfile="${session:+_${session}_}"

logdir="${outputdir}/logs/${subj}${sessionpath}"
tmpdir="${workdir}/${subj}${sessionpath}tmp"
mkdir -p "${logdir}" "${tmpdir}"


# define bind folders for apptainer 
bindcmd="${bidsdir},${workdir},${outputdir},${tmpdir}:/scratch"


# ============================================================
# Expected output files (used for skip-logic)
# ============================================================
preproc_nifti="${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz"
preproc_bvec="${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec"
preproc_qc="${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz"
preproc_anat="${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_template.nii.gz"
preproc_atlas="${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-300P7N_dseg.nii.gz"

tracto_file="${outputdir}/dwi-tracto/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck"
conn_file="${outputdir}/dwi-tracto/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-300P7N_desc-streams_connmatrix.csv"

# ============================================================
# Helper: submit a job, optionally depending on a previous job
#   $1 = job name
#   $2 = command string to run (the apptainer call)
#   $3 = job id to depend on (may be empty -> no dependency)
#   $4 = cpus-per-task
#   $5 = mem
#   $6 = time
# ============================================================
submit_job () {
    local jobname="$1"
    local cmd="$2"
    local dep_id="$3"
    local cpus="$4"
    local mem="$5"
    local time="$6"

    local dep_arg=()
    if [ -n "${dep_id}" ]; then
        dep_arg=(--kill-on-invalid-dep="yes" --dependency="afterok:${dep_id}")
    fi

    sbatch --parsable \
        --job-name="${jobname}" \
        --partition=defq \
        --time="${time}" \
        --mem="${mem}" \
        --cpus-per-task="${cpus}" \
        --output="${logdir}/${jobname}%j.log" \
        "${dep_arg[@]}" \
        --wrap="${cmd}"
}

# Per-stage resource settings
PREPROC_CPUS=8;   PREPROC_MEM=28G; PREPROC_TIME=07:00:00
TRACTO_CPUS=16;   TRACTO_MEM=12G;   TRACTO_TIME=04:00:00
QC_CPUS=1;        QC_MEM=2G;       QC_TIME=00:30:00

last_job=""

# ============================================================
# STAGE 1: dwi-preproc + eddy (skip if outputs already exist)
# ============================================================
if [ -f "${preproc_nifti}" ] && [ -f "${preproc_bvec}" ] && [ -f "${preproc_qc}" ]  && [ -f "${preproc_anat}" ]  && [ -f "${preproc_atlas}" ]; then
    echo "Preprocessing already done for ${subj} ${session:-}, skipping preproc job."
    job_id_preproc=""
else
    cmd_preproc="apptainer run --cleanenv --bind ${bindcmd} --env TMPDIR=${tmpdir} --env TMP=${tmpdir} --env TEMP=${tmpdir} ${containerpath} dwi-preproc ${subjspecjson}"
    job_id_preproc=$(submit_job "dwi-preproc_${subj}${sessionfile}" "${cmd_preproc}" "${last_job}" "${PREPROC_CPUS}" "${PREPROC_MEM}" "${PREPROC_TIME}")
    echo "Submitted dwi-preproc job: ${job_id_preproc}"
    last_job="${job_id_preproc}"
fi

# ============================================================
# STAGE 2: dwi-tracto (skip if outputs already exist)
# ============================================================
if [ -f "${tracto_file}" ] && [ -f "${conn_file}" ]; then
    echo "Tractography already done for ${subj} ${session:-}, skipping tracto job."
    job_id_tracto=""
else
   cmd_tracto="apptainer run --cleanenv --bind ${bindcmd} --env TMPDIR=${tmpdir} --env TMP=${tmpdir} --env TEMP=${tmpdir} ${containerpath} dwi-tracto ${subjspecjson}"
    job_id_tracto=$(submit_job "dwi-tracto_${subj}${sessionfile}" "${cmd_tracto}" "${last_job}" "${TRACTO_CPUS}" "${TRACTO_MEM}" "${TRACTO_TIME}")
    echo "Submitted dwi-tracto job: ${job_id_tracto}"
    last_job="${job_id_tracto}"
fi

# ============================================================
# STAGE 3: dwi-qc (always runs)
# ============================================================
cmd_qc="apptainer run --cleanenv --bind ${bindcmd} --env TMPDIR=${tmpdir} --env TMP=${tmpdir} --env TEMP=${tmpdir} ${containerpath} dwi-qc ${subjspecjson}"
job_id_qc=$(submit_job "dwi-qc_${subj}${sessionfile}" "${cmd_qc}" "${last_job}" "${QC_CPUS}" "${QC_MEM}" "${QC_TIME}")
echo "Submitted dwi-qc job: ${job_id_qc}"
last_job="${job_id_qc}"


echo "All jobs submitted for ${subj}. Final job in chain: ${last_job}"
