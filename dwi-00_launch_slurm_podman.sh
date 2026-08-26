#!/bin/bash

########################
# SLURM DIRECTIVES
########################

#SBATCH --job-name=dwipipeline
#SBATCH --partition=<your partition>
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --mem=20M
#SBATCH --time=00:10:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-1%1
# NOTE: --array above is a fallback only. Command-line flags override this.

set -euo pipefail

ml podman   # uncomment/adjust if podman is provided as an environment module
containerimage=docker.io/cvriend/tractoprep:v1.0.6

Usage() {
  echo "Usage: $0 <path_to_spec.json> <host_bidsdir> <host_outputdir> <host_workdir> [<host_freesurferdir>]"
  echo ""
  echo "  path_to_spec.json : template spec file. The bidsdir/outputdir/workdir"
  echo "                      keys inside this file are  CONTAINER-side paths"
  echo "                      (i.e. what the pipeline sees once inside the container)."
  echo "  host_bidsdir      : path to the BIDS dataset on the HOST filesystem"
  echo "  host_outputdir    : path to the output dir on the HOST filesystem"
  echo "  host_workdir      : path to the work dir on the HOST filesystem"
  echo "  host_freesurferdir : (optional) path to the existing freesurferdir on the HOST filesystem"

  exit 1
}

# 4 required args, 5th (freesurferdir) is optional
if [[ $# -lt 4 || $# -gt 5 ]]; then
  Usage
fi

templatejson=$1
host_bidsdir=$2
host_outputdir=$3
host_workdir=$4
# optional 5th arg -- must default with :- under `set -u`, otherwise
# referencing $5 when it wasn't passed is an unbound-variable error.
host_freesurferdir="${5:-}"


###########################################################
# Per-stage resource settings
PREPROC_CPUS=8;   PREPROC_MEM=28G; PREPROC_TIME=07:00:00
TRACTO_CPUS=16;   TRACTO_MEM=4G;   TRACTO_TIME=02:00:00
QC_CPUS=1;        QC_MEM=4G;       QC_TIME=00:10:00
###########################################################


# read spec.json file
# NOTE: bidsdir/outputdir/workdir declared here are CONTAINER-side paths,
# used only in the podman command line and internally by the pipeline.
for key in $(jq -r 'keys[]' "${templatejson}"); do
  value=$(jq -r --arg k "$key" '.[$k]' "${templatejson}")
  declare "$key"="$value"
done

# get rid of subj variable from spec.json, since we will override it with the SLURM_ARRAY_TASK_ID
unset subj

# check container-side path variables (from spec.json) are non-empty
for var in bidsdir outputdir workdir nstreamlines; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty in ${templatejson}."
    exit 1
  fi
done

# check host-side path arguments are non-empty
for var in host_bidsdir host_outputdir host_workdir; do
  if [[ -z "${!var}" ]]; then
    echo "Error: '$var' argument is not set or is empty."
    exit 1
  fi
done

# if a host freesurferdir was given, spec.json must also define the
# matching container-side path (freesurferdir key), and the host dir
# must actually exist. ${freesurferdir:-} guards against unbound-variable
# errors under `set -u` if the key is absent from spec.json.
if [[ -n "${host_freesurferdir}" ]]; then
  if [[ -z "${freesurferdir:-}" ]]; then
    echo "Error: host_freesurferdir was given, but '${templatejson}' has no 'freesurferdir' key (container-side path)."
    exit 1
  fi
  if [[ ! -d "${host_freesurferdir}" ]]; then
    echo "Error: host_freesurferdir '${host_freesurferdir}' does not exist or is not a directory."
    exit 1
  fi
fi

mkdir -p "${host_workdir}" "${host_outputdir}"

cd "${host_bidsdir}"

# FIX: build a proper bash array of subject dirs instead of piping
# `ls` output through sed on a variable (fragile, breaks on odd names).
mapfile -t subjects < <(find . -maxdepth 1 -mindepth 1 -type d -name 'sub-*' -printf '%f\n' | sort)

if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${#subjects[@]}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) exceeds number of subjects (${#subjects[@]})."
    exit 1
fi

# arrays are 0-indexed, SLURM_ARRAY_TASK_ID starts at 1
subj="${subjects[$((SLURM_ARRAY_TASK_ID - 1))]}"

# basedir is where we're writing spec_${subj}.json. Podman does NOT
# auto-mount the cwd like apptainer does, so we bind it explicitly.
basedir="$(pwd)"

jq --arg subj "$subj" '.subj = $subj' "${templatejson}" > "spec_${subj}.json"
subjspecjson="${basedir}/spec_${subj}.json"

echo "Submitting preprocessing and tracto jobs for ${subj} ${session:-}..."

sessionpath="/${session:+${session}/}"
sessionfile="_${session:+${session}_}"

# Logs are written directly by SLURM on the compute node, so they use
# the HOST output path.
logdir="${host_outputdir}/logs/${subj}${sessionpath}"
mkdir -p "${logdir}"

# bind mounts (podman -v flags): HOST path : CONTAINER path
mountargs="-v ${host_bidsdir}:${bidsdir} -v ${host_workdir}:${workdir} -v ${host_outputdir}:${outputdir} -v ${basedir}:${basedir}"
if [[ -n "${host_freesurferdir}" ]]; then
  mountargs="${mountargs} -v ${host_freesurferdir}:${freesurferdir}"
fi

# Rootless podman: --userns=keep-id maps your host uid/gid into the
# container 1:1, so files written by the pipeline stay owned by you
# on the host (equivalent purpose to docker's --user $(id -u):$(id -g),
# but the correct/idiomatic way to do it under rootless podman).
useropt="--userns=keep-id"

# ============================================================
# Expected output files (used for skip-logic)
# NOTE: checked against HOST paths, since that's where the bind-mounted
# files physically land. The pipeline itself writes to the container path.
# ============================================================
preproc_nifti="${host_outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz"
preproc_bvec="${host_outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec"
preproc_qc="${host_outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz"
preproc_anat="${host_outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_template.nii.gz"
preproc_atlas="${host_outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-400P17N_dseg.nii.gz"

tracto_file="${host_outputdir}/dwi-tracto/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck"
conn_file="${host_outputdir}/dwi-tracto/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-400P17N_desc-streams_connmatrix.csv"

# ============================================================
# Helper: submit a job, optionally depending on a previous job
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
        --partition=${SLURM_JOB_PARTITION} \
        --time="${time}" \
        --mem="${mem}" \
        --cpus-per-task="${cpus}" \
        --output="${logdir}/${jobname}%j.log" \
        "${dep_arg[@]}" \
        --wrap="${cmd}"
}

last_job=""

# ============================================================
# STAGE 1: dwi-preproc + eddy (skip if outputs already exist)
# ============================================================
if [ -f "${preproc_nifti}" ] && [ -f "${preproc_bvec}" ] && [ -f "${preproc_qc}" ] && [ -f "${preproc_anat}" ] && [ -f "${preproc_atlas}" ]; then
    echo "Preprocessing already done for ${subj} ${session:-}, skipping preproc job."
    job_id_preproc=""
else
    cmd_preproc="tmpdir_job=\"\${SLURM_TMPDIR:-${host_workdir}/tmp}/${subj}${sessionfile}\${SLURM_JOB_ID}\"; \
    mkdir -p \"\${tmpdir_job}\"; \
    podman run --rm ${useropt} \
      ${mountargs} -v \${tmpdir_job}:/scratch \
      -e TMPDIR=/scratch -e TMP=/scratch -e TEMP=/scratch \
      ${containerimage} dwi-preproc ${subjspecjson}"

    job_id_preproc=$(submit_job "dwi-preproc_${subj}${sessionfile}" "${cmd_preproc}" "${last_job}" "${PREPROC_CPUS}" "${PREPROC_MEM}" "${PREPROC_TIME}")
    echo "Submitted dwi-preproc job: ${job_id_preproc}"
    last_job="${job_id_preproc}"
fi

# ============================================================
# STAGE 2: dwi-qc (always runs)
# ============================================================
cmd_qc="podman run --rm ${useropt} ${mountargs} ${containerimage} dwi-qc ${subjspecjson}"
job_id_qc=$(submit_job "dwi-qc_${subj}${sessionfile}" "${cmd_qc}" "${last_job}" "${QC_CPUS}" "${QC_MEM}" "${QC_TIME}")
echo "Submitted dwi-qc job: ${job_id_qc}"
last_job="${job_id_qc}"

# ============================================================
# STAGE 3: dwi-tracto (skip if outputs already exist)
# ============================================================
if [ -f "${tracto_file}" ] && [ -f "${conn_file}" ]; then
    echo "Tractography already done for ${subj} ${session:-}, skipping tracto job."
    job_id_tracto=""
else
  cmd_tracto="tmpdir_job=\"\${SLURM_TMPDIR:-${host_workdir}/tmp}/${subj}${sessionfile}\${SLURM_JOB_ID}\"; \
  mkdir -p \"\${tmpdir_job}\"; \
  podman run --rm ${useropt} \
    ${mountargs} -v \${tmpdir_job}:/scratch \
    -e TMPDIR=/scratch -e TMP=/scratch -e TEMP=/scratch \
    ${containerimage} dwi-tracto ${subjspecjson}"
    job_id_tracto=$(submit_job "dwi-tracto_${subj}${sessionfile}" "${cmd_tracto}" "${last_job}" "${TRACTO_CPUS}" "${TRACTO_MEM}" "${TRACTO_TIME}")
    echo "Submitted dwi-tracto job: ${job_id_tracto}"
    last_job="${job_id_tracto}"
fi

# ============================================================
# STAGE 4: dwi-qc (always runs)
# ============================================================
cmd_qc="podman run --rm ${useropt} ${mountargs} ${containerimage} dwi-qc ${subjspecjson}"
job_id_qc=$(submit_job "dwi-qc_${subj}${sessionfile}" "${cmd_qc}" "${last_job}" "${QC_CPUS}" "${QC_MEM}" "${QC_TIME}")
echo "Submitted dwi-qc job: ${job_id_qc}"
last_job="${job_id_qc}"

echo "All jobs submitted for ${subj}. Final job in chain: ${last_job}"