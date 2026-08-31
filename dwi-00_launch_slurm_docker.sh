#!/bin/bash

########################
# SLURM DIRECTIVES
########################

#SBATCH --job-name=dwipipeline
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --mem=20M
#SBATCH --time=12:00:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-1%1
# NOTE: --array above is a fallback only. Command-line flags override this.

set -euo pipefail

# Set your Docker image reference here
containerimage="docker.io/cvriend/tractoprep:v1.0.6"

run_user="$(id -u):$(id -g)"

# Set to ":z" if your cluster enforces SELinux (Docker supports the same
# volume relabeling suffixes as Podman). Leave empty otherwise.
selinux_suffix=""

# # ensure that the container image is available locally (pull if not; default memory specs might not be sufficient)
# if ! docker image inspect "${containerimage}" >/dev/null 2>&1; then
#     echo "Container image ${containerimage} not found locally. Attempting to pull..."
#     docker pull "${containerimage}" || { echo "Failed to pull container image ${containerimage}. Exiting."; exit 1; }
# fi


Usage() {
  echo "Usage: $0 <path_to_spec.json> <host_bidsdir> <host_outputdir> <host_workdir> [<host_freesurferdir>] [--preproc-only]"
  echo ""
  echo "  path_to_spec.json  : template spec file. The bidsdir/outputdir/workdir"
  echo "                       keys inside this file are CONTAINER-side paths"
  echo "                       (i.e. what the pipeline sees once inside docker)."
  echo "  host_bidsdir       : path to the BIDS dataset on the HOST filesystem"
  echo "  host_outputdir     : path to the output dir on the HOST filesystem"
  echo "  host_workdir       : path to the work dir on the HOST filesystem"
  echo "  host_freesurferdir : (optional) path to the existing freesurferdir on the HOST filesystem"
  echo "  --preproc-only     : (optional flag) run only stages 02a, 02b, 03 and dwi-qc;"
  echo "                       skip tractography stages 04a and 04b"
  exit 1
}

# Parse arguments: up to 5 positional args + optional --preproc-only flag
templatejson=""
host_bidsdir=""
host_outputdir=""
host_workdir=""
host_freesurferdir=""

## preproc_only flag 0/1: if set, skip tractography stages
preproc_only=0

positional=()
for arg in "$@"; do
  case "${arg}" in
    --preproc-only)
      preproc_only=1
      ;;
    --*)
      echo "Error: Unknown option '${arg}'"
      Usage
      ;;
    *)
      positional+=("${arg}")
      ;;
  esac
done

if [[ ${#positional[@]} -ne 5 ]]; then
  Usage
fi

templatejson="${positional[0]}"
host_bidsdir="${positional[1]}"
host_outputdir="${positional[2]}"
host_workdir="${positional[3]}"
host_freesurferdir="${positional[4]}"

##########################################################
# Per-stage resource settings
PREPARE_CPUS=4;   PREPARE_MEM=8G; PREPARE_TIME=01:00:00
EDDY_CPUS=4;      EDDY_MEM=4G;    EDDY_TIME=06:00:00
ANAT_CPUS=4;      ANAT_MEM=28G;     ANAT_TIME=06:00:00
TRACTO_CPUS=16;   TRACTO_MEM=4G;   TRACTO_TIME=02:00:00
CONN_CPUS=1;      CONN_MEM=2G;     CONN_TIME=01:00:00
QC_CPUS=1;        QC_MEM=4G;       QC_TIME=00:10:00
##########################################################

# read spec.json file
# NOTE: bidsdir/outputdir/workdir declared here are CONTAINER-side paths,
# used only in the docker command line and internally by the pipeline.
for key in $(jq -r 'keys[]' "${templatejson}"); do
  value=$(jq -r --arg k "$key" '.[$k]' "${templatejson}")
  declare "$key"="$value"
done

# get rid of subj variable from spec.json, since we will override it with the SLURM_ARRAY_TASK_ID
unset subj

# check container-side path variables (from spec.json) are non-empty
for var in bidsdir outputdir workdir freesurferdir nstreamlines eddy_method; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty in ${templatejson}."
    exit 1
  fi
done

# check host-side path arguments are non-empty
for var in host_bidsdir host_outputdir host_workdir host_freesurferdir; do
  if [[ -z "${!var}" ]]; then
    echo "Error: '$var' argument is not set or is empty."
    exit 1
  fi
done

mkdir -p "${host_workdir}" "${host_outputdir}" "${host_freesurferdir}"

cd "${host_bidsdir}"

# Build a proper bash array of subject dirs
mapfile -t subjects < <(find . -maxdepth 1 -mindepth 1 -type d -name 'sub-*' -printf '%f\n' | sort)

if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${#subjects[@]}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) exceeds number of subjects (${#subjects[@]})."
    exit 1
fi

# arrays are 0-indexed, SLURM_ARRAY_TASK_ID starts at 1
subj="${subjects[$((SLURM_ARRAY_TASK_ID - 1))]}"

basedir="$(pwd)"

jq --arg subj "$subj" '.subj = $subj' "${templatejson}" > "spec_${subj}.json"
subjspecjson="${basedir}/spec_${subj}.json"

if [[ "${preproc_only}" -eq 1 ]]; then
    echo "Mode: preproc-only (stages denoise + topup + eddy + anat2dwi + dwi-qc; skipping tractography stages)"
else
    echo "Mode: full pipeline (stages denoise + topup + eddy + anat2dwi + tractography + connmatrix + dwi-qc)"
fi
echo "Submitting jobs for ${subj} ${session:-}..."

sessionpath="/${session:+${session}/}"
sessionfile="_${session:+${session}_}"

# Logs are written directly by SLURM on the compute node, so they use
# the HOST output path.
logdir="${host_outputdir}/logs/${subj}${sessionpath}"
mkdir -p "${logdir}"

# Volume mounts: HOST path:CONTAINER path (container path comes from spec.json)
# The optional SELinux relabel suffix (":z") is appended when selinux_suffix is set.
volcmd="${host_bidsdir}:${bidsdir}${selinux_suffix},${host_workdir}:${workdir}${selinux_suffix},${host_outputdir}:${outputdir}${selinux_suffix},${basedir}:${basedir}${selinux_suffix},${host_freesurferdir}:${freesurferdir}${selinux_suffix}"
# Build a single-line string of -v flags
volflags=$(echo "${volcmd}" | tr ',' '\n' | sed 's/^/-v /' | tr '\n' ' ')

# ============================================================
# Expected output files (used for skip-logic)
# NOTE: checked against HOST paths, since that's where the bind-mounted
# files physically land. The pipeline itself writes to the container path.
# ============================================================

# 02a: topup / denoising outputs
topup_nifti="${host_outputdir}/dwi-preproc/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-unwarped_epi.nii.gz"
degibbs_nifti="${host_outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz"

# 02b: eddy outputs
preproc_nifti="${host_outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz"
preproc_bvec="${host_outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec"
preproc_qc="${host_outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz"

# 03: anat-to-dwi outputs
preproc_anat="${host_outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_template.nii.gz"
preproc_atlas="${host_outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-400P17N_dseg.nii.gz"

# 04a: tractogram
tracto_file="${host_outputdir}/dwi-tracto/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck"

# 04b: connectivity matrix
conn_file="${host_outputdir}/dwi-tracto/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-400P17N_desc-streams_connmatrix.csv"

# ============================================================
# Helper: submit a job, optionally depending on one or more previous jobs
#   $1 = job name
#   $2 = command string to run (the docker call)
#   $3 = dependency string: colon-separated job IDs, or empty for no dependency
#   $4 = cpus-per-task
#   $5 = mem
#   $6 = time
# ============================================================
submit_job () {
    local jobname="$1"
    local cmd="$2"
    local dep_ids="$3"   # colon-separated list of job IDs, or empty
    local cpus="$4"
    local mem="$5"
    local time="$6"

    local dep_arg=()
    if [[ -n "${dep_ids}" ]]; then
        dep_arg=(--kill-on-invalid-dep="yes" --dependency="afterok:${dep_ids}")
    fi

    sbatch --parsable \
        --job-name="${jobname}" \
        --partition=${SLURM_JOB_PARTITION} \
        --time="${time}" \
        --mem="${mem}" \
        --cpus-per-task="${cpus}" \
        --output="${logdir}/${jobname}_%j.log" \
        "${dep_arg[@]}" \
        --wrap="${cmd}"
}

# Helper: build a colon-separated dependency string from an array of job IDs,
# filtering out any empty entries.
# Usage: dep_str=$(make_dep_string "${id1}" "${id2}" ...)
make_dep_string () {
    local ids=()
    for id in "$@"; do
        if [[ "${id}" =~ ^[0-9]+$ ]]; then
            ids+=("${id}")
        fi
    done
    local IFS=":"
    echo "${ids[*]}"
}

all_jobs=()

# ============================================================
# STAGE 02a: dwi-prepare (topup / denoising)
# ============================================================
if [[ -f "${topup_nifti}" && -f "${degibbs_nifti}" ]]; then
    echo "denoise and topup outputs already exist for ${subj} ${session:-}, skipping dwi-prepare."
    job_id_02a=""
else
    cmd_02a="tmpdir_job=\"\${SLURM_TMPDIR:-${host_workdir}/tmp}/${subj}${sessionfile}\${SLURM_JOB_ID}\"; \
    mkdir -p \"\${tmpdir_job}\"; \
    docker run --rm \
      --pull=always --user=\"${run_user}\" \
      ${volflags} \
      --env TMPDIR=/scratch --env TMP=/scratch --env TEMP=/scratch --env HOME=/tmp \
      -v \"\${tmpdir_job}:/scratch\" \
      ${containerimage} dwi-prepare ${subjspecjson}"

    job_id_02a=$(submit_job "dwi-prepare_${subj}${sessionfile}" "${cmd_02a}" "" \
        "${PREPARE_CPUS}" "${PREPARE_MEM}" "${PREPARE_TIME}")
    echo "Submitted dwi-prepare (02a) job: ${job_id_02a}"
    all_jobs+=("${job_id_02a}")
fi

# ============================================================
# STAGE 02b: dwi-eddy  — depends on 02a
# STAGE 03:  dwi-anat2dwi — depends on 02a  (runs in parallel with 02b)
# ============================================================
dep_02a=$(make_dep_string "${job_id_02a:-}")

# --- 02b ---
if [[ -f "${preproc_nifti}" && -f "${preproc_bvec}" && -f "${preproc_qc}" ]]; then
    echo "eddy outputs already exist for ${subj} ${session:-}, skipping dwi-eddy."
    job_id_02b=""
else
    cmd_02b="tmpdir_job=\"\${SLURM_TMPDIR:-${host_workdir}/tmp}/${subj}${sessionfile}\${SLURM_JOB_ID}\"; \
    mkdir -p \"\${tmpdir_job}\"; \
    docker run --rm \
      --pull=always --user=\"${run_user}\" \
      ${volflags} \
      --env TMPDIR=/scratch --env TMP=/scratch --env TEMP=/scratch --env HOME=/tmp \
      -v \"\${tmpdir_job}:/scratch\" \
      ${containerimage} dwi-eddy ${subjspecjson}"

    job_id_02b=$(submit_job "dwi-eddy_${subj}${sessionfile}" "${cmd_02b}" "${dep_02a}" \
        "${EDDY_CPUS}" "${EDDY_MEM}" "${EDDY_TIME}")
    echo "Submitted dwi-eddy (02b) job: ${job_id_02b}"
    all_jobs+=("${job_id_02b}")
fi

# --- 03 ---
if [[ -f "${preproc_anat}" && -f "${preproc_atlas}" ]]; then
    echo "anat2dwi outputs already exist for ${subj} ${session:-}, skipping dwi-anat2dwi."
    job_id_03=""
else
    cmd_03="tmpdir_job=\"\${SLURM_TMPDIR:-${host_workdir}/tmp}/${subj}${sessionfile}\${SLURM_JOB_ID}\"; \
    mkdir -p \"\${tmpdir_job}\"; \
    docker run --rm \
      --pull=always --user=\"${run_user}\" \
      ${volflags} \
      --env TMPDIR=/scratch --env TMP=/scratch --env TEMP=/scratch --env HOME=/tmp \
      -v \"\${tmpdir_job}:/scratch\" \
      ${containerimage} dwi-anat2dwi ${subjspecjson}"

    job_id_03=$(submit_job "dwi-anat2dwi_${subj}${sessionfile}" "${cmd_03}" "${dep_02a}" \
        "${ANAT_CPUS}" "${ANAT_MEM}" "${ANAT_TIME}")
    echo "Submitted dwi-anat2dwi (03) job: ${job_id_03}"
    all_jobs+=("${job_id_03}")
fi

# ============================================================
# QC after preproc stages (02a, 02b, 03)
# Depends on both 02b and 03 finishing successfully.
# ============================================================
dep_preproc=$(make_dep_string "${job_id_02b:-}" "${job_id_03:-}")

cmd_qc_preproc="docker run --rm \
  --pull=always --user=\"${run_user}\" \
  ${volflags} \
  --env HOME=/tmp \
  ${containerimage} dwi-qc ${subjspecjson}"
job_id_qc_preproc=$(submit_job "dwi-qc-preproc_${subj}${sessionfile}" "${cmd_qc_preproc}" "${dep_preproc}" \
    "${QC_CPUS}" "${QC_MEM}" "${QC_TIME}")
echo "Submitted dwi-qc (post-preproc) job: ${job_id_qc_preproc}"
all_jobs+=("${job_id_qc_preproc}")

# ============================================================
# Conditional: skip 04a/04b if --preproc-only
# ============================================================
if [[ "${preproc_only}" -eq 1 ]]; then
    echo "preproc-only mode: skipping tractography stages (04a, 04b)"
else

    # ============================================================
    # STAGE 04a: dwi-tractogram
    # Depends on 02b AND 03 finishing successfully.
    # ============================================================
    dep_04a=$(make_dep_string "${job_id_02b:-}" "${job_id_03:-}")

    if [[ -f "${tracto_file}" ]]; then
        echo "04a output already exists for ${subj} ${session:-}, skipping dwi-tractogram."
        job_id_04a=""
    else
        cmd_04a="tmpdir_job=\"\${SLURM_TMPDIR:-${host_workdir}/tmp}/${subj}${sessionfile}\${SLURM_JOB_ID}\"; \
        mkdir -p \"\${tmpdir_job}\"; \
        docker run --rm \
          --pull=always --user=\"${run_user}\" \
          ${volflags} \
          --env TMPDIR=/scratch --env TMP=/scratch --env TEMP=/scratch --env HOME=/tmp \
          -v \"\${tmpdir_job}:/scratch\" \
          ${containerimage} dwi-tractogram ${subjspecjson}"

        job_id_04a=$(submit_job "dwi-tractogram_${subj}${sessionfile}" "${cmd_04a}" "${dep_04a}" \
            "${TRACTO_CPUS}" "${TRACTO_MEM}" "${TRACTO_TIME}")
        echo "Submitted dwi-tractogram (04a) job: ${job_id_04a}"
        all_jobs+=("${job_id_04a}")
    fi

    # ============================================================
    # STAGE 04b: dwi-conn — depends on 04a
    # ============================================================
    dep_04b=$(make_dep_string "${job_id_04a:-}")

    if [[ -f "${conn_file}" ]]; then
        echo "04b output already exists for ${subj} ${session:-}, skipping dwi-conn."
        job_id_04b=""
    else
        cmd_04b="docker run --rm \
          --pull=always --user=\"${run_user}\" \
          ${volflags} \
          --env HOME=/tmp \
          ${containerimage} dwi-conn ${subjspecjson}"

        job_id_04b=$(submit_job "dwi-conn_${subj}${sessionfile}" "${cmd_04b}" "${dep_04b}" \
            "${CONN_CPUS}" "${CONN_MEM}" "${CONN_TIME}")
        echo "Submitted dwi-conn (04b) job: ${job_id_04b}"
        all_jobs+=("${job_id_04b}")
    fi

    # ============================================================
    # QC after tracto stages (04a, 04b)
    # ============================================================
    dep_tracto=$(make_dep_string "${job_id_04b:-}")

    cmd_qc_tracto="docker run --rm \
      --pull=always --user=\"${run_user}\" \
      ${volflags} \
      --env HOME=/tmp \
      ${containerimage} dwi-qc ${subjspecjson}"
    job_id_qc_tracto=$(submit_job "dwi-qc-tracto_${subj}${sessionfile}" "${cmd_qc_tracto}" "${dep_tracto}" \
        "${QC_CPUS}" "${QC_MEM}" "${QC_TIME}")
    echo "Submitted dwi-qc (post-tracto) job: ${job_id_qc_tracto}"
    all_jobs+=("${job_id_qc_tracto}")

fi  # end preproc_only

# ============================================================
# FINAL SENTINEL (waits for all submitted jobs)
# ============================================================

# Collect all valid job IDs (filter out empty strings from skipped stages)
valid_jobs=()
for job_id in "${all_jobs[@]}"; do
    if [[ "${job_id}" =~ ^[0-9]+$ ]]; then
        valid_jobs+=("${job_id}")
    fi
done

if [[ "${#valid_jobs[@]}" -gt 0 ]]; then
    dep_string=$(IFS=:; printf "%s" "${valid_jobs[*]}")
    dep_arg=(--dependency=afterok:${dep_string} --kill-on-invalid-dep=yes)
    echo "Final sentinel job will depend on jobs: ${dep_string}"
else
    dep_arg=()
    dep_string=""
fi

if [[ -z "${dep_string}" ]]; then
    echo "No jobs were submitted, skipping final sentinel job."
else
    final_job_id=$(sbatch --wait --parsable \
        "${dep_arg[@]}" \
        --job-name="dwi_Hodor_${subj}" \
        --time=00:01:00 -c 1 --mem=10M \
        --wrap "echo 'Pipeline finished for ${subj}'")

    echo "Pipeline completed for ${subj} (final job ${final_job_id})"
fi
