#!/bin/bash
###############################################################################
# eddyCPU.sh
# Author: C. Vriend - AUMC (modified: respect multi-run index.txt from 02a)
# Date: Nov 05 2025 (modified Aug 2026)
# Description: perform FSL eddy
#
# No changes are required here for the *scenario* itself -- dwi-02a now
# writes out identically-named DWImain/DWIbvecs/DWIbvals/DWIacqp/DWIjson
# files regardless of whether it took the single-run or multi-run
# (reversed phase-encode) branch. The only thing that differs between the
# two scenarios is the *content* of acqparams.tsv (1 line vs N lines) and
# index.txt (all-1s vs a block per run). dwi-02a already writes a correct
# multi-line index.txt for the multi-run case; this script must not
# clobber it. That's the only functional change below.
###############################################################################

set -euo pipefail


# usage instructions
Usage() {
    cat <<EOF
(C) C.Vriend - 9/7/2025 - dwi-02b-eddy.sh

Usage: ./dwi-02b-eddyCPU_container.sh -i <bidsdir> -o <outputdir> -c <scriptdir> -w <workdir> -s <subj> [-z <session>] -m <method> -t <nthreads>
EOF
    exit 1
}

# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

run_qc() {
    log "$YELLOW" "running QC"
    eddy_quad "$@"
}

# Define color variables
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[34m'
NC='\033[0m' # No Color

# Initialize variables
bidsdir=""
outputdir=""
scriptdir=""
workdir=""
subj=""
session=""
method=""
nthreads=""

# Parse command line arguments
while getopts ":i:o:c:w:s:z:m:t:" opt; do
    case $opt in
    i) bidsdir="$OPTARG" ;;
    o) outputdir="$OPTARG" ;;
    c) scriptdir="$OPTARG" ;;
    w) workdir="$OPTARG" ;;
    s) subj="$OPTARG" ;;
    z) session="$OPTARG" ;;
    m) method="$OPTARG" ;;
    t) nthreads="$OPTARG" ;;
    \?) log "$RED" "Invalid option: -$OPTARG"; exit 1 ;;
    :) log "$RED" "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

missing=0
for var in bidsdir outputdir workdir scriptdir subj method nthreads; do
    if [[ -z "${!var}" ]]; then
        log "$RED" "Error: $var is required."
        missing=1
    fi
done

if [[ $missing -eq 1 ]]; then
    Usage
fi

# define session-specific paths and filenames (if session is empty, these will be empty strings)
sessionpath="/${session:+${session}/}"
sessionfile="_${session:+${session}_}"

# Check if eddy already completed
if [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" ] &&
    [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" ] &&
    [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval" ] &&
    [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz" ]; then
    log "$GREEN" "Eddy already completed for ${subj} | ${session}"
    log "$GREEN" "skipping"
    echo
    exit 0
fi

log "$YELLOW" "----------------------"
log "$YELLOW" "running EDDY on dwi data"
log "$YELLOW" "w/ eddy method: ${method}"
log "$YELLOW" "& nthreads: ${nthreads}"
log "$YELLOW" "${subj}"
log "$YELLOW" "${session}"
log "$YELLOW" "----------------------"

# inputs
dwiworkdir="${workdir}/${subj}${sessionpath}dwi"
DWImain="${dwiworkdir}/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz"
DWImask="${dwiworkdir}/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz"
DWIacqp="${dwiworkdir}/${subj}${sessionfile}acq-dwi_desc-acqparams.tsv"
DWIbvecs="${dwiworkdir}/${subj}${sessionfile}dwi.bvec"
DWIbvals="${dwiworkdir}/${subj}${sessionfile}dwi.bval"
DWIjson="${dwiworkdir}/${subj}${sessionfile}dwi.json"
topup="${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-topup"
DWIout="${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc"

basedir="$(dirname "${DWImain}")"
cd ${basedir}

# Check required files
required_files=(
    "${DWImain}"
    "${DWImask}"
    "${DWIacqp}"
    "${DWIbvecs}"
    "${DWIbvals}"
    "${DWIjson}"
)

for f in "${required_files[@]}"; do
    if [[ ! -f "$f" ]]; then
        log "$RED" "Missing required file: $f"
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# create index.txt, UNLESS dwi-02a already wrote one (multi-run / reversed
# phase-encode branch: one block of indices per acqparams line, in the
# same volume order as DWImain/DWIbvecs/DWIbvals). If dwi-02a ran the
# classic single-run branch, no index.txt exists yet here and we generate
# the usual "every volume -> acqparams line 1" version.
# ---------------------------------------------------------------------------
nvols=$(fslnvols "${DWImain}")

if [[ -f "${basedir}/index.txt" ]]; then
    idxcount=$(wc -w <"${basedir}/index.txt")
    if [[ "$idxcount" -ne "$nvols" ]]; then
        log "$RED" "index.txt (${idxcount} entries) does not match number of volumes in ${DWImain} (${nvols})"
        exit 1
    fi
    log "$BLUE" "using existing index.txt from dwi-02a"
else
    printf '1 %.0s' $(seq 1 "$nvols") >"${basedir}/index.txt"
fi

# sanity check: every index must reference a line that actually exists in acqparams
acqp_lines=$(grep -c . "${DWIacqp}")
idx_max=$(tr ' ' '\n' <"${basedir}/index.txt" | sed '/^$/d' | sort -un | tail -1)
if [[ -n "$idx_max" && "$idx_max" -gt "$acqp_lines" ]]; then
    log "$RED" "index.txt references acqparams line ${idx_max} but ${DWIacqp} only has ${acqp_lines} line(s)"
    exit 1
fi

# json available with slice-timing?
if jq -e '.SliceTiming' "${DWIjson}" >/dev/null; then
    STavail=1
else
    STavail=0
fi


# run detect whole / half sphere sampling
log "$BLUE" "infering sphere sampling scheme from bvecs"

sampling=$(python ${scriptdir}/helpers/detect_sphere_sampling.py \
    "${DWIbvecs}")

if [[ "$sampling" == "HALF_SPHERE" ]]; then
    log "$YELLOW" "Detected half-sphere sampling"
    log "$YELLOW" "Setting eddy method to slmlinear"
    method="slmlinear"
elif [[ "$sampling" == "WHOLE_SPHERE" ]]; then
    log "$YELLOW" "Detected whole-sphere sampling"
    log "$YELLOW" "leaving eddy method as ${method}"
else
    log "$RED" "Could not detect sampling scheme from bvecs"
    log "$RED" "Please check your bvecs file: ${DWIbvecs}"
    exit 1
fi


case "$method" in
default)
    eddy diffusion \
        --imain="${DWImain}" \
        --mask="${DWImask}" \
        --acqp="${DWIacqp}" \
        --index=index.txt \
        --bvecs="${DWIbvecs}" \
        --bvals="${DWIbvals}" \
        --out="${DWIout}" \
        --topup="${topup}" \
        --repol --cnr_maps \
        --verbose \
        --nthr=${nthreads} >"${basedir}/eddy.log"

    run_qc "${DWIout}" \
        -idx index.txt \
        -par "${DWIacqp}" \
        -m "${DWImask}" \
        -b "${DWIbvals}" \
        -f "${topup}_fieldmap.nii.gz"
    ;;
slmlinear)
    eddy diffusion \
        --imain="${DWImain}" \
        --mask="${DWImask}" \
        --acqp="${DWIacqp}" \
        --index=index.txt \
        --bvecs="${DWIbvecs}" \
        --bvals="${DWIbvals}" \
        --out="${DWIout}" \
        --topup="${topup}" \
        --repol --cnr_maps \
        --slm=linear \
        --verbose \
        --nthr=${nthreads} >"${basedir}/eddy.log"

    run_qc "${DWIout}" \
        -idx index.txt \
        -par "${DWIacqp}" \
        -m "${DWImask}" \
        -b "${DWIbvals}" \
        -f "${topup}_fieldmap.nii.gz"
    ;;
slmquadratic)
    eddy diffusion \
        --imain="${DWImain}" \
        --mask="${DWImask}" \
        --acqp="${DWIacqp}" \
        --index=index.txt \
        --bvecs="${DWIbvecs}" \
        --bvals="${DWIbvals}" \
        --out="${DWIout}" \
        --topup="${topup}" \
        --repol --cnr_maps \
        --slm=quadratic \
        --verbose \
        --nthr=${nthreads} >"${basedir}/eddy.log"

    run_qc "${DWIout}" \
        -idx index.txt \
        -par "${DWIacqp}" \
        -m "${DWImask}" \
        -b "${DWIbvals}" \
        -f "${topup}_fieldmap.nii.gz"
    ;;
volcorr | volcorrnosdc)
    if ((STavail == 1)); then
        eddy_args=(
            --imain="${DWImain}"
            --mask="${DWImask}"
            --acqp="${DWIacqp}"
            --index=index.txt
            --json="${DWIjson}"
            --bvecs="${DWIbvecs}"
            --bvals="${DWIbvals}"
            --out="${DWIout}"
            --topup="${topup}"
            --repol --cnr_maps
            --slm=linear
            --mbs_niter=10 --mbs_lambda=10 --mbs_ksp=10
            --niter=6 --fwhm=15,10,4,2,0,0
            --mporder=8 --s2v_niter=8 --json="${DWIjson}"
            --s2v_lambda=1 --s2v_interp=trilinear
            --nthr=${nthreads}
        )
        [[ "$method" == "volcorr" ]] && eddy_args+=(--estimate_move_by_susceptibility)

        eddy diffusion "${eddy_args[@]}" --verbose >"${basedir}/eddy.log"

        run_qc "${DWIout}" \
            -idx index.txt \
            -par "${DWIacqp}" \
            -m "${DWImask}" \
            -b "${DWIbvals}" \
            -f "${topup}_fieldmap.nii.gz" \
            -g "${DWIout}.eddy_rotated_bvecs" \
            -j "${DWIjson}" \
            -v
    else
        log "$RED" "Slice to volume correction not possible without SliceTime information in json"
        exit 1
    fi
    ;;
nofmap)
    eddy diffusion \
        --imain="${DWImain}" \
        --mask="${DWImask}" \
        --acqp="${DWIacqp}" \
        --index="${basedir}/index.txt" \
        --bvecs="${DWIbvecs}" \
        --bvals="${DWIbvals}" \
        --out="${DWIout}" \
        --repol --cnr_maps \
        --slm=linear \
        --nthr=${nthreads} \
        --verbose >"${basedir}/eddy.log"

    run_qc "${DWIout}" \
        -idx "${basedir}/index.txt" \
        -par "${DWIacqp}" \
        -m "${DWImask}" \
        -b "${DWIbvals}"
    ;;
*)
    log "$RED" "Proper method for eddy not set"
    exit 1
    ;;
esac

cp ${basedir}/eddy.log "${outputdir}/dwi-preproc/${subj}/log/${subj}${sessionfile}eddy_$(date +"%Y-%m-%d_%H-%M").log"

# rename output
cd "${workdir}/${subj}${sessionpath}dwi"
cp "${subj}${sessionfile}space-dwi_desc-preproc.eddy_rotated_bvecs" \
    "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec"
cp "${subj}${sessionfile}space-dwi_desc-preproc.nii.gz" \
    "${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz"
cp "${subj}${sessionfile}space-dwi_desc-preproc.eddy_cnr_maps.nii.gz" \
    "${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz"
cp "${DWIbvals}" \
    "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval"

mv *.qc eddyqc

# needed for QC (first two might no longer be necessary when parameters are available in json file)
mv ${DWIout}.eddy_movement_rms ${DWIout}.eddy_outlier_report \
    "${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz" ./eddyqc

rsync -rltpDv eddyqc/* "${outputdir}/dwi-preproc/${subj}${sessionpath}qc"
rm -r eddyqc

rsync -rltpDv ${subj}${sessionfile}*acqparams.tsv "${subj}${sessionfile}space-dwi_desc-preproc.eddy.json" \
    ${subj}${sessionfile}space-dwi*_dwi.* "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" \
    "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi"

# also keep the run-key (multi-run branch) if present, for provenance
if [[ -f "${dwiworkdir}/${subj}${sessionfile}acq-dwi_desc-runkey.tsv" ]]; then
    rsync -rltpDv "${dwiworkdir}/${subj}${sessionfile}acq-dwi_desc-runkey.tsv" \
        "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi"
fi

# clean-up
if [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" ] &&
    [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" ] &&
    [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval" ] &&
    [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}qc/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz" ]; then
    rm ${workdir}/${subj}${sessionpath}dwi/*desc-preproc*
    rm "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"*meanb0* \
        "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"*dns+degibbs*
    log "$GREEN" "FINISHED preprocessing ${subj}${sessionpath}"
else
    log "$RED" "ERROR! not all output was created successfully"
fi
