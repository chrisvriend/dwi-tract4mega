#!/bin/bash

#SBATCH --job-name=eddy
#SBATCH --mem=4G
#SBATCH --partition=luna-cpu-short
#SBATCH --cpus-per-task=4
#SBATCH --time=00-8:00:00
#SBATCH --nice=2000
#SBATCH --qos=anw-cpu
#SBATCH --output=%x_%A.log

###############################################################################
# eddyCPU.sh
# Author: C. Vriend - AUMC
# Date: Nov 05 2025
# Description: perform FSL eddy
###############################################################################

set -euo pipefail

# usage instructions
Usage() {
    cat <<EOF

    (C) C.Vriend - 9/7/2025 - dwi-02b-eddy.sh
   
    Usage: ./dwi-02b-eddy.sh -i <bidsdir> -o <outputdir> -w <workdir> -s <subj> [-z <session>] -m <method>
  
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

module load fsl/6.0.7.6

# Initialize variables
bidsdir=""
outputdir=""
workdir=""
subj=""
session=""
method=""

# Parse command line arguments
while getopts ":i:o:w:s:z:m:" opt; do
    case $opt in
        i) bidsdir="$OPTARG" ;;
        o) outputdir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        z) session="$OPTARG" ;;
        m) method="$OPTARG" ;;
        \?) log "$RED" "Invalid option: -$OPTARG"; exit 1 ;;
        :) log "$RED" "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

missing=0
for var in bidsdir outputdir workdir subj method; do
    if [[ -z "${!var}" ]]; then
        log "$RED" "Error: -${var:0:1} ($var) is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    Usage
fi

log "$BLUE" "chosen method for eddy: ${method}"

# Set session path/file
if [[ -z "${session}" ]]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
fi

log "$YELLOW" "----------------------"
log "$YELLOW" "running EDDY on dwi data"
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
cd "${basedir}"

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

# create index.txt file
idx=$(fslnvols "${DWImain}")
printf '1 %.0s' $(seq 1 "$idx") >"${basedir}/index.txt"

# json available with slice-timing?
if jq -e '.SliceTiming' "${DWIjson}" >/dev/null; then
    STavail=1
else
    STavail=0
fi

case "$method" in
    default)
        eddy \
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
            --estimate_move_by_susceptibility --verbose

        run_qc "${DWIout}" \
            -idx index.txt \
            -par "${DWIacqp}" \
            -m "${DWImask}" \
            -b "${DWIbvals}" \
            -f "${topup}_fieldmap.nii.gz"
        ;;
    volcorr|volcorrnosdc)
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
            )
            [[ "$method" == "volcorr" ]] && eddy_args+=(--estimate_move_by_susceptibility)
            eddy "${eddy_args[@]}" --verbose >"${basedir}/eddy.log"

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
        eddy \
            --imain="${DWImain}" \
            --mask="${DWImask}" \
            --acqp="${DWIacqp}" \
            --index="${basedir}/index.txt" \
            --bvecs="${DWIbvecs}" \
            --bvals="${DWIbvals}" \
            --out="${DWIout}" \
            --repol --cnr_maps \
            --slm=linear \
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

cp eddy_*.log "${outputdir}/dwi-preproc/${subj}${sessionpath}logs/${subj}${sessionfile}eddy.log"

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

rsync -av "${subj}${sessionfile}space-dwi*_dwi.*" "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" eddyqc \
    "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi"

# clean-up
if [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" ] &&
   [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" ] &&
   [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz" ]; then
    rm -r "${workdir}/${subj}${sessionpath}"
    rm "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"*meanb0* \
       "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"*dns+degibbs*

    log "$GREEN" "FINISHED preprocessing ${subj}${sessionpath}"
else
    log "$RED" "ERROR! not all output was created successfully"
fi
