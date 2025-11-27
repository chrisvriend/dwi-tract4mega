#!/bin/bash

#SBATCH --job-name=NODDI
#SBATCH --mem=4G
#SBATCH --partition=luna-gpu-short
#SBATCH --cpus-per-task=1
#SBATCH --time=00-00:20:00
#SBATCH --nice=2000
#SBATCH --qos=anw
#SBATCH --output=%x_%A.log

###############################################################################
# noddi.sh
# Author: C. Vriend - AUMC
# Date: Nov 19 2025
# Description: fit NODDI model with CUDIMOT
###############################################################################

set -euo pipefail

# Define color variables
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

usage() {
    cat <<EOF

    (C) C.Vriend - Amsterdam UMC - Nov 2025
    fit NODDI model with CUDIMOT
    Usage: ./dwi-02d-noddi.sh -o <outputdir> -w <workdir> -s <subj> [-z <session>]

EOF
    exit 1
}

module load fsl/6.0.7.6
CUDIMOT=/scratch/anw/cvriend/cudimot/FSLDEV
export CUDIMOT

# Initialize variables
outputdir=""
workdir=""
subj=""
session=""
sessionID=""   


# Parse command line arguments
while getopts ":o:w:s:z:" opt; do
    case $opt in
        o) outputdir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        z) session="$OPTARG" ;;
        \?) log "$RED" "Invalid option: -$OPTARG"; usage ;;
        :) log "$RED" "Option -$OPTARG requires an argument."; usage ;;
    esac
done

missing=0
for var in outputdir workdir subj; do
    if [[ -z "${!var}" ]]; then
        log "$RED" "Error: -${var:0:1} ($var) is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    usage
fi

if [ -z "${session}" ]; then
    sessionpath="/"
    sessionfile="_"
    sessionID=""
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
    sessionID="-${session}"
fi

cd "${workdir}/${subj}${sessionpath}noddi"

log "$BLUE" "Running CUDIMOT NODDI pipeline for ${subj}${sessionID}..."
"${CUDIMOT}/bin/Pipeline_NODDI_Watson.sh" "${subj}${sessionID}"

if [ -f "${subj}${sessionID}.NODDI_Watson/mean_fiso.nii.gz" ]; then
    log "$GREEN" "NODDI output found, renaming and moving files..."

    mv "${subj}${sessionID}.NODDI_Watson/OD.nii.gz" \
        "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-odi_noddi.nii.gz"
    mv "${subj}${sessionID}.NODDI_Watson/mean_fintra.nii.gz" \
        "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-ndi_noddi.nii.gz"
    mv "${subj}${sessionID}.NODDI_Watson/mean_fiso.nii.gz" \
        "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-isovf_noddi.nii.gz"
else
    log "$RED" "ERROR! NODDI failed: output not found."
    exit 1
fi

log "$BLUE" "Copying NODDI outputs and related files to output directory..."
rsync -av "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi"*noddi.nii.gz \
    "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodif-brain_dwi.nii.gz" \
    "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" \
    "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"

if [ $? -eq 0 ]; then
    log "$GREEN" "NODDI processing completed successfully for ${subj}${sessionfile}"
    rm -r "${workdir}/${subj}${sessionpath}"
else
    log "$RED" "ERROR! Failed to copy NODDI outputs to ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"
    exit 1
fi
