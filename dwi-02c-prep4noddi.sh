#!/bin/bash

#SBATCH --job-name=prepNODDI
#SBATCH --mem=6G
#SBATCH --partition=luna-cpu-short
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=2
#SBATCH --time=00-00:30:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A.log

###############################################################################
# prep4noddi.sh
# Author: C. Vriend - AUMC
# Date: Nov 19 2025
# Description: prepare files for NODDI fitting with CUDIMOT
###############################################################################

set -euo pipefail

###############################################################################
# source software
###############################################################################
module load fsl/6.0.7.6
module load ANTs/2.5.1
module load Anaconda3/2023.03
synthstrippath=/scratch/anw/share-np/fmridenoiser/synthstrip.1.2.sif
conda activate /scratch/anw/share/python-env/mrtrix

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
    prepare DWI data for NODDI fitting with CUDIMOT
    Usage: ./dwi-02c-prep4noddi.sh -o <outputdir> -w <workdir> -s <subj> [-z <session>]

EOF
    exit 1
}

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
        \?) log "$RED" "Invalid option: -$OPTARG"; exit 1 ;;
        :) log "$RED" "Option -$OPTARG requires an argument."; exit 1 ;;
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

  if [ ! -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" ] &&
    [ ! -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" ] &&
   [ ! -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz" ]; then

    log "${RED}" "eddy output not available in output dir for ${subj} ${session:-nosession}, abort."
    exit 1
    fi

# determine whether it is single or multishell using bval file
bval_file="${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval"
if [[ ! -f "$bval_file" ]]; then
    log "$RED" "ERROR! bval file not found: $bval_file"
    exit 1
fi

read -r line <"$bval_file"
IFS=" " read -ra values <<<"$line"
unique_values=$(printf "%s\n" "${values[@]}" | awk '$1 > 0' | sort -u)
Nshells=$(echo "$unique_values" | wc -w)

if ((Nshells == 1)); then
    log "$YELLOW" "DWI ${subj}${sessionpath} is single shell | skipping NODDI"
    exit 1
elif ((Nshells > 1)); then
    if [[ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-isovf_noddi.nii.gz" \
       && -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-ndi_noddi.nii.gz" \
       && -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-odi_noddi.nii.gz" ]]; then 
        log "$GREEN" "${subj}${sessionfile} already has NODDI output"
        exit 1
    fi

    if [[ -f ${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}/data.nii.gz \
         && -f ${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}/bvecs \
         && -f ${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}/bvals ]]; then
          log "$GREEN" "NODDI preparation already done for ${subj}${sessionID}, skipping."
          exit 0
     fi
     

    log "$GREEN" "${subj}${sessionfile} is multishell. Preparing for NODDI"
    mkdir -p "${workdir}/${subj}${sessionpath}dwi"
    rsync -a "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc"* \
        "${workdir}/${subj}${sessionpath}dwi"
    cd "${workdir}/${subj}${sessionpath}dwi"

    if [[ ! -f "${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz" ]] ||
       [[ ! -f "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" ]]; then
        log "$BLUE" "skullstrip dwi and create mask"
        dwiextract -nthreads "${SLURM_CPUS_PER_TASK}" \
            "${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" - -bzero \
            -fslgrad "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" \
            "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval" | \
            mrmath - mean "${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz" -axis 3 -force
        # skullstrip mean b0 (nodif_brain)
        apptainer run --cleanenv "${synthstrippath}" \
            -i "${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz" \
            -o "${subj}${sessionfile}space-dwi_desc-nodif-brain_dwi.nii.gz" \
            --mask "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz"
    fi

    # transfer nodif images
    rsync -av ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz \
    ${subj}${sessionfile}space-dwi_desc-nodif-brain_dwi.nii.gz \
    ${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz \
    ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/
    
    # transfer and rename files to NODDI compatible filenames
    mkdir -p "${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}"
    cp "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" \
        "${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}/nodif_brain_mask.nii.gz"
    cp "${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" \
        "${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}/data.nii.gz"
    cp "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" \
        "${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}/bvecs"
    cp "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval" \
        "${workdir}/${subj}${sessionpath}noddi/${subj}${sessionID}/bvals"
    log "$GREEN" "Preparation for NODDI fitting completed for ${subj}${sessionID}"
else
    log "$RED" "ERROR! something went wrong with reading the bval file to determine number of shells for ${subj}${sessionpath}"
    exit 1
fi