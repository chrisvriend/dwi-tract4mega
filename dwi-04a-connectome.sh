#!/bin/bash

#SBATCH --job-name=dwi-fod-tck
#SBATCH --mem=12G
#SBATCH --partition=luna-cpu-short
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=16
#SBATCH --time=00-04:00:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A.log

###############################################################################
# fod-tck.sh
# Author: C. Vriend - AUMC
# Date: Nov 05 2025
# Description: FOD estimation and tractography
###############################################################################

set -euo pipefail

usage() {
    cat <<EOF

    (C) C.Vriend - Amsterdam UMC - May 24 2023 
    performs tractography
    Usage: ./dwi-04a_connectome.sh -i <bidsdir> -o <outputdir> -w <workdir> -s <subj> [-z <session>]

EOF
    exit 1
}

# Define color variables
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[34m'
NC='\033[0m' # No Color

# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

# Initialize variables
bidsdir=""
outputdir=""
workdir=""
subj=""
session=""

# Parse command line arguments
while getopts ":i:o:w:s:z:" opt; do
    case $opt in
        i) bidsdir="$OPTARG" ;;
        o) outputdir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        z) session="$OPTARG" ;;
        \?) log "$RED" "Invalid option: -$OPTARG"; exit 1 ;;
        :) log "$RED" "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

# Check required arguments
missing=0
for var in bidsdir outputdir workdir subj; do
    if [[ -z "${!var}" ]]; then
        log "$RED" "Error: $var is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    usage
fi

nstreamlines=50M
threads=${SLURM_CPUS_PER_TASK}

###############################################################################
# source software
###############################################################################
module load fsl/6.0.6.5
module load FreeSurfer/7.3.2-centos8_x86_64
module load ANTs/2.5.0
module load Anaconda3/2023.03
conda activate /scratch/anw/share/python-env/mrtrix
export PATH=${PATH}:/scratch/anw/share/python-env/mrtrix/MRtrix3Tissue/bin

###############################################################################

# Set session path/file
if [[ -z "${session}" ]]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
fi


if [ ! -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" ]; then 
    log "$RED" "ERROR!! no preprocessed dwi scan found for ${subj} - ${session}"
    exit 1

else 

    rsync -a "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc"* \
    "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" \
        "${workdir}/${subj}${sessionpath}dwi"
fi


if [ ! -f "${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" ] ||
    [ ! -f "${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-gmwm_probseg.nii.gz" ]; then 
    log "$RED" "ERROR!! no hsvs or gmwm anatomical files found for ${subj} - ${session}"
    exit 1
fi

# Check if output already available
if [ -f "${outputdir}/dwi-connectome/${subj}${sessionpath}dwi/${subj}${sessionfile}tissue-WM-norm_fod.nii.gz" ] &&
   [ -f "${outputdir}/dwi-connectome/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}_desc-sift_weights.txt" ] &&
   [ -f "${outputdir}/dwi-connectome/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck" ]; then
    log "$GREEN" "${subj}${sessionfile} already has tractogram and sift weights"
    log "$GREEN" "...skip..."
    echo
    exit 0
fi

mkdir -p "${workdir}/${subj}/${sessionpath}dwi"

for folder in dwi figures logs; do
    mkdir -p "${outputdir}/dwi-connectome/${subj}${sessionpath}/${folder}"
done
   
rsync -a "${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc"* \
    "${workdir}/${subj}${sessionpath}anat" 

# In case there is a previous run
if compgen -G "${outputdir}/dwi-connectome/${subj}${sessionpath}dwi/${subj}${sessionfile}*" > /dev/null; then
    rsync -a "${outputdir}/dwi-connectome/${subj}${sessionpath}dwi/${subj}${sessionfile}"* \
        "${outputdir}/dwi-connectome/${subj}${sessionpath}rpf/"* \
        "${workdir}/${subj}${sessionpath}dwi"
fi

#----------------------------------------------------------------------
# Import to MRtrix 
#----------------------------------------------------------------------

cd "${workdir}/${subj}${sessionpath}dwi"
mrconvert "${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" \
    -fslgrad "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" \
    "${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval" \
    "${subj}${sessionfile}space-dwi_desc-preproc_dwi.mif" -force

#----------------------------------------------------------------------
# Bias Correction 
#----------------------------------------------------------------------
if [ ! -f "${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif" ]; then
    dwibiascorrect ants "${subj}${sessionfile}space-dwi_desc-preproc_dwi.mif" \
        "${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif" -nthreads "${threads}" \
        -bias "${subj}${sessionfile}space-dwi_desc-biasest_dwi.mif" \
        -scratch "${workdir}/${subj}${sessionpath}tempbiascorrect" -force
fi

# Dilate/erode brain mask
for manipulation in dilate erode; do
    maskfilter -npass 2 "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" \
        "${manipulation}" "${subj}${sessionfile}space-dwi_desc-brain-${manipulation}d_mask.mif" \
        -nthreads "${threads}" -info -force
done

#----------------------------------------------------------------------
#       Response Function Estimation & Fiber Orientation Distribution 
#----------------------------------------------------------------------

# Determine whether it is single or multishell
dwishells=$(mrinfo "${subj}${sessionfile}space-dwi_desc-preproc_dwi.mif" -shell_bvalues | \
  tr ' ' '\n' | awk '$1 > 0' )
Nshells=$(echo "$dwishells" | wc -w)
shells=$(echo $dwishells | tr ' ' '\n' | paste -sd, -)


if (( Nshells == 1 )); then
    if [ ! -f "${subj}${sessionfile}space-dwi_tissue-WM_response.txt" ] ||
       [ ! -f "${subj}${sessionfile}space-dwi_tissue-GM_response.txt" ] ||
       [ ! -f "${subj}${sessionfile}space-dwi_tissue-CSF_response.txt" ]; then
    echo
    log "$BLUE" "Estimate response functions - dhollander"
    echo
    dwi2response dhollander "${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif" \
        "${subj}${sessionfile}space-dwi_tissue-WM_response.txt" \
        "${subj}${sessionfile}space-dwi_tissue-GM_response.txt" \
        "${subj}${sessionfile}space-dwi_tissue-CSF_response.txt" \
        -nthreads "${threads}" -scratch "${workdir}/${subj}${sessionpath}tempdwiresponse"
        rm -rf "${workdir}/${subj}${sessionpath}tempdwiresponse"
    fi
    if [ ! -f "${subj}${sessionfile}FOD-wm.mif" ] ||
       [ ! -f "${subj}${sessionfile}FOD-gm.mif" ] ||
       [ ! -f "${subj}${sessionfile}FOD-csf.mif" ]; then
    log "$BLUE" "Spherical Deconvolution - ss3t"
    ss3t_csd_beta1 "${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif" \
        "${subj}${sessionfile}space-dwi_tissue-WM_response.txt" "${subj}${sessionfile}FOD-wm.mif" \
        "${subj}${sessionfile}space-dwi_tissue-GM_response.txt" "${subj}${sessionfile}FOD-gm.mif" \
        "${subj}${sessionfile}space-dwi_tissue-CSF_response.txt" "${subj}${sessionfile}FOD-csf.mif" \
        -mask "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz"

    fi
elif (( Nshells > 1 )); then
    if [ ! -f "${subj}${sessionfile}space-dwi_tissue-WM_response.txt" ] ||
       [ ! -f "${subj}${sessionfile}space-dwi_tissue-GM_response.txt" ] ||
       [ ! -f "${subj}${sessionfile}space-dwi_tissue-CSF_response.txt" ]; then
        log "$BLUE" "estimate response functions - msmt"
        dwi2response msmt_5tt \
            "${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif" \
            "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" \
            "${subj}${sessionfile}space-dwi_tissue-WM_response.txt" \
            "${subj}${sessionfile}space-dwi_tissue-GM_response.txt" \
            "${subj}${sessionfile}space-dwi_tissue-CSF_response.txt" \
            -shell 0,"${shells}" -mask "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" \
            -nthreads "${threads}" -scratch "${workdir}/${subj}${sessionpath}tempdwiresponse"
        rm -rf "${workdir}/${subj}${sessionpath}tempdwiresponse"
    fi

    if [ ! -f "${subj}${sessionfile}FOD-wm.mif" ] ||
       [ ! -f "${subj}${sessionfile}FOD-gm.mif" ] ||
       [ ! -f "${subj}${sessionfile}FOD-csf.mif" ]; then
        log "$BLUE" "Spherical Deconvolution - msmt csd"
        dwi2fod msmt_csd "${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif" \
            "${subj}${sessionfile}space-dwi_tissue-WM_response.txt" "${subj}${sessionfile}FOD-wm.mif" \
            "${subj}${sessionfile}space-dwi_tissue-GM_response.txt" "${subj}${sessionfile}FOD-gm.mif" \
            "${subj}${sessionfile}space-dwi_tissue-CSF_response.txt" "${subj}${sessionfile}FOD-csf.mif" \
            -mask "${subj}${sessionfile}space-dwi_desc-brain-dilated_mask.mif" \
            -shell 0,"${shells}" -nthreads "${threads}"
    fi
fi

log "$BLUE" "Multi-tissue informed log-domain intensity normalisation"
mtnormalise "${subj}${sessionfile}FOD-wm.mif" "${subj}${sessionfile}FOD-wm-norm.mif" \
    "${subj}${sessionfile}FOD-gm.mif" "${subj}${sessionfile}FOD-gm-norm.mif" \
    "${subj}${sessionfile}FOD-csf.mif" "${subj}${sessionfile}FOD-csf-norm.mif" \
    -mask "${subj}${sessionfile}space-dwi_desc-brain-eroded_mask.mif" \
    -nthreads "${threads}" -force

mrconvert "${subj}${sessionfile}FOD-wm.mif" - -coord 3 0 | \
    mrcat "${subj}${sessionfile}FOD-csf.mif" "${subj}${sessionfile}FOD-gm.mif" - \
    "${subj}${sessionfile}space-dwi_tissue-RGB.mif" -axis 3 -force

#----------------------------------------------------------------------
#                       Generate Tractogram 
#----------------------------------------------------------------------

if [ ! -f "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck" ]; then
    log "$BLUE" "start tractography"
    echo
    tckgen "${subj}${sessionfile}FOD-wm-norm.mif" \
        "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck" \
        -seed_gmwmi "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-gmwm_probseg.nii.gz" \
        -act "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" \
        -maxlength 250 \
        -cutoff 0.1 \
        -seeds "${nstreamlines}" \
        -select 0 \
        -nthreads "${threads}" \
        -info
fi

tckedit "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck" \
    "${subj}${sessionfile}space-dwi_tracto-100k.tck" -number 100k -force


#----------------------------------------------------------------------
#  Spherical-deconvolution Informed Filtering of space-dwi_tractos (SIFT) 
#----------------------------------------------------------------------


if [ ! -f "${subj}${sessionfile}space-dwi_desc-sift-${nstreamlines}_stats.csv" ]; then
    log "$BLUE" "start tck2sift"
    tcksift2 "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck" \
        "${subj}${sessionfile}FOD-wm-norm.mif" \
        "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}_desc-sift_weights.txt" \
        -act "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" \
        -out_mu "${subj}${sessionfile}space-dwi_mu.txt" \
        -csv "${subj}${sessionfile}space-dwi_desc-sift-${nstreamlines}_stats.csv" \
        -force -nthreads "${threads}"
fi

tckmap "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck" \
    "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_dwi.nii.gz" \
    -template "${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz" \
    -tck_weights "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}_desc-sift_weights.txt" \
    -force -nthreads "${threads}"

# QC
mrconvert ${subj}${sessionfile}space-dwi_tissue-RGB.mif \
 ${subj}${sessionfile}space-dwi_tissue-RGB.nii.gz

fslmaths "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_dwi.nii.gz" \
    -bin "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_mask.nii.gz"
overlay 1 0 "${subj}${sessionfile}space-dwi_tissue-RGB.nii.gz" \
    -a "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_mask.nii.gz" 0 1 \
    "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_overlay.nii.gz"
slicer "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_overlay.nii.gz" \
    -i 0 1 -a "${outputdir}/dwi-connectome/${subj}${sessionpath}figures/${subj}${sessionfile}siftoverlay3D.png"
rm "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_mask.nii.gz" \
   "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}-sift_overlay.nii.gz"

rsync -a "${subj}${sessionfile}space-dwi_tracto-${nstreamlines}"* *sift* *mu* \
    "${outputdir}/dwi-connectome/${subj}${sessionpath}dwi"
rsync -a *response* "${outputdir}/dwi-connectome/${subj}${sessionpath}rpf"
mrconvert "${subj}${sessionfile}FOD-wm-norm.mif" \
    "${outputdir}/dwi-connectome/${subj}${sessionpath}dwi/${subj}${sessionfile}tissue-WM-norm_fod.nii.gz"

echo "--------------------------------------------------"
log "$GREEN" "finished fod and tractography for subject = ${subj}${sessionfile}"
echo "---------------------------------------------------"

