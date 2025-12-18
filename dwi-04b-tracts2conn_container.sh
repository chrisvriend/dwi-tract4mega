#!/bin/bash


###############################################################################
# tck2conn.sh
# Author: C. Vriend - AUMC
# Date: Nov 05 2025
# Description: connectome generation from tractography
###############################################################################

set -euo pipefail

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

Usage() {
    cat <<EOF

    (C) C.Vriend - Amsterdam UMC - Nov 6 2025
    performs tractography
    Usage: ./dwi-04b_tracts2connectome_container.sh -i <bidsdir> -o <outputdir> -w <workdir> -s <subj> [-z <session>] -x <nstreamlines> -t <nthreads>

EOF
    exit 1
}

[ _$1 = _ ] && Usage

# Initialize variables
bidsdir=""
outputdir=""
workdir=""
subj=""
session=""
nstreamlines=50M
nthreads=16


# Parse command line arguments
while getopts ":i:o:w:s:z:t:x:" opt; do
    case $opt in
        i) bidsdir="$OPTARG" ;;
        o) outputdir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        z) session="$OPTARG" ;;
         t) nthreads="$OPTARG" ;;
        x) nstreamlines="$OPTARG" ;;
        \?) log "$RED" "Invalid option: -$OPTARG"; exit 1 ;;
        :) log "$RED" "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

# Check required arguments
missing=0
for var in bidsdir outputdir workdir subj nstreamlines; do
    if [[ -z "${!var}" ]]; then
        log "$RED" "Error: $var is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    Usage
fi

# Set session path/file
if [[ -z "${session}" ]]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
fi

##############
# CHECK FILES
##############
files=$(echo "
    ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck
    ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_tracto-${nstreamlines}_desc-sift_weights.txt
    ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif")
for file in ${files}; do
    if [ ! -f ${file} ]; then
        log "$RED" "!!!ERROR!!!"
        log "$RED" "A scan was not found in the workdir"
        log "$RED" "${file}"
        log "$RED" "Cannot continue without this file"
        exit 1
    fi
done

if [ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-odi_noddi.nii.gz ]; then
    log "$YELLOW" "Found noddi output in output directory"
    log "$YELLOW" "...copying to workdir"
    rsync -a ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-odi_noddi.nii.gz \
        ${workdir}/${subj}${sessionpath}dwi/
fi

if compgen -G "${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi*atlas-*" > /dev/null; then
    log "$YELLOW" "Found atlas files in output directory"
    log "$YELLOW" "...copying to workdir"
    rsync -a ${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi*atlas-* \
        ${workdir}/${subj}${sessionpath}anat/
else 
    log "$RED" "!!!ERROR!!!"
    log "$RED" "No atlas files found in output directory"
    log "$RED" "Cannot continue without these files"
    exit 1
fi

log "$BLUE" "-- ------------------ --"
log "$BLUE" "tractogram -2- connectome"
log "$BLUE" "-- ------------------ --"

###############################################################################
# DWI CONNECTOME GENERATION
###############################################################################
mkdir -p "${outputdir}/dwi-tracto/${subj}${sessionpath}conn"
mkdir -p "${workdir}/${subj}${sessionpath}conn"

cd ${workdir}/${subj}${sessionpath}dwi

# Determine whether it is single or multishell
dwishells=$(mrinfo "${subj}${sessionfile}space-dwi_desc-preproc_dwi.mif" -shell_bvalues | \
tr ' ' '\n' | awk '$1 > 0' )
Nshells=$(echo "$dwishells" | wc -w)

if (( Nshells > 1 )); then
    lowshell=$(mrinfo "${subj}${sessionfile}space-dwi_desc-preproc_dwi.mif" -shell_bvalues |  tr ' ' '\n' | awk '$1 > 0' | head -n 1)
    log "$YELLOW" "Multishell dwi detected"
    log "$YELLOW" "Using b=${lowshell} for DTI fitting"
    dwiextract ${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif \
        b0b${lowshell}.mif -shells 0,${lowshell} -force
    mrconvert b0b${lowshell}.mif ${subj}${sessionfile}space-dwi_desc-lowbval_dwi.nii.gz \
        -export_grad_fsl b${lowshell}.bvec b${lowshell}.bval -force
    dtifit -k ${subj}${sessionfile}space-dwi_desc-lowbval_dwi.nii.gz \
        -m ${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz \
        -r b${lowshell}.bvec -b b${lowshell}.bval \
        -o ${subj}${sessionfile}space-dwi_desc-lowbval --sse
    rm b${lowshell}.bv* b0b${lowshell}.mif
else
    mrconvert ${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.mif ${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.nii.gz \
        -export_grad_fsl bvec bval -force
    dtifit -k ${subj}${sessionfile}space-dwi_desc-preproc-biascor_dwi.nii.gz \
        -m ${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz \
        -r bvec -b bval \
        -o ${subj}${sessionfile}space-dwi_desc-lowbval --sse
    rm bval bvec 
fi

#----------------------------------------------------------------------
#                  Connectome generation 
#----------------------------------------------------------------------
# https://mrtrix.readthedocs.io/en/latest/quantitative_structural_connectivity/structural_connectome.html

if [ ! -f ${subj}${sessionfile}space-dwi_desc-lengths_stats.csv ]; then
    # extract lengths
    tckstats -dump ${subj}${sessionfile}space-dwi_desc-lengths_stats.csv ${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck \
        -tck_weights_in ${subj}${sessionfile}space-dwi_tracto-${nstreamlines}_desc-sift_weights.txt -force -nthreads ${nthreads}
fi

# FA / ND
for diff in FA ndi; do
    if [[ ${diff} == FA ]]; then
        inputfile=${subj}${sessionfile}space-dwi_desc-lowbval_FA.nii.gz
    elif [[ ${diff} == ndi ]]; then
        inputfile=${subj}${sessionfile}space-dwi_desc-odi_noddi.nii.gz
    fi

    if [ ! -f ${subj}${sessionfile}space-dwi_desc-${diff}_stats.csv ] &&
        [ -f ${inputfile} ]; then
        tcksample ${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck \
            ${inputfile} \
            ${subj}${sessionfile}space-dwi_desc-${diff}_stats.csv -stat mean -nthreads ${nthreads}
    fi
done

for atlas in BNA 300P7N; do
    if [ ! -f ${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_dseg.nii.gz ]; then
        log "$YELLOW" "!!WARNING!! ${atlas} atlas not available"
        continue
    fi
    log "$BLUE" ""
    log "$BLUE" "atlas = ${atlas}"
    # streamline count
    # normalization: https://community.mrtrix.org/t/normalization-of-connectomes/4363
    if [ ! -f ${workdir}/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-${atlas}_desc-streams_connmatrix.csv ]; then
        log "$BLUE" "...streamlines..."
        tck2connectome ${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck \
            ${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_dseg.nii.gz \
            ${workdir}/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-${atlas}_desc-streams_connmatrix.csv \
            -zero_diagonal \
            -tck_weights_in ${subj}${sessionfile}space-dwi_tracto-${nstreamlines}_desc-sift_weights.txt \
            -nthreads ${nthreads} -force -symmetric \
            -assignment_radial_search 4 
            # -out_assignments ${workdir}/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-${atlas}_trackassign.txt
    fi
    for scalar in lengths FA ndi; do
        if [ ! -f ${workdir}/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-${atlas}_desc-${scalar}_connmatrix.csv ]; then
            if [ -f ${subj}${sessionfile}space-dwi_desc-${scalar}_stats.csv ]; then
                log "$BLUE" "...${scalar}..."
                tck2connectome ${subj}${sessionfile}space-dwi_tracto-${nstreamlines}.tck \
                    ${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_dseg.nii.gz \
                    ${workdir}/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-${atlas}_desc-${scalar}_connmatrix.csv \
                    -scale_file ${subj}${sessionfile}space-dwi_desc-${scalar}_stats.csv \
                    -zero_diagonal \
                    -tck_weights_in ${subj}${sessionfile}space-dwi_tracto-${nstreamlines}_desc-sift_weights.txt \
                    -stat_edge mean \
                    -assignment_radial_search 4 \
                    -nthreads ${nthreads} -force -symmetric
            else
                log "$RED" ""
                log "$RED" "!ERROR! ${scalar} scalar file does not exist"
                log "$RED" ""
                continue
            fi
        fi
    done
done

log "$GREEN" ""
log "$GREEN" "finished tck2connectome. Transfer files"
rsync -av ${workdir}/${subj}${sessionpath}conn/* ${outputdir}/dwi-tracto/${subj}${sessionpath}conn

if compgen -G "${outputdir}/dwi-tracto/${subj}${sessionpath}conn/${subj}${sessionfile}atlas-*_desc-*_connmatrix.csv" >/dev/null ; then
    log "$GREEN" "--------------------------------------------------"
    log "$GREEN" "finished connectome generation subject = ${subj}"
    log "$GREEN" "--------------------------------------------------"
    #clean-up
    rm ${workdir}/${subj}${sessionpath}dwi/*.mif
else
    log "$RED" ""
    log "$RED" "ERROR! connectome generation failed for subject = ${subj}${sessionfile}"
    log "$RED" "inspect the log file"
    exit 1
fi
