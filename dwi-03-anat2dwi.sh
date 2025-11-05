#!/bin/bash
#SBATCH --job-name=anat2dwi
#SBATCH --mem=8G
#SBATCH --partition=luna-cpu-short
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=2
#SBATCH --time=00-03:00:00
#SBATCH --nice=2000
#SBATCH --output=anat2dwi_%A.log

###############################################################################
# anat2dwi.sh
# Author: C. Vriend - AUMC
# Date: Nov 05 2025
# Description: Preprocessing pipeline for anatomical to DWI registration and atlas mapping.
###############################################################################

set -euo pipefail

# Load modules
module load fsl/6.0.6.5
module load FreeSurfer/7.3.2-centos8_x86_64
module load ANTs/2.4.1
module load art
module load Anaconda3/2023.03
conda activate /scratch/anw/share/python-env/mrtrix

MRTRIXapp="/opt/aumc-containers/apptainer/mrtrix3/MRtrix3-3.0.4.sif"
export APPTAINER_BINDPATH="/scratch,/data/anw/anw-work"

# Color variables
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[34m'
NC='\033[0m' # No Color

# Initialize variables
bidsdir=""
outputdir=""
workdir=""
freesurferdir=""
subj=""
session=""
scriptdir=""

# Parse command line arguments
usage() {
    echo "Usage: $0 -i <bidsdir> -o <outputdir> -w <workdir> -f <freesurferdir> -s <subj> -c <scriptdir> [-z <session>]"
    exit 1
}

while getopts ":i:o:w:f:z:s:c:" opt; do
    case $opt in
        i) bidsdir="$OPTARG" ;;
        o) outputdir="$OPTARG" ;;
        f) freesurferdir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        z) session="$OPTARG" ;;
        c) scriptdir="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check required arguments
missing=0
for var in bidsdir outputdir workdir subj scriptdir freesurferdir; do
    if [[ -z "${!var}" ]]; then
        echo "Error: $var is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    usage
fi

# Set session path/file
if [[ -z "${session}" ]]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
fi

synthstrippath="/scratch/anw/share-np/fmridenoiser/synthstrip.1.2.sif"
atlasdir="/data/anw/anw-work/NP/doorgeefluik/atlas4FreeSurfer"
threads="${SLURM_CPUS_PER_TASK:-8}"

# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}
modify_5tt_hsvs() {
    local tmpdir="$1"
    local freesurferdir="$2"
    local subj="$3"

    cd $tmpdir
    mrmath Left-Inf-Lat-Vent.mif 3rd-Ventricle.mif 4th-Ventricle.mif CSF.mif \
        Right-Inf-Lat-Vent.mif 5th-Ventricle.mif Left_LatVent_ChorPlex.mif Right_LatVent_ChorPlex.mif sum - | \
        mrcalc - 1.0 -min tissue3_init.mif -force

    mrcalc tissue3_init.mif tissue3_init.mif tissue4.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue3.mif -force
    mrmath tissue3.mif tissue4.mif sum tissuesum_34.mif -force
    mrcalc tissue1_init.mif tissue1_init.mif tissuesum_34.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue1.mif -force
    mrmath tissue1.mif tissue3.mif tissue4.mif sum tissuesum_134.mif -force
    mrcalc tissue2_init.mif tissue2_init.mif tissuesum_134.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue2.mif -force
    mrmath tissue1.mif tissue2.mif tissue3.mif tissue4.mif sum tissuesum_1234.mif -force
    mrcalc tissue0_init.mif tissue0_init.mif tissuesum_1234.mif -add 1.0 -sub 0.0 -max -sub 0.0 -max tissue0.mif -force
    mrmath tissue0.mif tissue1.mif tissue2.mif tissue3.mif tissue4.mif sum tissuesum_01234.mif -force

    mrcalc aparc.mif 6 -eq aparc.mif 7 -eq -add aparc.mif 8 -eq -add aparc.mif 45 -eq -add aparc.mif 46 -eq -add aparc.mif 47 -eq -add Cerebellum_volume.mif -force

    mrcalc T1.nii Cerebellum_volume.mif -mult T1_cerebellum_precrop.mif -force
    mrgrid T1_cerebellum_precrop.mif crop -mask Cerebellum_volume.mif T1_cerebellum.nii -force

    mrcalc Cerebellum_volume.mif tissuesum_01234.mif -add 0.5 -gt 1.0 tissuesum_01234.mif -sub 0.0 -if Cerebellar_multiplier.mif -force
    mrconvert tissue0.mif tissue0_fast.mif -force
    mrcalc tissue1.mif Cerebellar_multiplier.mif FAST_1.mif -mult -add tissue1_fast.mif -force
    mrcalc tissue2.mif Cerebellar_multiplier.mif FAST_2.mif -mult -add tissue2_fast.mif -force
    mrcalc tissue3.mif Cerebellar_multiplier.mif FAST_0.mif -mult -add tissue3_fast.mif -force
    mrconvert tissue4.mif tissue4_fast.mif -force
    mrmath tissue0_fast.mif tissue1_fast.mif tissue2_fast.mif tissue3_fast.mif tissue4_fast.mif sum tissuesum_01234_fast.mif -force

    mrcalc 1.0 tissuesum_01234_fast.mif -sub tissuesum_01234_fast.mif 0.0 -gt \
        "${freesurferdir}/${subj}/mri/brainmask.mgz" \
        -add 1.0 -min -mult 0.0 -max csf_fill.mif -force

    mrcalc tissue3_fast.mif csf_fill.mif -add tissue3_fast_filled.mif -force
    mrcat tissue0_fast.mif tissue1_fast.mif tissue2_fast.mif tissue3_fast_filled.mif tissue4_fast.mif - -axis 3 | \
        5ttedit - 5TT.mif -none brain_stem_crop.mif

    mv 5TT.mif result.mif
    mrconvert result.mif "${subj}_5TThsvs.nii.gz" -force
    5ttcheck result.mif
}


# Check if output already exists
if compgen -G "${outputdir}/dwi-preproc/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_atlas-*_dseg.nii.gz" > /dev/null; then
    log "$GREEN" "${subj}${sessionfile} already has atlases in dwi-space"
    log "$GREEN" "...skip..."
    exit 0

else 
    echo "-----------------------------------"
    log "$BLUE" "Processing subject: ${subj}${sessionfile}"
    echo "-----------------------------------"
fi

export SUBJECTS_DIR="${workdir}/${subj}/freesurfer"

# --- FastSurfer/FreeSurfer block ---
if [[ ! -d "${freesurferdir}/${subj}" || ! -f "${freesurferdir}/${subj}/surf/lh.pial" ]]; then
    log "$BLUE" "No pre-run FreeSurfer output available"
    log "$BLUE" "Initializing FastSurfer"

    mkdir -p "${workdir}/${subj}/anat"
    mkdir -p "${workdir}/${subj}${sessionpath}xfms/"
    mkdir -p "${workdir}/${subj}/freesurfer"


    #----------------------------------------------------------------------
    #                           Register T1w to dwi space 
    #----------------------------------------------------------------------
    if [[ ! -f "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-FS_T1w.nii.gz" ]]; then
        # Convert T1w to FreeSurfer compatible resolution
        mri_convert --conform "${bidsdir}/${subj}${sessionpath}anat/${subj}${sessionfile}T1w.nii.gz" \
            "${workdir}/${subj}/anat/${subj}${sessionfile}res-FS_T1w.nii.gz"
        #reorient to std?

        # Brainstrip T1w
        apptainer run --cleanenv "${synthstrippath}" \
            -i "${workdir}/${subj}/anat/${subj}${sessionfile}res-FS_T1w.nii.gz" \
            -o "${workdir}/${subj}/anat/${subj}${sessionfile}res-FS_desc-brain_T1w.nii.gz" \
            --mask "${workdir}/${subj}/anat/${subj}${sessionfile}res_FS_desc-brain_mask.nii.gz"

        log "$BLUE" "Register T1w to dwi space"
        flirt -in "${workdir}/${subj}/anat/${subj}${sessionfile}res-FS_desc-brain_T1w.nii.gz" \
            -ref "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodifbrain_epi.nii.gz" \
            -dof 6 -cost normmi -omat "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}T1w-2-dwi.mat"

        transformconvert "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}T1w-2-dwi.mat" \
            "${workdir}/${subj}/anat/${subj}${sessionfile}res-FS_desc-brain_T1w.nii.gz" \
            "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodifbrain_epi.nii.gz" \
            flirt_import \
            "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}desc-mrtrix_T1w-2-dwi.txt"

        mrtransform "${workdir}/${subj}/anat/${subj}${sessionfile}res-FS_T1w.nii.gz" \
            -linear "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}desc-mrtrix_T1w-2-dwi.txt" \
            "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-FS_T1w.nii.gz"
    fi
    #----------------------------------------------------------------------
    #                           FastSurfer 
    #----------------------------------------------------------------------

    log "$BLUE" "Start FastSurfer"
    apptainer exec -B "${workdir}:/data" \
        -B "${SUBJECTS_DIR}:/output" \
        -B "/opt/aumc-apps-eb/software/FreeSurfer:/fs_license" \
        /scratch/anw/cvriend/fastsurfer-2.4.2.sif \
        /fastsurfer/run_fastsurfer.sh \
        --sd /output --sid "${subj}" --t1 "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-FS_T1w.nii.gz" \
        --3T --threads ${threads} \
        --fs_license /fs_license/license.txt

    rsync -av "${workdir}/${subj}/freesurfer/${subj}" "${freesurferdir}"
fi

# --- 5TT generation and GM/WM boundary ---
if [[ -d "${freesurferdir}/${subj}" && -f "${freesurferdir}/${subj}/scripts/deep-seg.log" ]]; then
    if [[ ! -f "${workdir}/${subj}/freesurfer/${subj}/scripts/deep-seg.log" ]]; then
        mkdir -p "${workdir}/${subj}/freesurfer/"
        rsync -av "${freesurferdir}/${subj}" "${workdir}/${subj}/freesurfer/"
    fi

    if [[ ! -f "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" ]]; then
        log "$BLUE" "5ttgen"
        5ttgen hsvs "${workdir}/${subj}/freesurfer/${subj}" \
            "${subj}${sessionfile}5TThsvs.nii.gz" \
            -hippocampi aseg -thalami aseg -white_stem -nthreads "${threads}" \
            -nocrop -nocleanup -scratch "${workdir}/${subj}/temp_5ttgen" -force
        rm "${subj}${sessionfile}5TThsvs.nii.gz"
        modify_5tt_hsvs "${workdir}/${subj}/temp_5ttgen" "${workdir}/${subj}/freesurfer/" "${subj}"
        mv "${workdir}/${subj}/temp_5ttgen/${subj}${sessionfile}5TThsvs.nii.gz" \
            "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz"
        rm -r "${workdir}/${subj}/temp_5ttgen"
    fi

    log "$BLUE" "5tt GM/WM boundary estimation"
    5tt2gmwmi "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" \
        "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-gmwm_probseg.nii.gz" \
        -nthreads "${threads}" -info -force

    for label in 5tt-hsvs gmwm; do
        mri_convert "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-${label}_probseg.nii.gz" \
        --out_orientation RAS \
            "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-${label}_probseg.nii.gz"
        echo "{
        \"Resolution\": \"based on T1w used as input for FastSurfer\",
        \"Orientation\": \"RAS\",
        \"Space\":\"dwi\"}" > "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-${label}_probseg.json"
    done

    rsync -a ${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc*.* "${outputdir}/dwi-preproc/${subj}/anat"
fi
  
#----------------------------------------------------------------------
#                           Existing FreeSurfer run
#----------------------------------------------------------------------

if [[ -d "${freesurferdir}/${subj}" && ! -f "${freesurferdir}/${subj}/scripts/deep-seg.log" ]]; then
    log "$YELLOW" "Relying on existing FreeSurfer run"
    mkdir -p "${workdir}/${subj}/freesurfer/"
    rsync -av "${freesurferdir}/${subj}" "${workdir}/${subj}/freesurfer/"
    cd "${workdir}/${subj}${sessionpath}"

    # Check if required output files exist
    if [[ ! -f "${outputdir}/dwi-preproc/${subj}/anat/${subj}_res-FS_desc-5tt-hsvs_probseg.nii.gz" ]] ||
       [[ ! -f "${outputdir}/dwi-preproc/${subj}/anat/${subj}_res-FS_desc-gmwm_probseg.nii.gz" ]] ||
       [[ ! -f "${outputdir}/dwi-preproc/${subj}/anat/${subj}_res-FS_desc-wm_probseg.nii.gz" ]]; then

        # Copy T1 and brain images from FreeSurfer directory if needed
        if [[ ! -f "${workdir}/${subj}/anat/${subj}_res-FS_desc-preproc_T1w.nii.gz" ]] ||
           [[ ! -f "${workdir}/${subj}/anat/${subj}_res-FS_desc-brain_T1w.nii.gz" ]]; then

            rsync -azv --ignore-existing \
                "${workdir}/${subj}/freesurfer/${subj}/mri/T1.mgz" \
                "${workdir}/${subj}/freesurfer/${subj}/mri/brain.mgz" \
                "${workdir}/${subj}/anat"

            cd "${workdir}/${subj}/anat"
            # Convert to nii.gz
            mri_convert --in_type mgz --out_type nii \
                --out_orientation RAS brain.mgz "${subj}_res-FS_desc-brain_T1w.nii.gz"
            mri_convert --in_type mgz --out_type nii \
                --out_orientation RAS T1.mgz "${subj}_res-FS_desc-preproc_T1w.nii.gz"
            # Binarize brain mask
            fslmaths "${subj}_res-FS_desc-brain_T1w.nii.gz" -bin "${subj}_res-FS_desc-brain_mask.nii.gz"
            rm T1.mgz brain.mgz

            # QC overlay
            mkdir -p "${workdir}/${subj}/figures"
            slicer "${subj}_res-FS_desc-brain_T1w.nii.gz" "${subj}_res-FS_desc-brain_T1w.nii.gz" \
                -a "${workdir}/${subj}/figures/${subj}_res-FS_BETQC.png"
        fi

        # 5TT estimation
        if [[ ! -f "${workdir}/${subj}/anat/${subj}_res-FS_desc-5tt-hsvs_probseg.nii.gz" ]]; then
            log "$YELLOW" "Prepare 5TT estimation"
            5ttgen hsvs "${workdir}/${subj}/freesurfer/${subj}" \
                "${subj}_5TThsvs.nii.gz" \
                -hippocampi aseg -thalami aseg -white_stem -nthreads "${threads}" \
                -nocrop -nocleanup -scratch "${workdir}/${subj}/temp_5ttgen" -force
            rm "${subj}_5TThsvs.nii.gz"
            modify_5tt_hsvs "${workdir}/${subj}/temp_5ttgen" "${workdir}/${subj}/freesurfer/" "${subj}"
            mv "${workdir}/${subj}/temp_5ttgen/${subj}_5TThsvs.nii.gz" \
                "${workdir}/${subj}/anat/${subj}_res-FS_desc-5tt-hsvs_probseg.nii.gz"
            rm -r "${workdir}/${subj}/temp_5ttgen"
        fi

        if [[ ! -f "${workdir}/${subj}/anat/${subj}_res-FS_desc-gmwm_probseg.nii.gz" ]]; then
            5tt2gmwmi "${workdir}/${subj}/anat/${subj}_res-FS_desc-5tt-hsvs_probseg.nii.gz" \
             "${workdir}/${subj}/anat/${subj}_res-FS_desc-gmwm_probseg.nii.gz" \
                -nthreads "${threads}" -info -force
        fi

        # Reorient to FSL RAS
        for label in 5tt-hsvs gmwm; do
            mri_convert "${workdir}/${subj}/anat/${subj}_res-FS_desc-${label}_probseg.nii.gz" \
             --out_orientation RAS \
             "${workdir}/${subj}/anat/${subj}_res-FS_desc-${label}_probseg.nii.gz"
                echo "{
        \"Resolution\": \"based on FreeSurfer seg\",
        \"Orientation\": \"RAS\",
        \"Space\":\"FreeSurfer\"
        }" > "${workdir}/${subj}/anat/${subj}_res-FS_desc-${label}_probseg.json"
        done     

        # Used in BBR to speed up and prevent NaN errors
        fslroi "${workdir}/${subj}/anat/${subj}_res-FS_desc-5tt-hsvs_probseg.nii.gz" \
            "${workdir}/${subj}/anat/${subj}_res-FS_desc-wm_probseg.nii.gz" 2 1

        # Transfer from work directory to output directory
        rsync -a ${workdir}/${subj}/anat/*.* "${outputdir}/dwi-preproc/${subj}/anat"
    fi

    ###########################
    # T1 to DWI registration
    ###########################

    if [[ -d "${outputdir}/dwi-preproc/${subj}${sessionpath}xfms" ]]; then
        rsync -a "${outputdir}/dwi-preproc/${subj}${sessionpath}xfms" "${workdir}/${subj}${sessionpath}"
    fi
    mkdir -p "${workdir}/${subj}${sessionpath}xfms"
    cd "${workdir}/${subj}${sessionpath}dwi"

    log "$BLUE" "Register T1w to dwi space"

    flirt -in "${workdir}/${subj}/anat/${subj}_res-FS_desc-brain_T1w.nii.gz" \
        -ref "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodifbrain_epi.nii.gz" \
        -dof 6 -cost normmi -omat "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}T1w-2-dwi.mat"

    transformconvert "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}T1w-2-dwi.mat" \
        "${workdir}/${subj}/anat/${subj}_res-FS_desc-brain_T1w.nii.gz" \
        "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodifbrain_epi.nii.gz" \
        flirt_import \
        "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}desc-mrtrix_T1w-2-dwi.txt" -force

    # Apply linear transformation to T1w image:
    if [[ ! -f "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" ]] ||
       [[ ! -f "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-gmwm_probseg.nii.gz" ]]; then

        for label in 5tt-hsvs gmwm wm; do
            log "$BLUE" "Register ${label} to DWI-space"
            mrtransform "${workdir}/${subj}/anat/${subj}_res-FS_desc-${label}_probseg.nii.gz" \
                -linear "${workdir}/${subj}${sessionpath}xfms/${subj}${sessionfile}desc-mrtrix_T1w-2-dwi.txt" \
                "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-${label}_probseg.nii.gz"

            echo "{
            \"Resolution\": \"based on T1w used as input for FastSurfer\",
            \"Orientation\": \"RAS\",
            \"Space\":\"dwi\"
            }" > "${workdir}/${subj}/anat/${subj}${sessionfile}space-dwi_res-high_desc-${label}_probseg.json"
        done
    fi

    # Transfer to output directory
    mkdir -p "${outputdir}/dwi-preproc/${subj}/anat"
    mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"
    rsync -a "${workdir}/${subj}/anat/"* "${outputdir}/dwi-preproc/${subj}/anat"
    rsync -a "${workdir}/${subj}${sessionpath}xfms" "${outputdir}/dwi-preproc/${subj}${sessionpath}"

    rsync -a "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-gmwm_probseg.nii.gz" \
        "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz" \
        "${outputdir}/dwi-preproc/${subj}${sessionpath}anat/"

    ##################################
    # FREESURFER to DWI registration
    ##################################
    if [[ ! -f "${SUBJECTS_DIR}/${subj}/dwi/${subj}${sessionfile}register.dat" ]]; then
        mkdir -p "${SUBJECTS_DIR}/${subj}/dwi/"
        bbregister --s "${subj}" \
            --mov "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz" \
            --init-best --reg "${SUBJECTS_DIR}/${subj}/dwi/${subj}${sessionfile}register.dat" --dti
    fi
fi

# --- Atlas warping and registration to FreeSurfer space---
if [[ ! -d "${freesurferdir}/${subj}/mri" ]]; then
    log "$RED" "FreeSurfer output not available"
    log "$RED" "Processing stopped for ${subj}"
    sleep 1
    exit 1
fi

echo
log "$BLUE" "Warp atlases to FreeSurfer output"

# BRAINNETOME ATLAS
log "$BLUE" "---BRAINNETOME"
if [[ ! -f "${SUBJECTS_DIR}/${subj}/label/lh.BN_Atlas.annot" ]] ||
   [[ ! -f "${SUBJECTS_DIR}/${subj}/label/rh.BN_Atlas.annot" ]]; then
    log "$BLUE" "Warping cortical BNA to individual FreeSurfer space"
    for hemi in lh rh; do
        mris_ca_label -seed 1234 -l "${SUBJECTS_DIR}/${subj}/label/${hemi}.cortex.label" \
            "${subj}" "${hemi}" \
            "${SUBJECTS_DIR}/${subj}/surf/${hemi}.sphere.reg" \
            "${atlasdir}/BNA/${hemi}.BN_Atlas.gcs" "${SUBJECTS_DIR}/${subj}/label/${hemi}.BN_Atlas.annot"

        mris_anatomical_stats -mgz -cortex "${SUBJECTS_DIR}/${subj}/label/${hemi}.cortex.label" \
            -f "${SUBJECTS_DIR}/${subj}/stats/${hemi}.BN_Atlas.stats" -b -a "${SUBJECTS_DIR}/${subj}/label/${hemi}.BN_Atlas.annot" \
            -c "${atlasdir}/BNA/BNA_labels_orig_wsubcortex.txt" "${subj}" "${hemi}" white
    done
fi

# if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/BN_Atlas_subcortex.mgz" ]]; then
#     mri_ca_label -threads 2 "${SUBJECTS_DIR}/${subj}/mri/brain.mgz" \
#         "${SUBJECTS_DIR}/${subj}/mri/transforms/talairach.m3z" "${atlasdir}/BNA/BN_Atlas_subcortex.gca" \
#         "${SUBJECTS_DIR}/${subj}/mri/BN_Atlas_subcortex.mgz"
#     mri_segstats --seg "${SUBJECTS_DIR}/${subj}/mri/BN_Atlas_subcortex.mgz" \
#         --ctab "${atlasdir}/BNA/BNA_labels_orig_wsubcortex.txt" --excludeid 0 \
#         --sum "${SUBJECTS_DIR}/${subj}/stats/BN_Atlas_subcortex.stats"
# fi

if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.nii.gz" ]]; then
    mri_aparc2aseg --threads ${threads} --s "${subj}" --annot BN_Atlas --o "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.mgz"
    mrconvert "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.mgz" "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.nii.gz" -force
fi

# Schaefer Atlas
log "$BLUE" "---Schaefer"
rsync -av --ignore-existing "${FREESURFER_HOME}/subjects/fsaverage" "${SUBJECTS_DIR}"

for parcel in 300P7N ; do
    log "$BLUE" "Parcellation = ${parcel}"
    case "${parcel}" in
        300P7N) ID="300Parcels_7Networks" ;;
        300P17N) ID="300Parcels_17Networks" ;;
        200P7N) ID="200Parcels_7Networks" ;;
        100P7N) ID="100Parcels_7Networks" ;;
        400P7N) ID="400Parcels_7Networks" ;;
        400P17N) ID="400Parcels_17Networks" ;;
        *) log "$RED" "Error: atlas not found"; exit 1 ;;
    esac

    if [[ ! -f "${SUBJECTS_DIR}/${subj}/label/lh.${parcel}.annot" ]] ||
       [[ ! -f "${SUBJECTS_DIR}/${subj}/label/rh.${parcel}.annot" ]]; then
        for hemi in lh rh; do
            mri_surf2surf --srcsubject fsaverage --trgsubject "${subj}" --hemi "${hemi}" \
                --sval-annot "${atlasdir}/Schaefer/fsaverage/label/${hemi}.Schaefer2018_${ID}_order.annot" \
                --tval "${SUBJECTS_DIR}/${subj}/label/${hemi}.${parcel}.annot"
        done
    fi

    if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/${parcel}+aseg.mgz" ]]; then
        mri_aparc2aseg --s "${subj}" --o "${SUBJECTS_DIR}/${subj}/mri/${parcel}+aseg.mgz" --annot "${parcel}"
    fi
done

# --- Atlas to DWI space ---
for atlas in BNA 300P7N; do
    if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.mgz" ]]; then
        log "$YELLOW" "WARNING! atlas: ${atlas} - not available in FreeSurfer directory of ${subj}"
        continue
    fi

    # Declare atlas Ids
    if [[ "${atlas}" =~ ^(100P7N|200P7N|300P7N|300P17N|400P7N|400P17N)$ ]]; then
        ID="Schaefer_${atlas}"
    else
        declare -A map=(
            [aparc500]="aparc500_labels"
            [BNA]="BNA_labels"
            ["BNA+cerebellum"]="BNA+CER_labels"
        )
        ID="${map[$atlas]}"
        [[ -z "${ID}" ]] && { log "$RED" "Atlas not found!"; exit 1; }
    fi

    if [[ ! -f "${SUBJECTS_DIR}/${subj}/dwi/${subj}${sessionfile}register.dat" && -f "${SUBJECTS_DIR}/${subj}/scripts/deep-seg.log" ]]; then
        mri_convert --in_type mgz --out_type nii \
            --out_orientation RAS "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.mgz" \
            "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_temp.nii.gz"
    else
        mri_vol2vol --mov "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz" \
            --targ "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.mgz" \
            --o "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_temp.nii.gz" \
            --reg "${SUBJECTS_DIR}/${subj}/dwi/${subj}${sessionfile}register.dat" --inv --no-save-reg --interp nearest \
            --no-resample

        fslreorient2std "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_temp.nii.gz" \
            "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_temp.nii.gz"
    fi

    if [[ "${ID}" != *"Schaefer"* ]]; then
        atlaspath="${atlasdir}/${atlas}"
    else
        atlaspath="${atlasdir}/Schaefer"
    fi

    labelconvert "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_temp.nii.gz" \
        "${atlaspath}/${ID}_orig.txt" \
        "${atlaspath}/${ID}_modified.txt" \
        "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_dseg.nii.gz" -force

    rm "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_temp.nii.gz"

    # QC registration of atlas to dwi
    "${scriptdir}/check_atlasreg.py" \
        --subjid "${subj}${sessionfile}" \
        --atlas "${atlas}" \
        --atlas_image "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_dseg.nii.gz" \
        --nodif "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz" \
        --output "${outputdir}/dwi-preproc/${subj}${sessionpath}figures"

       
    # fslstats -K "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_dseg.nii.gz" \
    #     "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz" \
    #     -V > roivols.txt

    # labelfile="${atlaspath}/${ID}_modified.txt"
    # grep -v "Unknown" "${labelfile}" | awk '{ print $1,$2 }' > labels.txt
    # paste labels.txt roivols.txt > "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas-${atlas}_roivols.tsv"
    # rm roivols.txt labels.txt
done

# Transfer files
rsync -av "${workdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-dwi_res-high_atlas"* \
    "${outputdir}/dwi-preproc/${subj}${sessionpath}anat/"
rsync -av --ignore-existing "${SUBJECTS_DIR}/${subj}" "${freesurferdir}"

# Clean up
chmod -R u+w "${workdir}/${subj}/freesurfer/fsaverage"
# rm -rf "${workdir}/${subj}"

echo "-----------------------------------"
echo "finished anat2dwi subject = ${subj}"
echo "-----------------------------------"
