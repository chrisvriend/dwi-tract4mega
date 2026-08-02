#!/bin/bash
###############################################################################
# preproc.sh
# Author: C. Vriend - AUMC (modified: added generic multi-run reversed-PE branch)
# Date: Nov 05 2025 (modified Aug 2026)
# Description: prepare dwi scan for eddy and perform topup
#
# Two scenarios are now supported, auto-detected from what's on disk:
#
#  (A) SINGLE-RUN: one dwi.nii.gz for the subject/session, optionally with
#      short b0-only fieldmap(s) in fmap/ used purely to drive topup.
#      -> original pipeline, unchanged.
#
#  (B) MULTI-RUN / REVERSED PHASE-ENCODE: two or more FULL dwi runs for the
#      subject/session, named with the BIDS `dir-<label>` entity
#      (e.g. dir-AP_dwi.nii.gz / dir-PA_dwi.nii.gz, or dir-LR/dir-RL,
#      dir-IS/dir-SI, or any custom label). All runs are denoised,
#      deringed, used together to estimate topup's field, and then
#      concatenated (via dwicat, with b0 intensity matching) into a single
#      dataset with a matching multi-line acqparams/index.txt for eddy.
#      This is direction-LABEL agnostic: only the JSON
#      PhaseEncodingDirection (i/i-/j/j-/k/k-) and TotalReadoutTime fields
#      are used for the actual math; the `dir-<label>` string is only used
#      for naming intermediate files.
###############################################################################

set -euo pipefail
#export PATH=/opt/c3d/bin:$PATH
export FSLOUTPUTTYPE=NIFTI_GZ

on_error() {
    local exit_code=$?
    echo -e "\033[0;31m[dwi-02a] FAILED (exit ${exit_code}) at line ${BASH_LINENO[0]}: ${BASH_COMMAND}\033[0m" >&2
    echo -e "\033[0;31m[dwi-02a] subject=${subj:-?} session=${session:-?}\033[0m" >&2
}
trap on_error ERR

Usage() {
    echo "Usage: $0 -i <bidsdir> -o <outputdir> -w <workdir> -s <subj> -c <scriptdir> -t <nthreads> [-z <session>]"
    echo ""
    echo "Automatically detects and handles either:"
    echo "  - a single dwi run (dwi.nii.gz [+ fmap/ b0 fieldmaps]), or"
    echo "  - two or more full dwi runs with a BIDS dir-<label> entity"
    echo "    (e.g. dir-AP/dir-PA, dir-LR/dir-RL, dir-IS/dir-SI, or any"
    echo "    custom label) representing reversed phase-encode acquisitions,"
    echo "    which are combined for topup + concatenated for eddy."
    exit 1
}

# Define color variables
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[34m'
NC='\033[0m' # No Color

# initialize functions
# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

# Function to get the "opposite" of a PE direction
get_opposite_PE() {
    local PEdir="$1"
    if [[ "$PEdir" == *"-" ]]; then
        echo "${PEdir/-/}"
    else
        echo "${PEdir}-"
    fi
}

# Convert a BIDS PhaseEncodingDirection code (i, i-, j, j-, k, k-) into the
# FSL-style 3-vector used in acqparams files. This is agnostic to which
# anatomical axis (AP/PA, LR/RL, IS/SI) the code actually corresponds to --
# it only depends on the JSON field, never on filename labels.
pe_to_fsl() {
    local pe="$1"
    case "$pe" in
        i)  echo "1 0 0" ;;
        i-) echo "-1 0 0" ;;
        j)  echo "0 1 0" ;;
        j-) echo "0 -1 0" ;;
        k)  echo "0 0 1" ;;
        k-) echo "0 0 -1" ;;
        *)  echo "" ;;
    esac
}

# Human-readable direction label purely for logging (not used for any math)
pe_to_readable() {
    local pe="$1"
    case "$pe" in
        i)  echo "LR" ;;
        i-) echo "RL" ;;
        j)  echo "AP" ;;
        j-) echo "PA" ;;
        k)  echo "IS" ;;
        k-) echo "SI" ;;
        *)  echo "unknown" ;;
    esac
}

# Extract the dir-<label> BIDS entity from a dwi filename, whatever the
# label happens to be (AP, PA, LR, RL, IS, SI, or any custom string).
# Used purely for naming intermediate/output files.
get_dir_label() {
    local f
    f=$(basename "$1")
    echo "$f" | sed -n 's/.*_dir-\([A-Za-z0-9]*\)_dwi.*/\1/p'
}

# Initialize variables
bidsdir=""
outputdir=""
workdir=""
subj=""
session=""
scriptdir=""
nthreads=""

# input variables
# Parse command line arguments
while getopts ":i:o:w:s:c:z:t:" opt; do
    case $opt in
    i) bidsdir="$OPTARG" ;;
    o) outputdir="$OPTARG" ;;
    w) workdir="$OPTARG" ;;
    s) subj="$OPTARG" ;;
    z) session="$OPTARG" ;;
    c) scriptdir="$OPTARG" ;;
    t) nthreads="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done

missing=0
for var in bidsdir outputdir workdir subj scriptdir nthreads; do
    if [[ -z "${!var}" ]]; then
        echo "Error: $var is required."
        missing=1
    fi
done

if [[ $missing -eq 1 ]]; then
    Usage
fi

# Check if directories exist
for dir in "$bidsdir" "$scriptdir"; do
    if [[ ! -d "$dir" ]]; then
        echo "Error: Directory $dir does not exist."
        exit 1
    fi
done

# Preflight: fail fast (with a clear message) rather than deep inside the
# pipeline if a required tool is missing from the container/environment.
# dwicat in particular is only available in MRtrix3 >=3.0.4.
required_cmds=(
    jq fslnvols fslinfo fslmerge slicer
    mrconvert mrinfo mrmath dwidenoise mrdegibbs dwiextract dwicat
    topup mri_synthstrip antsRegistrationSyN.sh antsApplyTransforms N4BiasFieldCorrection
)
missing_cmds=()
for cmd in "${required_cmds[@]}"; do
    command -v "$cmd" >/dev/null 2>&1 || missing_cmds+=("$cmd")
done
if [[ ${#missing_cmds[@]} -gt 0 ]]; then
    echo -e "${RED}missing required command(s): ${missing_cmds[*]}${NC}"
    echo -e "${RED}(dwicat requires MRtrix3 >=3.0.4 -- check the container/module version)${NC}"
    exit 1
fi

mkdir -p ${workdir}
mkdir -p "${outputdir}/dwi-preproc"

# Set session path/file
if [[ -z "${session}" ]]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
fi

dwi_bids_dir="${bidsdir}/${subj}${sessionpath}dwi"

#----------------------------------------------------------------------
# detect scenario: single dwi run vs multiple full dwi runs (dir-<label>)
#----------------------------------------------------------------------
mapfile -t dwi_runs < <(ls ${dwi_bids_dir}/${subj}${sessionfile}dir-*_dwi.nii.gz 2>/dev/null | sort)

multirun=false
if [[ ${#dwi_runs[@]} -ge 2 ]]; then
    multirun=true
elif [[ ${#dwi_runs[@]} -eq 1 ]]; then
    echo -e "${RED}found exactly one dir-<label> dwi run ($(basename "${dwi_runs[0]}")) for ${subj} - ${session}${NC}"
    echo -e "${RED}the reversed phase-encode branch needs >=2 full dwi runs; a single dir-labeled run is not supported${NC}"
    echo -e "${RED}either provide a second dir-<label> run, or rename this one to the plain '${subj}${sessionfile}dwi.nii.gz' convention (with matching .bval/.bvec/.json) to use the single-run branch${NC}"
    exit 1
elif [[ ! -f ${dwi_bids_dir}/${subj}${sessionfile}dwi.nii.gz || ! -f ${dwi_bids_dir}/${subj}${sessionfile}dwi.bvec ]]; then
    echo -e "${RED}no dwi scan/bvec found for ${subj} - ${session}${NC}"
    exit 1
fi

echo -e ${YELLOW}----------------------${NC}
echo -e ${YELLOW}Preprocessing dwi data${NC}
echo -e ${YELLOW}${subj}${NC}
echo -e ${YELLOW}${session}${NC}
if [[ "$multirun" == true ]]; then
    echo -e ${YELLOW}${#dwi_runs[@]} full dwi runs detected -- reversed phase-encode combination${NC}
fi
echo -e ${YELLOW}----------------------${NC}

fmap_samePE=()
fmap_otherPE=()

if [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz ]] &&
    [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval ]] &&
    [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec ]]; then
    log "${GREEN}" "${subj}${sessionfile} already preprocessed with eddy"
    exit 0
fi

mkdir -p "${workdir}/${subj}${sessionpath}dwi"
mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi"
mkdir -p "${outputdir}/dwi-preproc/${subj}/log"
mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}fmap"
mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}figures"
mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}qc"

###############################################################################
# MULTI-RUN (reversed phase-encode) BRANCH
###############################################################################
if [[ "$multirun" == true ]]; then

    dwiworkdir="${workdir}/${subj}${sessionpath}dwi"
    fmapworkdir="${workdir}/${subj}${sessionpath}fmap"
    mkdir -p "${dwiworkdir}" "${fmapworkdir}"

    log "${BLUE}" "runs found:"
    for r in "${dwi_runs[@]}"; do
        log "${BLUE}" "  $(basename "$r")"
    done

    run_labels=()
    run_pe=()
    run_pefsl=()
    run_trt=()
    run_mif=()
    run_nvols=()
    run_json=()

    for nii in "${dwi_runs[@]}"; do
        base="${nii%.nii.gz}"
        json="${base}.json"
        bvec="${base}.bvec"
        bval="${base}.bval"
        label=$(get_dir_label "$nii")

        if [[ -z "$label" ]]; then
            log "${RED}" "could not parse a dir-<label> entity from $(basename "$nii")"
            exit 1
        fi
        for f in "$json" "$bvec" "$bval"; do
            if [[ ! -f "$f" ]]; then
                log "${RED}" "missing $f for run dir-${label}"
                exit 1
            fi
        done

        pe=$(jq -r '.PhaseEncodingDirection' "$json")
        trt=$(jq -r '.TotalReadoutTime' "$json")
        pefsl=$(pe_to_fsl "$pe")

        if [[ -z "$pe" || "$pe" == "null" || -z "$trt" || "$trt" == "null" ]]; then
            log "${RED}" "no PhaseEncodingDirection or TotalReadoutTime in ${json}"
            exit 1
        fi
        if [[ -z "$pefsl" ]]; then
            log "${RED}" "unrecognized PhaseEncodingDirection '${pe}' in ${json}"
            exit 1
        fi

        log "${BLUE}" "  dir-${label}: PE=${pe} ($(pe_to_readable "$pe"), ${pefsl})  TRT=${trt}"

        # round bvals up front so they're embedded correctly from the start
        rounded_bval="${dwiworkdir}/${subj}${sessionfile}dir-${label}_dwi.bval"
        cp "$bval" "$rounded_bval"
        chmod u+rw "$rounded_bval"
        ${scriptdir}/round_bvals.py "$rounded_bval"

        mif_raw="${dwiworkdir}/${subj}${sessionfile}dir-${label}_dwi.mif"
        mif_dns="${dwiworkdir}/${subj}${sessionfile}dir-${label}_desc-dns_dwi.mif"
        mif_deg="${dwiworkdir}/${subj}${sessionfile}dir-${label}_desc-dns+degibbs_dwi.mif"
        noise_nii="${dwiworkdir}/${subj}${sessionfile}dir-${label}_desc-noise_dwi.nii.gz"

        if [[ ! -f "$mif_deg" ]]; then
            mrconvert "$nii" "$mif_raw" -fslgrad "$bvec" "$rounded_bval" \
                -json_import "$json" -force

            dwidenoise "$mif_raw" "$mif_dns" \
                -noise "$noise_nii" -nthreads ${nthreads} -force

            mrdegibbs "$mif_dns" "$mif_deg" -nthreads ${nthreads} -force

            rm -f "$mif_raw" "$mif_dns"
        fi

        nvols=$(mrinfo -size "$mif_deg" | awk '{print $NF}')

        run_labels+=("$label")
        run_pe+=("$pe")
        run_pefsl+=("$pefsl")
        run_trt+=("$trt")
        run_mif+=("$mif_deg")
        run_nvols+=("$nvols")
        run_json+=("$json")
    done

    # sanity check: need at least two distinct PE directions to pair for topup
    uniq_pe=$(printf '%s\n' "${run_pe[@]}" | sort -u | wc -l)
    if [[ "$uniq_pe" -lt 2 ]]; then
        log "${RED}" "all detected dwi runs share the same PhaseEncodingDirection -- nothing to pair for topup"
        exit 1
    fi

    # sanity check: all runs must be on the same voxel grid, or fslmerge/
    # dwicat will error (or silently misbehave) further down the pipeline
    ref_size=$(mrinfo -size "${run_mif[0]}" | awk '{print $1, $2, $3}')
    ref_vox=$(mrinfo -vox "${run_mif[0]}" | awk '{printf "%.3f %.3f %.3f", $1, $2, $3}')
    for i in "${!run_labels[@]}"; do
        this_size=$(mrinfo -size "${run_mif[$i]}" | awk '{print $1, $2, $3}')
        this_vox=$(mrinfo -vox "${run_mif[$i]}" | awk '{printf "%.3f %.3f %.3f", $1, $2, $3}')
        if [[ "$this_size" != "$ref_size" ]]; then
            log "${RED}" "run dir-${run_labels[$i]} has matrix size (${this_size}) that does not match dir-${run_labels[0]} (${ref_size})"
            exit 1
        fi
        if [[ "$this_vox" != "$ref_vox" ]]; then
            log "${RED}" "run dir-${run_labels[$i]} has voxel size (${this_vox}) that does not match dir-${run_labels[0]} (${ref_vox})"
            exit 1
        fi
    done

    # sanity check: warn (don't fail) if TotalReadoutTime differs a lot
    # between runs -- usually indicates a metadata/protocol mismatch
    ref_trt="${run_trt[0]}"
    for i in "${!run_labels[@]}"; do
        pct=$(awk -v a="${run_trt[$i]}" -v b="${ref_trt}" 'BEGIN { if (b==0) print "nan"; else printf "%.1f", (100*(a-b)/b < 0 ? -(100*(a-b)/b) : 100*(a-b)/b) }')
        if awk -v p="$pct" 'BEGIN { exit !(p != "nan" && p+0 > 20) }'; then
            log "${YELLOW}" "TotalReadoutTime for dir-${run_labels[$i]} (${run_trt[$i]}) differs by ${pct}% from dir-${run_labels[0]} (${ref_trt}) -- check the protocol/metadata"
        fi
    done

    # informational: are the runs sampling different gradient directions
    # (more angular coverage when combined) or repeats of the same scheme
    # (mainly a distortion-correction/SNR benefit)?
    if [[ ${#run_labels[@]} -eq 2 ]]; then
        b0_1="${dwiworkdir}/${subj}${sessionfile}dir-${run_labels[0]}_dwi.bval"
        b0_2="${dwiworkdir}/${subj}${sessionfile}dir-${run_labels[1]}_dwi.bval"
        if [[ -f "$b0_1" && -f "$b0_2" ]] && diff -q "$b0_1" "$b0_2" >/dev/null 2>&1; then
            log "${BLUE}" "dir-${run_labels[0]} and dir-${run_labels[1]} share identical b-values -- likely a repeated scheme (distortion-correction + SNR benefit; angular coverage unchanged)"
        else
            log "${BLUE}" "dir-${run_labels[0]} and dir-${run_labels[1]} sample different b-values/directions -- combining increases angular coverage"
        fi
    fi

    combined_label=$(IFS=+; echo "${run_labels[*]}")

    #----------------------------------------------------------------------
    # build merged b0 series + refparams for topup (any number of runs,
    # any PE axes, any direction labels)
    #----------------------------------------------------------------------
    refparams="${fmapworkdir}/${subj}${sessionfile}dir-${combined_label}_desc-refparams.tsv"
    rm -f "$refparams"

    b0_niis=()
    for i in "${!run_labels[@]}"; do
        label="${run_labels[$i]}"
        b0_nii="${fmapworkdir}/${subj}${sessionfile}dir-${label}_desc-b0s_epi.nii.gz"
        dwiextract -nthreads ${nthreads} "${run_mif[$i]}" - -bzero |
            mrconvert - "$b0_nii" -force
        nb0=$(fslnvols "$b0_nii")
        if [[ "$nb0" -eq 0 ]]; then
            log "${RED}" "run dir-${label} has no b=0 volumes -- cannot use it for topup"
            exit 1
        fi
        for ((v = 0; v < nb0; v++)); do
            echo "${run_pefsl[$i]} ${run_trt[$i]}" >>"$refparams"
        done
        b0_niis+=("$b0_nii")
    done

    # quick QC: mean b0 per run, pre-topup
    figdir="${outputdir}/dwi-preproc/${subj}${sessionpath}figures"
    for i in "${!run_labels[@]}"; do
        mean_b0="${fmapworkdir}/${subj}${sessionfile}dir-${run_labels[$i]}_desc-meanb0_epi.nii.gz"
        mrmath "${b0_niis[$i]}" mean "$mean_b0" -axis 3 -force
        slicer "$mean_b0" -a "${figdir}/${subj}${sessionfile}dir-${run_labels[$i]}_desc-pretopup_epi.png" >/dev/null 2>&1 || true
        rm -f "$mean_b0"
    done

    topup_input="${fmapworkdir}/${subj}${sessionfile}dir-${combined_label}_space-dwi_desc-4topup_epi.nii.gz"
    fslmerge -t "$topup_input" "${b0_niis[@]}"
    rm -f "${b0_niis[@]}"

    rsync -a "$topup_input" "$refparams" "${outputdir}/dwi-preproc/${subj}${sessionpath}fmap/"

    #----------------------------------------------------------------------
    # topup
    #----------------------------------------------------------------------
    unwarped="${fmapworkdir}/${subj}${sessionfile}space-dwi_desc-unwarped_epi.nii.gz"
    topup_out="${fmapworkdir}/${subj}${sessionfile}space-dwi_desc-topup"

    if [[ ! -f "$unwarped" || ! -f "${topup_out}_fieldcoef.nii.gz" ]]; then
        cd "$fmapworkdir"

        dim1=$(fslinfo "$topup_input" | grep -w dim1 | awk '{ print $2 }' | awk '{print int($0)}')
        dim2=$(fslinfo "$topup_input" | grep -w dim2 | awk '{ print $2 }' | awk '{print int($0)}')
        dim3=$(fslinfo "$topup_input" | grep -w dim3 | awk '{ print $2 }' | awk '{print int($0)}')

        if ((dim1 % 4 == 0 && dim2 % 4 == 0 && dim3 % 4 == 0)); then
            log "${BLUE}" "All dimensions are integer multiples of 4; using b02b0_4.cnf for topup"
            configfile=b02b0_4.cnf
        elif ((dim1 % 2 == 0 && dim2 % 2 == 0 && dim3 % 2 == 0)); then
            log "${BLUE}" "All dimensions are integer multiples of 2; using b02b0_2.cnf for topup"
            configfile=b02b0_2.cnf
        else
            log "${BLUE}" "At least one dimension is odd; using b02b0_1.cnf as config file for topup"
            configfile=b02b0_1.cnf
        fi

        echo
        echo -e "${BLUE}running topup${NC}"
        echo

        topup --imain=$(basename "$topup_input") \
            --datain=$(basename "$refparams") \
            --config=${configfile} \
            --out=${subj}${sessionfile}space-dwi_desc-topup \
            --iout=${subj}${sessionfile}space-dwi_desc-unwarped_epi \
            --fout=${subj}${sessionfile}space-dwi_desc-topup_fieldmap \
            --verbose >${subj}${sessionfile}topup_$(date +"%Y-%m-%d_%H-%M").log

        cp ${subj}${sessionfile}topup_*.log ${outputdir}/dwi-preproc/${subj}/log
    fi

    # quick QC: mean unwarped b0, post-topup
    mean_unwarped="${fmapworkdir}/${subj}${sessionfile}desc-meanunwarped_epi.nii.gz"
    mrmath "$unwarped" mean "$mean_unwarped" -axis 3 -force
    slicer "$mean_unwarped" -a "${figdir}/${subj}${sessionfile}space-dwi_desc-posttopup_epi.png" >/dev/null 2>&1 || true
    rm -f "$mean_unwarped"

    #----------------------------------------------------------------------
    # concatenate the full denoised/deringed series with dwicat (matches
    # b0 intensity scaling across runs), then export back to nii+bval+bvec
    #----------------------------------------------------------------------
    cat_mif="${dwiworkdir}/${subj}${sessionfile}dir-${combined_label}_desc-dns+degibbs_dwi.mif"

    if [[ ! -f "$cat_mif" ]]; then
        # rough brain mask (from the first run's mean b0, native/distorted
        # space) so dwicat's b0 intensity matching isn't skewed by
        # background/non-brain voxels. This is a throwaway mask, not the
        # analysis mask used later for eddy.
        dwicat_mask="${dwiworkdir}/${subj}${sessionfile}desc-dwicattmp_mask.nii.gz"
        first_meanb0="${dwiworkdir}/${subj}${sessionfile}desc-dwicattmp_meanb0.nii.gz"
        dwiextract -nthreads ${nthreads} "${run_mif[0]}" - -bzero |
            mrmath - mean "$first_meanb0" -axis 3 -force
        mri_synthstrip -i "$first_meanb0" --mask "$dwicat_mask" >/dev/null 2>&1 || true

        if [[ -f "$dwicat_mask" ]]; then
            dwicat "${run_mif[@]}" "$cat_mif" -mask "$dwicat_mask" -nthreads ${nthreads} -force
        else
            log "${YELLOW}" "could not build a rough mask for dwicat -- running without -mask"
            dwicat "${run_mif[@]}" "$cat_mif" -nthreads ${nthreads} -force
        fi
        rm -f "$dwicat_mask" "$first_meanb0"
    fi

    cat_nii="${dwiworkdir}/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz"
    cat_bvec="${dwiworkdir}/${subj}${sessionfile}dwi.bvec"
    cat_bval="${dwiworkdir}/${subj}${sessionfile}dwi.bval"

    mrconvert "$cat_mif" "$cat_nii" \
        -export_grad_fsl "$cat_bvec" "$cat_bval" -force

    # acqparams: one line per run (same PE/TRT used for that run's b0s above)
    acqp="${dwiworkdir}/${subj}${sessionfile}acq-dwi_desc-acqparams.tsv"
    rm -f "$acqp"
    for i in "${!run_labels[@]}"; do
        echo "${run_pefsl[$i]} ${run_trt[$i]}" >>"$acqp"
    done

    # index.txt: maps every concatenated volume to the correct acqparams
    # line, in the same run order used when calling dwicat above
    index="${dwiworkdir}/index.txt"
    rm -f "$index"
    for i in "${!run_labels[@]}"; do
        line=$((i + 1))
        printf "${line} %.0s" $(seq 1 "${run_nvols[$i]}") >>"$index"
    done
    echo >>"$index"

    # bookkeeping: which acqparams line / run each block of volumes came from
    runkey="${dwiworkdir}/${subj}${sessionfile}acq-dwi_desc-runkey.tsv"
    {
        for i in "${!run_labels[@]}"; do
            echo -e "line$((i + 1))\tdir-${run_labels[$i]}\tPE=${run_pe[$i]}\tTRT=${run_trt[$i]}\tnvols=${run_nvols[$i]}"
        done
    } >"$runkey"

    # combined json for eddy: only carry SliceTiming through (needed for
    # method=volcorr/volcorrnosdc in dwi-02b) if it's identical across runs
    combined_json="${dwiworkdir}/${subj}${sessionfile}dwi.json"
    first_st=$(jq -c '.SliceTiming // empty' "${run_json[0]}")
    same_st=true
    for j in "${run_json[@]}"; do
        st=$(jq -c '.SliceTiming // empty' "$j")
        if [[ "$st" != "$first_st" ]]; then
            same_st=false
        fi
    done
    if [[ "$same_st" == true && -n "$first_st" ]]; then
        cp "${run_json[0]}" "$combined_json"
    else
        log "${YELLOW}" "SliceTiming differs or is absent across runs -- slice-to-volume eddy methods (volcorr/volcorrnosdc) will not be usable"
        jq 'del(.SliceTiming)' "${run_json[0]}" >"$combined_json"
    fi

    #######################
    ## create brain mask ##
    #######################
    cd "$dwiworkdir"

    # mean of the topup-corrected (undistorted, mid-space) b0s
    mrmath "${unwarped}" mean \
        ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz -axis 3 -force

    # mean of the uncorrected b0s, in the native (distorted) space that
    # eddy's imain sits in
    dwiextract -nthreads ${nthreads} "$cat_nii" - -bzero \
        -fslgrad "$cat_bvec" "$cat_bval" |
        mrmath - mean ${subj}${sessionfile}space-dwi_desc-meanb0-uncorrected_dwi.nii.gz -axis 3 -force

    # rigid-register the topup mid-space mean b0 onto the native distorted
    # mean b0, so the brain mask ends up in the space eddy expects
    antsRegistrationSyN.sh -d 3 -m ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz \
        -f ${subj}${sessionfile}space-dwi_desc-meanb0-uncorrected_dwi.nii.gz \
        -o ${subj}${sessionfile}rigidreg -t r -n ${nthreads} -p d
    mv ${subj}${sessionfile}rigidregWarped.nii.gz ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz
    rm -f *rigidreg*

    if [[ ! -f ${dwiworkdir}/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz ]]; then
        mri_synthstrip \
            -i ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz \
            -o ${subj}${sessionfile}space-dwi_desc-nodif-brain_dwi.nii.gz \
            --mask ${dwiworkdir}/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz
    fi

    cd "$fmapworkdir"
    rsync -av ${subj}${sessionfile}space-dwi_desc-unwarped_epi* \
        ${subj}${sessionfile}space-dwi_desc-topup_fieldmap* \
        ${subj}${sessionfile}space-dwi_desc-topup* \
        ${outputdir}/dwi-preproc/${subj}${sessionpath}fmap

    # for QC
    for i in "${!run_labels[@]}"; do
        rsync -a "${dwiworkdir}/${subj}${sessionfile}dir-${run_labels[$i]}_desc-noise_dwi.nii.gz" \
            "${outputdir}/dwi-preproc/${subj}${sessionpath}qc/" 2>/dev/null || true
    done
    rsync -a "$runkey" "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/"

    log "${GREEN}" "Preprocessing complete (multi-run / reversed phase-encode branch)."
    exit 0
fi
###############################################################################
# END MULTI-RUN BRANCH -- everything below is the original single-run
# pipeline, unchanged.
###############################################################################

# Specify the path to the DWI JSON sidecar
dwi_json_path=$(ls ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}dwi.json)

# extract TotalReadoutTime
dwi_trt=$(cat ${dwi_json_path} | jq -r '.TotalReadoutTime')
dwi_PE=$(cat ${dwi_json_path} | jq -r '.PhaseEncodingDirection')

if [ -z ${dwi_trt} ] || [ -z ${dwi_PE} ]; then
    echo -e "${RED}no TotalReadOutTime or PhaseEncodingDirection found in dwi json file${NC}"
    echo
    exit 1
fi

# determine settings for topup
if [ ${dwi_PE} == "i" ]; then PE_dwi_FSL="1 0 0"
elif [ ${dwi_PE} == "i-" ]; then PE_dwi_FSL="-1 0 0"
elif [ ${dwi_PE} == "j" ]; then PE_dwi_FSL="0 1 0"
elif [ ${dwi_PE} == "j-" ]; then PE_dwi_FSL="0 -1 0"
elif [ ${dwi_PE} == "k" ]; then PE_dwi_FSL="0 0 1"
elif [ ${dwi_PE} == "k-" ]; then PE_dwi_FSL="0 0 -1"
fi

#----------------------------------------------------------------------
# MP-PCA denoising & deringing of dwi scan
#----------------------------------------------------------------------
if [ ! -f ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz ]; then
    dwidenoise ${bidsdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.nii.gz \
        ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns_dwi.mif \
        -noise ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-noise_dwi.nii.gz \
        -nthreads ${nthreads} -force

    #Remove Gibbs Ringing Artifacts
    mrdegibbs ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns_dwi.mif \
        ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz \
        -nthreads ${nthreads} -force

    rm ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns_dwi.mif
fi

# write dwi acqparams
echo "${PE_dwi_FSL} ${dwi_trt}" >${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}acq-dwi_desc-acqparams.tsv

#----------------------------------------------------------------------
# Check for fieldmaps
#----------------------------------------------------------------------
# Specify the path to the fieldmap folder
fieldmap_folder=${bidsdir}/${subj}${sessionpath}/fmap

# Check if the fieldmap folder exists
if [ -d ${fieldmap_folder} ]; then
    echo -e "${YELLOW} fieldmap folder found${NC}"
    mkdir -p ${workdir}/${subj}${sessionpath}fmap

    for fmap_json in ${fieldmap_folder}/*acq-dwi*dir*epi.json; do
        if [[ $(jq -r '.IntendedFor' "${fmap_json}") == *dwi* ]]; then
            fmap_nii=${fmap_json%%.json}.nii.gz
            fmap_PE=$(cat ${fmap_json} | jq -r '.PhaseEncodingDirection')
            fmap_trt=$(jq -r '.TotalReadoutTime' "$fmap_json")

            echo "${fmap_json}"
            echo -e "${BLUE}PhaseEncodingDirection: $fmap_PE${NC}"
            echo -e "${BLUE}TotalReadoutTime: $fmap_trt${NC}"
            echo

            if [[ -f ${fmap_nii} ]]; then
                if [ "${fmap_PE}" == "${dwi_PE}" ]; then
                    fmap_samePE+=("${fmap_nii}")
                else
                    fmap_otherPE+=("${fmap_nii}")
                fi
            fi
        fi
    done

    for fmap in "${fmap_samePE:-empty}" "${fmap_otherPE:-empty}"; do
        if [[ "$fmap" == "empty" ]]; then
            continue
        fi
        if [ ! -z ${fmap} ]; then
            fmap_json=${fmap%%.nii.gz}.json
            fmap_PE=$(cat ${fmap_json} | jq -r '.PhaseEncodingDirection')

            if [[ "$fmap_PE" == "j" ]]; then dir=AP
            elif [[ "$fmap_PE" == "j-" ]]; then dir=PA
            elif [[ "$fmap_PE" == "i" ]]; then dir=LR
            elif [[ "$fmap_PE" == "i-" ]]; then dir=RL
            elif [[ "$fmap_PE" == "k" ]]; then dir=IS
            elif [[ "$fmap_PE" == "k-" ]]; then dir=SI
            else echo "Unknown Phase Encoding Direction"
            fi

            # if fmap has multiple volumes
            if [[ $(fslnvols ${fieldmap_folder}/${fmap}) -gt 1 ]]; then
                if [[ ! -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dir}_space-dwi_desc-degibbs_epi.nii.gz ]]; then
                    dwidenoise ${fmap} \
                        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dir}_space-dwi_desc-dns_epi.mif -force

                    #Remove Gibbs Ringing Artifacts
                    mrdegibbs ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dir}_space-dwi_desc-dns_epi.mif \
                        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dir}_space-dwi_desc-degibbs_epi.nii.gz -force

                    rm ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dir}_space-dwi_desc-dns_epi.mif
                fi
            else
                mrdegibbs ${fmap} \
                    ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dir}_space-dwi_desc-degibbs_epi.nii.gz -force
            fi
        fi
    done
fi

#----------------------------------------------------------------------
# One Fieldmap available
#----------------------------------------------------------------------
# in case only fmap with opposite but not same PE available
if [ ${#fmap_samePE[@]} -eq 0 ] && [ ${#fmap_otherPE[@]} -ne 0 ]; then
    echo -e "${BLUE} one fieldmap available in fmap folder${NC}"
    echo

    fmap_otherjson=${fmap_otherPE%%.nii.gz}.json
    PE_other=$(cat ${fmap_otherjson} | jq -r '.PhaseEncodingDirection')
    other_trt=$(cat ${fmap_otherjson} | jq -r '.TotalReadoutTime')
    opposite_pe1=$(get_opposite_PE "$PE_other")

    # Check if the second direction matches the opposite of the first
    if [[ "$dwi_PE" == "$opposite_pe1" ]]; then
        echo -e "${GREEN} PE directions of dwi and fmap are opposites.${NC}"
    else
        echo -e "${RED}PE directions of dwi and fmap are NOT opposites.${NC}"
        exit 1
    fi

    # determine letter for opposite PE
    if [[ "$dwi_PE" == "j" ]]; then dwidir=AP
    elif [[ "$dwi_PE" == "j-" ]]; then dwidir=PA
    elif [[ "$dwi_PE" == "i" ]]; then dwidir=LR
    elif [[ "$dwi_PE" == "i-" ]]; then dwidir=RL
    elif [[ "$dwi_PE" == "k" ]]; then dwidir=IS
    elif [[ "$dwi_PE" == "k-" ]]; then dwidir=SI
    else echo "Unknown Phase Encoding Direction"
    fi

    if [[ "$PE_other" == "j" ]]; then otherdir=AP
    elif [[ "$PE_other" == "j-" ]]; then otherdir=PA
    elif [[ "$PE_other" == "i" ]]; then otherdir=LR
    elif [[ "$PE_other" == "i-" ]]; then otherdir=RL
    elif [[ "$PE_other" == "k" ]]; then otherdir=IS
    elif [[ "$PE_other" == "k-" ]]; then otherdir=SI
    else echo "Unknown Phase Encoding Direction"
    fi

    # extract mean b0 from dwi
    dwiextract -nthreads ${nthreads} \
        ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz - -bzero \
        -fslgrad ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}*dwi.bvec ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}*dwi.bval |
        mrmath - mean ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-temp_epi.nii.gz -axis 3

    # create mean b0 from PA fieldmap
    if [[ $(fslnvols ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz) -gt 1 ]]; then
        mrmath ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz \
            mean \
            ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-temp_epi.nii.gz -axis 3
    else
        ln -s ${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz \
            ${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-temp_epi.nii.gz
    fi

    # rigid registration of dwi mean b0 and other-PE fieldmap
    antsRegistrationSyN.sh -d 3 -m ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-temp_epi.nii.gz \
        -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-temp_epi.nii.gz \
        -o ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}rigidreg -t r -n ${nthreads} -p d

    # apply to multi-volume fieldmap
    antsApplyTransforms -d 3 -e 3 -i ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz \
        -r ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-temp_epi.nii.gz \
        -t ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}rigidreg0GenericAffine.mat \
        -o ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-warped-degibbs_epi.nii.gz -v -u int

    rm -f ${workdir}/${subj}${sessionpath}fmap/*rigidreg*

    # merge dwi mean b0 and registered fieldmap
    fslmerge -t ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-temp_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-warped-degibbs_epi.nii.gz

    rm ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-warped-degibbs_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}*temp*.nii.gz

    if [ ${PE_other} == "i" ]; then PE_other_FSL="1 0 0"
    elif [ ${PE_other} == "i-" ]; then PE_other_FSL="-1 0 0"
    elif [ ${PE_other} == "j" ]; then PE_other_FSL="0 1 0"
    elif [ ${PE_other} == "j-" ]; then PE_other_FSL="0 -1 0"
    elif [ ${PE_other} == "k" ]; then PE_other_FSL="0 0 1"
    elif [ ${PE_other} == "k-" ]; then PE_other_FSL="0 0 -1"
    fi

    # write TRT to refparams file
    cd ${workdir}/${subj}${sessionpath}fmap
    echo "${PE_dwi_FSL} ${dwi_trt}" >${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv
    for ((i = 0; i < $(fslnvols ${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz); i++)); do
        echo "${PE_other_FSL} ${other_trt}" >>"${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv"
    done

    # set dwidir to samedir for later steps
    samedir=${dwidir}

elif [ ${#fmap_samePE[@]} -ne 0 ] && [ ${#fmap_otherPE[@]} -ne 0 ]; then
    #----------------------------------------------------------------------
    # two Fieldmaps available
    #----------------------------------------------------------------------
    fmap_otherjson=${fmap_otherPE%%.nii.gz}.json
    fmap_samejson=${fmap_samePE%%.nii.gz}.json

    PE_other=$(cat ${fmap_otherjson} | jq -r '.PhaseEncodingDirection')
    PE_same=$(cat ${fmap_samejson} | jq -r '.PhaseEncodingDirection')
    other_trt=$(cat ${fmap_otherjson} | jq -r '.TotalReadoutTime')
    same_trt=$(cat ${fmap_samejson} | jq -r '.TotalReadoutTime')

    # Get the opposite of the first direction
    opposite_pe1=$(get_opposite_PE "$PE_other")

    if [[ "$PE_same" == "$opposite_pe1" ]]; then
        echo -e "${GREEN}The fmap PE directions are opposites.${NC}"
    else
        echo -e "${RED}The fmap PE directions are NOT opposites.${NC}"
        exit 1
    fi

    if [[ "$dwi_PE" == "$PE_same" ]]; then
        echo -e "${GREEN} PE directions of dwi and 'same' fmap are consistent.${NC}"
    else
        echo -e "${RED}PE directions of dwi and fmap are NOT opposites.${NC}"
        exit 1
    fi

    opposite_pe1=$(get_opposite_PE "$PE_other")
    if [[ "$dwi_PE" == "$opposite_pe1" ]]; then
        echo -e "${GREEN} PE directions of dwi and 'opposite' fmap are consistent.${NC}"
    else
        echo -e "${RED}PE directions of dwi and fmap are NOT opposites.${NC}"
        exit 1
    fi

    # get directions
    if [[ "$PE_same" == "j" ]]; then samedir=AP
    elif [[ "$PE_same" == "j-" ]]; then samedir=PA
    elif [[ "$PE_same" == "i" ]]; then samedir=LR
    elif [[ "$PE_same" == "i-" ]]; then samedir=RL
    elif [[ "$PE_same" == "k" ]]; then samedir=IS
    elif [[ "$PE_same" == "k-" ]]; then samedir=SI
    else echo "Unknown Phase Encoding Direction"
    fi

    if [[ "$PE_other" == "j" ]]; then otherdir=AP
    elif [[ "$PE_other" == "j-" ]]; then otherdir=PA
    elif [[ "$PE_other" == "i" ]]; then otherdir=LR
    elif [[ "$PE_other" == "i-" ]]; then otherdir=RL
    elif [[ "$PE_other" == "k" ]]; then otherdir=IS
    elif [[ "$PE_other" == "k-" ]]; then otherdir=SI
    else echo "Unknown Phase Encoding Direction"
    fi

    if [ ${PE_same} == "i" ]; then PE_same_FSL="1 0 0"
    elif [ ${PE_same} == "i-" ]; then PE_same_FSL="-1 0 0"
    elif [ ${PE_same} == "j" ]; then PE_same_FSL="0 1 0"
    elif [ ${PE_same} == "j-" ]; then PE_same_FSL="0 -1 0"
    elif [ ${PE_same} == "k" ]; then PE_same_FSL="0 0 1"
    elif [ ${PE_same} == "k-" ]; then PE_same_FSL="0 0 -1"
    fi

    if [ ${PE_other} == "i" ]; then PE_other_FSL="1 0 0"
    elif [ ${PE_other} == "i-" ]; then PE_other_FSL="-1 0 0"
    elif [ ${PE_other} == "j" ]; then PE_other_FSL="0 1 0"
    elif [ ${PE_other} == "j-" ]; then PE_other_FSL="0 -1 0"
    elif [ ${PE_other} == "k" ]; then PE_other_FSL="0 0 1"
    elif [ ${PE_other} == "k-" ]; then PE_other_FSL="0 0 -1"
    fi

    echo -e "${BLUE}merge blip up/down scans for topup${NC}"
    cd ${workdir}/${subj}${sessionpath}fmap
    rm -f ${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv

    for ((i = 0; i < $(fslnvols ${subj}${sessionfile}dir-${samedir}_space-dwi_desc-degibbs_epi.nii.gz); i++)); do
        echo "${PE_same_FSL} ${same_trt}" >>"${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv"
    done
    for ((i = 0; i < $(fslnvols ${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz); i++)); do
        echo "${PE_other_FSL} ${other_trt}" >>"${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv"
    done

    fslmerge -t ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${samedir}_space-dwi_desc-degibbs_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz

    refnvols=$(fslnvols ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz)

    rsync -a ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz \
        ${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv \
        ${outputdir}/dwi-preproc/${subj}${sessionpath}fmap

elif [ ${#fmap_samePE[@]} -eq 0 ] && [ ${#fmap_otherPE[@]} -eq 0 ]; then
    #----------------------------------------------------------------------
    # Syn b0 (no fieldmaps available)
    #----------------------------------------------------------------------
    echo -e "${BLUE}no fmaps found - creating syn b0 for topup${NC}"
    mkdir -p "${workdir}/${subj}${sessionpath}fmap/"

    if [[ "$dwi_PE" == "j" ]]; then dwidir=AP
    elif [[ "$dwi_PE" == "j-" ]]; then dwidir=PA
    elif [[ "$dwi_PE" == "i" ]]; then dwidir=LR
    elif [[ "$dwi_PE" == "i-" ]]; then dwidir=RL
    elif [[ "$dwi_PE" == "k" ]]; then dwidir=IS
    elif [[ "$dwi_PE" == "k-" ]]; then dwidir=SI
    else echo "Unknown Phase Encoding Direction"
    fi

    PE_other=$(get_opposite_PE "$dwi_PE")

    if [[ "$PE_other" == "j" ]]; then otherdir=AP
    elif [[ "$PE_other" == "j-" ]]; then otherdir=PA
    elif [[ "$PE_other" == "i" ]]; then otherdir=LR
    elif [[ "$PE_other" == "i-" ]]; then otherdir=RL
    elif [[ "$PE_other" == "k" ]]; then otherdir=IS
    elif [[ "$PE_other" == "k-" ]]; then otherdir=SI
    else echo "Unknown Phase Encoding Direction"
    fi

    # write TRT to refparams file
    echo "${PE_dwi_FSL} ${dwi_trt}" >${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv
    echo "${PE_dwi_FSL} 0.00" >>${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv

    if [[ ! -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz ]]; then
        mkdir -p "${workdir}/${subj}${sessionpath}fmap/synb0/tmp" \
            "${workdir}/${subj}${sessionpath}fmap/synb0/input" \
            "${workdir}/${subj}${sessionpath}fmap/synb0/output"

        cp ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv \
            ${workdir}/${subj}${sessionpath}fmap/synb0/input/

        # extract first b0 vol from dwi
        dwiextract -nthreads ${nthreads} \
            ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz - -bzero \
            -fslgrad ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}*dwi.bvec ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}*dwi.bval |
            mrconvert - -coord 3 0 \
                ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-b0_epi.nii.gz -force

        rsync -a ${bidsdir}/${subj}${sessionpath}anat/${subj}${sessionfile}T1w.nii.gz \
            ${workdir}/${subj}${sessionpath}fmap/synb0/input

        N4BiasFieldCorrection -d 3 \
            -i ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}T1w.nii.gz \
            -o ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}desc-n4_T1w.nii.gz

        mri_synthstrip \
            -i ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}desc-n4_T1w.nii.gz \
            -o ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}desc-brain_T1w.nii.gz \
            --mask ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}space-T1w_desc-brain_mask.nii.gz

        cd ${workdir}/${subj}${sessionpath}fmap/synb0/input
        if [ -L T1.nii.gz ]; then
            unlink T1.nii.gz
            unlink BRAIN.nii.gz
            unlink acqparams.txt
            unlink b0.nii.gz
        fi
        ln -s ${subj}${sessionfile}desc-n4_T1w.nii.gz T1.nii.gz
        ln -s ${subj}${sessionfile}desc-brain_T1w.nii.gz BRAIN.nii.gz
        ln -s ${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv acqparams.txt
        ln -s ${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-b0_epi.nii.gz b0.nii.gz

        #Run Synb0-DISCO for fieldmap-free distortion correction
        if [[ ! -f ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_d_smooth.nii.gz ]] ||
            [[ ! -f ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_u.nii.gz ]]; then
            if [[ "${lowmem}" -eq 1 ]]; then
                synb0 --input ${workdir}/${subj}${sessionpath}fmap/synb0/input \
                    --output ${workdir}/${subj}${sessionpath}fmap/synb0/output --notopup --lowmem
            else
                synb0 --input ${workdir}/${subj}${sessionpath}fmap/synb0/input \
                    --output ${workdir}/${subj}${sessionpath}fmap/synb0/output --notopup
            fi
        fi

        fslmerge -t ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_all.nii.gz \
            ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_d_smooth.nii.gz \
            ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_u.nii.gz &&
            mv ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_all.nii.gz \
                ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz

        mv ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv \
            ${workdir}/${subj}${sessionpath}fmap/

        if [[ -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz ]]; then
            rm -r ${workdir}/${subj}${sessionpath}fmap/synb0/
        else
            echo -e "${RED}something went wrong with synb0${NC}"
            exit 1
        fi
    fi

    samedir=${dwidir}
    cd ${workdir}/${subj}${sessionpath}fmap
    rsync -a ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz \
        ${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv \
        ${outputdir}/dwi-preproc/${subj}${sessionpath}fmap
fi

#----------------------------------------------------------------------
# topup
#----------------------------------------------------------------------
if [[ ! -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-unwarped_epi.nii.gz ]] ||
    [[ ! -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-topup_fieldcoef.nii.gz ]]; then
    cd ${workdir}/${subj}${sessionpath}fmap

    dim1=$(fslinfo ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz | grep -w dim1 | awk '{ print $2 }' | awk '{print int($0)}')
    dim2=$(fslinfo ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz | grep -w dim2 | awk '{ print $2 }' | awk '{print int($0)}')
    dim3=$(fslinfo ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz | grep -w dim3 | awk '{ print $2 }' | awk '{print int($0)}')

    if ((dim1 % 4 == 0 && dim2 % 4 == 0 && dim3 % 4 == 0)); then
        log "${BLUE}" "All dimensions are integer multiples of 4; using b02b0_4.cnf for topup"
        configfile=b02b0_4.cnf
    elif ((dim1 % 2 == 0 && dim2 % 2 == 0 && dim3 % 2 == 0)); then
        log "${BLUE}" "All dimensions are integer multiples of 2; using b02b0_2.cnf for topup"
        configfile=b02b0_2.cnf
    else
        log "${BLUE}" "At least one dimension is odd; using b02b0_1.cnf as config file for topup"
        configfile=b02b0_1.cnf
    fi

    echo
    echo -e "${BLUE}running topup${NC}"
    echo

    topup --imain=${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz \
        --datain=${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv \
        --config=${configfile} \
        --out=${subj}${sessionfile}space-dwi_desc-topup \
        --iout=${subj}${sessionfile}space-dwi_desc-unwarped_epi \
        --fout=${subj}${sessionfile}space-dwi_desc-topup_fieldmap --verbose >${subj}${sessionfile}topup_$(date +"%Y-%m-%d_%H-%M").log

    cp ${subj}${sessionfile}topup_*.log ${outputdir}/dwi-preproc/${subj}/log
fi

###################
### round bvals ###
###################
cd ${workdir}/${subj}${sessionpath}dwi
rsync -av ${bidsdir}/${subj}${sessionpath}/dwi/{${subj}${sessionfile}dwi.bv*,${subj}${sessionfile}dwi.json} ${workdir}/${subj}${sessionpath}dwi/
chmod u+rw ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.bval
${scriptdir}/round_bvals.py ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.bval

#######################
## create brain mask ##
#######################
# mean of unwarped image to allow registration
mrmath ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-unwarped_epi.nii.gz mean \
    ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz -axis 3 -force

# Get the mean b-zero (un-corrected)
dwiextract -nthreads ${nthreads} \
    ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz - -bzero \
    -fslgrad ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}*dwi.bvec ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.bval |
    mrmath - mean ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-meanb0-uncorrected_dwi.nii.gz -axis 3 -force

if [[ ! -f ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz ]]; then
    # rigid registration of nodif_dwi to b0
    antsRegistrationSyN.sh -d 3 -m ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz \
        -f "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-meanb0-uncorrected_dwi.nii.gz" \
        -o ${subj}${sessionfile}rigidreg -t r -n ${nthreads} -p d
    mv ${subj}${sessionfile}rigidregWarped.nii.gz ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz
    rm *rigidreg*
fi

if [[ ! -f ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz ]]; then
    mri_synthstrip \
        -i ${subj}${sessionfile}space-dwi_desc-nodif_dwi.nii.gz \
        -o ${subj}${sessionfile}space-dwi_desc-nodif-brain_dwi.nii.gz \
        --mask ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz
fi

cd ${workdir}/${subj}${sessionpath}fmap
rsync -av ${subj}${sessionfile}space-dwi_desc-unwarped_epi* \
    ${subj}${sessionfile}space-dwi_desc-topup_fieldmap* \
    ${subj}${sessionfile}space-dwi_desc-topup* \
    ${outputdir}/dwi-preproc/${subj}${sessionpath}fmap

# for QC
rsync -av ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-noise_dwi.nii.gz \
    ${outputdir}/dwi-preproc/${subj}${sessionpath}qc

echo "Preprocessing complete."
exit 0
