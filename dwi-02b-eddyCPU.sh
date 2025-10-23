#!/bin/bash

#SBATCH --job-name=eddy
#SBATCH --mem=4G
#SBATCH --partition=luna-cpu-short
#SBATCH --cpus-per-task=4
#SBATCH --time=00-8:00:00
#SBATCH --nice=2000
#SBATCH --qos=anw-cpu
#SBATCH --output %x_%A.log

# Written by C. Vriend - AmsUMC Jun 2024
# c.vriend@amsterdamumc.nl

# usage instructions
Usage() {
    cat <<EOF

    (C) C.Vriend - 9/7/2025 - dwi-02b-eddy.sh
   
    Usage: ./dwi-02b-eddy.sh <5 inputs>
  

EOF
    exit 1
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
method=""
# input variables

# Parse command line arguments
while getopts ":i:o:w:s:m:" opt; do
    case $opt in
        i) bidsdir="$OPTARG" ;;
        o) outputdir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        m) method="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done
missing=0
for var in bidsdir outputdir workdir subj method; do
    if [[ -z "${!var}" ]]; then
        echo "Error: -${var:0:1} ($var) is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    Usage
fi

echo
echo -e "${BLUE}chosen method for eddy: ${method}${NC}"
echo

for dwidir in ${bidsdir}/${subj}/{,ses*/}dwi; do
    if [ ! -d ${dwidir} ]; then
        continue
    fi
    sessiondir=$(dirname ${dwidir})

    # if [[ $(ls ${sessiondir}/dwi/*dwi.nii.gz | wc -l) -gt 1 ]]; then
    #     echo -e "${RED}ERROR! this script cannot handle >1 dwi scan per session${NC}"
    #     echo -e "${RED}exiting script${NC}"
    #     exit
    # fi

    session=$(echo "${sessiondir}" | grep -oP "(?<=${subj}/).*")
    if [ -z ${session} ]; then
        sessionpath=/
        sessionfile=_
    else
        sessionpath=/${session}/
        sessionfile=_${session}_

    fi
    echo -e ${YELLOW}----------------------${NC}
    echo -e ${YELLOW}running EDDY on dwi data${NC}
    echo -e ${YELLOW}${subj}${NC}
    echo -e ${YELLOW}${session}${NC}
    echo -e ${YELLOW}----------------------${NC}

    # inputs
    dwiworkdir=${workdir}/${subj}${sessionpath}dwi
    DWImain=${dwiworkdir}/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz
    DWImask=${dwiworkdir}/${subj}${sessionfile}space-dwi_desc-brain-uncorrected_mask.nii.gz
    DWIacqp=${dwiworkdir}/${subj}${sessionfile}acq-dwi_desc-acqparams.tsv
    DWIbvecs=${dwidir}/${subj}${sessionfile}dwi.bvec
    DWIbvals=${dwidir}/${subj}${sessionfile}dwi.bval
    DWIjson=${dwidir}/${subj}${sessionfile}dwi.json
    topup=${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-topup
    DWIout=${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc

    # choose method of eddy correction (default or volcorr)

    basedir=$(dirname DWImain)
    cd ${basedir}

    # create index.txt file
    idx=$(fslnvols ${DWImain})
    printf '1 %.0s' $(seq 1 "$idx") >${basedir}/index.txt

        # json available with slice-timing?
    if jq -e '.SliceTiming' "${DWIjson}" >/dev/null; then
        STavail=1
    else
        STavail=0
    fi

    # default
    if [[ ${method} == "default" ]]; then
    # adding json input will let the cmd fail if json does not contain ST info
        eddy \
            --imain=${DWImain} \
            --mask=${DWImask} \
            --acqp=${DWIacqp} \
            --index=index.txt \
            --bvecs=${DWIbvecs} \
            --bvals=${DWIbvals} \
            --out=${DWIout} \
            --topup=${topup} \
            --repol --cnr_maps \
            --slm=linear \
            --estimate_move_by_susceptibility --verbose

        # run QC
        echo
        echo "running QC"
        eddy_quad ${DWIout} \
            -idx index.txt \
            -par ${DWIacqp} \
            -m ${DWImask} \
            -b ${DWIbvals} \
            -f ${topup}_fieldmap.nii.gz

    elif [[ ${method} == "volcorr" ]]; then

        if ((STavail == 1)); then
            # w/ slice-to-vol correction
            eddy \
                --imain=${DWImain} \
                --mask=${DWImask} \
                --acqp=${DWIacqp} \
                --index=index.txt \
                --json=${DWIjson} \
                --bvecs=${DWIbvecs} \
                --bvals=${DWIbvals} \
                --out=${DWIout} \
                --topup=${topup} \
                --repol --cnr_maps \
                --slm=linear \
                --estimate_move_by_susceptibility --verbose \
                --mbs_niter=10 --mbs_lambda=10 --mbs_ksp=10 \
                --niter=6 --fwhm=15,10,4,2,0,0 \
                --mporder=8 --s2v_niter=8 --json=${DWIjson} \
                --s2v_lambda=1 --s2v_interp=trilinear >${basedir}/eddy.log

            echo
            echo "running QC"
            eddy_quad ${DWIout} \
                -idx index.txt \
                -par ${DWIacqp} \
                -m ${DWImask} \
                -b ${DWIbvals} \
                -f ${topup}_fieldmap.nii.gz \
        -g ${DWIout}.eddy_rotated_bvecs \
        -j ${DWIjson} \
        -v

        else 
        echo -e "${RED}Slice to volume correction not possible without SliceTime information in json${NC}"
        continue
        fi


    elif [[ ${method} == "volcorrnosdc" ]]; then

        if ((STavail == 1)); then
            # w/ slice-to-vol correction nomove by suscep.
            eddy \
                --imain=${DWImain} \
                --mask=${DWImask} \
                --acqp=${DWIacqp} \
                --index=index.txt \
                --json=${DWIjson} \
                --bvecs=${DWIbvecs} \
                --bvals=${DWIbvals} \
                --out=${DWIout} \
                --topup=${topup} \
                --repol --cnr_maps \
                --slm=linear \
                 --verbose \
                --mbs_niter=10 --mbs_lambda=10 --mbs_ksp=10 \
                --niter=6 --fwhm=15,10,4,2,0,0 \
                --mporder=8 --s2v_niter=8 --json=${DWIjson} \
                --s2v_lambda=1 --s2v_interp=trilinear >${basedir}/eddy.log

  	    echo
        echo "running QC"
        eddy_quad ${DWIout} \
            -idx index.txt \
            -par ${DWIacqp} \
            -m ${DWImask} \
            -b ${DWIbvals} \
            -f ${topup}_fieldmap.nii.gz \
	-g ${DWIout}.eddy_rotated_bvecs \
	-j ${DWIjson} \
	-v
        else 
        echo -e "${RED}Slice to volume correction not possible without SliceTime information in json${NC}"
        continue
        fi


    elif [[ ${method} == "nofmap" ]]; then

        eddy \
            --imain=${DWImain} \
            --mask=${DWImask} \
            --acqp=${DWIacqp} \
            --index=${basedir}/index.txt \
            --bvecs=${DWIbvecs} \
            --bvals=${DWIbvals} \
            --out=${DWIout} \
            --repol --cnr_maps \
            --slm=linear \
            --verbose >${basedir}/eddy.log

        # run QC
        echo
        echo "running QC"
        echo
        eddy_quad ${DWIout} \
            -idx ${basedir}/index.txt \
            -par ${DWIacqp} \
            -m ${DWImask} \
            -b ${DWIbvals}


    else

        echo "proper method for eddy not set"
        echo "exiting script"
        exit
    fi

    cp eddy_*.log ${outputdir}/dwi-preproc/${subj}${sessionpath}logs/${subj}${sessionfile}eddy.log

    # rename output
    cd ${workdir}/${subj}${sessionpath}dwi

    cp ${subj}${sessionfile}space-dwi_desc-preproc.eddy_rotated_bvecs \
        ${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec
    cp ${subj}${sessionfile}space-dwi_desc-preproc.nii.gz \
        ${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz
    cp ${subj}${sessionfile}space-dwi_desc-preproc.eddy_cnr_maps.nii.gz \
        ${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz
    cp ${DWIbvals} \
        ${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval
    mv *.qc eddyqc

    cd ${workdir}/${subj}${sessionpath}dwi

    rsync -av ${subj}${sessionfile}space-dwi*_dwi.* ${subj}${sessionfile}space-dwi_desc-brain-uncorrected_mask.nii.gz eddyqc \
        ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi

    # clean-up
    if [ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz ] &&
        [ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec ] &&
        [ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz ]; then
        rm -r ${workdir}/${subj}${sessionpath}
        rm ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/*meanb0* \
        ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/*dns+degibbs*

        echo
        echo -e ${GREEN}FINISHED preprocessing ${subj}${sessionpath}${NC}
        echo

    else
        echo -e "${RED}ERROR! not all output was created successfully${NC}"

    fi

done
