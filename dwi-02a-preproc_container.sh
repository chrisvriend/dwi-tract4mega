#!/bin/bash

###############################################################################
# preproc.sh
# Author: C. Vriend - AUMC
# Date: Nov 05 2025
# Description: prepare dwi scan for eddy and perform topup
###############################################################################
set -euo pipefail


usage() {
    echo "Usage: $0 -i <bidsdir> -o <outputdir> -w <workdir> -s <subj> -c <scriptdir>"
    exit 1
}

threads=8


# Define color variables
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[34m'
NC='\033[0m' # No Color


# initialize function
# Function to get the "opposite" of a PE direction
get_opposite_PE() {
    local PEdir="$1"
    if [[ "$PEdir" == *"-" ]]; then
        echo "${PEdir/-/}"
    else
        echo "${PEdir}-"
    fi
}

# Initialize variables
bidsdir=""
outputdir=""
workdir=""
subj=""
session=""
scriptdir=""
# input variables

# Parse command line arguments
while getopts ":i:o:w:s:c:z:" opt; do
    case $opt in
        i) bidsdir="$OPTARG" ;;
        o) outputdir="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        z) session="$OPTARG" ;;
        c) scriptdir="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done
missing=0
for var in bidsdir outputdir workdir subj scriptdir; do
    if [[ -z "${!var}" ]]; then
        echo "Error: -${var:0:1} ($var) is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    echo "Usage: $0 -i <bidsdir> -o <outputdir> -w <workdir>  -c <scriptdir> -s <subj> [-z <session>]"
    exit 1
fi

# Check if directories exist
for dir in "$bidsdir" "$scriptdir"; do
    if [[ ! -d "$dir" ]]; then
        echo "Error: Directory $dir does not exist."
        exit 1
    fi
done

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

    
    if [[ ! -f ${bidsdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.nii.gz || ! -f ${bidsdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.bvec ]]; then
        echo -e "${RED}no dwi scan/bvec found for ${subj} - ${session}${NC}"
        exit 1
      
    fi
    
    echo -e ${YELLOW}----------------------${NC}
    echo -e ${YELLOW}Preprocessing dwi data${NC}
    echo -e ${YELLOW}${subj}${NC}
    echo -e ${YELLOW}${session}${NC}
    echo -e ${YELLOW}----------------------${NC}

    fmap_samePE=()
    fmap_otherPE=()

    if [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz ]] &&
    [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval ]] &&
    [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec ]]; then
        echo -e "${GREEN}${subj}${sessionfile} already preprocessed with eddy${NC}"
        exit 1
    fi
    
    
    mkdir -p "${workdir}/${subj}${sessionpath}dwi"
    mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi"
    mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}logs"
    mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}fmap"
    mkdir -p "${outputdir}/dwi-preproc/${subj}${sessionpath}figures"
    
    # Specify the path to the DWI JSON sidecar
    dwi_json_path=$(ls ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}dwi.json)
    # extract TotalReadoutTime
    dwi_trt=$(cat ${dwi_json_path} | jq -r '.TotalReadoutTime')
    dwi_PE=$(cat ${dwi_json_path} | jq -r '.PhaseEncodingDirection')
    
    if [ -z ${dwi_trt} ] || [ -z ${dwi_PE} ]; then
    echo -e "${RED}no TotalReadOutTime or PhaseEncodingDirection found in dwi json file${NC}"
    echo 
    continue
    fi

    # determine settings for topup
    if [ ${dwi_PE} == "i" ];then PE_dwi_FSL="1 0 0"
    elif [ ${dwi_PE} == "i-" ];then PE_dwi_FSL="-1 0 0"
    elif [ ${dwi_PE} == "j" ];then PE_dwi_FSL="0 1 0"
    elif [ ${dwi_PE} == "j-" ];then PE_dwi_FSL="0 -1 0"
    elif [ ${dwi_PE} == "k" ];then PE_dwi_FSL="0 0 1"
    elif [ ${dwi_PE} == "k-" ];then PE_dwi_FSL="0 0 -1"
    fi
    

    #----------------------------------------------------------------------
    #                           MP-PCA denoising & deringing of dwi scan
    #----------------------------------------------------------------------
    if [ ! -f ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz ]; then
        dwidenoise ${bidsdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.nii.gz \
        ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns_dwi.mif \
        -nthreads ${threads} -force
        #Remove Gibbs Ringing Artifacts
        mrdegibbs ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns_dwi.mif \
        ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz \
        -nthreads ${threads} -force
        rm ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns_dwi.mif
    fi

    # write dwi acqparams
    echo "${PE_dwi_FSL} ${dwi_trt}" >${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}acq-dwi_desc-acqparams.tsv
    
    #----------------------------------------------------------------------
    #                           Check for fieldmaps
    #----------------------------------------------------------------------
    
    # Specify the path to the fieldmap folder
    fieldmap_folder=${bidsdir}/${subj}${sessionpath}/fmap

    # Check if the fieldmap folder exists
    if [ -d ${fieldmap_folder} ]; then

            echo -e "${YELLOW} fieldmap folder found${NC}"

            mkdir -p ${workdir}/${subj}${sessionpath}fmap

            # if [[ $(ls ${fieldmap_folder}/*acq-dwi*dir*epi.json | wc -l) -gt 2 ]]; then
            #     echo -e "${RED}ERROR! more than 2 fieldmaps for acq-dwi found in ${fieldmap_folder}${NC}"
            #     exit 0
            # fi


          
            for fmap_json in ${fieldmap_folder}/*acq-dwi*dir*epi.json; do
                    if [[ $(jq -r '.IntendedFor' "${fmap_json}") == *dwi* ]]; then
                        fmap_nii=${fmap_json%%.json}.nii.gz
                        fmap_PE=$(cat ${fmap_json} | jq -r '.PhaseEncodingDirection')
                        fmap_trt=$(jq -r '.TotalReadoutTime' "$fmap_json")
                    echo "${fmap_json}"
                    echo -e "${BLUE}PhaseEncodingDirection: $fmap_PE${NC}"
                    echo -e "${BLUE}TotalReadoutTime: $fmap_trt${NC}"
                    echo
                        if [[ -f ${fmap_nii} ]];then
                                if [ "${fmap_PE}" == "${dwi_PE}" ];then
                                        fmap_samePE+=("${fmap_nii}")
                                else
                                        fmap_otherPE+=("${fmap_nii}")
                                fi
                        fi
                fi
            done

        for fmap in "${fmap_samePE}" "${fmap_otherPE}"; do

            if [ ! -z ${fmap} ]; then 

                fmap_json=${fmap%%.nii.gz}.json
                fmap_PE=$(cat ${fmap_json} | jq -r '.PhaseEncodingDirection')

                if [[ "$fmap_PE" == "j" ]]; then
                    dir=AP
                elif [[ "$fmap_PE" == "j-" ]]; then
                    dir=PA
                elif [[ "$fmap_PE" == "i" ]]; then
                    dir=LR
                elif [[ "$fmap_PE" == "i-" ]]; then
                    dir=RL
                elif [[ "$fmap_PE" == "k" ]]; then
                    dir=IS
                elif [[ "$fmap_PE" == "k-" ]]; then
                    dir=SI
                else
                    echo "Unknown Phase Encoding Direction"
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
    #                           One Fieldmap available 
    #----------------------------------------------------------------------
    # in case only fmap with opposite but not same PE available 
    #if [ -z ${fmap_samePE} ] && [ ! -z ${fmap_otherPE} ]; then
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
            continue
        fi

            
        # options:
        # --> 1) extract mean b0 from dwi and merge with fieldmap but rigid reg with interpolation necessary
        # 2) take first b0 of dwi map as fieldmap is acquired right before dwi (on Siemens Vida); then uneven number of nvols
            
        # determine letter for opposite PE

                if [[ "$dwi_PE" == "j" ]]; then
                    dwidir=AP
                elif [[ "$dwi_PE" == "j-" ]]; then
                    dwidir=PA
                elif [[ "$dwi_PE" == "i" ]]; then
                    dwidir=LR
                elif [[ "$dwi_PE" == "i-" ]]; then
                    dwidir=RL
                elif [[ "$dwi_PE" == "k" ]]; then
                    dwidir=IS
                elif [[ "$dwi_PE" == "k-" ]]; then
                    dwidir=SI
                else
                    echo "Unknown Phase Encoding Direction"
                fi
                if [[ "$PE_other" == "j" ]]; then
                    otherdir=AP
                elif [[ "$PE_other" == "j-" ]]; then
                    otherdir=PA
                elif [[ "$PE_other" == "i" ]]; then
                    otherdir=LR
                elif [[ "$PE_other" == "i-" ]]; then
                    otherdir=RL
                elif [[ "$PE_other" == "k" ]]; then
                    otherdir=IS
                elif [[ "$PE_other" == "k-" ]]; then
                    otherdir=SI
                else
                    echo "Unknown Phase Encoding Direction"
                fi


        # extract mean b0 from dwi
        dwiextract -nthreads ${threads} \
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
        
        # rigid registration of dwi mean b0 (AP) and PA fieldmap
        antsRegistrationSyN.sh -d 3 -m ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-temp_epi.nii.gz \
        -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-temp_epi.nii.gz \
        -o ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}rigidreg -t r -n ${threads} -p d
        # apply to multi-volume PA fieldmap
        antsApplyTransforms -d 3 -e 3 -i ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz \
        -r ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-temp_epi.nii.gz \
        -t ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}rigidreg0GenericAffine.mat \
        -o ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-warped-degibbs_epi.nii.gz -v -u int
        rm -f ${workdir}/${subj}${sessionpath}fmap/*rigidreg*
        fslmerge -t ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-4topup_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-temp_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-warped-degibbs_epi.nii.gz
        
        rm ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-warped-degibbs_epi.nii.gz \
        ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}*temp*.nii.gz
        


        if [ ${PE_other} == "i" ];then PE_other_FSL="1 0 0"
        elif [ ${PE_other} == "i-" ];then PE_other_FSL="-1 0 0"
        elif [ ${PE_other} == "j" ];then PE_other_FSL="0 1 0"
        elif [ ${PE_other} == "j-" ];then PE_other_FSL="0 -1 0"
        elif [ ${PE_other} == "k" ];then PE_other_FSL="0 0 1"
        elif [ ${PE_other} == "k-" ];then PE_other_FSL="0 0 -1"
        fi

                # write TRT to refparams file
                cd ${workdir}/${subj}${sessionpath}fmap
                echo "${PE_dwi_FSL} ${dwi_trt}" >${subj}${sessionfile}dir-${dwidir}_desc-refparams.tsv
                for ((i = 0; i < $(fslnvols ${subj}${sessionfile}dir-${otherdir}_space-dwi_desc-degibbs_epi.nii.gz); i++)); do
                    echo "${PE_other_FSL} ${other_trt}" >>"${subj}${sessionfile}dir-${dwidir}_desc-refparams.tsv"
                done

        # set dwidir to samedir for later steps
        samedir=${bidsdir}/${subj}${sessionpath}

    #elif [ ! -z ${fmap_samePE} ] && [ ! -z ${fmap_otherPE} ]; then
    elif [ ${#fmap_samePE[@]} -ne 0 ] && [ ${#fmap_otherPE[@]} -ne 0 ]; then
        #----------------------------------------------------------------------
        #                           two Fieldmaps available 
        #----------------------------------------------------------------------

        fmap_otherjson=${fmap_otherPE%%.nii.gz}.json
        fmap_samejson=${fmap_samePE%%.nii.gz}.json

        PE_other=$(cat ${fmap_otherjson} | jq -r '.PhaseEncodingDirection')
        PE_same=$(cat ${fmap_samejson} | jq -r '.PhaseEncodingDirection')
        other_trt=$(cat ${fmap_otherjson} | jq -r '.TotalReadoutTime')
        same_trt=$(cat ${fmap_samejson} | jq -r '.TotalReadoutTime')



        # Get the opposite of the first direction
        opposite_pe1=$(get_opposite_PE "$PE_other")

        # Check if the second direction matches the opposite of the first
        if [[ "$PE_same" == "$opposite_pe1" ]]; then
            echo -e "${GREEN}The fmap PE directions are opposites.${NC}"
        else
            echo -e "${RED}The fmap PE directions are NOT opposites.${NC}"
            continue
        fi

        # check that the PE directions are consistent with dwi scan
        
        # Check if the second direction matches the opposite of the first
        if [[ "$dwi_PE" == "$PE_same" ]]; then
            echo -e "${GREEN} PE directions of dwi and 'same' fmap are consistent.${NC}"
        else
            echo -e "${RED}PE directions of dwi and fmap are NOT opposites.${NC}"
            continue
        fi

        opposite_pe1=$(get_opposite_PE "$PE_other")

        # Check if the second direction matches the opposite of the first
        if [[ "$dwi_PE" == "$opposite_pe1" ]]; then
            echo -e "${GREEN} PE directions of dwi and 'opposite' fmap are consistent.${NC}"
        else
            echo -e "${RED}PE directions of dwi and fmap are NOT opposites.${NC}"
            continue
        fi

        # get directions 
                if [[ "$PE_same" == "j" ]]; then
                    samedir=AP
                elif [[ "$PE_same" == "j-" ]]; then
                    samedir=PA
                elif [[ "$PE_same" == "i" ]]; then
                    samedir=LR
                elif [[ "$PE_same" == "i-" ]]; then
                    samedir=RL
                elif [[ "$PE_same" == "k" ]]; then
                    samedir=IS
                elif [[ "$PE_same" == "k-" ]]; then
                    samedir=SI
                else
                    echo "Unknown Phase Encoding Direction"
                fi
                if [[ "$PE_other" == "j" ]]; then
                    otherdir=AP
                elif [[ "$PE_other" == "j-" ]]; then
                    otherdir=PA
                elif [[ "$PE_other" == "i" ]]; then
                    otherdir=LR
                elif [[ "$PE_other" == "i-" ]]; then
                    otherdir=RL
                elif [[ "$PE_other" == "k" ]]; then
                    otherdir=IS
                elif [[ "$PE_other" == "k-" ]]; then
                    otherdir=SI
                else
                    echo "Unknown Phase Encoding Direction"
                fi


            # determine settings for topup
            if [ ${PE_same} == "i" ];then PE_same_FSL="1 0 0"
            elif [ ${PE_same} == "i-" ];then PE_same_FSL="-1 0 0"
            elif [ ${PE_same} == "j" ];then PE_same_FSL="0 1 0"
            elif [ ${PE_same} == "j-" ];then PE_same_FSL="0 -1 0"
            elif [ ${PE_same} == "k" ];then PE_same_FSL="0 0 1"
            elif [ ${PE_same} == "k-" ];then PE_same_FSL="0 0 -1"
            fi
            if [ ${PE_other} == "i" ];then PE_other_FSL="1 0 0"
            elif [ ${PE_other} == "i-" ];then PE_other_FSL="-1 0 0"
            elif [ ${PE_other} == "j" ];then PE_other_FSL="0 1 0"
            elif [ ${PE_other} == "j-" ];then PE_other_FSL="0 -1 0"
            elif [ ${PE_other} == "k" ];then PE_other_FSL="0 0 1"
            elif [ ${PE_other} == "k-" ];then PE_other_FSL="0 0 -1"
            fi
                
                # merge blip up/down scans for topup (if available)
                echo -e "${BLUE}merge blip up/down scans for topup${NC}"

                cd ${workdir}/${subj}${sessionpath}fmap

                rm -f ${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv

                # same PE
                for ((i = 0; i < $(fslnvols ${subj}${sessionfile}dir-${samedir}_space-dwi_desc-degibbs_epi.nii.gz); i++)); do
                    echo "${PE_same_FSL} ${same_trt}" >>"${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv"
                done
                # other PE
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
   # elif [ -z ${fmap_samePE} ] && [ -z ${fmap_otherPE} ]; then
        #----------------------------------------------------------------------
        #                           Syn b0 (no fieldmaps available) 
        #----------------------------------------------------------------------
            
    
            echo -e "${BLUE}no fmaps found - creating syn b0 for topup${NC}"
            
            mkdir -p "${workdir}/${subj}${sessionpath}fmap/"

            

                if [[ "$dwi_PE" == "j" ]]; then
                    dwidir=AP
                elif [[ "$dwi_PE" == "j-" ]]; then
                    dwidir=PA
                elif [[ "$dwi_PE" == "i" ]]; then
                    dwidir=LR
                elif [[ "$dwi_PE" == "i-" ]]; then
                    dwidir=RL
                elif [[ "$dwi_PE" == "k" ]]; then
                    dwidir=IS
                elif [[ "$dwi_PE" == "k-" ]]; then
                    dwidir=SI
                else
                    echo "Unknown Phase Encoding Direction"
                fi

                PE_other=$(get_opposite_PE "$dwi_PE")

                if [[ "$PE_other" == "j" ]]; then
                    otherdir=AP
                elif [[ "$PE_other" == "j-" ]]; then
                    otherdir=PA
                elif [[ "$PE_other" == "i" ]]; then
                    otherdir=LR
                elif [[ "$PE_other" == "i-" ]]; then
                    otherdir=RL
                elif [[ "$PE_other" == "k" ]]; then
                    otherdir=IS
                elif [[ "$PE_other" == "k-" ]]; then
                    otherdir=SI
                else
                    echo "Unknown Phase Encoding Direction"
                fi


            # write TRT to refparams file
            echo "${PE_dwi_FSL} ${dwi_trt}" > ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv
            echo "${PE_dwi_FSL} 0.00" >> ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv
                

            if [[ ! -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz  ]]; then

             
                mkdir -p "${workdir}/${subj}${sessionpath}fmap/synb0/tmp" \
                "${workdir}/${subj}${sessionpath}fmap/synb0/input" \
                "${workdir}/${subj}${sessionpath}fmap/synb0/output"

                cp ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv \
                    ${workdir}/${subj}${sessionpath}fmap/synb0/input/ 

                # extract first b0 vol from dwi
                dwiextract -nthreads ${threads} \
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
                
   
            
                cd  ${workdir}/${subj}${sessionpath}fmap/synb0/input
                if [ -L T1.nii.gz ]; then 
                    unlink T1.nii.gz
                    unlink BRAIN.nii.gz
                    unlink acqparams.txt
                    unlink b0.nii.gz
                fi
                ln -s  ${subj}${sessionfile}desc-n4_T1w.nii.gz T1.nii.gz
                ln -s  ${subj}${sessionfile}desc-brain_T1w.nii.gz BRAIN.nii.gz
                ln -s  ${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv acqparams.txt
                ln -s  ${subj}${sessionfile}dir-${dwidir}_space-dwi_desc-b0_epi.nii.gz b0.nii.gz
                
                
                #Run Synb0-DISCO for fieldmap-free distortion correction
                if [[ ! -f ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_d_smooth.nii.gz ]] || \
                [[ ! -f ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_u.nii.gz ]]; then
                    echo " test" 
                    synb0-disco --input ${workdir}/${subj}${sessionpath}fmap/synb0/input \
                    --output ${workdir}/${subj}${sessionpath}fmap/synb0/output --notopup
                     
                fi
                
                fslmerge -t ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_all.nii.gz \
                ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_d_smooth.nii.gz \
                ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_u.nii.gz &&
                
                mv ${workdir}/${subj}${sessionpath}fmap/synb0/output/b0_all.nii.gz \
                ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz
                mv ${workdir}/${subj}${sessionpath}fmap/synb0/input/${subj}${sessionfile}dir-${dwidir}${otherdir}_desc-refparams.tsv \
                ${workdir}/${subj}${sessionpath}fmap/
                

                if [[ -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}dir-${dwidir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz ]]; then
                    #clean-up
                    rm -r ${workdir}/${subj}${sessionpath}fmap/synb0/
                else
                    echo -e "${RED}something went wrong with synb0${NC}"
                    exit 1
                fi
            fi 

            samedir=${dwidir}
            
            cd  ${workdir}/${subj}${sessionpath}fmap
              rsync -a ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz \
                ${subj}${sessionfile}dir-${samedir}${otherdir}_desc-refparams.tsv \
                ${outputdir}/dwi-preproc/${subj}${sessionpath}fmap

    fi
        
    #----------------------------------------------------------------------
    #                           topup
    #----------------------------------------------------------------------
    
        if [[ ! -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-unwarped_epi.nii.gz ]] ||
        [[ ! -f ${workdir}/${subj}${sessionpath}fmap/${subj}${sessionfile}space-dwi_desc-topup_fieldcoef.nii.gz ]]; then
            
            cd ${workdir}/${subj}${sessionpath}fmap
            #    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;6c4c9591.2002
            #      b02b0_4.cnf  -- Recommended when the data matrix is an integer multiple of 4 in all direction
            #      b02b0_2.cnf  -- Recommended when the data matrix is an integer multiple of 2 in all direction
            #      b02b0_1.cnf  -- Recommended when the data matrix is odd in one or more directions
            
            dim3=$(fslinfo ${subj}${sessionfile}dir-${samedir}${otherdir}_space-dwi_desc-4topup_epi.nii.gz | grep -w dim3 | awk '{ print $2 }' | awk '{print int($0)}')
            if ((dim3 % 4 == 0)); then
                echo -e "${BLUE}slices are integer multiple of 4; using b02b0_4.cnf for topup${NC}"
                configfile=b02b0_4.cnf
                elif ((dim3 % 2 == 0)); then
                echo -e "${BLUE}slices are integer multiple of 2; using b02b0_2.cnf for topup${NC}"
                configfile=b02b0_2.cnf
            else
                echo -e "${BLUE}odd number of slices; using b02b0_1.cnf as config file for topup${NC}"
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
            --fout=${subj}${sessionfile}space-dwi_desc-topup_fieldmap --verbose >${subj}${sessionfile}topup.log
            cp ${subj}${sessionfile}topup.log ${outputdir}/dwi-preproc/${subj}${sessionpath}logs
            
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
        ${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz -axis 3 -force
        
        # Get the mean b-zero (un-corrected)
        dwiextract -nthreads ${threads} \
        ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-dns+degibbs_dwi.nii.gz - -bzero \
        -fslgrad ${bidsdir}/${subj}${sessionpath}/dwi/${subj}${sessionfile}*dwi.bvec ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}dwi.bval |
        mrmath - mean ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-meanb0-uncorrected_dwi.nii.gz -axis 3 -force
        
        if [[ ! -f ${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz ]]; then
            # rigid registration of nodif_epi to b0
            antsRegistrationSyN.sh -d 3 -m ${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz \
            -f "${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-meanb0-uncorrected_dwi.nii.gz" \
            -o ${subj}${sessionfile}rigidreg -t r -n ${threads} -p d
            mv ${subj}${sessionfile}rigidregWarped.nii.gz ${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz
            rm *rigidreg*
        fi
        if [[ ! -f ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz ]]; then
            mri_synthstrip \
            -i ${subj}${sessionfile}space-dwi_desc-nodif_epi.nii.gz \
            --mask ${workdir}/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-brain_mask.nii.gz
        fi
    
        cd ${workdir}/${subj}${sessionpath}fmap

        rsync -av ${subj}${sessionfile}space-dwi_desc-unwarped_epi* \
        ${subj}${sessionfile}space-dwi_desc-topup_fieldmap* \
        ${subj}${sessionfile}space-dwi_desc-topup* \
        ${outputdir}/dwi-preproc/${subj}${sessionpath}fmap
        
echo "Preprocessing complete."
exit 0

