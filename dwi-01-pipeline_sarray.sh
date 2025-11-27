#!/bin/bash

# (C) C. Vriend - Aumc ANW/Psy - June '23
# c.vriend@amsterdamumc.nl


#################
# tree of outputs
#################
# derivatives
# │
# ├── dwi-connectome
# │   └── sub-[subjID]_[sessionID]
# │       ├── conn
# │       │   ├── sub-[subjID]_[sessionID]_atlas-300P7N_desc-FA_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-300P7N_desc-lengths_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-300P7N_desc-ndi_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-300P7N_desc-streams_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-300P7N_trackassign.txt
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P17N_desc-FA_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P17N_desc-lengths_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P17N_desc-ndi_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P17N_desc-streams_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P17N_trackassign.txt
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P7N_desc-FA_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P7N_desc-lengths_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P7N_desc-ndi_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P7N_desc-streams_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-400P7N_trackassign.txt
# │       │   ├── sub-[subjID]_[sessionID]_atlas-aparc500_desc-FA_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-aparc500_desc-lengths_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-aparc500_desc-ndi_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-aparc500_desc-streams_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-aparc500_trackassign.txt
# │       │   ├── sub-[subjID]_[sessionID]_atlas-BNA_desc-FA_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-BNA_desc-lengths_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-BNA_desc-ndi_connmatrix.csv
# │       │   ├── sub-[subjID]_[sessionID]_atlas-BNA_desc-streams_connmatrix.csv
# │       │   └── sub-[subjID]_[sessionID]_atlas-BNA_trackassign.txt
# │       ├── dwi
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-sift-50M_stats.csv
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_mu.txt
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_tracto-50M_desc-sift_weights.txt
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_tracto-50M-sift_dwi.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_tracto-50M.tck
# │       │   └── sub-[subjID]_[sessionID]_tissue-WM-norm_fod.nii.gz
# │       ├── figures
# │       │   ├── sub-[subjID]_[sessionID]_siftoverlay3D.png
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_tissue-overlay.png
# │       │   └── sub-[subjID]_[sessionID]_space-dwi_tractoverlay.png
# │       ├── logs
# │       │   └── sub-[subjID]_[sessionID]_dwi-tckconn.log
# │       └── rpf
# │           ├── sub-[subjID]_[sessionID]_space-dwi_tissue-CSF_response.txt
# │           ├── sub-[subjID]_[sessionID]_space-dwi_tissue-GM_response.txt
# │           └── sub-[subjID]_[sessionID]_space-dwi_tissue-WM_response.txt
# ├── dwi-preproc
# │   └── sub-[subjID]_[sessionID]
# │       ├── anat
# │       │   ├── sub-[subjID]_[sessionID]_desc-5tt-hsvs_probseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_desc-brain_mask.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_desc-brain_T1w.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_desc-gmwm_probseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_desc-preproc_T1w.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_desc-wm_probseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-300P7N_dseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-300P7N_roivols.txt
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-400P17N_dseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-400P17N_roivols.txt
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-400P7N_dseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-400P7N_roivols.txt
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-aparc500_dseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-aparc500_roivols.txt
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_atlas-BNA_dseg.nii.gz
# │       │   └── sub-[subjID]_[sessionID]_space-dwi_atlas-BNA_roivols.txt
# │       ├── dwi
# │       │   ├── eddyqc
# │       │   │   ├── avg_b0_pe0.png
# │       │   │   ├── avg_b0.png
# │       │   │   ├── avg_b1000.png
# │       │   │   ├── avg_b2000.png
# │       │   │   ├── avg_b3000.png
# │       │   │   ├── cnr0000.nii.gz.png
# │       │   │   ├── cnr0001.nii.gz.png
# │       │   │   ├── cnr0002.nii.gz.png
# │       │   │   ├── cnr0003.nii.gz.png
# │       │   │   ├── qc.json
# │       │   │   ├── qc.pdf
# │       │   │   ├── ref_list.png
# │       │   │   ├── ref.txt
# │       │   │   └── vdm.png
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-5tt-hsvs_probseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-brain_mask.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-brain-uncorrected_mask.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-gmwm_probseg.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-isovf_noddi.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-ndi_noddi.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-nodif-brain_dwi.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-nodif_dwi.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-odi_noddi.nii.gz
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-preproc_dwi.bval
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-preproc_dwi.bvec
# │       │   ├── sub-[subjID]_[sessionID]_space-dwi_desc-preproc_dwi.nii.gz
# │       │   └── sub-[subjID]_[sessionID]_space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz
# │       ├── figures
# │       │   ├── sub-[subjID]_[sessionID]_label-300P7N_overlay.png
# │       │   ├── sub-[subjID]_[sessionID]_label-400P17N_overlay.png
# │       │   ├── sub-[subjID]_[sessionID]_label-400P7N_overlay.png
# │       │   ├── sub-[subjID]_[sessionID]_label-aparc500_overlay.png
# │       │   ├── sub-[subjID]_[sessionID]_label-BNA_overlay.png
# │       │   └── sub-[subjID]_[sessionID]_maskQC.png
# │       ├── fmap
# │       │   ├── sub-[subjID]_[sessionID]_acq-APPA_desc-refparams.tsv
# │       │   └── sub-[subjID]_[sessionID]_acq-APPA_space-dwi_desc-4topup_epi.nii.gz
# │       ├── logs
# │       │   ├── sub-[subjID]_[sessionID]_dwi-anat2dwi.log
# │       │   ├── sub-[subjID]_[sessionID]_dwi-noddi.log
# │       │   ├── sub-[subjID]_[sessionID]_dwi-preproc.log
# │       │   ├── sub-[subjID]_[sessionID]_eddy.log
# │       │   └── sub-[subjID]_[sessionID]_topup.log
# │       └── xfms
# │           ├── diff-2-T1w.mat
# │           ├── sub-[subjID]_[sessionID]_epireg_fast_wmedge.nii.gz
# │           ├── sub-[subjID]_[sessionID]_epireg_fast_wmseg.nii.gz
# │           ├── sub-[subjID]_[sessionID]_epireg_init.mat
# │           ├── sub-[subjID]_[sessionID]_epireg_inversed.mat
# │           ├── sub-[subjID]_[sessionID]_epireg.mat
# │           ├── sub-[subjID]_[sessionID]_epireg.nii.gz
# │           └── T1w-2-diff.mat


# SLURM INPUTS
#SBATCH --job-name=dwipipeline
#SBATCH --mem=20M
#SBATCH --partition=luna-cpu-short
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=1
#SBATCH --time=00-4:00:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A_%a.log

# Load subject for this array task
subj=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${subjects})


# random delay
duration=$((RANDOM % 10 + 2))
echo "INITIALIZING...(wait a sec)"
sleep ${duration}

echo "Processing subject: ${subj}"

mkdir -p ${scriptdir}/${subj}
cd ${scriptdir}/${subj}
###########################
##  DWI-PREPROCESSING    ##
###########################
total_sessions=0
done_sessions=0

for dwidir in ${bidsdir}/${subj}/{,ses*/}dwi; do
    total_sessions=$((total_sessions+1))

    if [ ! -d ${dwidir} ]; then
        done_sessions=$((done_sessions+1))
        continue
        
    fi
    sessiondir=$(dirname ${dwidir})
    session=$(echo "${sessiondir}" | grep -oP "(?<=${subj}/).*")
    if [ -z "${session}" ]; then
        sessionpath=/
        sessionfile=_
    else
        sessionpath=/${session}/
        sessionfile=_${session}_
        
    fi
    
    if [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz ]] &&
    [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bval ]] &&
    [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec ]]; then
        echo -e "${GREEN}${subj}${sessionfile} already preprocessed with eddy${NC}"
        
        done_sessions=$((done_sessions+1))
        continue
    fi
    
done

if [ "$total_sessions" -gt 0 ] && [ "$done_sessions" -eq "$total_sessions" ]; then
    echo "All sessions already have eddy output for ${subj}."
    echo "or no dwi present"
else
    
    echo "Submitting preprocessing job for ${subj}"
    job_id_preproc=$(sbatch --parsable ${scriptdir}/dwi-02a-preproc.sh -i ${bidsdir} -o ${outputdir} -w ${workdir} -s ${subj} -c ${scriptdir})
    echo "Submitting eddy job for ${subj}"
    job_id_eddy=$(sbatch --gres=gpu:1g.10gb:1 --parsable --kill-on-invalid-dep=yes --dependency=afterok:$job_id_preproc ${scriptdir}/dwi-02b-eddyGPU.sh -w ${workdir} -o ${outputdir} -s ${subj} -z ses-t0 -m default)
fi

###########################
##      DWI - NODDI      ##
###########################
if [[ ${noddi} == 1 ]]; then
    total_sessions=0
    done_sessions=0
    
    for dwidir in ${bidsdir}/${subj}/{,ses*/}dwi; do
        if [ ! -d ${dwidir} ]; then
            
            continue
            
        fi
        total_sessions=$((total_sessions+1))
        sessiondir=$(dirname ${dwidir})
        session=$(echo "${sessiondir}" | grep -oP "(?<=${subj}/).*")
        if [ -z "${session}" ]; then
            sessionpath=/
            sessionfile=_
        else
            sessionpath=/${session}/
            sessionfile=_${session}_
            
        fi
        
        if [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-ndi_noddi.nii.gz ]] &&
        [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-odi_noddi.nii.gz ]] &&
        [[ -f ${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-isovf_noddi.nii.gz ]]; then
            echo -e "${GREEN}${subj}${sessionfile} already preprocessed with NODDI${NC}"
            
            done_sessions=$((done_sessions+1))
            continue
        fi
        
    done
    
    if [ "$total_sessions" -gt 0 ] && [ "$done_sessions" -eq "$total_sessions" ]; then
        echo "All sessions already have NODDI output for ${subj}."
    else
        
        echo "Submitting NODDI job for ${subj}"
        if [ -z "${job_id_eddy}" ]; then 
        job_id_noddi=$(sbatch --parsable ${scriptdir}/dwi-02c-prep4noddi.sh -w ${workdir} -o ${outputdir} -s ${subj} -c ${scriptdir})

        else 
        job_id_noddi=$(sbatch --parsable --dependency=afterok:$job_id_eddy ${scriptdir}/dwi-02c-prep4noddi.sh -w ${workdir} -o ${outputdir} -s ${subj} -c ${scriptdir})
        echo "NODDI job ID: $job_id_noddi"
        fi


    fi
    
else
    job_id_noddi=$job_id_eddy
fi

###########################
## DWI-2-T1 registration ##
###########################
# echo "Submitting Anat-2-DWI registration job for ${subj}"
# job_id_anat2dwi=$(sbatch --parsable --dependency=afterok:$job_id_noddi ${scriptdir}/dwi-03-anat2dwi.sh -i ${bidsdir} -o ${outputdir} -f ${freesurferdir} -w ${workdir} -s ${subj} -c ${scriptdir})
# echo "Anat-2-DWI registration job ID: $job_id_anat2dwi"

# ###########################
# ##   DWI-TRACTOGRAPHY    ##
# ###########################
# echo "Submitting FOD + Tractogram job for ${subj}"
# job_id_fodtck=$(sbatch --parsable --dependency=afterok:$job_id_anat2dwi ${scriptdir}/dwi-04a-connectome.sh -i ${bidsdir} -o ${outputdir} -w ${workdir} -s ${subj} -c ${scriptdir})
# echo "FOD + Tractogram job ID: $job_id_fodtck"

# echo "Submitting Tract-to-Connectome job for ${subj}"
# job_id_tck2conn=$(sbatch --parsable --dependency=afterok:$job_id_fodtck ${scriptdir}/dwi-04b-tracts2conn_v2.sh -i ${bidsdir} -o ${outputdir} -w ${workdir} -s ${subj} -n 50M)
# echo "Tract-to-Connectome job ID: $job_id_tck2conn"


final_job_id=$(sbatch --wait --parsable --kill-on-invalid-dep=yes --dependency=afterok:$job_id_eddy --mem=20M -c 1 --time=00-00:00:10 --wrap "echo 'All jobs completed for ${subj}'")


#final_job_id=$(sbatch --parsable --dependency=afterok:$job_id_tck2conn --wrap "echo 'All jobs completed for ${subj}'")
echo "Final dummy job ID: $final_job_id"
