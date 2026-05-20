#!/bin/bash

############################################
#  DWI PIPELINE – SLURM-CORRECT VERSION
#  (C) C. Vriend – revised & aangepast
############################################

########################
# SLURM DIRECTIVES
########################

#SBATCH --job-name=dwipipeline
#SBATCH --partition=luna-cpu-short
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=20M
#SBATCH --time=04:00:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A_%a.log

set -euo pipefail

########################
# INPUTS
########################
subj=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${subjects}")
method=slmlinearsdc
nstreamlines=50M
preproc_only=1

echo "Processing subject: ${subj}"

mkdir -p ${scriptdir}/${subj}
cd ${scriptdir}/${subj}

########################
# RANDOM DELAY (ARRAY SAFETY)
########################
sleep $((RANDOM % 10 + 2))

########################
# DISCOVER SESSIONS
########################
sessions=()
for dwidir in ${bidsdir}/${subj}/{,ses*/}dwi; do
    [ -d "$dwidir" ] || continue
    sessiondir=$(dirname "$dwidir")
    session=$(echo "$sessiondir" | grep -oP "(?<=${subj}/).*" || true)
    sessions+=("$session")
done

if [ "${#sessions[@]}" -eq 0 ]; then
    echo "No DWI data found for ${subj}"
    exit 1
fi


########################
# PREPROCESSING
########################
echo "Submitting preprocessing and eddy jobs for ${subj}"

preproc_jobs=()
eddy_jobs=()
anat_jobs=()
tck_jobs=()
tck2conn_jobs=()    
eddy_submitted=()
anat_submitted=()
tck_submitted=()
tck2conn_submitted=()
noddi_jobs=()

for session in "${sessions[@]}"; do
    echo "Processing session: ${subj} ${session:-}"

    # Set session path/file
    if [[ -n "${session}" ]]; then
        sessionpath="/${session}/"
        sessionfile="_${session}_"
    else
        sessionpath="/"
        sessionfile="_"
    fi

    # Build session argument if session is not empty
    session_arg=()
    if [[ -n "${session}" ]]; then
        session_arg=(-z "${session}")
    fi

    # Check if preprocessing is already done
    if [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.nii.gz" ] &&
       [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-preproc_dwi.bvec" ] &&
       [ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz" ]; then

        echo "Preprocessing already done for ${subj} ${session:-}, skipping preproc and eddy."
        preproc_jobs+=("")
        eddy_jobs+=("")
        eddy_submitted+=(0)
        # Set eddy dependency to nothing, since it's already done
        eddy_dep=""
    else
        # check presence of fmap to determine memory usage
        if [ ! -d "${bidsdir}/${subj}/${session}/fmap" ]; then 
            echo "No fieldmap found for ${subj} ${session:-}, using higher memory"
            mem_request="24G"
        else 
            echo "Fieldmap found for ${subj} ${session:-}, using lower memory"
            mem_request="8G"
        fi 

        # Submit preproc job
        job_id_preproc=$(sbatch --mem=${mem_request} --parsable \
            "${scriptdir}/dwi-02a-preproc.sh" \
            -i "${bidsdir}" \
            -o "${outputdir}" \
            -w "${workdir}" \
            -s "${subj}" \
            -c "${scriptdir}" \
            "${session_arg[@]}"
        )
        preproc_jobs+=("${job_id_preproc}")

        echo "Submitting eddy job for ${subj} ${session:-nosession}"

        # Submit eddy job with dependency on preproc
        job_id_eddy=$(sbatch --parsable \
            --gres=gpu:1g.10gb:1 \
            --dependency=afterok:${job_id_preproc} \
            --kill-on-invalid-dep=yes \
            "${scriptdir}/dwi-02b-eddyGPU.sh" \
            -w "${workdir}" \
            -o "${outputdir}" \
            -s "${subj}" \
            -m "${method}" \
            "${session_arg[@]}"
        )
        eddy_jobs+=("${job_id_eddy}")
        eddy_submitted+=(1)
        eddy_dep="--dependency=afterok:${job_id_eddy} --kill-on-invalid-dep=yes"
    fi

    # NODDI-specific code, inside the session loop, after eddy
    if [[ "${noddi:-0}" == 1 ]]; then

        if [[ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-ndi_noddi.nii.gz" ]] &&
            [[ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-isovf_noddi.nii.gz" ]] &&
            [[ -f "${outputdir}/dwi-preproc/${subj}${sessionpath}dwi/${subj}${sessionfile}space-dwi_desc-odi_noddi.nii.gz" ]]; then
                        echo "Preprocessing already done for ${subj} ${session:-}, skipping NODDI jobs."
                        noddi_jobs+=("")
        else

            echo "Submitting NODDI jobs for ${subj} ${session:-nosession}"
            job_id_noddi_prep=$(sbatch --parsable \
                ${eddy_dep} \
                "${scriptdir}/dwi-02c-prep4noddi.sh" \
                -w "${workdir}" \
                -o "${outputdir}" \
                -s "${subj}" \
                "${session_arg[@]}"
            )
            job_id_noddi_gpu=$(sbatch --parsable \
                --gres=gpu:1g.10gb:1 \
                --dependency=afterok:${job_id_noddi_prep} \
                --kill-on-invalid-dep=yes \
                "${scriptdir}/dwi-02d-noddi.sh" \
                -w "${workdir}" \
                -o "${outputdir}" \
                -s "${subj}" \
                "${session_arg[@]}"
            )
            noddi_jobs+=("${job_id_noddi_gpu}")

        fi
    else
        noddi_jobs+=("")
    fi

    # Continue with the rest of the pipeline
    job_id_anat2dwi=$(sbatch --parsable \
        ${eddy_dep} \
        "${scriptdir}/dwi-03-anat2dwi.sh" \
        -i "${bidsdir}" \
        -f "${freesurferdir}" \
        -w "${workdir}" \
        -o "${outputdir}" \
        -s "${subj}" \
        -c "${scriptdir}" \
        "${session_arg[@]}"
    )
    anat_jobs+=("${job_id_anat2dwi}")
    anat_submitted+=(1)
    
    # perform tractography and connectome steps only if preproc_only is not set
    if [ "${preproc_only:-0}" == 1 ]; then
        echo "Preprocessing only flag is set, skipping tractography and connectome steps for ${subj} ${session:-nosession}"
        tck_jobs+=("")
        tck2conn_jobs+=("")
        tck_submitted+=(0)
        tck2conn_submitted+=(0)
        continue
    fi

    job_id_fodtck=$(sbatch --parsable \
        --dependency=afterok:${job_id_anat2dwi} \
        --kill-on-invalid-dep=yes \
        "${scriptdir}/dwi-04a-connectome.sh" \
        -i "${bidsdir}" \
        -o "${outputdir}" \
        -w "${workdir}" \
        -s "${subj}" \
        "${session_arg[@]}"
    )
    tck_jobs+=("${job_id_fodtck}")
    tck_submitted+=(1)

    job_id_tck2conn=$(sbatch --parsable \
        --dependency=afterok:${job_id_fodtck} \
        --kill-on-invalid-dep=yes \
        "${scriptdir}/dwi-04b-tracts2conn.sh" \
        -i "${bidsdir}" \
        -w "${workdir}" \
        -o "${outputdir}" \
        -s "${subj}" \
        -n "${nstreamlines}" \
        "${session_arg[@]}"
    )
    tck2conn_jobs+=("${job_id_tck2conn}")
    tck2conn_submitted+=(1)

   
done

########################
# FINAL SENTINEL (WAIT HERE)
########################

all_jobs=()

collect_jobs() {
    local array_name="$1"
    local -n arr_ref="$array_name"   # nameref to the actual array
    for job_id in "${arr_ref[@]}"; do
        if [[ "$job_id" =~ ^[0-9]+$ ]]; then
            all_jobs+=("$job_id")
        fi
    done
}

if [ "${preproc_only:-0}" -eq 1 ]; then
    collect_jobs eddy_jobs
    collect_jobs noddi_jobs
    collect_jobs anat_jobs
else
    collect_jobs eddy_jobs
    collect_jobs noddi_jobs
    collect_jobs anat_jobs
    collect_jobs tck_jobs
    collect_jobs tck2conn_jobs
fi





if [[ "${#all_jobs[@]}" -gt 0 ]]; then
    dep_string=$(IFS=:; printf "%s" "${all_jobs[*]// /:}")
    dep_arg=(--dependency=afterok:${dep_string} --kill-on-invalid-dep=yes)
else
    dep_arg=()
    dep_string=""
fi

if [ -z "$dep_string" ]; then
    echo "No jobs were submitted, skipping final sentinel job."
else
    echo "Final sentinel job will depend on jobs: ${dep_string}"
fi

target_dir="${workdir}/${subj}/freesurfer"
final_job_id=$(sbatch --wait --parsable \
    "${dep_arg[@]}" \
    --job-name="dwi_Hodor" \
    --time=00:01:00 -c 1 --mem=10M \
    --wrap "echo 'Pipeline finished for ${subj}'; rm -rf \"$target_dir\"")


echo "Pipeline completed for ${subj} (final job ${final_job_id})"
