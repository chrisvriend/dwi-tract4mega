#!/usr/bin/env bash

########################################
# CONFIGURATION
########################################

jsonfile=$1

# Colors
NC='\033[0m'
BLUE='\033[0;34m'
GREEN='\033[0;32m'
RED='\033[0;31m'

Usage() {
  echo "Usage: $0 <path_to_spec.json>"
  echo "prepare datapackage to send to for analyses"
  exit 1
}

# Helper function for colored output
log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

# Check input json
if [[ -z "$jsonfile" || ! -f "$jsonfile" ]]; then
  log "$RED" "Error: spec.json file not provided or does not exist."
  Usage
fi

# read spec.json file into variables
for key in $(jq -r 'keys[]' "${jsonfile}"); do
  value=$(jq -r --arg k "$key" '.[$k]' "${jsonfile}")
  declare "$key"="$value"
done

# check if required variables are non-empty
for var in outputdir sitename ; do
  if [[ -z "${!var}" ]]; then
    log "$RED" "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  

mkdir -p "${outputdir}/tosend"


### DWI parameters 
if [[ ! -d ${outputdir}/params && ! -f ${outputdir}/params/parameter_report.html ]]; then

  log "${YELLOW}" "extract dwi / T1w parameters (to write methods section)"

  /tracto/dwi-extract_parameters.sh ${jsonfile}

else 

  cp -r ${outputdir}/params ${outputdir}/tosend/ 
fi 

### output files


# connectivity matrices 
cd "${outputdir}/dwi-tracto"

# List of required files (relative to conn/, use ${subj} placeholder)
# NOTE: single-quoted so ${subj} is NOT expanded here — it's expanded later via eval,
# once ${subj} is actually set inside the loop.
required_files=(
    '${subj}_atlas-300P7N_desc-streams_connmatrix.csv'
    '${subj}_atlas-300P7N_desc-FA_connmatrix.csv'
    '${subj}_atlas-400P7N_desc-streams_connmatrix.csv'
    '${subj}_atlas-400P7N_desc-FA_connmatrix.csv'
    '${subj}_atlas-BNA_desc-streams_connmatrix.csv'
    '${subj}_atlas-BNA_desc-FA_connmatrix.csv'

    # '${subj}_another_required_file.ext'
    # '${subj}_yet_another_file.ext'
)

# Usage: have_all_files <basepath>
# (relies on $subj being set in the calling scope)
have_all_files() {
    local basepath="$1"
    local f fname
    for f in "${required_files[@]}"; do
        fname=$(eval echo "$f")
        [ -f "${basepath}/${fname}" ] || return 1
    done
    return 0
}

for subj in sub-*; do
    [ -d "${subj}" ] || continue

    src_conn=""
    ses=""

    # 1) Try without session: sub-XXX/conn/
    if have_all_files "${subj}/conn"; then
        src_conn="${subj}/conn"
    else
        # 2) Try with any session: sub-XXX/ses-YYY/conn/
        for sesdir in "${subj}"/ses-*; do
            [ -d "${sesdir}" ] || continue
            if have_all_files "${sesdir}/conn"; then
                src_conn="${sesdir}/conn"
                ses="$(basename "${sesdir}")"
                break
            fi
        done
    fi

    if [ -z "${src_conn}" ]; then
        log "$RED" "ERROR! ${subj} does not have all required output files"
        log "$RED" "rerun tractography for ${subj} or delete the folder"
        exit 1
    fi

    # Build destination, preserving session subdir if present
    if [ -n "${ses}" ]; then
        dest_conn="${outputdir}/tosend/${subj}/${ses}/conn"
    else
        dest_conn="${outputdir}/tosend/${subj}/conn"
    fi
    mkdir -p "${dest_conn}"

    for f in "${required_files[@]}"; do
        fname=$(eval echo "$f")
        rsync -a "${src_conn}/${fname}" "${dest_conn}/"
    done
done


# qc html files
cd "${outputdir}/dwi-preproc"

for subjhtml in sub-*.html; do

  rsync -a ${subjhtml} ${outputdir}/tosend/

done


cd "${outputdir}"
mv tosend ${sitename}_dwi_output
tar -jxcvf ${sitename}_dwi_output.tar.bz2 ${sitename}_dwi_output
rm -r ${sitename}_dwi_output

log "${GREEN}" "----------------"
log "${GREEN}" "Zipping complete"
log "${GREEN}" "----------------"
echo
log "${GREEN}" "Please upload this tar file along with the Covariates file using the provided Surfdrive link"
echo
log "${GREEN}" "Many thanks for your contribution to this ENIGMA project!!"




 

