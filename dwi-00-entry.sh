#!/bin/bash 
# entrypoint script 

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

# Parse command line arguments
Usage() {
    echo "Usage: dwi-preproc|dwi-tracto|dtitk <specfile>"
    echo "select one of these pipelines and supply the path to the specfile (spec.json)"
    exit 1
}

[ "$2" = "" ] && Usage


pipeline=$1
specfile=$2
scriptdir=/tracto

# read spec.json file
for key in $(jq -r 'keys[]' ${specfile}); do
  value=$(jq -r --arg k "$key" '.[$k]' ${specfile})
  declare "$key"="$value"
done


# check if all variables are non-empty
for var in pipeline specfile scriptdir; do
  if [[ -z "${!var}" ]]; then
    echo "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  

echo
log "$BLUE" "------------------------"
log "$BLUE" "------------------------"
log "$BLUE" "--- INPUT VARIABLES  ---"
log "$BLUE" "pipeline: ${pipeline}"
log "$BLUE" "BIDS path: ${bidsdir}"
log "$BLUE" "Output path: ${outputdir}"
log "$BLUE" "Work path: ${workdir}"
log "$BLUE" "Freesurferdir: ${freesurferdir}"
log "$BLUE" "------------------------"
log "$BLUE" "------------------------"


case "$pipeline" in
    dwi-preproc)
        ${scriptdir}/dwi-preproc_wrapper.sh ${specfile}
    ;;
    dwi-tracto)
        ${scriptdir}/dwi-tracto_wrapper.sh ${specfile}
    ;;
    dtitk)
    :
    ;;
esac





