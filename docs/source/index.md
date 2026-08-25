# TractoPrep Processing Manual

TractoPrep packages a set of tools (denoising, topup/eddy, FreeSurfer,
anatomical-to-diffusion registration, tractography, and connectome
construction) into a single container image. It draws on tools from
several open-source packages — MRtrix3, FSL, ANTs, FreeSurfer, and
Synb0-DISCO.

The pipeline is agnostic to the parameters of the diffusion MRI scan and
can work with single- and multishell diffusion MRI, single- and
multiband acquisitions, acquisitions with multiple phase-encoding
directions, and with or without PEPolar fieldmaps for susceptibility-
induced distortion correction.

## Quick start

```bash
# Docker
docker run --rm -v <bids>:/bids -v <out>:/output -v <work>:/work \
  -v <fs>:/freesurfer -v <spec.json>:/spec/spec.json \
  cvriend/tractoprep dwi-preproc /spec/spec.json

# Podman
podman run --rm -v <bids>:/bids:Z -v <out>:/output:Z -v <work>:/work:Z \
  -v <fs>:/freesurfer:Z -v <spec.json>:/spec/spec.json:Z \
  docker://cvriend/tractoprep dwi-preproc /spec/spec.json

# Apptainer/Singularity (local .sif)
apptainer exec --bind <bids>:/bids --bind <out>:/output --bind <work>:/work \
  --bind <fs>:/freesurfer --bind <spec.json>:/spec/spec.json \
  tractoprep.sif /tracto/dwi-00-entry.sh dwi-preproc /spec/spec.json
```

See {doc}`installation` and {doc}`spec-json` before running this for real —
you need a container engine, the image, and a filled-in `spec.json`.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
data-organization
running
spec-json
output-structure
workflow
troubleshooting
Software
```
