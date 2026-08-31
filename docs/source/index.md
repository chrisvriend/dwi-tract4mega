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

::::{tab-set}
:::{tab-item} Docker
````bash
docker run --rm -v <bidsdir>:/bids -v <outputdir>:/output -v <workdir>:/work \
  -v <freesurferdir>:/freesurfer -v <path/to/spec.json>:/spec/spec.json \
  cvriend/tractoprep:{{RELEASE_TAG}} //pipeline// /spec/spec.json
````
:::

:::{tab-item} Podman
`````bash
podman run --rm -v <bidsdir>:/bids -v <outputdir>:/output -v <workdir>:/work \
  -v <freesurferdir>:/freesurfer -v <path/to/spec.json>:/spec/spec.json \
  cvriend/tractoprep:{{RELEASE_TAG}} //pipeline// /spec/spec.json
`````
:::

:::{tab-item} Apptainer
`````bash
apptainer run --cleanenv --bind \
      /host/path/to/bids:/bids, \
      /host/path/to/output:/derivatives, \
      /host/path/to/work:/work, \
      /host/path/to/freesurfer:/freesurfer \
      /host/path/to/spec.json:/spec/spec.json
      tractoprep_{{RELEASE_TAG}}.sif //pipeline// /spec/spec.json
`````
:::
::::


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
QC
troubleshooting
Acknowledgement
```
