# Running the Container

Regardless of engine (Docker, Podman, or Apptainer), the invocation
pattern is always:

```text
<container-run-command> <container-image> <pipeline> <path-to-spec.json>
```

where `<pipeline>` is one of: `dwi-preproc`, `dwi-tracto`, `dwi-qc`,
`dwi-params`, `dwi-send`.

The `spec.json` file contains the full paths to the input and output
directories, and the settings for the pipeline — see {doc}`spec-json`.

## dwi-preproc

1. Preprocessing of diffusion MRI scans, including denoising,
   susceptibility-induced distortion correction (with PEPolar field
   maps or using synthetic field maps), and eddy and motion correction.
   By default (parameter in `spec.json`) eddy is run with the
   `--repol` flag (outlier replacement), but without intra-volume
   motion or susceptibility-by-motion correction.

   The pipeline checks whether sampling was performed on a whole or
   half sphere and adjusts the eddy settings accordingly — see the
   [FSL eddy documentation](https://fsl.fmrib.ox.ac.uk/fsl/docs/diffusion/eddy/index.html)
   for more information.

2. Processing of the T1-weighted anatomical image and registration to
   the diffusion MRI scan. This includes running FreeSurfer (v8.2.0) —
   if the participant doesn't already have FreeSurfer output available
   in `freesurferdir` — creating a "5 tissue type" (5TT) segmentation
   file for anatomically constrained tractography (ACT), and warping
   atlases to diffusion space.

```{note}
Eddy correction and FreeSurfer + anat2dwi registration run in parallel
if `nthreads` > 2.
```

**Outputs (in dwi space):**

- eddy-corrected diffusion MRI
  (`/dwi/<subj>[_<session>]_space-dwi_desc-preproc_dwi.nii.gz`)
- 5TT image
  (`/anat/<subj>[_<session>]_space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz`)
- gray/white matter boundary image
  (`/anat/<subj>[_<session>]_space-dwi_res-high_desc-gmwm_probseg.nii.gz`)
- atlas segmentations
  (`/anat/<subj>[_<session>]_space-dwi_res-high_atlas-<atlas>_dseg.nii.gz`)

## dwi-tracto

Performs tractography (fiber orientation estimation and streamline
generation) on the preprocessed DWI data from the previous step. Throws
an error if no preprocessed diffusion MRI or 5TT segmentation file is
available. The default for the ENIGMA OCD project is 20 million seeds
(set in `spec.json`). After producing a tractogram, it determines the
number of streamlines between the different atlas regions and produces
structural connectivity matrices.


**Outputs:**

- Tractogram
  (`/dwi/<subj>[_<session>]_space-dwi_tracto-<nstreamlines>.tck`)
- SIFT weights
  (`/dwi/<subj>[_<session>]_space-dwi_tracto-<nstreamlines>_desc-sift_weights.txt`)
- Connectivity matrix
  (`/conn/<subj>[_<session>]_atlas-<atlas>_desc-streams_connmatrix.csv`)


## dwi-qc

Produces a QC HTML page per participant to quickly assess the quality
of the diffusion MRI preprocessing, the registration between the T1w
and diffusion scans, and the tractography. This step can be run right
after `dwi-preproc` (recommended) or after `dwi-tracto`. When run after
`dwi-preproc`, the `dwi-tracto` outputs will (naturally) not be shown.

**Outputs:**
`dwi-preproc/<subj>[_<session>].html`


## dwi-params

Collects the scanning parameters from the T1w structural MRI, diffusion
MRI, and field maps (if available) and summarizes them in an HTML
report. This lets users quickly spot discrepancies in imaging
parameters between participants in a sample, and is used to write up
the methods for the ENIGMA OCD paper.

## dwi-send

Copies the completed derivatives and makes them ready to send to the
project lead, along with the clinical covariates file.

```{tip}
Run these **in dependency order** for a fresh subject:
`dwi-preproc` → `dwi-qc` → `dwi-tracto` → `dwi-qc` (updates previous qc html) → `dwi-send`, using `dwi-params`
at any point you just need parameter extraction.
```

## Helper scripts

Because the entry point reads *all* its locations from the `spec.json`
file (BIDS dir, output dir, work dir, FreeSurfer dir), you must:

1. Bind-mount every host directory referenced in `spec.json` into the
   container.
2. Make sure the **paths inside `spec.json` are the container-side
   paths**, not the host paths — e.g. mount `/data/bids` on the host to
   `/bids` in the container, and set `"bidsdir": "/bids"` in
   `spec.json`, not the host path.

## Engine Specific Commands
::::{tab-set}

:::{tab-item} Docker
````bash
docker run --rm \
  -v /host/path/to/bids:/bids \
  -v /host/path/to/output:/derivatives \
  -v /host/path/to/work:/work \
  -v /host/path/to/freesurfer:/freesurfer \
  -v /host/path/to/spec.json:/spec/spec.json \
  cvriend/tractoprep{{RELEASE_TAG}} \
  dwi-preproc /spec/spec.json
````
:::

:::{tab-item} Podman
`````bash
podman run --rm \
  -v /host/path/to/bids:/bids \
  -v /host/path/to/output:/derivatives \
  -v /host/path/to/work:/work \
  -v /host/path/to/freesurfer:/freesurfer \
  -v /host/path/to/spec.json:/spec/spec.json \
  cvriend/tractoprep{{RELEASE_TAG}} \
  dwi-preproc /spec/spec.json
`````
You can add a `:Z` suffix to relabel bind mounts for
SELinux-enforcing hosts — omit it if not applicable.
:::

:::{tab-item} Apptainer
`````bash
 apptainer run --cleanenv \
      --bind /host/path/to/bids:/bids,/host/path/to/output:/derivatives,/host/path/to/work:/work,/host/path/to/freesurfer:/freesurfer \
      tractoprep_{{RELEASE_TAG}}.sif dwi-preproc /spec/spec.json
`````
:::
::::

- Make sure you have sufficient disk space for a workdir (scratch
  space, deleted automatically after successful completion of a
  subject) and an output dir.


```{note}
**Resource recommendations:** on SLURM clusters, ...
*(fill in your site's SLURM resource recommendations here.)*
```
