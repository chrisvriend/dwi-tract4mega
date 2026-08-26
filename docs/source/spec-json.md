# Spec.json file

```json
{
  "Software_version": "1.0.4",
  "subj": "sub-0001",
  "session": "",
  "eddy_method": "default",
  "nstreamlines": "20M",
  "bidsdir": "/scratch/user/bids",
  "outputdir": "/net/beegfs/user/derivatives",
  "workdir": "/scratch/users/work",
  "freesurferdir": "/net/beegfs/user/derivatives/freesurfer",
  "nthreads": 8,
  "atlas": "400P17N",
  "sitename": "OCD_site1"
}
```

| Key | Meaning |
|---|---|
| `subj` | Subject ID of the scan to run `dwi-preproc`, `dwi-tracto` or `dwi-qc` on |
| `session` | Session ID, e.g. `ses-T0` — leave blank if there are no sessions |
| `eddy_method` | Method used for eddy correction (`default`) |
| `nstreamlines` | Number of streamlines for tractography (default `20M`) |
| `bidsdir` | Path to the BIDS directory |
| `outputdir` | Path to the output directory (derivatives) |
| `workdir` | Path to the working directory (subject directories are deleted after processing) |
| `freesurferdir` | Path to the FreeSurfer output directory |
| `nthreads` | Number of threads to use per subject — see {doc}`running` for parallelization notes |
| `atlas` | Atlas used to produce the connectivity matrix in the QC HTML report |
| `sitename` | Name of your site/sample (only used when sending data) |

## Generating a template

A helper script, `create_spec_template.sh`, is provided to create an
empty `spec.json` for you rather than writing one by hand:


::::{tab-set}

:::{tab-item} Docker
`````bash
docker run --rm -v /host/path/to/specjsonlocation:/spec \
docker://cvriend/tractoprep:{{RELEASE_TAG}} \
/tracto/helpers/create_spec_template.sh > /spec/spec.json
`````
:::

:::{tab-item} Podman
`````bash
podman run --rm -v /host/path/to/specjsonlocation:/spec \
docker://cvriend/tractoprep:{{RELEASE_TAG}} \
/tracto/helpers/create_spec_template.sh > /spec/spec.json
`````
:::
:::{tab-item} Apptainer
`````bash
apptainer exec tractoprep.{{RELEASE_TAg}} /tracto/helpers/create_spec_template.sh > spec.json
`````
:::
::::

## Low-memory systems

The preprocessing script supports an additional `lowmem` flag in `spec.json` for
running FreeSurfer 8.2.0 on memory-constrained nodes. Set it to enable
this behavior (expect longer run time)

