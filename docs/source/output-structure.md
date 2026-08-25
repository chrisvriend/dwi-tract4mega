# Output Structure

`dwi-preproc` writes BIDS-derivatives-style output under
`outputdir/dwi-preproc/<subj>/[<session>/]`, including the following
subfolders (names as produced by the pipeline):

```text
outputdir/
└── dwi-preproc/
    └── sub-01/
              ses-Tx / (optional)
              ├── dwi/       # preprocessed DWI (denoised, degibbs'd,
              │              # eddy/topup-corrected), bval/bvec
              ├── fmap/      # topup fieldmap outputs, unwarped reference images
              ├── anat/      # anatomical-to-DWI registration, parcellations
              │              # (e.g., atlas-300P17N, atlas-400P17N), 5TT/GM-WM
              │              # interface images
              ├── qc/        # CNR maps, denoising residuals, and other
              │              # QC-relevant intermediates
              └── log/       # per-step timestamped log files (preproc, eddy,
                            # anat2dwi)
```

The wrapper verifies that a specific set of expected output files
exist before declaring success and cleaning up the workdir for that
subject. If any are missing, it leaves the workdir intact for
inspection and reports an error rather than deleting your scratch
data.

## dwi-preproc key final outputs 

- `*_space-dwi_desc-preproc_dwi.nii.gz` / `.bval` / `.bvec`
- `*_space-dwi_label-cnr-maps_desc-preproc_dwi.nii.gz`
- `*_space-dwi_res-high_atlas-300P17N_dseg.nii.gz` and the
  `atlas-400P17N` variant
- `*_space-dwi_res-high_desc-5tt-hsvs_probseg.nii.gz`
- `*_space-dwi_res-high_desc-gmwm_probseg.nii.gz`

`dwi-tracto` 
`outputdir/dwi-tracto/<subj>/[<session>/]`, including the following
subfolders (names as produced by the pipeline):

```text
outputdir/
└── dwi-tracto/
    └── sub-01/
              ses-Tx / (optional)
              ├── dwi/       # WM FOD image, tractogram and SIFT weights,
              ├── conn/      # structural connectivity matrices 
              │              # (e.g., atlas-300P17N, atlas-400P17N)
              ├── rpf/      # response functions  
              │              # (i.e. for CSF, WM and GM)
              ├── qc/        # tissue RGB, downsampled tractogram
              └── log/       # per-step timestamped log files (preproc, eddy,
                            # anat2dwi)
```

## dwi-tracto key final outputs 

- `*_space-dwi_tracto-<nstreamlines>.tck`
- `*_space-dwi_tracto-<nstreamlines_desc-sift_weights.txt`
- `*_atlas-300P17N_desc-streams_connmatrix.csv`
- `*_atlas-400P17N_desc-streams_connmatrix.csv`
