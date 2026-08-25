# Imaging Data Organization Requirements

- A **BIDS-organized** dataset containing at minimum `anat/` (T1w) and
  `dwi/` data with valid JSON sidecars. The pipeline reads
  `PhaseEncodingDirection` and `TotalReadoutTime` directly from the
  diffusion MRI JSON file — make sure these fields are populated.
- Optionally, field maps with `PhaseEncodingDirection`,
  `TotalReadoutTime`, and a correct `IntendedFor` entry (pointing to the
  diffusion MRI image(s)). If no field maps are present, the pipeline
  automatically falls back to **Synb0-DISCO** (synthetic b0) for
  susceptibility distortion correction — this needs the anatomical T1w
  to be available.
- (Optional) A **FreeSurfer output directory** (`freesurferdir`) if you
  already have FreeSurfer output available for these participants.

The files can be organized with or without a session label, as long as
the organization is BIDS-compliant. Either of these is correct:

```text
sub-0001/dwi/sub-0001_dwi.[nii.gz|json|bval|bvec]
sub-0001/ses-T0/dwi/sub-0001_ses-T0_dwi.[nii.gz|json|bval|bvec]
```

If you use a session label, indicate it in the `spec.json` file — see
{doc}`spec-json`.
