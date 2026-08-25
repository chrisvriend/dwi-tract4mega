# Typical End-to-End Workflow

1. **Organize data** as BIDS and populate `PhaseEncodingDirection` and `TotalReadoutTime` in the JSON sidecar for the dwi + `IntendedFor` in fieldmaps (if available)
2. **Pull the container** — see {doc}`installation`.
3. **Generate a spec.json** with `helpers/create_spec_template.sh` and
   fill in subject, session (if applicable), and all paths — see {doc}`spec-json`.
4. **Run `dwi-preproc`** for each subject/session. Loop this per
   subject on a cluster using the repo's SLURM array launchers, or
   script your own loop calling Docker/Podman/Apptainer per subject. 
5. **Inspect QC outputs** (`dwi-qc`) before proceeding — check
   registration quality (`check_atlasreg.py`) and CNR/denoising maps.
6. **Run `dwi-tracto`** on subjects that pass QC.
7. **Run `dwi-send`** to export the completed derivatives.
8. Use **`dwi-params`** at any stage to pull out acquisition/processing
   parameters for logging or a shared summary spreadsheet.

```{seealso}
{doc}`running` for details on each pipeline step, and
{doc}`troubleshooting` if a step fails.
```
## SLURM array launcher

```{literalinclude} ../../dwi-00_launch_slurm_container.sh 
:language: bash
:linenos:
```
