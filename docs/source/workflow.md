# Typical Workflow

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
7. **Run `dwi-qc`** (again) to update qc html pages.
8. **Run `dwi-send`** to export the completed derivatives.
9. Use **`dwi-params`** at any stage to pull out acquisition/processing
   parameters for logging or a shared summary spreadsheet.

```{seealso}
{doc}`running` for details on each pipeline step, and
{doc}`troubleshooting` if a step fails.
```


## Pipeline launcher
Below are engine specific scripts to launch the entire workflow for your entire sample (in the bids directory) using a SLURM array or without. 
At the top of the script you can adjust how many participants are processed in parallel.
Based on testing we provide best guess settings under the **Per-stage resource settings** section — adjust where necessary.

```{note}
Although it is recommended to first run **dwi-qc** before continuing to **dwi-tracto** 
in these pipeline launchers dwi-preproc and dwi-tracto are run back to back. 
```

::::::{tab-set}

:::::{tab-item} SLURM

::::{tab-set}

:::{tab-item} Docker
```{literalinclude} ../../dwi-00_launch_slurm_docker.sh
:language: bash
:linenos:
```
:::

:::{tab-item} Podman
```{literalinclude} ../../dwi-00_launch_slurm_podman.sh
:language: bash
:linenos:
```
:::

:::{tab-item} Apptainer
```{literalinclude} ../../dwi-00_launch_slurm_apptainer.sh
:language: bash
:linenos:
```
:::

::::

:::::

:::::{tab-item} NO-SLURM

::::{tab-set}

:::{tab-item} Apptainer
```{literalinclude} ../../dwi-00_launch_noslurm.sh
:language: bash
:linenos:
```
:::

:::{tab-item} Docker

code here

:::

:::{tab-item} Podman

code here

:::

::::

:::::

::::::


