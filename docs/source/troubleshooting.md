# Troubleshooting

| Symptom | Likely cause + solution |
|---|---|
| `Usage: dwi-preproc\|dwi-tracto <specfile>` and immediate exit | You didn't supply both the pipeline name and the specfile path as the two arguments to the entry point. |
| `FreeSurfer throws an error related to /scratch/.tmp` | We noticed that on some systems that have a /scratch directory, FreeSurfer tries to write files here on the host system but gets a permission denied error. the fix is to define /scratch in the container as a tmpdir and bind to a location on the host that does have write permissions and define several TMPdir environmental variables. apptainer example: apptainer run --bind ${bidsdir}:/bids,${workdir}:/work,${outputdir}:/derivatives,${tmpdir_job}:/scratch --env TMPDIR=/scratch --env TMP=/scratch --env TEMP=/scratch tractoprep.sif ... |
| `Error: Variable '<x>' is not set or is empty.` | A required key is missing or empty in `spec.json` — regenerate from `create_spec_template.sh` and check every field is filled in, and that the key names exactly match (case-sensitive). |
| `no dwi scan/bvec found for <subj> - <session>` | BIDS naming mismatch, or the bind-mounted `bidsdir` inside the container doesn't point to the path you think it does — double check your `-v`/`--bind` mounts against the paths written in `spec.json`. |
| `no TotalReadOutTime or PhaseEncodingDirection found in dwi json file` | The DWI `.json` sidecar is missing required BIDS metadata fields — fix at the source/BIDS-conversion stage. |
| `PE directions of dwi and fmap are NOT opposites` | Field map(s) present in `fmap/` don't have a phase-encoding direction opposite (or matching, for the "same-PE" case) the DWI acquisition, or the `IntendedFor` linking is wrong. |
| Pipeline falls back to synb0 unexpectedly | No usable field maps were found under `fmap/` (check `IntendedFor` tagging) — this is expected fallback behavior if you don't have field maps, but confirm it's intentional if you do. Note that B0 scans are not supported; Phase Encoding POLARity (PEPOLAR) techniques (also called blip-up/blip-down) are what the pipeline expects. |
| Something went wrong with synb0 | The Synb0-DISCO step failed — check that the T1w N4-corrected/skull-stripped inputs were generated successfully and that enough memory was available (see the low-memory note in {doc}`spec-json`). |
| Output verification fails, workdir not cleaned up | One or more expected output files (see {doc}`output-structure`) weren't produced — check the per-step log files under `outputdir/dwi-preproc/<subj>/log/` for the actual tool error, then re-run once fixed. Already-completed subjects/steps are skipped on re-run. |
| SELinux permission denied on bind mounts (Podman on RHEL/Fedora-family hosts) | Add the `:Z` (or `:z` for shared mounts) suffix to your `-v` mount specs. |
| Container can't write to output directories | Run with `--user $(id -u):$(id -g)` (Docker/Podman), or check that the host directories are writable by whichever UID Apptainer runs as (typically your own, since Apptainer runs unprivileged by default). |

## Idempotency

At each stage, the pipeline checks for existing expected output files
before re-running a step, so re-launching the container on a subject
that's already (partially) processed will skip completed steps rather
than redoing them. This makes it safe to re-run after fixing an error
partway through a subject.
