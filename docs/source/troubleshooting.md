# Troubleshooting

---

````{dropdown} Test DropDown
This should collapse




```{admonition} `Usage: dwi-preproc|dwi-tracto <specfile>` and immediate exit
:class: dropdown

You didn't supply both the pipeline name and the specfile path as the two
arguments to the entry point. The correct call is:

```bash
dwi-preproc /spec/spec.json
```
```



````{dropdown} `Error: Variable '<x>' is not set or is empty.`
A required key is missing or empty in `spec.json`. Regenerate the file
from `helpers/create_spec_template.sh` and check that every field is
filled in and that key names exactly match (case-sensitive).
````


````{dropdown} `no dwi scan/bvec found for <subj> - <session>`
BIDS naming mismatch, or the bind-mounted `bidsdir` inside the container
doesn't point to the path you think it does. Double-check your
`-v`/`--bind` mounts against the paths written in `spec.json`.
````

````{dropdown} `no TotalReadOutTime or PhaseEncodingDirection found in dwi json file`
The DWI `.json` sidecar is missing required BIDS metadata fields. Fix
this at the source/BIDS-conversion stage before re-running.
````

````{dropdown} `PE directions of dwi and fmap are NOT opposites`
The field map(s) in `fmap/` don't have a phase-encoding direction
opposite the DWI acquisition, or the `IntendedFor` linking is wrong.
Check both the PE direction metadata in the fmap sidecars and the
`IntendedFor` field pointing to the DWI file.
````

````{dropdown} Pipeline falls back to synb0 unexpectedly
No usable field maps were found under `fmap/` (check `IntendedFor`
tagging). This is expected fallback behavior if you don't have field
maps, but confirm it's intentional if you do.

```{note}
B0 scans are **not** supported. The pipeline expects Phase Encoding
POLARity (PEPOLAR / blip-up blip-down) field maps.
```
````


````{dropdown} `FreeSurfer throws an error related to /scratch/.tmp`
On some systems FreeSurfer tries to write to `/scratch` on the host and
gets a permission-denied error. The fix is to bind a writable host
directory to `/scratch` inside the container and set the relevant
`TMPDIR` environment variables.

Apptainer example:

```bash
apptainer run \
  --bind ${bidsdir}:/bids,${workdir}:/work,${outputdir}:/derivatives,${tmpdir_job}:/scratch \
  --env TMPDIR=/scratch \
  --env TMP=/scratch \
  --env TEMP=/scratch \
  tractoprep.sif ...
```
```{note}
see for further examples the launch scripts under {doc}`workflow`
```

````

````{dropdown} Output verification fails, workdir not cleaned up
One or more expected output files (see {doc}`output-structure`) weren't
produced. Check the per-step log files under
`<outputdir>/dwi-preproc/<subj>/log/` for the actual tool error, then
re-run once fixed.

Already-completed subjects and steps are skipped on re-run — see
[Idempotency](#idempotency) below.
````

````{dropdown} SELinux permission denied on bind mounts (Podman on RHEL/Fedora)
Add the `:Z` suffix (or `:z` for shared mounts) to your `-v` mount
specs:

```bash
-v /host/path:/container/path:Z
```
````

````{dropdown} Container can't write to output directories
- **Docker / Podman:** run with `--user $(id -u):$(id -g)`.
- **Apptainer:** check that the host directories are writable by the UID
  Apptainer runs as (typically your own, since Apptainer runs
  unprivileged by default).
````

---

## Idempotency

At each stage, the pipeline checks for existing expected output files
before re-running a step, so re-launching the container on a subject
that's already (partially) processed will skip completed steps rather
than redoing them. This makes it safe to re-run after fixing an error
partway through a subject.
