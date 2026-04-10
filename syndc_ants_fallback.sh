#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------
# ants_syn_susc_fallback.sh
# ANTs SyN-based susceptibility distortion fallback for DWI (BIDS)
#
# Usage:
#   ./ants_syn_susc_fallback.sh <sub-id> <${bidsdir}> [DERIV_DIR]
#
# Example:
#   ./ants_syn_susc_fallback.sh 01 /data/myproject/BIDS /data/myproject/derivatives
#
# Requirements:
#   - ANTs (antsRegistration, antsApplyTransforms)
#   - FSL (fslroi, fslmaths, flirt)
#   - ANTs + Python (nibabel & numpy) for PE-projection step (optional)
#
# Behavior:
#   * Extract first b0 from DWI (or mean of all b0s)
#   * N4 + brain-extract T1 -> T1_brain
#   * Brain-extract b0 -> b0_brain
#   * antsRegistration: fixed=T1_brain, moving=b0_brain, do Rigid->Affine->SyN
#   * Optionally project warp to PE axis and apply to DWI
#   * Output: derivatives/ants_syn_susc/sub-<ID>_dwi_dc.nii.gz
# ------------------------------------------------------------------


# Initialize variables
workdir=""
subj=""
session=
# input variables

# Parse command line arguments
while getopts "w:s:z:" opt; do
    case $opt in
        w) workdir="$OPTARG" ;;
        s) subj="$OPTARG" ;;
        z) session="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done
missing=0
for var in  workdir subj ; do
    if [[ -z "${!var}" ]]; then
        echo "Error: -${var:0:1} ($var) is required."
        missing=1
    fi
done
if [[ $missing -eq 1 ]]; then
    Usage
fi

# Set session path/file
if [[ -z "${session}" ]]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${session}/"
    sessionfile="_${session}_"
fi

PE_AXIS="${PE_AXIS:-y}"         # set environment var PE_AXIS=x|y|z if different
USE_PE_PROJECTION="${USE_PE_PROJECTION:-true}"  # set to false to skip projection

workdir="${DERIV_DIR}/${subj}"
mkdir -p "${workdir}"

echo "Subject: ${SUB}"
echo "BIDS: ${bidsdir}"
echo "Output: ${workdir}"
echo "Phase-encoding axis (PE_AXIS): ${PE_AXIS}"
echo "PE projection enabled: ${USE_PE_PROJECTION}"


# -------------------------
# 4) Run ANTs registration (fixed = T1_brain, moving = b0_brain)
#    Rigid -> Affine -> SyN
# -------------------------
echo "[4] Running antsRegistration (SyN)"
ANTS_PREFIX="${workdir}/ants_b0_to_T1_"
antsRegistration \
  --dimensionality 3 \
  --float 0 \
  --output [${ANTS_PREFIX},${ANTS_PREFIX}Warped.nii.gz] \
  --interpolation Linear \
  --winsorize-image-intensities [0.005,0.995] \
  --use-histogram-matching 0 \
  --initial-moving-transform [${T1_brain},${B0_brain},1] \
  \
  --transform Rigid[0.1] \
  --metric MI[${T1_brain},${B0_brain},1,32,Regular,0.25] \
  --convergence [1000x500x250,1e-6,10] \
  --shrink-factors 8x4x2 \
  --smoothing-sigmas 3x2x1vox \
  \
  --transform Affine[0.1] \
  --metric MI[${T1_brain},${B0_brain},1,32,Regular,0.25] \
  --convergence [1000x500,1e-6,10] \
  --shrink-factors 4x2 \
  --smoothing-sigmas 2x1vox \
  \
  --transform SyN[0.1,3,0] \
  --metric CC[${T1_brain},${B0_brain},1,4] \
  --convergence [200x200x50,1e-6,10] \
  --shrink-factors 4x2x1 \
  --smoothing-sigmas 2x1x0vox

# ants outputs:
# ${ANTS_PREFIX}0GenericAffine.mat   (affine)
# ${ANTS_PREFIX}1Warp.nii.gz         (displacement field image - mapping moving->fixed)
# (${ANTS_PREFIX}Warped.nii.gz is warped moving -> fixed)

# -------------------------
# 5) Optionally project warp to PE axis (reduce deformation to PE direction)
#    This uses a small Python step to zero other vector components of 1Warp.nii.gz.
# -------------------------
WARP="${ANTS_PREFIX}1Warp.nii.gz"
WARP_PROJECTED="${workdir}/1Warp_proj_${PE_AXIS}.nii.gz"

if [ "${USE_PE_PROJECTION}" = "true" ] ; then
  echo "[5] Projecting warp to PE axis (${PE_AXIS})"
  python3 - <<PY
import nibabel as nb, numpy as np, sys
wfile = "${WARP}"
out = "${WARP_PROJECTED}"
pe = "${PE_AXIS}".lower()
img = nb.load(wfile)
data = img.get_fdata()
# data shape: (X,Y,Z,3) -- vector displacement
if data.ndim != 4 or data.shape[3] < 3:
    print("Expected vector displacement field with 3 components; aborting projection.")
    sys.exit(1)
# index of axis
axmap = {'x':0,'y':1,'z':2}
if pe not in axmap:
    print("Invalid PE axis. Use x,y or z.")
    sys.exit(1)
keep = axmap[pe]
proj = np.zeros_like(data)
proj[..., keep] = data[..., keep]   # keep only PE component
nb.Nifti1Image(proj, img.affine, img.header).to_filename(out)
print("Wrote projected warp:", out)
PY
  # set warp to projected version
  WARP_USED="${WARP_PROJECTED}"
else
  echo "[5] Skipping warp projection"
  WARP_USED="${WARP}"
fi

# -------------------------
# 6) Apply transforms to full DWI (transform maps moving -> fixed (b0 -> T1))
#    We first apply affine+warp to the DWI, resampling to T1 space (undistorted).
# -------------------------
echo "[6] Applying transforms to full DWI (this will resample DWI into T1 space)"

# Compose transforms in the order antsApplyTransforms expects: last transform first on cmdline.
# Use the (possibly projected) warp, then the affine
affine="${ANTS_PREFIX}0GenericAffine.mat"
OUT_DWI="${workdir}/${subj}_dwi_dc_ants_syn.nii.gz"

antsApplyTransforms -d 3 \
  -e 0 \
  -i "${DWI}" \
  -r "${T1}" \
  -o "${OUT_DWI}" \
  -n Linear \
  -t "${WARP_USED}" \
  -t "${affine}"

echo "Wrote corrected DWI: ${OUT_DWI}"

# -------------------------
# 7) Rotate bvecs using affine component (approximate)
#    Warning: Nonlinear warp rotations are not applied to bvecs here.
# -------------------------
echo "[7] Rotating bvecs by the affine part (approximate; nonlinear rotations ignored)"

if [ -f "${BVEC}" ] ; then
  # convert affine to text matrix via ants to FSL? Use ants to get matrix is fine; we'll parse it.
  # Use c3d_affine_tool if available, otherwise apply a simple rotation extraction (works for rigid+affine if no shear)
  if command -v c3d_affine_tool >/dev/null 2>&1 ; then
    c3d_affine_tool -ref "${T1}" -src "${DWI}" -itk "${affine}" -o "${workdir}/affine_for_bvecs.txt"
    # then use MRtrix or in-house rotation; here we'll fallback to a simple approach using flirt to get rotation matrix
  fi

  # Fallback method: use flirt to get matrix in FSL format from the affine
  flirt -in "${DWI}" -ref "${T1}" -applyxfm -init "${affine}" -out /dev/null -omat "${workdir}/affine_fsl.mat" || true

  # Python rotate bvecs by rotation component of affine (approximate)
  python3 - <<PY
import numpy as np, sys
bv_in="${BVEC}"
bv_out="${workdir}/${subj}_dwi_dc.bvec"
matfile="${workdir}/affine_fsl.mat"
# load bvecs (3 x N or N x 3)
bv = np.loadtxt(bv_in)
if bv.shape[0] == 3:
    bv = bv.copy()
else:
    bv = bv.T
# load affine matrix
M = np.loadtxt(matfile)
if M.shape != (4,4):
    # try reading as 4x4 with header; fallback to identity
    print("Couldn't read 4x4 affine; writing original bvecs unchanged")
    np.savetxt(bv_out, bv, fmt='%.6f')
    sys.exit(0)
# extract rotation+scaling from upper-left 3x3
R = M[:3,:3]
# remove scaling by normalizing column vectors (approx)
rot = np.zeros((3,3))
for i in range(3):
    col = R[:,i]
    norm = np.linalg.norm(col)
    if norm == 0:
        rot[:,i] = col
    else:
        rot[:,i] = col / norm
# rotate bvecs
bv_rot = rot.dot(bv)
# normalize each vector
for j in range(bv_rot.shape[1]):
    n = np.linalg.norm(bv_rot[:,j])
    if n > 0:
        bv_rot[:,j] /= n
# write in FSL-style (3 x N)
np.savetxt(bv_out, bv_rot, fmt='%.6f')
print("Wrote rotated bvecs:", bv_out)
PY
else
  echo "No bvec file found; skipping rotation step."
fi

echo "Done. Outputs are in ${workdir}"
echo "Caveats: see script header and notes below."

