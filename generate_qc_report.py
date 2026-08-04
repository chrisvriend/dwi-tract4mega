#!/usr/bin/env python3
"""
Generate QC HTML report for DWI preprocessing outputs.

Currently implements:
  - Noise map section (dwidenoise output)
  - Eddy QC section (eddy_quad)
  - Topup / susceptibility distortion correction QC
  - Brainmask QC
  - Response function voxel selection QC (dwi2response -voxels)
  - T1–DWI coregistration QC with optional 5tt2vis overlay
  - Tractogram QC (whole-brain .tck overlaid on T1, colored by fibre orientation)
  - Atlas connectivity matrix QC (normalized tck2connectome CSV heatmap)

"""

import argparse
import base64
import json
import re
from io import BytesIO
from pathlib import Path

import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import LineCollection
from PIL import Image
from scipy.ndimage import median_filter

DARK_BG = "#1b1f24"
MUTED = "#9aa4af"
GRID = "#2a2f36"
ACCENT = "#4fa3ff"
ACCENT2 = "#ff9f4f"

# --------------------------------------------------------------------------
# Image helpers
# --------------------------------------------------------------------------

def load_volume(path):
    img = nib.load(str(path))
    data = img.get_fdata()
    if data.ndim == 4:
        data = data[..., 0]
    return data, img.header


def load_volume_no_header(path):
    """Original behavior where header is not needed."""
    img = nib.load(str(path))
    data = img.get_fdata()
    if data.ndim == 4:
        data = data[..., 0]
    return data


def get_voxel_sizes_from_header(header):
    """Return (dx, dy, dz) from NIfTI header, fall back to 1.0 if missing."""
    zooms = header.get_zooms()
    if len(zooms) < 3:
        return 1.0, 1.0, 1.0
    dx, dy, dz = zooms[:3]
    dx = float(dx) if dx > 0 else 1.0
    dy = float(dy) if dy > 0 else 1.0
    dz = float(dz) if dz > 0 else 1.0
    return dx, dy, dz


def _row_col_mm(axis, header):
    """
    Physical mm-per-pixel spacing of (row, col) as slices are ACTUALLY
    plotted, i.e. after the np.rot90() that every slice/MIP extraction
    function applies before imshow. np.rot90 transposes rows/cols, so
    this is the post-rotation mapping, not the raw array-slicing one.

    axis: 0=sagittal, 1=coronal, 2=axial
      Axial   : raw slice is (X, Y) -> after rot90: rows=Y (dy), cols=X (dx)
      Coronal : raw slice is (X, Z) -> after rot90: rows=Z (dz), cols=X (dx)
      Sagittal: raw slice is (Y, Z) -> after rot90: rows=Z (dz), cols=Y (dy)
    """
    dx, dy, dz = get_voxel_sizes_from_header(header)
    if axis == 2:
        return dy, dx
    elif axis == 1:
        return dz, dx
    else:
        return dz, dy


def slice_display_dims(shape3d, header, axis):
    """
    Physical (width_mm, height_mm) of a slice from `axis`, as displayed
    after np.rot90(). axis: 0=sagittal, 1=coronal, 2=axial.
    """
    nx, ny, nz = shape3d[:3]
    row_mm, col_mm = _row_col_mm(axis, header)
    if axis == 2:
        rows, cols = ny, nx
    elif axis == 1:
        rows, cols = nz, nx
    else:
        rows, cols = nz, ny
    return cols * col_mm, rows * row_mm  # (width_mm, height_mm)


def make_triplanar_axes(shape3d, header, facecolor="black", target_max_dim=8.0,
                         gap=0.15, right_margin=0.0, top_margin=0.45,
                         bottom_margin=0.0):
    """
    Build a figure with 3 axes (axial, coronal, sagittal), laid out left to
    right, all sharing ONE physical mm-to-inch scale. This is important:
    if the acquisition FOV is much shorter in Z than in X/Y (e.g. a DWI
    slab covering 30-60 slices vs a 256x256 in-plane matrix), the coronal
    and sagittal views are genuinely shorter in physical extent. Sizing
    each panel independently to a shared *height* (as an earlier version
    of this function did) forces those short-FOV panels to become very
    WIDE to preserve their aspect ratio, which then dominates the overall
    figure and makes the axial panel look tiny once everything is
    downscaled to a fixed output width. Using one shared scale for all
    three panels instead keeps every panel's size proportional to its
    true physical size relative to the others -- axial stays full-size,
    and coronal/sagittal come out correctly smaller (not wider).

    Each axes is sized exactly to its own panel (no numeric 'aspect' is
    used for the box itself) and vertically centered in the row, so
    imshow(..., aspect='auto') fills the box exactly with no stretching
    and no cropping. Any leftover space around a shorter panel is just
    blank figure background -- which is now the *correct*, minimal
    amount, reflecting real anatomy rather than a layout bug.

    Returns (fig, [ax_axial, ax_coronal, ax_sagittal]).
    """
    views = [2, 1, 0]  # axial, coronal, sagittal -- order used throughout
    dims = [slice_display_dims(shape3d, header, axis) for axis in views]
    max_extent = max(max(w, h) for w, h in dims) or 1.0
    scale = target_max_dim / max_extent
    sizes_in = [(w * scale, h * scale) for w, h in dims]

    max_h = max(h for _, h in sizes_in)
    fig_h = max_h + top_margin + bottom_margin
    n_gaps = len(sizes_in) - 1
    total_w = sum(w for w, _ in sizes_in) + n_gaps * gap + right_margin

    fig = plt.figure(figsize=(total_w, fig_h), facecolor=facecolor)
    axes = []
    x_cursor = 0.0
    for w, h in sizes_in:
        left = x_cursor / total_w
        bottom = (bottom_margin + (max_h - h) / 2.0) / fig_h
        width_frac = w / total_w
        height_frac = h / fig_h
        ax = fig.add_axes([left, bottom, width_frac, height_frac])
        ax.set_facecolor(facecolor)
        axes.append(ax)
        x_cursor += w + gap

    return fig, axes



def panel_figsize_from_slice(slice2d, axis, header, target_height=4.0):
    """
    Figure size (width, height) in inches for a SINGLE already-extracted
    (post-rot90) 2D slice, so its displayed aspect matches physical mm,
    for use with aspect='auto'.
    """
    row_mm, col_mm = _row_col_mm(axis, header)
    rows, cols = slice2d.shape[:2]
    height_mm = rows * row_mm
    width_mm = cols * col_mm
    width = target_height * (width_mm / height_mm) if height_mm else target_height
    return width, target_height


def make_mosaic(data, header, axis, n_slices=7, cmap="gray", vmin=None, vmax=None):
    """Mosaic of evenly spaced slices along `axis` (0=sagittal, 1=coronal, 2=axial)."""
    n = data.shape[axis]
    idxs = np.linspace(int(n * 0.15), int(n * 0.85), n_slices).astype(int)

    width_mm, height_mm = slice_display_dims(data.shape, header, axis)
    panel_height = 4.0
    panel_width = panel_height * (width_mm / height_mm) if height_mm else panel_height

    fig, axes = plt.subplots(1, n_slices, figsize=(panel_width * n_slices, panel_height),
                              facecolor="black")
    if n_slices == 1:
        axes = [axes]
    for ax_, idx in zip(axes, idxs):
        if axis == 0:
            sl = data[idx, :, :]
        elif axis == 1:
            sl = data[:, idx, :]
        else:
            sl = data[:, :, idx]
        sl = np.rot90(sl)
        ax_.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax_.axis("off")
    fig.subplots_adjust(wspace=0.02, hspace=0, left=0, right=1, top=1, bottom=0)
    return fig


def optimize_png_bytes(png_bytes, max_width=1000):
    """Downscale (if wider than max_width) and PNG-optimize image bytes."""
    im = Image.open(BytesIO(png_bytes))
    if im.mode not in ("RGB", "RGBA", "L"):
        im = im.convert("RGBA")
    if im.width > max_width:
        new_h = int(im.height * (max_width / im.width))
        im = im.resize((max_width, new_h), Image.Resampling.LANCZOS)
    out = BytesIO()
    im.save(out, format="PNG", optimize=True, compress_level=9)
    return out.getvalue()


def fig_to_base64(fig, max_width=1000):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0.05,
                facecolor=fig.get_facecolor(), dpi=100)
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    return base64.b64encode(optimized).decode("utf-8")

# --------------------------------------------------------------------------
# QC sections
# --------------------------------------------------------------------------

def noise_section(noise_path):
    data, header = load_volume(noise_path)
    positive = data[data > 0]
    vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

    imgs = {}
    for axis, name in zip([2, 1, 0], ["axial", "coronal", "sagittal"]):
        fig = make_mosaic(data, header, axis, cmap="viridis", vmin=vmin, vmax=vmax)
        imgs[name] = fig_to_base64(fig)

    return f"""
    <section id="noise" class="qc-section">
      <h2>Denoising &mdash; Noise Map</h2>
      <p class="qc-desc">Spatial distribution of the noise level estimated by
      <code>dwidenoise</code> (MP-PCA).</p>
      <div class="slice-block">
        <h3>Axial</h3>
        <img src="data:image/png;base64,{imgs['axial']}" class="mosaic" alt="axial noise map"/>
      </div>
      <div class="slice-block">
        <h3>Coronal</h3>
        <img src="data:image/png;base64,{imgs['coronal']}" class="mosaic" alt="coronal noise map"/>
      </div>
      <div class="slice-block">
        <h3>Sagittal</h3>
        <img src="data:image/png;base64,{imgs['sagittal']}" class="mosaic" alt="sagittal noise map"/>
      </div>
    </section>
    """

def load_eddy_qc_json(path):
    with open(path) as f:
        return json.load(f)

def load_movement_rms(path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]

_OUTLIER_RE = re.compile(
    r"Slice (\d+) in scan (\d+) is an outlier with mean ([-\d.]+) standard "
    r"deviations off, and mean squared ([-\d.]+) standard deviations off\."
)

def parse_outlier_report(path):
    outliers = []
    if not path or not Path(path).exists():
        return outliers
    text = Path(path).read_text().strip()
    for line in text.splitlines():
        m = _OUTLIER_RE.match(line.strip())
        if m:
            slice_idx, scan_idx, mean_dev, mean_sq_dev = m.groups()
            outliers.append({
                "slice": int(slice_idx),
                "scan": int(scan_idx),
                "mean_dev": float(mean_dev),
                "mean_sq_dev": float(mean_sq_dev),
            })
    return outliers

def _find_qc_image(qc_dir, filename):
    p = Path(qc_dir) / filename
    return p if p.exists() else None

def _img_file_to_base64(path, max_width=1000):
    optimized = optimize_png_bytes(Path(path).read_bytes(), max_width=max_width)
    return base64.b64encode(optimized).decode("utf-8")

def collect_eddyqc_images(qc_dir, bvals, skip_b0_snr=False):
    qc_dir = Path(qc_dir)
    items = []

    p = _find_qc_image(qc_dir, "avg_b0.png")
    if p:
        items.append(("Average b0 (all PE directions)", p,
            "Average of all b=0 volumes after eddy/topup correction."))

    p = _find_qc_image(qc_dir, "avg_b0_pe0.png")
    if p:
        items.append(("Average b0 (phase-encoding direction 0)", p,
            "Average b=0 volume for the primary PE direction."))

    for bval in bvals:
        p = _find_qc_image(qc_dir, f"avg_b{int(bval)}.png")
        if p:
            items.append((f"Average b={int(bval)} shell", p,
                f"Average of all diffusion-weighted volumes in the b={int(bval)} shell."))

    cnr_labels = ["b0 SNR"] + [f"b={int(b)} CNR" for b in bvals]
    for i, label in enumerate(cnr_labels):
        if i == 0 and skip_b0_snr:
            continue
        p = _find_qc_image(qc_dir, f"cnr{i:04d}.nii.gz.png")
        if p:
            if i == 0:
                desc = "Voxel-wise b0 SNR map."
            else:
                desc = f"Voxel-wise CNR map for b={int(bvals[i-1])}."
            items.append((f"{label} map", p, desc))

    return items

def eddyqc_images_html(items, extra_prefix_html=""):
    if not items and not extra_prefix_html:
        return ""
    blocks = "\n".join(f"""
      <div class="slice-block">
        <h3>{title}</h3>
        <p class="qc-desc">{desc}</p>
        <img src="data:image/png;base64,{_img_file_to_base64(path)}" class="mosaic" alt="{title}"/>
      </div>
    """ for title, path, desc in items)
    return f"""
    <h3 class="subsection-title">eddy_quad Summary Images</h3>
    {extra_prefix_html}
    {blocks}
    """

def eddyqc_acknowledgment(eddy_input=None):
    repol_note = ""
    if eddy_input and str(eddy_input.get("repol", "False")).lower() == "true":
        repol_note = """
        <p>Outlier replacement (<code>--repol</code>) was used.</p>
        """
    return f"""
    <div class="ack-box">
      <strong>Acknowledgment</strong>
      <p>QC metrics and summary images were generated with FSL's
      <code>eddy_quad</code> (EDDY QC toolbox).</p>
      {repol_note}
    </div>
    """

def _style_axes(ax):
    ax.set_facecolor(DARK_BG)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.grid(alpha=0.15, color=MUTED)

def plot_motion_rms(abs_rms, rel_rms):
    n = len(abs_rms)
    fig, ax = plt.subplots(figsize=(10, 3), facecolor=DARK_BG)
    x = np.arange(n)
    ax.plot(x, abs_rms, color=ACCENT, lw=1.5, label="Absolute RMS (mm)")
    ax.plot(x, rel_rms, color=ACCENT2, lw=1.5, label="Relative RMS (mm)")
    ax.set_xlabel("Volume")
    ax.set_ylabel("Displacement (mm)")
    ax.set_xlim(0, max(n - 1, 1))
    _style_axes(ax)
    leg = ax.legend(facecolor=DARK_BG, edgecolor=GRID, fontsize=8, loc="upper left")
    for text in leg.get_texts():
        text.set_color("#e6e6e6")
    fig.tight_layout()
    return fig

def plot_outlier_scatter(outliers, n_vols):
    fig, ax = plt.subplots(figsize=(10, 3), facecolor=DARK_BG)
    if outliers:
        scans = [o["scan"] for o in outliers]
        slices = [o["slice"] for o in outliers]
        severity = np.abs([o["mean_sq_dev"] for o in outliers])
        sizes = np.clip(severity * 12, 25, 300)
        sc = ax.scatter(scans, slices, c=severity, cmap="Reds", s=sizes,
                         edgecolors="white", linewidths=0.4, vmin=0)
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label("Mean sq. stdev off", color=MUTED)
        cb.ax.yaxis.set_tick_params(color=MUTED)
        plt.setp(cb.ax.get_yticklabels(), color=MUTED)
        cb.outline.set_edgecolor(GRID)
    else:
        ax.text(0.5, 0.5, "No outlier slices detected", color=MUTED,
                 ha="center", va="center", transform=ax.transAxes)
    ax.set_xlim(-1, max(n_vols, 1))
    ax.set_xlabel("Volume (scan index)")
    ax.set_ylabel("Slice index")
    _style_axes(ax)
    fig.tight_layout()
    return fig

def load_4d_volume(path):
    img = nib.load(str(path))
    return img.get_fdata(), img.header

def load_4d_volume_no_header(path):
    return nib.load(str(path)).get_fdata()

def get_vol(data, idx):
    return data[..., idx] if data.ndim == 4 else data

def extract_slice(vol3d, axis, frac=0.5):
    """axis: 0=sagittal, 1=coronal, 2=axial."""
    n = vol3d.shape[axis]
    idx = int(n * frac)
    if axis == 0:
        sl = vol3d[idx, :, :]
    elif axis == 1:
        sl = vol3d[:, idx, :]
    else:
        sl = vol3d[:, :, idx]
    return np.rot90(sl)

def mip_slice(vol3d, axis):
    sl = vol3d.max(axis=axis)
    return np.rot90(sl)


def make_triplanar_mip_multimask_data_uri(underlay_vol, underlay_header,
                                           masks, colors, vmin, vmax,
                                           alpha=0.6, max_width=1000):
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]
    fig, axes = make_triplanar_axes(underlay_vol.shape, underlay_header)
    for ax_, (label, axis) in zip(axes, views):
        base_sl = mip_slice(underlay_vol, axis)
        ax_.imshow(base_sl, cmap="gray", vmin=vmin, vmax=vmax, aspect="auto")
        for mask_vol, color in zip(masks, colors):
            msl = mip_slice(mask_vol, axis)
            msl = np.ma.masked_less_equal(msl, 0.5)
            ax_.imshow(msl, cmap=ListedColormap([color]), vmin=0, vmax=1,
                       alpha=alpha, aspect="auto")
        ax_.axis("off")
    buf = BytesIO()

    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def slice_to_data_uri(slice2d, axis, header, vmin=None, vmax=None,
                      cmap="gray", max_width=500):
    fig_w, fig_h = panel_figsize_from_slice(slice2d, axis, header, target_height=4.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="black")
    ax.imshow(slice2d, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def slice_overlay_to_data_uri(base_slice, overlay_slice, axis, header,
                               vmin, vmax, overlay_absmax,
                               overlay_cmap="bwr", alpha=0.55,
                               max_width=500):
    fig_w, fig_h = panel_figsize_from_slice(base_slice, axis, header, target_height=4.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="black")
    ax.imshow(base_slice, cmap="gray", vmin=vmin, vmax=vmax, aspect="auto")
    ax.imshow(overlay_slice, cmap=overlay_cmap, vmin=-overlay_absmax,
              vmax=overlay_absmax, alpha=alpha, aspect="auto")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def make_triplanar_data_uri(vol3d, header, vmin, vmax, overlay_vol=None,
                             overlay_absmax=None, cmap="gray",
                             overlay_cmap="bwr", alpha=0.55,
                             max_width=1000, overlay_vmin=None,
                             overlay_mask_zero=False):
    """Axial, coronal, sagittal mid-slices side by side."""
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]
    fig, axes = make_triplanar_axes(vol3d.shape, header)
    for ax_, (label, axis) in zip(axes, views):
        sl = extract_slice(vol3d, axis)
        ax_.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        if overlay_vol is not None:
            osl = extract_slice(overlay_vol, axis)
            if overlay_mask_zero:
                osl = np.ma.masked_less_equal(osl, 0)
            ov_vmin = -overlay_absmax if overlay_vmin is None else overlay_vmin
            ax_.imshow(osl, cmap=overlay_cmap, vmin=ov_vmin,
                       vmax=overlay_absmax, alpha=alpha, aspect="auto")
        ax_.axis("off")
        ax_.set_title(label, color=MUTED, fontsize=11)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def make_triplanar_contour_data_uri(base_vol, base_header, mask_vol,
                                     vmin, vmax, contour_color="#ff4f4f",
                                     max_width=1000):
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]
    fig, axes = make_triplanar_axes(base_vol.shape, base_header)
    for ax_, (label, axis) in zip(axes, views):
        base_sl = extract_slice(base_vol, axis)
        mask_sl = extract_slice(mask_vol, axis)
        ax_.imshow(base_sl, cmap="gray", vmin=vmin, vmax=vmax, aspect="auto")
        if np.any(mask_sl > 0.5):
            ax_.contour(mask_sl, levels=[0.5], colors=[contour_color], linewidths=1.3)
        ax_.axis("off")
        ax_.set_title(label, color=MUTED, fontsize=11)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def brainmask_section(nodif_path, mask_path):
    nodif_img = nib.load(str(nodif_path))
    nodif = nodif_img.get_fdata()
    if nodif.ndim == 4:
        nodif = nodif.mean(axis=-1)
    mask = nib.load(str(mask_path)).get_fdata()

    positive = nodif[nodif > 0]
    vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

    img_uri = make_triplanar_contour_data_uri(
        nodif, nodif_img.header, mask, vmin, vmax
    )

    mask_voxels = int(np.sum(mask > 0.5))
    total_voxels = int(mask.size)
    coverage_pct = 100 * mask_voxels / total_voxels if total_voxels else 0

    return f"""
    <section id="brainmask" class="qc-section">
      <h2>Brain Mask</h2>
      <p class="qc-desc">Brain mask boundary overlaid on nodif (b0).</p>
      <div class="stat-grid">
        <div class="stat-card">
          <div class="stat-label">Mask Volume</div>
          <div class="stat-value">{mask_voxels:,} vox</div>
        </div>
        <div class="stat-card">
          <div class="stat-label">Mask Coverage</div>
          <div class="stat-value">{coverage_pct:.1f}%</div>
        </div>
      </div>
      <div class="slice-block">
        <img src="{img_uri}" class="mosaic" alt="brain mask overlay on nodif, axial/coronal/sagittal"/>
      </div>
    </section>
    """

# --------------------------------------------------------------------------
# Response voxels QC
# --------------------------------------------------------------------------

RESPONSE_TISSUE_LABELS = ["WM", "GM", "CSF"]
RESPONSE_TISSUE_COLORS = ["#4fa3ff", "#ff9f4f", "#e14fff"]

def response_voxels_section(voxels_path, underlay_path, underlay_label="nodif (b0)"):
    voxels_img = nib.load(str(voxels_path))
    voxels_data = voxels_img.get_fdata()

    if voxels_data.ndim != 4 or voxels_data.shape[-1] != 3:
        return f"""
        <section id="responsevoxels" class="qc-section">
          <h2>Response Function Voxel Selection (dwi2response)</h2>
          <p class="qc-desc">Could not render: expected shape (..., 3) but got {voxels_data.shape}.</p>
        </section>
        """

    underlay_img = nib.load(str(underlay_path))
    underlay = underlay_img.get_fdata()
    if underlay.ndim == 4:
        underlay = underlay.mean(axis=-1)

    if underlay.shape != voxels_data.shape[:3]:
        return f"""
        <section id="responsevoxels" class="qc-section">
          <h2>Response Function Voxel Selection (dwi2response)</h2>
          <p class="qc-desc">Voxel selection mask shape {voxels_data.shape[:3]} "
          and underlay shape {underlay.shape} differ.</p>
        </section>
        """

    positive = underlay[underlay > 0]
    vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

    masks = [voxels_data[..., i] for i in range(3)]
    counts = [int(np.sum(m > 0.5)) for m in masks]

    draw_order = [2, 1, 0]
    img_uri = make_triplanar_mip_multimask_data_uri(
        underlay, underlay_img.header,
        [masks[i] for i in draw_order],
        [RESPONSE_TISSUE_COLORS[i] for i in draw_order],
        vmin, vmax,
    )

    legend_items = "".join(
        f'<div class="legend-item"><span class="legend-swatch" '
        f'style="background:{color};"></span>{label} response voxels '
        f'&mdash; {count:,} vox</div>'
        for label, color, count in zip(RESPONSE_TISSUE_LABELS, RESPONSE_TISSUE_COLORS, counts)
    )

    return f"""
    <section id="responsevoxels" class="qc-section">
      <h2>Response Function Voxel Selection (dwi2response)</h2>
      <p class="qc-desc">Voxels selected by <code>dwi2response ... -voxels</code> for each tissue.</p>
      <div class="legend-row">{legend_items}</div>
      <div class="slice-block">
        <img src="{img_uri}" class="mosaic"
             alt="response voxels overlaid on {underlay_label}, axial/coronal/sagittal MIP"/>
      </div>
    </section>
    """

# --------------------------------------------------------------------------
# T1–DWI coregistration QC
# --------------------------------------------------------------------------

def coreg_section(t1w_dwi_path, nodif_path, five_tt_vis_path=None):
    t1_img = nib.load(str(t1w_dwi_path))
    t1 = t1_img.get_fdata()
    nodif_img = nib.load(str(nodif_path))
    nodif = nodif_img.get_fdata()

    if t1.ndim == 4:
        t1 = t1.mean(axis=-1)
    if nodif.ndim == 4:
        nodif = nodif.mean(axis=-1)

    if t1.shape != nodif.shape:
        return f"""
        <section id="coreg" class="qc-section">
          <h2>T1–DWI Coregistration</h2>
          <p class="qc-desc">T1w and nodif shapes differ: {t1.shape} vs {nodif.shape}.</p>
        </section>
        """

    t1_positive = t1[t1 > 0]
    nodif_positive = nodif[nodif > 0]

    t1_vmin, t1_vmax = (
        np.percentile(t1_positive, [2, 98]) if t1_positive.size else (None, None)
    )
    nodif_vmin, nodif_vmax = (
        np.percentile(nodif_positive, [1, 99]) if nodif_positive.size else (None, None)
    )

    t1_uri = make_triplanar_data_uri(t1, t1_img.header, t1_vmin, t1_vmax, cmap="gray")
    nodif_uri = make_triplanar_data_uri(nodif, nodif_img.header, nodif_vmin, nodif_vmax, cmap="magma")

    has_overlay = five_tt_vis_path is not None

    if has_overlay:
        five_tt_img = nib.load(str(five_tt_vis_path))
        five_tt = five_tt_img.get_fdata()
        if five_tt.ndim == 4:
            five_tt = five_tt.mean(axis=-1)

        if five_tt.shape != t1.shape:
            has_overlay = False
            overlay_expl = f"""
            <p class="qc-desc">
              5tt2vis overlay shape ({five_tt.shape}) does not match T1/nodif ({t1.shape}), overlay disabled.
            </p>
            """
            t1_overlay_uri = t1_uri
            nodif_overlay_uri = nodif_uri
        else:
            pos = five_tt[five_tt > 0]
            ov_absmax = float(np.percentile(pos, 99)) if pos.size else 1.0

            t1_overlay_uri = make_triplanar_data_uri(
                t1, t1_img.header,
                t1_vmin, t1_vmax,
                overlay_vol=five_tt,
                overlay_absmax=ov_absmax,
                overlay_vmin=0,
                overlay_mask_zero=True,
                cmap="gray",
                overlay_cmap="Set1",
                alpha=1.0,
            )
            nodif_overlay_uri = make_triplanar_data_uri(
                nodif, nodif_img.header,
                nodif_vmin, nodif_vmax,
                overlay_vol=five_tt,
                overlay_absmax=ov_absmax,
                overlay_vmin=0,
                overlay_mask_zero=True,
                cmap="magma",
                overlay_cmap="tab10",
                alpha=0.6,
            )

            overlay_expl = """
            <p class="qc-desc">
              5tt2vis tissue-type overlay is drawn on top of the current base image.
            </p>
            """
    else:
        overlay_expl = ""
        t1_overlay_uri = t1_uri
        nodif_overlay_uri = nodif_uri

    view_data = {
        "t1_plain": t1_uri,
        "nodif_plain": nodif_uri,
        "t1_overlay": t1_overlay_uri,
        "nodif_overlay": nodif_overlay_uri,
        "has_overlay": has_overlay,
        "t1_label": "T1w in DWI space",
        "nodif_label": "Nodif (b0)",
    }

    overlay_control = ""
    if has_overlay:
        overlay_control = """
          <label class="overlay-toggle">
            <input type="checkbox" id="coreg-overlay" onchange="updateCoregView()">
            Show 5tt2vis tissue-type overlay
          </label>
        """

    block = f"""
    <div class="slice-block">
      <div class="coreg-stack">
        <img id="coreg-t1-img" class="mosaic slider-img coreg-base" src="{t1_uri}"
             alt="T1w in DWI space, axial/coronal/sagittal"/>
        <img id="coreg-nodif-img" class="mosaic slider-img coreg-overlay-img" src="{nodif_uri}"
             alt="Nodif (b0, DWI space), axial/coronal/sagittal" style="opacity: 0;">
      </div>
      <div class="topup-slider-row">
        <span class="topup-slider-endlabel">T1w</span>
        <input type="range" min="0" max="1" step="0.1" value="0" class="topup-range"
               id="coreg-range" oninput="updateCoregView()">
        <span class="topup-slider-endlabel">Nodif</span>
      </div>
      <div class="slider-label" id="coreg-label">T1w in DWI space</div>
    </div>
    """

    script = f"""
    <script>
      window.coregData = {json.dumps(view_data)};
      function updateCoregView() {{
        var d = window.coregData;
        var alpha = parseFloat(document.getElementById('coreg-range').value);
        var overlayCb = document.getElementById('coreg-overlay');
        var showOverlay = overlayCb && overlayCb.checked && d.has_overlay;

        var t1Img = document.getElementById('coreg-t1-img');
        var nodifImg = document.getElementById('coreg-nodif-img');
        var labelEl = document.getElementById('coreg-label');

        if (showOverlay) {{
          t1Img.src = d.t1_overlay;
          nodifImg.src = d.nodif_overlay;
        }} else {{
          t1Img.src = d.t1_plain;
          nodifImg.src = d.nodif_plain;
        }}

        nodifImg.style.opacity = alpha.toString();
        t1Img.style.opacity = (1.0 - alpha).toString();

        var baseLabel = alpha < 0.5 ? d.t1_label : d.nodif_label;
        if (showOverlay) {{
          labelEl.innerHTML = baseLabel + " + 5tt2vis overlay";
        }} else {{
          labelEl.innerHTML = baseLabel;
        }}
      }}
    </script>
    """

    return f"""
    <section id="coreg" class="qc-section">
      <h2>T1–DWI Coregistration</h2>
      <p class="qc-desc">
        Slider transitions between T1w and nodif, across axial/coronal/sagittal views.
      </p>
      {overlay_expl}
      {overlay_control}
      {block}
      {script}
    </section>
    """

# --------------------------------------------------------------------------
# Tractogram QC
# --------------------------------------------------------------------------

def voxel_to_plot_xy(vox_pts, axis, shape):
    if axis == 0:      # sagittal: row=j, col=k
        r, c, C = vox_pts[:, 1], vox_pts[:, 2], shape[2]
    elif axis == 1:     # coronal: row=i, col=k
        r, c, C = vox_pts[:, 0], vox_pts[:, 2], shape[2]
    else:               # axial: row=i, col=j
        r, c, C = vox_pts[:, 0], vox_pts[:, 1], shape[1]
    return r, (C - 1) - c

def build_tract_view_segments(streamlines, inv_affine, axis, shape, max_points_per_line=40):
    seg_list, color_list = [], []
    for sl in streamlines:
        pts = np.asarray(sl)
        if len(pts) < 2:
            continue
        if max_points_per_line and len(pts) > max_points_per_line:
            keep = np.linspace(0, len(pts) - 1, max_points_per_line).astype(int)
            pts = pts[keep]

        diffs = np.diff(pts, axis=0)
        norms = np.linalg.norm(diffs, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        dirs = np.abs(diffs / norms)

        vox = nib.affines.apply_affine(inv_affine, pts)
        x, y = voxel_to_plot_xy(vox, axis, shape)
        segs = np.stack(
            [np.column_stack([x[:-1], y[:-1]]), np.column_stack([x[1:], y[1:]])], axis=1
        )
        seg_list.append(segs)
        color_list.append(dirs)

    if not seg_list:
        return np.zeros((0, 2, 2)), np.zeros((0, 3))
    return np.concatenate(seg_list, axis=0), np.concatenate(color_list, axis=0)

def make_tractography_triplanar_data_uri(t1_data, t1_header, streamlines,
                                          inv_affine, max_points_per_line=40,
                                          linewidth=0.35, alpha=0.75,
                                          max_width=1000):
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]
    positive = t1_data[t1_data > 0]
    vmin, vmax = np.percentile(positive, [2, 98]) if positive.size else (None, None)
    shape = t1_data.shape

    fig, axes = make_triplanar_axes(t1_data.shape, t1_header)
    for ax_, (label, axis) in zip(axes, views):
        base_sl = mip_slice(t1_data, axis)
        ax_.imshow(base_sl, cmap="gray", vmin=vmin, vmax=vmax, aspect="auto")
        xlim, ylim = ax_.get_xlim(), ax_.get_ylim()

        segs, colors = build_tract_view_segments(
            streamlines, inv_affine, axis, shape, max_points_per_line=max_points_per_line
        )
        if len(segs):
            ax_.add_collection(LineCollection(segs, colors=colors,
                                              linewidths=linewidth, alpha=alpha))

        ax_.set_xlim(xlim)
        ax_.set_ylim(ylim)
        ax_.axis("off")
        ax_.set_title(label, color=MUTED, fontsize=11)

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def tractography_section(tck_path, t1_path, max_streamlines=6000,
                         max_points_per_line=40, seed=0):
    t1_img = nib.load(str(t1_path))
    t1_data = t1_img.get_fdata()
    if t1_data.ndim == 4:
        t1_data = t1_data.mean(axis=-1)

    try:
        tfile = nib.streamlines.load(str(tck_path))
    except Exception as exc:
        return f"""
        <section id="tractography" class="qc-section">
          <h2>Tractogram</h2>
          <p class="qc-desc">Could not read tractogram: {exc}</p>
        </section>
        """

    streamlines_full = tfile.streamlines
    n_total = len(streamlines_full)

    rng = np.random.default_rng(seed)
    if n_total > max_streamlines:
        keep_idx = rng.choice(n_total, size=max_streamlines, replace=False)
        streamlines = [streamlines_full[i] for i in keep_idx]
    else:
        streamlines = list(streamlines_full)

    img_uri = make_tractography_triplanar_data_uri(
        t1_data, t1_img.header, streamlines, np.linalg.inv(t1_img.affine),
        max_points_per_line=max_points_per_line,
    )

    subsample_note = ""
    if n_total > max_streamlines:
        subsample_note = f"""
        <p class="qc-desc">Showing a random subsample of {max_streamlines:,} of
        {n_total:,} streamlines for rendering performance.</p>
        """

    return f"""
    <section id="tractography" class="qc-section">
      <h2>Tractogram</h2>
      <p class="qc-desc">Whole-brain tractogram overlaid on T1w in DWI space.</p>
      {subsample_note}
      <div class="slice-block">
        <img src="{img_uri}" class="mosaic"
             alt="tractogram overlaid on T1, axial/coronal/sagittal MIP"/>
      </div>
    </section>
    """

# --------------------------------------------------------------------------
# Topup QC
# --------------------------------------------------------------------------

def read_acqparams(path):
    rows = []
    for line in Path(path).read_text().strip().splitlines():
        parts = line.split()
        if len(parts) >= 3:
            rows.append(tuple(float(x) for x in parts[:3]))
    return rows

def pick_pe_volume_indices(acq_rows):
    seen = {}
    for i, pe in enumerate(acq_rows):
        if pe not in seen:
            seen[pe] = i
    return sorted(seen.items(), key=lambda kv: kv[1])

def find_matching_pe_index(topup_pe_items, dwi_acqparams_path):
    if not dwi_acqparams_path:
        vec, idx = topup_pe_items[0]
        return vec, idx, False
    dwi_rows = read_acqparams(dwi_acqparams_path)
    if not dwi_rows:
        vec, idx = topup_pe_items[0]
        return vec, idx, False
    dwi_pe = dwi_rows[0]
    for vec, idx in topup_pe_items:
        if vec == dwi_pe:
            return vec, idx, True
    vec, idx = topup_pe_items[0]
    return vec, idx, False

def topup_section(before_path, after_path, topup_acqparams_path,
                   fieldmap_path=None, dwi_acqparams_path=None,
                   matched_pe_index=None):
    topup_acq_rows = read_acqparams(topup_acqparams_path)
    pe_items = pick_pe_volume_indices(topup_acq_rows)

    if not pe_items:
        return f"""
        <section id="topup" class="qc-section">
          <h2>Susceptibility Distortion Correction (topup)</h2>
          <p class="qc-desc">No phase-encode rows parsed from {Path(topup_acqparams_path).name}.</p>
        </section>
        """

    if matched_pe_index is not None:
        pe_vec = topup_acq_rows[matched_pe_index]
        pe_idx = matched_pe_index
        matched = True
    else:
        pe_vec, pe_idx, matched = find_matching_pe_index(pe_items, dwi_acqparams_path)

    def pe_label(vec):
        return f"({vec[0]:g}, {vec[1]:g}, {vec[2]:g})"

    match_note = ""
    if not matched:
        match_note = f"""
        <p class="qc-desc" style="color:#e0a54e;">DWI phase-encode direction could not be confirmed;
        showing first PE direction {pe_label(pe_vec)}.</p>
        """

    before_data, before_header = load_4d_volume(before_path)
    after_data, after_header = load_4d_volume(after_path)
    before_vol = get_vol(before_data, pe_idx)
    after_vol = get_vol(after_data, pe_idx)

    both_vals = np.concatenate([before_vol.flatten(), after_vol.flatten()])
    positive = both_vals[both_vals > 0]
    vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

    fieldmap_data = None
    overlay_absmax = None
    if fieldmap_path and Path(fieldmap_path).exists():
        fm_img = nib.load(str(fieldmap_path))
        fieldmap_data = fm_img.get_fdata()
        overlay_absmax = np.percentile(np.abs(fieldmap_data), 99)

    before_label = f"Before correction &mdash; PE {pe_label(pe_vec)}"
    after_label = f"After correction &mdash; PE {pe_label(pe_vec)}"

    before_uri = make_triplanar_data_uri(before_vol, before_header, vmin, vmax)
    after_plain_uri = make_triplanar_data_uri(after_vol, after_header, vmin, vmax)

    has_overlay = fieldmap_data is not None
    if has_overlay:
        # assume same grid/header as after_vol
        after_overlay_uri = make_triplanar_data_uri(
            after_vol, after_header, vmin, vmax,
            overlay_vol=fieldmap_data, overlay_absmax=overlay_absmax
        )
    else:
        after_overlay_uri = after_plain_uri

    view_data = {
        "before": before_uri,
        "after_plain": after_plain_uri,
        "after_overlay": after_overlay_uri,
        "before_label": before_label,
        "after_label": after_label,
        "has_overlay": has_overlay,
    }

    global_overlay_control = f"""
      <label class="overlay-toggle">
        <input type="checkbox" id="topup-overlay" onchange="updateTopupView()">
        Show distortion overlay (applies to corrected image only)
      </label>
    """ if has_overlay else ""

    block = f"""
    <div class="slice-block">
      <img id="topup-img" class="mosaic slider-img" src="{before_uri}"
           alt="axial, coronal, sagittal topup before/after comparison"/>
      <div class="vtoggle-row">
        <span class="vtoggle-label-top">Before</span>
        <label class="vtoggle">
          <input type="checkbox" id="topup-range" onchange="updateTopupView()">
          <span class="vtoggle-track"><span class="vtoggle-thumb"></span></span>
        </label>
        <span class="vtoggle-label-bottom">After</span>
      </div>
      <div class="slider-label" id="topup-label">{before_label}</div>
    </div>
    """

    script = f"""
    <script>
      window.topupData = {json.dumps(view_data)};
      function updateTopupView() {{
        var d = window.topupData;
        var mode = document.getElementById('topup-range').checked ? 'after' : 'before';
        var overlayCb = document.getElementById('topup-overlay');
        var img = document.getElementById('topup-img');
        var label = document.getElementById('topup-label');
        if (mode === 'before') {{
          img.src = d.before;
          label.innerHTML = d.before_label;
        }} else {{
          var showOverlay = overlayCb && overlayCb.checked && d.has_overlay;
          img.src = showOverlay ? d.after_overlay : d.after_plain;
          label.innerHTML = d.after_label + (showOverlay ? ' + distortion overlay' : '');
        }}
      }}
    </script>
    """

    overlay_expl = ""
    if has_overlay:
        overlay_expl = """
        <p class="qc-desc">Distortion overlay shows topup's off-resonance field (Hz)
        on the corrected image.</p>
        """

    return f"""
    <section id="topup" class="qc-section">
      <h2>Susceptibility Distortion Correction (topup)</h2>
      <p class="qc-desc">Compare DWI phase-encode direction before and after topup.</p>
      {match_note}
      {overlay_expl}
      {global_overlay_control}
      {block}
      {script}
    </section>
    """

# --------------------------------------------------------------------------
# SNR map helpers
# --------------------------------------------------------------------------

def compute_snr_map(cnr_maps_path):
    img = nib.load(str(cnr_maps_path))
    data = img.get_fdata()
    snr_vol = data[..., 0] if data.ndim == 4 else data
    snr_vol = median_filter(snr_vol, size=3)
    return snr_vol, img.header

def make_triplanar_colorbar_data_uri(vol3d, header, cmap="viridis",
                                      vmin=0, vmax=None,
                                      cbar_label="", max_width=1000):
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]

    # Reserve extra figure width to the right of the three panels for a
    # colorbar that's a comfortable, fixed size regardless of how wide/
    # narrow the brain panels end up being. Layout of that reserved strip
    # (left to right): gap -> colorbar -> tick labels/axis label -> edge.
    gap_to_cbar = 0.35
    cbar_width_in = 0.40
    label_space_in = 1.45
    right_margin = gap_to_cbar + cbar_width_in + label_space_in

    fig, axes = make_triplanar_axes(vol3d.shape, header, right_margin=right_margin)
    fig_w, fig_h = fig.get_size_inches()

    im = None
    for ax_, (label, axis) in zip(axes, views):
        sl = extract_slice(vol3d, axis)
        im = ax_.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax_.axis("off")
        ax_.set_title(label, color=MUTED, fontsize=11)

    panels_right_in = fig_w - right_margin
    cbar_left_in = panels_right_in + gap_to_cbar
    cbar_height_frac = 0.75
    cbar_bottom_frac = (1 - cbar_height_frac) / 2

    cbar_ax = fig.add_axes([cbar_left_in / fig_w, cbar_bottom_frac,
                             cbar_width_in / fig_w, cbar_height_frac])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(cbar_label, color=MUTED)
    cbar.ax.yaxis.set_tick_params(color=MUTED)
    plt.setp(cbar.ax.get_yticklabels(), color=MUTED)
    cbar.outline.set_edgecolor(GRID)

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def snr_map_block(cnr_maps_path):
    snr_vol, header = compute_snr_map(cnr_maps_path)
    positive = snr_vol[snr_vol > 0]
    vmax = float(np.percentile(positive, 99)) if positive.size else None

    img_uri = make_triplanar_colorbar_data_uri(
        snr_vol, header, cmap="viridis", vmin=0, vmax=vmax,
        cbar_label="SNR (mean / std)"
    )

    return f"""
    <div class="slice-block">
      <h3>SNR Map (b0)</h3>
      <p class="qc-desc">Voxel-wise b0 SNR.</p>
      <img src="{img_uri}" class="mosaic" alt="SNR map, axial/coronal/sagittal, with colorbar"/>
    </div>
    """

def eddyqc_section(json_path, rms_path, outlier_path=None, qc_dir=None,
                    raw_dwi_path=None, preproc_dwi_path=None,
                    cnr_maps_path=None, mot_abs_thresh=1.0,
                    mot_rel_thresh=0.5, outlier_pct_thresh=5.0):
    qc = load_eddy_qc_json(json_path)
    abs_rms, rel_rms = load_movement_rms(rms_path)
    outliers = parse_outlier_report(outlier_path)
    n_vols = len(abs_rms)
    qc_dir = qc_dir or Path(json_path).parent

    motion_b64 = fig_to_base64(plot_motion_rms(abs_rms, rel_rms))
    outlier_b64 = fig_to_base64(plot_outlier_scatter(outliers, n_vols))

    mot_abs = qc.get("qc_mot_abs")
    mot_rel = qc.get("qc_mot_rel")
    outlier_pct = qc.get("qc_outliers_tot", 0) * 100
    bvals = qc.get("data_unique_bvals", [])
    cnr_avg = qc.get("qc_cnr_avg", [])

    def flag(value, threshold):
        if value is None:
            return ""
        return "stat-bad" if value > threshold else "stat-ok"

    stat_cards = f"""
      <div class="stat-card {flag(mot_abs, mot_abs_thresh)}">
        <div class="stat-label">Mean Abs. Motion</div>
        <div class="stat-value">{mot_abs:.2f} mm</div>
      </div>
      <div class="stat-card {flag(mot_rel, mot_rel_thresh)}">
        <div class="stat-label">Mean Rel. Motion</div>
        <div class="stat-value">{mot_rel:.2f} mm</div>
      </div>
      <div class="stat-card {flag(outlier_pct, outlier_pct_thresh)}">
        <div class="stat-label">Outlier Slices</div>
        <div class="stat-value">{outlier_pct:.1f}%</div>
      </div>
      <div class="stat-card">
        <div class="stat-label">Volumes (b0 / dwi)</div>
        <div class="stat-value">{qc.get('data_no_b0_vols', '?')} / {qc.get('data_no_dw_vols', '?')}</div>
      </div>
    """

    cnr_cards = ""
    if cnr_avg:
        labels = ["b0"] + [f"b={b:g}" for b in bvals]
        cnr_cards = "".join(
            f'<div class="stat-card"><div class="stat-label">{lbl} SNR/CNR</div>'
            f'<div class="stat-value">{val:.2f}</div></div>'
            for lbl, val in zip(labels, cnr_avg)
        )

    outlier_raw = ""
    if outlier_path and Path(outlier_path).exists():
        raw_text = Path(outlier_path).read_text().strip()
        n_lines = len(raw_text.splitlines()) if raw_text else 0
        outlier_raw = f"""
        <details class="raw-details">
          <summary>Raw eddy_outlier_report ({n_lines} entries)</summary>
          <pre>{raw_text if raw_text else "No outliers reported."}</pre>
        </details>
        """

    snr_html = snr_map_block(cnr_maps_path) if cnr_maps_path else ""
    image_items = collect_eddyqc_images(qc_dir, bvals, skip_b0_snr=bool(cnr_maps_path))
    images_html = eddyqc_images_html(image_items, extra_prefix_html=snr_html)
    outlier_volumes_html = outlier_volumes_block(raw_dwi_path, preproc_dwi_path, outliers)
    ack_html = eddyqc_acknowledgment(qc.get("eddy_input"))

    return f"""
    <section id="eddyqc" class="qc-section">
      <h2>Eddy Current &amp; Motion Correction (eddy_quad)</h2>
      <p class="qc-desc">Summary statistics and motion/outlier diagnostics from <code>eddy_quad</code>.</p>

      <div class="stat-grid">
        {stat_cards}
        {cnr_cards}
      </div>

      <div class="slice-block">
        <h3>Volume-to-volume Motion</h3>
        <img src="data:image/png;base64,{motion_b64}" class="mosaic" alt="motion rms plot"/>
      </div>

      <div class="slice-block">
        <h3>Outlier Slices</h3>
        <img src="data:image/png;base64,{outlier_b64}" class="mosaic" alt="outlier scatter plot"/>
      </div>

      {outlier_raw}

      {outlier_volumes_html}

      {images_html}

      {ack_html}
    </section>
    """

def outlier_volumes_block(raw_dwi_path, preproc_dwi_path, outliers):
    if not outliers:
        return ""

    if not raw_dwi_path or not preproc_dwi_path:
        return """
        <div id="outliervol" class="subsection-anchor">
        <h3 class="subsection-title">Outlier Volume Inspection</h3>
        <p class="qc-desc">Outlier slices were reported, but raw/preprocessed DWI volumes were not provided.</p>
        </div>
        """

    raw_data, raw_header = load_4d_volume(raw_dwi_path)
    proc_data, proc_header = load_4d_volume(preproc_dwi_path)

    by_scan = {}
    for o in outliers:
        by_scan.setdefault(o["scan"], []).append(o)

    blocks = []
    all_view_data = {}
    for scan_idx, entries in sorted(by_scan.items()):
        raw_vol = get_vol(raw_data, scan_idx)
        proc_vol = get_vol(proc_data, scan_idx)

        both_vals = np.concatenate([raw_vol.flatten(), proc_vol.flatten()])
        positive = both_vals[both_vals > 0]
        vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

        before_uri = make_triplanar_data_uri(raw_vol, raw_header, vmin, vmax)
        after_uri = make_triplanar_data_uri(proc_vol, proc_header, vmin, vmax)

        view_id = f"outliervol-{scan_idx}"
        all_view_data[view_id] = {"before": before_uri, "after": after_uri}

        slice_list = ", ".join(str(e["slice"]) for e in sorted(entries, key=lambda e: e["slice"]))
        worst = max(abs(e["mean_sq_dev"]) for e in entries)

        blocks.append(f"""
        <div class="slice-block">
          <h3>Volume {scan_idx}</h3>
          <p class="qc-desc">Flagged slice(s): {slice_list} &middot; worst mean-sq. deviation: {worst:.2f}</p>
          <img id="{view_id}-img" class="mosaic slider-img" src="{before_uri}"
               alt="volume {scan_idx} raw vs eddy-processed comparison"/>
          <div class="vtoggle-row">
            <span class="vtoggle-label-top">Raw</span>
            <label class="vtoggle">
              <input type="checkbox" id="{view_id}-range" onchange="updateOutlierVol('{view_id}')">
              <span class="vtoggle-track"><span class="vtoggle-thumb"></span></span>
            </label>
            <span class="vtoggle-label-bottom">Eddy-processed</span>
          </div>
        </div>
        """)

    script = f"""
    <script>
      window.outlierVolData = Object.assign(window.outlierVolData || {{}}, {json.dumps(all_view_data)});
      function updateOutlierVol(id) {{
        var d = window.outlierVolData[id];
        var mode = document.getElementById(id + '-range').checked ? 'after' : 'before';
        document.getElementById(id + '-img').src = d[mode];
      }}
    </script>
    """

    return f"""
    <div id="outliervol" class="subsection-anchor">
    <h3 class="subsection-title">Outlier Volume Inspection</h3>
    <p class="qc-desc">Compare raw and eddy-processed volumes for each outlier scan.</p>
    {"".join(blocks)}
    {script}
    </div>
    """

# --------------------------------------------------------------------------
# Connectivity matrix
# --------------------------------------------------------------------------

def connectivity_matrix_section(csv_path, atlas_name=None):
    if not csv_path:
        return ""

    csv_path = Path(csv_path)
    if not csv_path.exists():
        return f"""
        <section id="connectivity" class="qc-section">
          <h2>Connectivity Matrix</h2>
          <p class="qc-desc">Connectivity matrix file not found: <code>{csv_path}</code></p>
        </section>
        """

    try:
        raw = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    except Exception:
        raw = np.genfromtxt(csv_path, delimiter=",", skip_header=1, filling_values=0)

    if raw.ndim != 2:
        raw = np.atleast_2d(raw)

    if raw.shape[1] == raw.shape[0] + 1:
        mat = raw[:, 1:]
    else:
        mat = raw

    mat = mat.astype(float)
    mat[mat < 0] = 0

    vals = mat[mat > 0]
    if vals.size:
        p95 = float(np.percentile(vals, 95))
        scale = p95 if p95 > 0 else float(vals.max())
        mat_norm = mat / scale
        mat_norm[mat_norm > 1] = 1.0
    else:
        mat_norm = mat

    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(mat_norm, cmap="viridis", interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Normalized streamline count (≈95th percentile = 1.0)", color=MUTED)
    cbar.ax.yaxis.set_tick_params(color=MUTED)
    plt.setp(cbar.ax.get_yticklabels(), color=MUTED)
    cbar.outline.set_edgecolor(GRID)

    ax.set_xlabel("Target node", color=MUTED)
    ax.set_ylabel("Source node", color=MUTED)

    title = f"Connectivity matrix ({atlas_name})" if atlas_name else "Connectivity matrix"
    ax.set_title(title, color="white")

    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_BG)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.tick_params(colors=MUTED, which="both")

    fig.tight_layout()

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=900)
    b64 = base64.b64encode(optimized).decode("utf-8")
    uri = f"data:image/png;base64,{b64}"

    atlas_html = f" for atlas <code>{atlas_name}</code>" if atlas_name else ""

    return f"""
    <section id="connectivity" class="qc-section">
      <h2>Connectivity Matrix{atlas_html}</h2>
      <p class="qc-desc">Normalized connectivity matrix derived from <code>tck2connectome</code>.</p>
      <div class="slice-block">
        <img src="{uri}" class="mosaic"
             alt="Normalized connectivity matrix heatmap"/>
      </div>
    </section>
    """

# --------------------------------------------------------------------------
# Page assembly
# --------------------------------------------------------------------------

PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>DWI QC Report{subject_suffix}</title>
<style>
  :root {{
    --bg: #111417;
    --panel: #1b1f24;
    --accent: #4fa3ff;
    --text: #e6e6e6;
    --muted: #9aa4af;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: var(--bg);
    color: var(--text);
  }}
  nav {{
    position: sticky;
    top: 0;
    z-index: 100;
    background: var(--panel);
    border-bottom: 1px solid #2a2f36;
    padding: 0.75rem 1.5rem;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    flex-wrap: wrap;
  }}
  nav .brand {{
    font-weight: 600;
    margin-right: 1.5rem;
    color: var(--accent);
  }}
  nav a {{
    color: var(--text);
    text-decoration: none;
    padding: 0.4rem 0.9rem;
    border-radius: 6px;
    font-size: 0.9rem;
    border: 1px solid transparent;
  }}
  nav a:hover {{
    background: #262b32;
    border-color: var(--accent);
  }}
  nav a.nav-sublink {{
    font-size: 0.8rem;
    color: var(--muted);
    padding: 0.4rem 0.7rem;
    margin-left: -0.3rem;
  }}
  nav a.nav-sublink::before {{
    content: "\\2192  ";
  }}
  main {{
    max-width: 1200px;
    margin: 0 auto;
    padding: 2rem 1.5rem 4rem;
  }}
  .qc-section {{
    background: var(--panel);
    border-radius: 10px;
    padding: 1.5rem 1.75rem;
    margin-bottom: 2rem;
    border: 1px solid #262b32;
    scroll-margin-top: 70px;
  }}
  .qc-section h2 {{
    margin-top: 0;
    color: var(--accent);
    border-bottom: 1px solid #262b32;
    padding-bottom: 0.5rem;
  }}
  .qc-desc {{
    color: var(--muted);
    font-size: 0.92rem;
    max-width: 900px;
  }}
  .slice-block h3 {{
    font-size: 0.85rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: var(--muted);
    margin-bottom: 0.4rem;
  }}
  .mosaic {{
    max-width: 100%;
    width: auto;
    max-height: 480px;
    display: block;
    margin-left: auto;
    margin-right: auto;
    border-radius: 6px;
    background: black;
    margin-bottom: 1.2rem;
  }}
  code {{
    background: #262b32;
    padding: 0.1rem 0.4rem;
    border-radius: 4px;
    font-size: 0.85rem;
  }}
  .stat-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 0.75rem;
    margin: 1rem 0 1.5rem;
  }}
  .stat-card {{
    background: #20242b;
    border: 1px solid #2a2f36;
    border-radius: 8px;
    padding: 0.9rem 1rem;
  }}
  .stat-label {{
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--muted);
    margin-bottom: 0.3rem;
  }}
  .stat-value {{
    font-size: 1.4rem;
    font-weight: 600;
  }}
  .stat-ok {{ border-left: 3px solid #4caf50; }}
  .stat-bad {{ border-left: 3px solid #e05a4e; }}
  .subsection-title {{
    margin-top: 2rem;
    padding-top: 1rem;
    border-top: 1px solid #262b32;
    font-size: 1rem;
    color: var(--text);
  }}
  .subsection-anchor {{
    scroll-margin-top: 70px;
  }}
  .ack-box {{
    margin-top: 1.75rem;
    background: #191d22;
    border-left: 3px solid var(--accent);
    border-radius: 6px;
    padding: 0.9rem 1.1rem;
    font-size: 0.8rem;
    color: var(--muted);
  }}
  .ack-box strong {{
    color: var(--text);
    display: block;
    margin-bottom: 0.4rem;
    font-size: 0.85rem;
  }}
  .ack-box p {{
    margin: 0.35rem 0;
  }}
  .raw-details {{
    margin-top: 1rem;
    color: var(--muted);
  }}
  .raw-details summary {{
    cursor: pointer;
    color: var(--accent);
    font-size: 0.85rem;
  }}
  .raw-details pre {{
    background: #14171b;
    padding: 0.75rem 1rem;
    border-radius: 6px;
    font-size: 0.78rem;
    overflow-x: auto;
    margin-top: 0.5rem;
    color: #c8ccd1;
    white-space: pre-wrap;
  }}
  .topup-range {{
    width: 100%;
    margin: 0.5rem 0 0.3rem;
    accent-color: var(--accent);
  }}
  .slider-label {{
    font-size: 0.85rem;
    color: var(--text);
    text-align: center;
    background: #20242b;
    border-radius: 6px;
    padding: 0.35rem 0.5rem;
    margin-bottom: 1.2rem;
  }}
  .slider-img {{
    margin-bottom: 0.4rem;
  }}
  .topup-slider-row {{
    display: flex;
    align-items: center;
    gap: 0.6rem;
    margin-bottom: 0.5rem;
  }}
  .topup-slider-endlabel {{
    font-size: 0.75rem;
    color: var(--muted);
    white-space: nowrap;
  }}
  .vtoggle-row {{
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0.4rem;
    margin-bottom: 0.6rem;
  }}
  .vtoggle {{
    position: relative;
    display: inline-block;
    width: 28px;
    height: 52px;
    cursor: pointer;
  }}
  .vtoggle input {{
    position: absolute;
    opacity: 0;
    width: 0;
    height: 0;
  }}
  .vtoggle-track {{
    position: absolute;
    inset: 0;
    background: #20242b;
    border: 1px solid #2a2f36;
    border-radius: 14px;
    transition: background 0.2s ease;
  }}
  .vtoggle-thumb {{
    position: absolute;
    left: 3px;
    top: 3px;
    width: 20px;
    height: 20px;
    border-radius: 50%;
    background: var(--muted);
    transition: transform 0.2s ease, background 0.2s ease;
  }}
  .vtoggle input:checked + .vtoggle-track {{
    background: #28405c;
  }}
  .vtoggle input:checked + .vtoggle-track .vtoggle-thumb {{
    transform: translateY(24px);
    background: var(--accent);
  }}
  .vtoggle input:focus-visible + .vtoggle-track {{
    outline: 2px solid var(--accent);
    outline-offset: 2px;
  }}
  .vtoggle-label-top,
  .vtoggle-label-bottom {{
    font-size: 0.75rem;
    color: var(--muted);
    white-space: nowrap;
  }}
  .legend-row {{
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
    margin: 0.75rem 0 1.25rem;
  }}
  .legend-item {{
    display: flex;
    align-items: center;
    gap: 0.4rem;
    font-size: 0.85rem;
    color: var(--text);
  }}
  .legend-swatch {{
    width: 14px;
    height: 14px;
    border-radius: 3px;
    display: inline-block;
    border: 1px solid rgba(255, 255, 255, 0.25);
  }}
  .overlay-toggle {{
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.4rem;
    font-size: 0.85rem;
    color: var(--text);
    margin-bottom: 1.2rem;
    cursor: pointer;
    background: #20242b;
    border: 1px solid #2a2f36;
    border-radius: 6px;
    padding: 0.6rem 0.9rem;
  }}
  .coreg-stack {{
    position: relative;
    display: inline-block;
  }}
  .coreg-stack .coreg-overlay-img {{
    position: absolute;
    top: 0;
    left: 0;
    margin-bottom: 0;
  }}
</style>
</head>
<body>
<nav>
  <span class="brand">DWI QC{subject_suffix}</span>
  {nav_links}
</nav>
<main>
  {sections}
</main>
</body>
</html>
"""

def build_report(sections, output_path, subject=None, extra_nav=None):
    subject_suffix = f" &mdash; {subject}" if subject else ""
    extra_nav = extra_nav or []
    nav_parts = []
    for sid, label, _ in sections:
        nav_parts.append(f'<a href="#{sid}">{label}</a>')
        for after_id, anchor_id, sub_label in extra_nav:
            if after_id == sid:
                nav_parts.append(f'<a href="#{anchor_id}" class="nav-sublink">{sub_label}</a>')
    nav_links = "\n  ".join(nav_parts)
    body = "\n".join(html for _, _, html in sections)
    page = PAGE_TEMPLATE.format(subject_suffix=subject_suffix,
                                nav_links=nav_links, sections=body)
    Path(output_path).write_text(page)

def main():
    p = argparse.ArgumentParser(description="Generate DWI preprocessing QC HTML report")
    p.add_argument("--noise", default=None, help="Path to dwidenoise noise map")
    p.add_argument("--eddy-json", default=None, help="Path to eddy_quad qc.json")
    p.add_argument("--eddy-rms", default=None, help="Path to *.eddy_movement_rms")
    p.add_argument("--eddy-outliers", default=None, help="Path to *.eddy_outlier_report")
    p.add_argument("--eddy-qc-dir", default=None,
                    help="Directory containing eddy_quad summary PNGs")
    p.add_argument("--eddy-raw-dwi", default=None,
                    help="Raw (pre-eddy) 4D DWI volume")
    p.add_argument("--eddy-preproc-dwi", default=None,
                    help="Eddy-processed 4D DWI volume")
    p.add_argument("--topup-before", default=None,
                    help="Pre-correction 4D volume fed into topup")
    p.add_argument("--topup-after", default=None,
                    help="Post-correction volume from topup")
    p.add_argument("--topup-acqparams", default=None,
                    help="PE-direction table for topup input volumes")
    p.add_argument("--topup-fieldmap", default=None,
                    help="topup's --fout field map in Hz")
    p.add_argument("--topup-dwi-acqparams", default=None,
                    help="DWI acqparams.tsv to match PE direction")
    p.add_argument("--topup-match-pe-index", type=int, default=None,
                    help="Manual topup input volume index")
    p.add_argument("--brainmask-nodif", default=None,
                    help="Nodif (b0) reference image")
    p.add_argument("--brainmask-mask", default=None,
                    help="Brain mask NIfTI")
    p.add_argument("--eddy-cnr-maps", default=None,
                    help="eddy's CNR-maps 4D volume")

    p.add_argument("--reg-t1w-dwi", default=None,
                    help="T1w image in DWI space")
    p.add_argument("--reg-nodif", default=None,
                    help="Nodif resampled to same grid as --reg-t1w-dwi")
    p.add_argument("--reg-5ttvis", default=None,
                    help="5tt2vis image in same space as --reg-t1w-dwi")

    p.add_argument("--response-voxels", default=None,
                    help="4D voxel-selection image from dwi2response -voxels")
    p.add_argument("--response-underlay", default=None,
                    help="Underlay image for response voxels")

    p.add_argument("--tract-tck", default=None,
                    help="Tractogram (.tck) to visualize")
    p.add_argument("--tract-t1", default=None,
                    help="T1 image to overlay tractogram on")
    p.add_argument("--tract-max-streamlines", type=int, default=6000,
                    help="Max streamlines to render")

    p.add_argument("--connectivity-matrix", default=None,
                    help="tck2connectome CSV connectivity matrix")
    p.add_argument("--connectivity-atlas-name", default=None,
                    help="Atlas label for connectivity section")

    p.add_argument("--output", default="qc_report.html", help="Output HTML path")
    p.add_argument("--subject", default=None, help="Subject label, e.g. sub-01")
    args = p.parse_args()

    sections=[]
    
    if args.noise:

        sections = [
        ("noise", "Noise Map", noise_section(args.noise)),
        ]

    extra_nav = []

    if args.eddy_json and args.eddy_rms:
        eddyqc_html = eddyqc_section(
            args.eddy_json, args.eddy_rms, args.eddy_outliers,
            qc_dir=args.eddy_qc_dir,
            raw_dwi_path=args.eddy_raw_dwi,
            preproc_dwi_path=args.eddy_preproc_dwi,
            cnr_maps_path=args.eddy_cnr_maps,
        )
        sections.append(("eddyqc", "Eddy QC", eddyqc_html))
        if 'id="outliervol"' in eddyqc_html:
            extra_nav.append(("eddyqc", "outliervol", "Outlier Volumes"))

    if args.topup_before and args.topup_after and args.topup_acqparams:
        sections.append(
            ("topup", "Topup", topup_section(
                args.topup_before, args.topup_after, args.topup_acqparams,
                fieldmap_path=args.topup_fieldmap,
                dwi_acqparams_path=args.topup_dwi_acqparams,
                matched_pe_index=args.topup_match_pe_index,
            ))
        )

    if args.brainmask_nodif and args.brainmask_mask:
        sections.append(
            ("brainmask", "Brain Mask",
             brainmask_section(args.brainmask_nodif, args.brainmask_mask))
        )

    if args.response_voxels:
        underlay_path = args.response_underlay
        underlay_label = "response-underlay"
        if not underlay_path:
            if args.reg_nodif:
                underlay_path, underlay_label = args.reg_nodif, "nodif (b0)"
            elif args.brainmask_nodif:
                underlay_path, underlay_label = args.brainmask_nodif, "nodif (b0)"
            elif args.reg_t1w_dwi:
                underlay_path, underlay_label = args.reg_t1w_dwi, "T1w in DWI space"
        if underlay_path:
            sections.append(
                (
                    "responsevoxels",
                    "Response Voxels",
                    response_voxels_section(
                        args.response_voxels, underlay_path,
                        underlay_label=underlay_label
                    ),
                )
            )

    if args.reg_t1w_dwi and args.reg_nodif:
        sections.append(
            (
                "coreg",
                "T1–DWI Coregistration",
                coreg_section(args.reg_t1w_dwi,
                              args.reg_nodif,
                              five_tt_vis_path=args.reg_5ttvis),
            )
        )

    if args.tract_tck:
        tract_t1 = args.tract_t1 or args.reg_t1w_dwi
        if tract_t1:
            sections.append(
                (
                    "tractography",
                    "Tractogram",
                    tractography_section(
                        args.tract_tck, tract_t1,
                        max_streamlines=args.tract_max_streamlines
                    ),
                )
            )

    if args.connectivity_matrix:
        sections.append(
            (
                "connectivity",
                "Connectivity",
                connectivity_matrix_section(
                    args.connectivity_matrix,
                    atlas_name=args.connectivity_atlas_name,
                ),
            )
        )

    build_report(sections, args.output,
                 subject=args.subject, extra_nav=extra_nav)
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()