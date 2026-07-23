#!/usr/bin/env python3
"""
Generate an fMRIPrep-style QC HTML report for DWI preprocessing outputs.

Currently implements:
  - Noise map section (dwidenoise output)
  - Eddy QC section (eddy_quad)
  - Topup / susceptibility distortion correction QC
  - Brainmask QC
  - T1–DWI coregistration QC with optional 5tt2vis overlay

Designed to be extended: add a new `..._section(path)` function that
returns an HTML string, and append it to the `sections` list in main().
No external network/CDN dependency, no display server required
(matplotlib runs headless via the Agg backend) -- safe to run inside
a container with no internet access.

Usage:
    python generate_qc_report.py --noise noise.nii.gz --output qc_report.html --subject sub-01
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
    return data

def make_mosaic(data, axis, n_slices=7, cmap="gray", vmin=None, vmax=None):
    """Mosaic of evenly spaced slices along `axis` (0=sagittal, 1=coronal, 2=axial)."""
    n = data.shape[axis]
    idxs = np.linspace(int(n * 0.15), int(n * 0.85), n_slices).astype(int)

    fig, axes = plt.subplots(1, n_slices, figsize=(2.2 * n_slices, 2.6), facecolor="black")
    for ax_, idx in zip(axes, idxs):
        if axis == 0:
            sl = data[idx, :, :]
        elif axis == 1:
            sl = data[:, idx, :]
        else:
            sl = data[:, :, idx]
        sl = np.rot90(sl)
        ax_.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax)
        ax_.axis("off")
    fig.subplots_adjust(wspace=0.02, hspace=0, left=0, right=1, top=1, bottom=0)
    return fig

def optimize_png_bytes(png_bytes, max_width=1000):
    """Downscale (if wider than max_width) and PNG-optimize image bytes.
    This is what keeps the self-contained HTML from ballooning in size --
    full-resolution QC images are overkill for on-screen review."""
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
# QC sections -- add new ones here following the same pattern
# --------------------------------------------------------------------------

def noise_section(noise_path):
    data = load_volume(noise_path)
    positive = data[data > 0]
    vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

    imgs = {}
    for axis, name in zip([2, 1, 0], ["axial", "coronal", "sagittal"]):
        fig = make_mosaic(data, axis, cmap="viridis", vmin=vmin, vmax=vmax)
        imgs[name] = fig_to_base64(fig)

    return f"""
    <section id="noise" class="qc-section">
      <h2>Denoising &mdash; Noise Map</h2>
      <p class="qc-desc">Spatial distribution of the noise level estimated by
      <code>dwidenoise</code> (MP-PCA). Structure here that resembles anatomy
      (rather than uniform background noise) can indicate loss of signal
      during denoising and is worth a closer look.</p>
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
    """temp.eddy_movement_rms: two unlabeled columns, one row per volume:
    absolute RMS displacement, relative RMS displacement (mm)."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]

_OUTLIER_RE = re.compile(
    r"Slice (\d+) in scan (\d+) is an outlier with mean ([-\d.]+) standard "
    r"deviations off, and mean squared ([-\d.]+) standard deviations off\."
)

def parse_outlier_report(path):
    """temp.eddy_outlier_report: one free-text line per outlier slice."""
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
    """Locate and describe the eddy_quad summary PNGs, if present in qc_dir.
    Returns a list of (title, filepath, description) tuples.
    skip_b0_snr: omit cnr0000.nii.gz.png (the static b0 SNR rendering) --
    used when a colorbar'd SNR map is derived directly from the CNR-maps
    NIfTI instead (see snr_map_block)."""
    qc_dir = Path(qc_dir)
    items = []

    p = _find_qc_image(qc_dir, "avg_b0.png")
    if p:
        items.append(("Average b0 (all PE directions)", p,
            "Average of all b=0 volumes after eddy/topup correction, combined across "
            "phase-encoding directions. This is the overall reference image for judging "
            "motion and distortion correction quality by eye."))

    p = _find_qc_image(qc_dir, "avg_b0_pe0.png")
    if p:
        items.append(("Average b0 (phase-encoding direction 0)", p,
            "Average b=0 volume restricted to the primary phase-encoding direction. Useful "
            "for checking that topup's susceptibility distortion correction is well aligned "
            "specifically for this PE direction."))

    for bval in bvals:
        p = _find_qc_image(qc_dir, f"avg_b{int(bval)}.png")
        if p:
            items.append((f"Average b={int(bval)} shell", p,
                f"Average of all diffusion-weighted volumes in the b={int(bval)} s/mm² "
                "shell after eddy current and motion correction."))

    # CNR maps follow eddy_quad's convention: index 0 is the b0 SNR map,
    # subsequent indices are the CNR map for each shell, in the order of data_unique_bvals.
    cnr_labels = ["b0 SNR"] + [f"b={int(b)} CNR" for b in bvals]
    for i, label in enumerate(cnr_labels):
        if i == 0 and skip_b0_snr:
            continue
        p = _find_qc_image(qc_dir, f"cnr{i:04d}.nii.gz.png")
        if p:
            if i == 0:
                desc = ("Voxel-wise signal-to-noise ratio map computed from the b0 volumes, "
                        "derived from eddy's predicted signal and the model residuals.")
            else:
                desc = (f"Voxel-wise contrast-to-noise ratio map for the b={int(bvals[i-1])} "
                        "s/mm² shell. Low CNR in white matter can indicate poor angular "
                        "contrast for downstream tractography or microstructure modelling.")
            items.append((f"{label} map", p, desc))

    # Note: vdm.png (voxel displacement map) intentionally omitted here -- the
    # topup section covers the same distortion information with an interactive
    # before/after comparison and overlay.

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
        <p>Outlier replacement (<code>--repol</code>) was used: Andersson, J.L.R., Graham, M.S.,
        Zsoldos, E., and Sotiropoulos, S.N. Incorporating outlier detection and replacement into
        a non-parametric framework for movement and distortion correction of diffusion MR images.
        NeuroImage, 141:556-572, 2016.</p>
        """
    return f"""
    <div class="ack-box">
      <strong>Acknowledgment</strong>
      <p>QC metrics and summary images in this section were generated with FSL's
      <code>eddy_quad</code> (part of the EDDY QC toolbox), developed by Matteo Bastiani and
      colleagues at FMRIB, Oxford. Please cite: Bastiani, M., Cottaar, M., Fitzgibbon, S.P.,
      Suri, S., Alfaro-Almagro, F., Sotiropoulos, S.N., Jbabdi, S., and Andersson, J.L.R.
      Automated quality control for within and between studies diffusion MRI data using a
      non-parametric framework for movement and distortion correction. NeuroImage,
      184:801-812, 2019.</p>
      <p>Underlying motion/distortion correction: Andersson, J.L.R. and Sotiropoulos, S.N.
      An integrated approach to correction for off-resonance effects and subject movement in
      diffusion MR imaging. NeuroImage, 125:1063-1078, 2016.</p>
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

def outlier_volumes_block(raw_dwi_path, preproc_dwi_path, outliers):
    """For each volume flagged in the outlier report, build a raw (pre-eddy) vs
    eddy-processed slider comparison. Uses the same triplanar (axial/coronal/
    sagittal) layout as the topup section -- coronal/sagittal panels are where
    slice-to-slice artifacts like 'Venetian blind' banding are most visible,
    since those planes cut across the full stack of axial slices."""
    if not outliers:
        return ""

    if not raw_dwi_path or not preproc_dwi_path:
        return """
        <div id="outliervol" class="subsection-anchor">
        <h3 class="subsection-title">Outlier Volume Inspection</h3>
        <p class="qc-desc">Outlier slices were reported, but the raw and eddy-processed
        DWI volumes were not provided (<code>--eddy-raw-dwi</code> /
        <code>--eddy-preproc-dwi</code>), so a visual before/after comparison could not
        be generated.</p>
        </div>
        """

    raw_data = load_4d_volume(raw_dwi_path)
    proc_data = load_4d_volume(preproc_dwi_path)

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

        before_uri = make_triplanar_data_uri(raw_vol, vmin, vmax)
        after_uri = make_triplanar_data_uri(proc_vol, vmin, vmax)

        view_id = f"outliervol-{scan_idx}"
        all_view_data[view_id] = {"before": before_uri, "after": after_uri}

        slice_list = ", ".join(str(e["slice"]) for e in sorted(entries, key=lambda e: e["slice"]))
        worst = max(abs(e["mean_sq_dev"]) for e in entries)

        blocks.append(f"""
        <div class="slice-block">
          <h3>Volume {scan_idx}</h3>
          <p class="qc-desc">Flagged slice(s): {slice_list} &middot; worst mean-sq.
          deviation: {worst:.2f}</p>
          <img id="{view_id}-img" class="mosaic slider-img" src="{before_uri}"
               alt="volume {scan_idx} raw vs eddy-processed comparison"/>
          <div class="topup-slider-row">
            <span class="topup-slider-endlabel">Raw</span>
            <input type="range" min="0" max="1" step="1" value="0" class="topup-range"
                   id="{view_id}-range" oninput="updateOutlierVol('{view_id}')">
            <span class="topup-slider-endlabel">Eddy-processed</span>
          </div>
        </div>
        """)

    script = f"""
    <script>
      window.outlierVolData = Object.assign(window.outlierVolData || {{}}, {json.dumps(all_view_data)});
      function updateOutlierVol(id) {{
        var d = window.outlierVolData[id];
        var mode = document.getElementById(id + '-range').value === '1' ? 'after' : 'before';
        document.getElementById(id + '-img').src = d[mode];
      }}
    </script>
    """

    return f"""
    <div id="outliervol" class="subsection-anchor">
    <h3 class="subsection-title">Outlier Volume Inspection</h3>
    <p class="qc-desc">For each volume flagged in the outlier report, drag the slider
    to compare the raw (pre-eddy) and eddy-processed image. Slice-to-volume or
    interpolation artifacts &mdash; such as "Venetian blind" banding across slices
    &mdash; are usually most visible in the coronal/sagittal panels, where alternating
    bright/dark stripes indicate inconsistent per-slice correction that outlier
    replacement may not have fully resolved.</p>
    {"".join(blocks)}
    {script}
    </div>
    """

def compute_snr_map(cnr_maps_path, median_filter_size=3):
    """Volume 0 of eddy's CNR-maps output is the b0 SNR map (mean / std across
    the b0 volumes, per eddy_quad's convention). A light median filter reduces
    the influence of random noise and small misalignment, per Tahedl,
    Tournier & Smith (2025)."""
    data = nib.load(str(cnr_maps_path)).get_fdata()
    snr_vol = data[..., 0] if data.ndim == 4 else data
    if median_filter_size and median_filter_size > 1:
        snr_vol = median_filter(snr_vol, size=median_filter_size)
    return snr_vol

def make_triplanar_colorbar_data_uri(vol3d, cmap="viridis", vmin=0, vmax=None,
                                      cbar_label="", max_width=1000):
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]
    fig, axes = plt.subplots(1, 3, figsize=(3 * 4, 4.4), facecolor="black")
    im = None
    for ax_, (label, axis) in zip(axes, views):
        sl = extract_slice(vol3d, axis)
        im = ax_.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax)
        ax_.axis("off")
        ax_.set_title(label, color=MUTED, fontsize=11)
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label(cbar_label, color=MUTED)
    cbar.ax.yaxis.set_tick_params(color=MUTED)
    plt.setp(cbar.ax.get_yticklabels(), color=MUTED)
    cbar.outline.set_edgecolor(GRID)
    fig.subplots_adjust(wspace=0.02, left=0.005, right=0.92, top=0.9, bottom=0.005)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def snr_map_block(cnr_maps_path):
    """Derives and renders the b0 SNR map (with colorbar) directly from eddy's
    CNR-maps NIfTI, in place of eddy_quad's static, colorbar-less PNG."""
    snr_vol = compute_snr_map(cnr_maps_path)
    positive = snr_vol[snr_vol > 0]
    vmax = float(np.percentile(positive, 99)) if positive.size else None

    img_uri = make_triplanar_colorbar_data_uri(
        snr_vol, cmap="viridis", vmin=0, vmax=vmax, cbar_label="SNR (mean / std)"
    )

    return f"""
    <div class="slice-block">
      <h3>SNR Map (b0)</h3>
      <p class="qc-desc">Voxel-wise signal-to-noise ratio (the mean divided by the standard
      deviation across the b0 volumes), taken from the first volume of eddy's CNR-maps
      output and median-filtered to reduce the influence of random noise and small
      misalignment. As a rule of thumb from Tahedl, Tournier &amp; Smith (2025,
      <em>Nature Protocols</em> 20(9):2652&ndash;2684), it's worth checking the SNR in
      regions that are typically low, such as the temporal lobes: if it's still
      reasonably high there, it implies it's high enough everywhere else. In their
      experience, an SNR of roughly 15 in such a region is acceptable for most
      analyses.</p>
      <img src="{img_uri}" class="mosaic" alt="SNR map, axial/coronal/sagittal, with colorbar"/>
    </div>
    """

def eddyqc_section(json_path, rms_path, outlier_path=None, qc_dir=None,
                    raw_dwi_path=None, preproc_dwi_path=None, cnr_maps_path=None,
                    mot_abs_thresh=1.0, mot_rel_thresh=0.5, outlier_pct_thresh=5.0):
    """Section for FSL eddy_quad outputs: qc.json + *.eddy_movement_rms +
    *.eddy_outlier_report + the eddy_quad summary PNGs (avg_b0*.png,
    avg_b<value>.png, cnr####.nii.gz.png), looked up in qc_dir
    (defaults to the directory containing json_path). If raw_dwi_path and
    preproc_dwi_path are given, also builds a raw-vs-processed comparison for
    every volume flagged in the outlier report.
    Thresholds are conventional soft guidelines, not diagnostic standards --
    pass your own to override."""
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
      <p class="qc-desc">Summary statistics and motion/outlier diagnostics from FSL's
      <code>eddy_quad</code>. The highlighted thresholds (abs &gt; {mot_abs_thresh} mm,
      rel &gt; {mot_rel_thresh} mm, outliers &gt; {outlier_pct_thresh}%) are conventional
      soft guidelines, not diagnostic standards &mdash; adjust to your own QC criteria.</p>

      <div class="stat-grid">
        {stat_cards}
        {cnr_cards}
      </div>

      <div class="slice-block">
        <h3>Volume-to-volume Motion (eddy_movement_rms)</h3>
        <img src="data:image/png;base64,{motion_b64}" class="mosaic" alt="motion rms plot"/>
      </div>

      <div class="slice-block">
        <h3>Outlier Slices (eddy_outlier_report)</h3>
        <img src="data:image/png;base64,{outlier_b64}" class="mosaic" alt="outlier scatter plot"/>
      </div>

      {outlier_raw}

      {outlier_volumes_html}

      {images_html}

      {ack_html}
    </section>
    """

# --------------------------------------------------------------------------
# Topup / susceptibility distortion correction QC
# --------------------------------------------------------------------------

def read_acqparams(path):
    """Parse acqparams.tsv: one row per volume fed into topup, columns are
    [PE_x, PE_y, PE_z, total_readout_time, ...]. Returns a list of
    (pe_x, pe_y, pe_z) tuples, one per row/volume."""
    rows = []
    for line in Path(path).read_text().strip().splitlines():
        parts = line.split()
        if len(parts) >= 3:
            rows.append(tuple(float(x) for x in parts[:3]))
    return rows

def pick_pe_volume_indices(acq_rows):
    """First volume index for each distinct phase-encode direction, in the
    order first encountered. Returns list of (pe_vector, index) tuples."""
    seen = {}
    for i, pe in enumerate(acq_rows):
        if pe not in seen:
            seen[pe] = i
    return sorted(seen.items(), key=lambda kv: kv[1])

def load_4d_volume(path):
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

def slice_to_data_uri(slice2d, vmin=None, vmax=None, cmap="gray", max_width=500):
    fig, ax = plt.subplots(figsize=(4, 4), facecolor="black")
    ax.imshow(slice2d, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def slice_overlay_to_data_uri(base_slice, overlay_slice, vmin, vmax, overlay_absmax,
                               overlay_cmap="bwr", alpha=0.55, max_width=500):
    """Grayscale base with a translucent diverging-colormap overlay on top
    (blue = negative field, red = positive field -- same convention as the
    fMRIPrep/QSIPrep SDC boundary plots)."""
    fig, ax = plt.subplots(figsize=(4, 4), facecolor="black")
    ax.imshow(base_slice, cmap="gray", vmin=vmin, vmax=vmax)
    ax.imshow(overlay_slice, cmap=overlay_cmap, vmin=-overlay_absmax, vmax=overlay_absmax, alpha=alpha)
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def make_triplanar_data_uri(vol3d, vmin, vmax, overlay_vol=None, overlay_absmax=None,
                             cmap="gray", overlay_cmap="bwr", alpha=0.55, max_width=1000):
    """Axial, coronal, and sagittal mid-slices side by side in one figure --
    same sizing convention (max_width) as the noise-map mosaics elsewhere in
    the report, so panels read consistently across sections."""
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]
    fig, axes = plt.subplots(1, 3, figsize=(3 * 4, 4.4), facecolor="black")
    for ax_, (label, axis) in zip(axes, views):
        sl = extract_slice(vol3d, axis)
        ax_.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax)
        if overlay_vol is not None:
            osl = extract_slice(overlay_vol, axis)
            ax_.imshow(osl, cmap=overlay_cmap, vmin=-overlay_absmax, vmax=overlay_absmax, alpha=alpha)
        ax_.axis("off")
        ax_.set_title(label, color=MUTED, fontsize=11)
    fig.subplots_adjust(wspace=0.02, left=0.005, right=0.995, top=0.9, bottom=0.005)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def make_triplanar_contour_data_uri(base_vol, mask_vol, vmin, vmax, contour_color="#ff4f4f",
                                     max_width=1000):
    """Axial/coronal/sagittal mid-slices with the mask boundary drawn as a
    contour outline on top -- keeps the underlying anatomy fully visible,
    unlike a filled translucent overlay."""
    views = [("Axial", 2), ("Coronal", 1), ("Sagittal", 0)]
    fig, axes = plt.subplots(1, 3, figsize=(3 * 4, 4.4), facecolor="black")
    for ax_, (label, axis) in zip(axes, views):
        base_sl = extract_slice(base_vol, axis)
        mask_sl = extract_slice(mask_vol, axis)
        ax_.imshow(base_sl, cmap="gray", vmin=vmin, vmax=vmax)
        if np.any(mask_sl > 0.5):
            ax_.contour(mask_sl, levels=[0.5], colors=[contour_color], linewidths=1.3)
        ax_.axis("off")
        ax_.set_title(label, color=MUTED, fontsize=11)
    fig.subplots_adjust(wspace=0.02, left=0.005, right=0.995, top=0.9, bottom=0.005)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    optimized = optimize_png_bytes(buf.read(), max_width=max_width)
    b64 = base64.b64encode(optimized).decode("utf-8")
    return f"data:image/png;base64,{b64}"

def brainmask_section(nodif_path, mask_path):
    """QC for the brain mask: its boundary overlaid as a contour on the nodif
    (b0) reference image, across all three views."""
    nodif = nib.load(str(nodif_path)).get_fdata()
    if nodif.ndim == 4:
        nodif = nodif.mean(axis=-1)
    mask = nib.load(str(mask_path)).get_fdata()

    positive = nodif[nodif > 0]
    vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

    img_uri = make_triplanar_contour_data_uri(nodif, mask, vmin, vmax)

    mask_voxels = int(np.sum(mask > 0.5))
    total_voxels = int(mask.size)
    coverage_pct = 100 * mask_voxels / total_voxels if total_voxels else 0

    return f"""
    <section id="brainmask" class="qc-section">
      <h2>Brain Mask</h2>
      <p class="qc-desc">The brain mask boundary (red outline) overlaid on the nodif (b0)
      reference image. Check that the outline hugs the brain surface without clipping
      cortex (especially at the vertex and cerebellum) and without including obvious
      skull, dura, or eyes.</p>
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
# New: T1–DWI coregistration QC section
# --------------------------------------------------------------------------

def coreg_section(t1w_dwi_path, nodif_path, five_tt_vis_path=None):
    """
    Coregistration QC between a T1w image in DWI space and a regridded nodif (b0)
    in the same space, with optional 5tt2vis overlay.

    - Shows axial, coronal, sagittal mid-slices (triplanar).
    - Slider smoothly cross-fades between T1w and nodif (nodif on top).
    - Optional checkbox toggles 5tt2vis overlay on whichever image is being shown.
    """
    # Load volumes
    t1 = nib.load(str(t1w_dwi_path)).get_fdata()
    nodif = nib.load(str(nodif_path)).get_fdata()

    # Collapse 4D to mean across time if needed
    if t1.ndim == 4:
        t1 = t1.mean(axis=-1)
    if nodif.ndim == 4:
        nodif = nodif.mean(axis=-1)

    if t1.shape != nodif.shape:
        return f"""
        <section id="coreg" class="qc-section">
          <h2>T1–DWI Coregistration</h2>
          <p class="qc-desc">
            Could not render coregistration QC: T1w-in-DWI-space image
            (shape {t1.shape}) and regridded nodif (shape {nodif.shape})
            have different dimensions.
          </p>
        </section>
        """

    # Intensity scaling: separate percentiles for T1 and nodif
    t1_positive = t1[t1 > 0]
    nodif_positive = nodif[nodif > 0]

    t1_vmin, t1_vmax = (
        np.percentile(t1_positive, [2, 98]) if t1_positive.size else (None, None)
    )
    nodif_vmin, nodif_vmax = (
        np.percentile(nodif_positive, [1, 99]) if nodif_positive.size else (None, None)
    )

    # Triplanar images:
    # - T1w: grayscale
    # - nodif: non-grayscale (magma) to highlight hyperintense CSF
    t1_uri = make_triplanar_data_uri(t1, t1_vmin, t1_vmax, cmap="gray")
    nodif_uri = make_triplanar_data_uri(nodif, nodif_vmin, nodif_vmax, cmap="magma")

    has_overlay = five_tt_vis_path is not None

    if has_overlay:
        five_tt = nib.load(str(five_tt_vis_path)).get_fdata()
        if five_tt.ndim == 4:
            five_tt = five_tt.mean(axis=-1)

        if five_tt.shape != t1.shape:
            # Overlay shape mismatch: disable overlay but explain why
            has_overlay = False
            overlay_expl = f"""
            <p class="qc-desc">
              Note: 5tt2vis overlay was not drawn because its shape ({five_tt.shape})
              does not match the T1/nodif grid ({t1.shape}).
            </p>
            """
            t1_overlay_uri = t1_uri
            nodif_overlay_uri = nodif_uri
        else:
            pos = five_tt[five_tt > 0]
            ov_absmax = float(np.percentile(pos, 99)) if pos.size else 1.0

            t1_overlay_uri = make_triplanar_data_uri(
                t1,
                t1_vmin,
                t1_vmax,
                overlay_vol=five_tt,
                overlay_absmax=ov_absmax,
                cmap="gray",
                overlay_cmap="tab10",
                alpha=0.6,
            )
            nodif_overlay_uri = make_triplanar_data_uri(
                nodif,
                nodif_vmin,
                nodif_vmax,
                overlay_vol=five_tt,
                overlay_absmax=ov_absmax,
                cmap="magma",
                overlay_cmap="tab10",
                alpha=0.6,
            )

            overlay_expl = """
            <p class="qc-desc">
              The tissue-type overlay comes from <code>5tt2vis</code>, drawn on top of
              the current base image (T1 or nodif) using a qualitative colormap.
              This makes it easier to judge whether GM/WM boundaries from the
              segmentation follow the anatomy in both contrasts.
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
        "nodif_label": "Nodif (b0, regridded to DWI space)",
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
        <input type="range" min="0" max="1" step="0.02" value="0" class="topup-range"
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

        // Cross-fade: T1 bottom, nodif on top
        nodifImg.style.opacity = alpha.toString();
        t1Img.style.opacity = (1.0 - alpha).toString();

        // Label reflects which image dominates visually
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
        Drag the slider to smoothly transition between the T1w image in DWI space and the
        regridded nodif (b0) in the same space, across axial, coronal and sagittal
        mid-slices. Borders between CSF and GM/WM, especially around the ventricles
        and cortical ribbon, should align well in both contrasts if coregistration
        is satisfactory. The nodif is shown with a non-grayscale colormap to make the
        hyperintense CSF stand out.
      </p>
      {overlay_expl}
      {overlay_control}
      {block}
      {script}
    </section>
    """

# --------------------------------------------------------------------------
# Topup helpers continued
# --------------------------------------------------------------------------

def find_matching_pe_index(topup_pe_items, dwi_acqparams_path):
    """Match the DWI's own phase-encode direction (first row of its
    acqparams.tsv) to one of the topup input volumes. Returns
    (vec, index, matched: bool). Falls back to the first topup PE direction
    with matched=False if no exact match is found or no dwi acqparams given."""
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
                   fieldmap_path=None, dwi_acqparams_path=None, matched_pe_index=None):
    """QC for susceptibility distortion correction (topup), regardless of
    whether the fieldmap came from a real reversed-PE acquisition or a
    synthetic SynB0-DisCo one -- both produce the same before/after file
    pair, so no branching is needed here.

    before_path: pre-correction 4D volume fed into topup (e.g. *_desc-4topup_epi.nii.gz)
    after_path:  corrected volume, same volume order (e.g. *_desc-unwarped_epi.nii.gz)
    topup_acqparams_path: PE-direction table for the topup input volumes (e.g. refparams.tsv)
    fieldmap_path: topup's --fout field map in Hz (e.g. *_desc-topup_fieldmap.nii.gz).
        NOTE: this is defined in the same space as `after_path`, not `before_path` --
        the overlay is only drawn on the after-correction frame for that reason.
        Do NOT pass the *_fieldcoef.nii.gz file here: it's topup's internal spline
        representation, not a directly plottable field.
    dwi_acqparams_path: the DWI's own acqparams.tsv, used to figure out which of the
        topup input volumes shares the DWI's phase-encode direction.
    matched_pe_index: manual override -- skip the matching logic and use this volume
        index directly.
    """
    topup_acq_rows = read_acqparams(topup_acqparams_path)
    pe_items = pick_pe_volume_indices(topup_acq_rows)

    if not pe_items:
        return f"""
        <section id="topup" class="qc-section">
          <h2>Susceptibility Distortion Correction (topup)</h2>
          <p class="qc-desc">No phase-encode rows could be parsed from
          <code>{Path(topup_acqparams_path).name}</code>.</p>
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
        <p class="qc-desc" style="color:#e0a54e;">Note: the DWI's phase-encode direction
        could not be confirmed against <code>{Path(topup_acqparams_path).name}</code>
        (no <code>--topup-dwi-acqparams</code> given, or no exact match found), so the
        first phase-encode direction PE {pe_label(pe_vec)} is shown below by default.
        Pass <code>--topup-match-pe-index</code> to select a specific volume explicitly.</p>
        """

    before_data = load_4d_volume(before_path)
    after_data = load_4d_volume(after_path)
    before_vol = get_vol(before_data, pe_idx)
    after_vol = get_vol(after_data, pe_idx)

    both_vals = np.concatenate([before_vol.flatten(), after_vol.flatten()])
    positive = both_vals[both_vals > 0]
    vmin, vmax = np.percentile(positive, [1, 99]) if positive.size else (None, None)

    fieldmap_data = None
    overlay_absmax = None
    if fieldmap_path and Path(fieldmap_path).exists():
        fieldmap_data = nib.load(str(fieldmap_path)).get_fdata()
        overlay_absmax = np.percentile(np.abs(fieldmap_data), 99)

    before_label = f"Before correction &mdash; PE {pe_label(pe_vec)} (matches DWI)"
    after_label = f"After correction &mdash; PE {pe_label(pe_vec)} (matches DWI)"

    before_uri = make_triplanar_data_uri(before_vol, vmin, vmax)
    after_plain_uri = make_triplanar_data_uri(after_vol, vmin, vmax)

    has_overlay = fieldmap_data is not None
    if has_overlay:
        after_overlay_uri = make_triplanar_data_uri(
            after_vol, vmin, vmax, overlay_vol=fieldmap_data, overlay_absmax=overlay_absmax
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
        Show distortion overlay (blue/red = field, applies to corrected image only)
      </label>
    """ if has_overlay else ""

    block = f"""
    <div class="slice-block">
      <img id="topup-img" class="mosaic slider-img" src="{before_uri}"
           alt="axial, coronal, sagittal topup before/after comparison"/>
      <div class="topup-slider-row">
        <span class="topup-slider-endlabel">Before</span>
        <input type="range" min="0" max="1" step="1" value="0" class="topup-range"
               id="topup-range" oninput="updateTopupView()">
        <span class="topup-slider-endlabel">After</span>
      </div>
      <div class="slider-label" id="topup-label">{before_label}</div>
    </div>
    """

    script = f"""
    <script>
      window.topupData = {json.dumps(view_data)};
      function updateTopupView() {{
        var d = window.topupData;
        var mode = document.getElementById('topup-range').value === '1' ? 'after' : 'before';
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
        <p class="qc-desc">The distortion overlay shows topup's estimated off-resonance
        field (in Hz) on top of the corrected image, using a blue/red diverging colormap
        the same way fMRIPrep and QSIPrep display SDC results: blue and red indicate
        opposite-signed field values, which translate to opposite-direction voxel
        displacement along the phase-encode axis. Areas with the strongest color
        (regardless of sign) are where topup found &mdash; and corrected for &mdash;
        the most distortion, typically near air/tissue boundaries such as the sinuses,
        ear canals, and orbitofrontal cortex. The overlay is only meaningful on the
        after-correction image, since the field map shares that image's (undistorted)
        geometry, not the before-correction image's.</p>
        """

    return f"""
    <section id="topup" class="qc-section">
      <h2>Susceptibility Distortion Correction (topup)</h2>
      <p class="qc-desc">Drag the slider to compare the DWI's own phase-encode direction
      before and after topup's distortion correction, across all three views at once.
      Anatomical boundaries (ventricles, brainstem, orbitofrontal and temporal cortex)
      that appear shifted or blurred before correction should look sharper and better
      aligned with expected anatomy after correction.</p>
      {match_note}
      {overlay_expl}
      {global_overlay_control}
      {block}
      {script}
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
    width: 100%;
    display: block;
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
  /* Coregistration stack: T1 base, nodif overlay */
  .coreg-stack {{
    position: relative;
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
    """extra_nav: list of (after_section_id, anchor_id, label) tuples for nav
    buttons that point to a subsection anchor rather than a full <section>,
    inserted right after their parent section's nav link."""
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
    page = PAGE_TEMPLATE.format(subject_suffix=subject_suffix, nav_links=nav_links, sections=body)
    Path(output_path).write_text(page)

def main():
    p = argparse.ArgumentParser(description="Generate DWI preprocessing QC HTML report")
    p.add_argument("--noise", required=True, help="Path to dwidenoise noise map (e.g. noise.nii.gz)")
    p.add_argument("--eddy-json", default=None, help="Path to eddy_quad qc.json")
    p.add_argument("--eddy-rms", default=None, help="Path to *.eddy_movement_rms")
    p.add_argument("--eddy-outliers", default=None, help="Path to *.eddy_outlier_report")
    p.add_argument("--eddy-qc-dir", default=None,
                    help="Directory containing eddy_quad summary PNGs (avg_b0.png, cnr*.png, "
                         "etc). Defaults to the directory containing --eddy-json.")
    p.add_argument("--eddy-raw-dwi", default=None,
                    help="Raw (pre-eddy) 4D DWI volume, e.g. *_desc-dns+degibbs_dwi.nii.gz "
                         "-- used to build raw-vs-processed comparisons for outlier volumes.")
    p.add_argument("--eddy-preproc-dwi", default=None,
                    help="Eddy-processed 4D DWI volume, e.g. *_desc-preproc_dwi.nii.gz, same "
                         "volume order as --eddy-raw-dwi.")
    p.add_argument("--topup-before", default=None,
                    help="Pre-correction 4D volume fed into topup (e.g. *_desc-4topup_epi.nii.gz)")
    p.add_argument("--topup-after", default=None,
                    help="Post-correction volume, same volume order (e.g. *_desc-unwarped_epi.nii.gz)")
    p.add_argument("--topup-acqparams", default=None,
                    help="PE-direction table for the topup input volumes (e.g. refparams.tsv)")
    p.add_argument("--topup-fieldmap", default=None,
                    help="topup's --fout field map in Hz (e.g. *_desc-topup_fieldmap.nii.gz). "
                         "Do NOT pass *_fieldcoef.nii.gz here -- that's the internal spline "
                         "representation, not a plottable field.")
    p.add_argument("--topup-dwi-acqparams", default=None,
                    help="The DWI's own acqparams.tsv, used to match which topup input volume "
                         "shares the DWI's phase-encode direction.")
    p.add_argument("--topup-match-pe-index", type=int, default=None,
                    help="Manual override: skip PE matching and use this topup input volume index.")
    p.add_argument("--brainmask-nodif", default=None,
                    help="Nodif (b0) reference image, e.g. *_desc-nodif_dwi.nii.gz or "
                         "*_desc-nodif-brain_dwi.nii.gz")
    p.add_argument("--brainmask-mask", default=None,
                    help="Brain mask, e.g. *_desc-brain_mask.nii.gz")
    p.add_argument("--eddy-cnr-maps", default=None,
                    help="eddy's CNR-maps 4D volume, e.g. *_label-cnr-maps_desc-preproc_dwi.nii.gz. "
                         "Volume 0 is the b0 SNR map (mean/std across b0 volumes) -- used to render "
                         "an SNR map with a colorbar in place of eddy_quad's static b0 SNR PNG.")

    # New: T1–DWI coregistration inputs
    p.add_argument(
        "--reg-t1w-dwi",
        default=None,
        help=(
            "T1w image already in DWI space "
            "(e.g. *desc-preproc_T1w_space-dwi.nii.gz), used for T1–DWI coregistration QC."
        ),
    )
    p.add_argument(
        "--reg-nodif",
        default=None,
        help=(
            "Nodif (b0) reference image resampled to the same grid as --reg-t1w-dwi "
            "(e.g. *desc-nodif_space-dwi.nii.gz)."
        ),
    )
    p.add_argument(
        "--reg-5ttvis",
        default=None,
        help=(
            "Optional 5tt2vis image in the same space as --reg-t1w-dwi, used as a "
            "tissue-type overlay in the coregistration QC section."
        ),
    )

    p.add_argument("--output", default="qc_report.html", help="Output HTML path")
    p.add_argument("--subject", default=None, help="Subject label, e.g. sub-01")
    args = p.parse_args()

    sections = [
        ("noise", "Noise Map", noise_section(args.noise)),
    ]

    extra_nav = []

    if args.eddy_json and args.eddy_rms:
        eddyqc_html = eddyqc_section(
            args.eddy_json, args.eddy_rms, args.eddy_outliers, qc_dir=args.eddy_qc_dir,
            raw_dwi_path=args.eddy_raw_dwi, preproc_dwi_path=args.eddy_preproc_dwi,
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
            ("brainmask", "Brain Mask", brainmask_section(args.brainmask_nodif, args.brainmask_mask))
        )

    # New: T1–DWI coregistration section
    if args.reg_t1w_dwi and args.reg_nodif:
        sections.append(
            (
                "coreg",
                "T1–DWI Coregistration",
                coreg_section(args.reg_t1w_dwi, args.reg_nodif, five_tt_vis_path=args.reg_5ttvis),
            )
        )

    build_report(sections, args.output, subject=args.subject, extra_nav=extra_nav)
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()