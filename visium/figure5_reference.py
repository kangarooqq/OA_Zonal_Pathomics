"""Visium spatial transcriptomics visualization: E67 developing human knee. 

Requires VISIUM_H5AD containing aligned spatial coordinates and expression. The
historical marker normalization and geometric ROI criteria are retained; this script
requires a source-data run to verify the resulting visualization.
"""
# =============================================================================
# FIGURE 5a & 5b | 10x VISIUM HUMAN KNEE SPATIAL TRANSCRIPTOMICS REFERENCE
#
# Compact visualization with fixed figure dimensions
#
# Output: PDF / SVG / TIFF (350 dpi, LZW)
# Output directory: OA_OUTPUT_DIR
# =============================================================================

from pathlib import Path
import os
import scanpy as sc

adata = sc.read_h5ad(os.environ["VISIUM_H5AD"])
import warnings
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")

# =============================================================================
# 1. OUTPUT DIRECTORY & COMPACT FIGURE SPECIFICATIONS
# =============================================================================

OUT_DIR = Path(os.environ.get("OA_OUTPUT_DIR", "outputs/visium"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

FIG_NAME = "Visium_Spatial_Reference_Compact"
DPI = 350

# Historical annotation removed during English-language code curation.
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "text.color": "#FFFFFF",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.default": "regular",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "black",
    "axes.facecolor": "black",
    "savefig.facecolor": "black",
})

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

def get_gene_vec(adata_obj, gene_name):
    if gene_name not in adata_obj.var_names:
        raise KeyError(f"Gene '{gene_name}' not found in adata.var_names")
    expr = adata_obj[:, gene_name].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    return np.asarray(expr).ravel()

def get_norm_channel(adata_obj, gene_input, mask):
    if isinstance(gene_input, (list, tuple)):
        mat = np.column_stack([get_gene_vec(adata_obj, g) for g in gene_input])
        vec = np.max(mat, axis=1)
    else:
        vec = get_gene_vec(adata_obj, gene_input)

    vec_masked = np.where(mask, vec, 0.0)
    pos_vals = vec_masked[vec_masked > 0]
    if len(pos_vals) == 0:
        return np.zeros_like(vec_masked)

    p97 = np.percentile(pos_vals, 97)
    return np.clip(vec_masked / (p97 + 1e-6), 0.0, 1.0)

# =============================================================================
# 3. GEOMETRY & -45° MECHANICAL ALIGNMENT
# =============================================================================

def rotate_coords(coords, angle_deg=-45):
    theta = np.radians(angle_deg)
    center = np.mean(coords, axis=0)
    shifted = coords - center
    rot_mat = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])
    return shifted @ rot_mat.T + center

coords_rot45 = rotate_coords(adata.obsm["spatial"][:, :2], angle_deg=-45)

pad_full = 0.03
fx_all_min, fx_all_max = coords_rot45[:, 0].min(), coords_rot45[:, 0].max()
fy_all_min, fy_all_max = coords_rot45[:, 1].min(), coords_rot45[:, 1].max()
w_all = fx_all_max - fx_all_min
h_all = fy_all_max - fy_all_min

full_xlim = (fx_all_min - w_all * pad_full, fx_all_max + w_all * pad_full)
full_ylim = (fy_all_min - h_all * pad_full, fy_all_max + h_all * pad_full)
full_aspect_ratio = (full_ylim[1] - full_ylim[0]) / (full_xlim[1] - full_xlim[0])

# Historical annotation removed during English-language code curation.
matn1_expr = get_gene_vec(adata, "MATN1")
core_pts = coords_rot45[matn1_expr > 1.5]
y_split = np.percentile(core_pts[:, 1], 45)
femur_core = core_pts[core_pts[:, 1] <= y_split]

fx_min, fx_max = np.percentile(femur_core[:, 0], [2, 98])
fy_min, fy_max = np.percentile(femur_core[:, 1], [2, 98])
w_span, h_span = fx_max - fx_min, fy_max - fy_min

femur_roi_base = (
    fx_min - w_span * 0.25,
    fx_max + w_span * 0.25,
    fy_min - h_span * 0.35,
    fy_max + h_span * 0.55,
)

x0, x1, y0, y1 = femur_roi_base
W, H = x1 - x0, y1 - y0
rx0 = x0 + W * 0.08
rx1 = x0 + W * 0.98
ry0 = y0 + H * 0.18
ry1 = y0 + H * 0.70

roi_h = ry1 - ry0
roi_w_matched = roi_h / full_aspect_ratio
rx_center = (rx0 + rx1) / 2.0
shift_left = roi_w_matched * 0.15

femur_roi_tight = (
    rx_center - roi_w_matched / 2.0 - shift_left,
    rx_center + roi_w_matched / 2.0 - shift_left,
    ry0,
    ry1,
)

if "cartilage_mask" not in globals():
    col2a1 = get_gene_vec(adata, "COL2A1")
    prg4 = get_gene_vec(adata, "PRG4")
    cartilage_mask = (col2a1 > 0.5) | (prg4 > 1.0)

# =============================================================================
# 4. COMPACT PLOTTING ENGINE
# =============================================================================

def plot_compact_visium_panel(
    coords,
    c1_in,
    c2_in,
    c3_in,
    mask,
    ax,
    panel_label="",
    sub_title="",
    c1_label="",
    c2_label="",
    c3_label="",
    crop_box=None,
    draw_roi=None,
):
    v1 = get_norm_channel(adata, c1_in, mask)
    v2 = get_norm_channel(adata, c2_in, mask)
    v3 = get_norm_channel(adata, c3_in, mask)

    hex1, hex2, hex3 = "#FF2A85", "#FFD700", "#00F0FF"

    def hex_to_rgb(h):
        h = h.lstrip("#")
        return np.array([int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)])

    rgb1, rgb2, rgb3 = hex_to_rgb(hex1), hex_to_rgb(hex2), hex_to_rgb(hex3)
    blended_rgb = np.clip(
        np.outer(v1, rgb1) + np.outer(v2, rgb2) + np.outer(v3, rgb3), 0.0, 1.0
    )
    active = (v1 + v2 + v3) > 0.05

    ax.set_facecolor("black")

    # Historical annotation removed during English-language code curation.
    ax.scatter(
        coords[~active, 0], coords[~active, 1],
        c="#181818",
        s=1.0 if crop_box is None else 7.0,
        marker="o", linewidths=0
    )
    ax.scatter(
        coords[active, 0], coords[active, 1],
        c=blended_rgb[active],
        s=3.0 if crop_box is None else 28.0,
        marker="o", linewidths=0
    )

    if crop_box is None:
        ax.set_xlim(full_xlim[0], full_xlim[1])
        ax.set_ylim(full_ylim[1], full_ylim[0])
    else:
        ax.set_xlim(crop_box[0], crop_box[1])
        ax.set_ylim(crop_box[3], crop_box[2])

    ax.set_aspect("equal", adjustable="box")
    ax.set_box_aspect(full_aspect_ratio)

    if draw_roi is not None:
        rect = patches.Rectangle(
            (draw_roi[0], draw_roi[2]),
            draw_roi[1] - draw_roi[0],
            draw_roi[3] - draw_roi[2],
            linewidth=1.3,
            edgecolor="#FFFFFF",
            facecolor="none",
            linestyle="--",
            zorder=8
        )
        ax.add_patch(rect)
        
        ax.text(
            draw_roi[0] + (draw_roi[1] - draw_roi[0]) * 0.05,
            draw_roi[2] - (draw_roi[3] - draw_roi[2]) * 0.035,
            "Femoral condyle ROI",
            color="#FFFFFF",
            fontsize=6.8,
            fontweight="bold",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="#000000")],
            zorder=9
        )

    # Historical annotation removed during English-language code curation.
    if panel_label:
        ax.text(
            -0.03, 1.05,
            panel_label,
            transform=ax.transAxes,
            color="#FFFFFF",
            fontsize=15.0,
            fontweight="bold",
            va="bottom", ha="left"
        )

    # Historical annotation removed during English-language code curation.
    ax.set_title(
        sub_title,
        fontsize=9,
        fontweight="bold",
        color="#FFFFFF",
        pad=7
    )

    # Historical annotation removed during English-language code curation.
    legend_elements = [
        patches.Patch(facecolor=hex1, edgecolor="#FFFFFF", label=c1_label),
        patches.Patch(facecolor=hex2, edgecolor="#FFFFFF", label=c2_label),
        patches.Patch(facecolor=hex3, edgecolor="#FFFFFF", label=c3_label),
    ]
    leg = ax.legend(
        handles=legend_elements,
        loc="lower right",
        frameon=True,
        facecolor="#0A0A0A",
        edgecolor="#333333",
        framealpha=0.90,
        fontsize=6.5,
        borderpad=0.4,
        labelspacing=0.35,
        handlelength=1.0,
        handleheight=0.7
    )
    for text in leg.get_texts():
        text.set_color("#FFFFFF")
        text.set_fontweight("bold")

    ax.set_xticks([])
    ax.set_yticks([])

# =============================================================================
# 5. RENDER COMPACT 2x2 FIGURE (7.08 in × 9.20 in)
# =============================================================================
# Historical annotation removed during English-language code curation.
fig, axes = plt.subplots(2, 2, figsize=(7.08, 9.20), facecolor="black")

# ROW 1 | FIGURE 5a
plot_compact_visium_panel(
    coords_rot45,
    ["PRG4", "CYTL1"], "COMP", "MATN1",
    cartilage_mask,
    axes[0, 0],
    panel_label="a",
    sub_title="3-layer sandwich (Overview)",
    c1_label="AS/SZ: PRG4 + CYTL1",
    c2_label="MZ: COMP",
    c3_label="DZ: MATN1",
    crop_box=None,
    draw_roi=femur_roi_tight,
)

plot_compact_visium_panel(
    coords_rot45,
    ["PRG4", "CYTL1"], "COMP", "MATN1",
    cartilage_mask,
    axes[0, 1],
    panel_label="",
    sub_title="Femoral condyle zoom",
    c1_label="AS/SZ (0–15%): PRG4 + CYTL1",
    c2_label="MZ (15–50%): COMP",
    c3_label="DZ (50–100%): MATN1",
    crop_box=femur_roi_tight,
    draw_roi=None,
)

# ROW 2 | FIGURE 5b
plot_compact_visium_panel(
    coords_rot45,
    "PRG4", "CYTL1", "MATN1",
    cartilage_mask,
    axes[1, 0],
    panel_label="b",
    sub_title="Sub-zonal dissection (Overview)",
    c1_label="AS front: PRG4",
    c2_label="SZ sub-surface: CYTL1",
    c3_label="MZ/DZ core: MATN1",
    crop_box=None,
    draw_roi=femur_roi_tight,
)

plot_compact_visium_panel(
    coords_rot45,
    "PRG4", "CYTL1", "MATN1",
    cartilage_mask,
    axes[1, 1],
    panel_label="",
    sub_title="Femoral condyle zoom",
    c1_label="AS (Surface): PRG4",
    c2_label="SZ (Sub-surface): CYTL1",
    c3_label="MZ/DZ (Deep core): MATN1",
    crop_box=femur_roi_tight,
    draw_roi=None,
)

# Historical annotation removed during English-language code curation.
plt.subplots_adjust(
    left=0.035, right=0.965,
    top=0.950, bottom=0.020,
    wspace=0.06, hspace=0.10
)

# =============================================================================
# 6. EXPORT PDF / SVG / TIFF (350 DPI)
# =============================================================================

for ext in ("pdf", "svg", "tiff"):
    out_file = OUT_DIR / f"{FIG_NAME}.{ext}"
    if ext == "tiff":
        fig.savefig(
            out_file, dpi=DPI, bbox_inches="tight", facecolor="black",
            pil_kwargs={"compression": "tiff_lzw"}
        )
    else:
        fig.savefig(out_file, dpi=DPI, bbox_inches="tight", facecolor="black")
    print(f"Saved: {out_file.name}")

plt.show()
plt.close(fig)

print("\n" + "=" * 70)
print(f"COMPACT FIGURE SAVED TO: {OUT_DIR}")
print("VISUALIZATION EXPORT COMPLETE.")
print("=" * 70)
