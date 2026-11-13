"""Visium B148: E67 developing-knee spatial reference. 

Single-gene displayed expression: natural log1p(CP10k).
Requires original raw H5 and processed coordinate H5AD via environment variables.
"""

# =============================================================================
# CELL 1 | PARAMETERS & JUPYTER LAB PLOTTING ENVIRONMENT
# =============================================================================
import gc
import os
from pathlib import Path
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

warnings.filterwarnings("ignore")

RAW_H5 = Path(os.environ["VISIUM_RAW_H5"])
SPATIAL_H5AD = Path(os.environ["VISIUM_H5AD"])
OUT_DIR = Path(os.environ.get("OA_OUTPUT_DIR", "outputs/visium"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

ROUND_1 = [
    "PRG4", "CRTAC1", "CLU", "ACAN", "CILP", "COMP", "COL10A1", "EPYC", "IBSP",
]
ROUND_2 = [
    "THBS4", "TNC", "POSTN", "UCMA", "MATN3", "COL9A1", "MATN1", "CYTL1", "CHRDL2",
]
GENES = ROUND_1 + ROUND_2
SELECTED_GENES = {"PRG4", "CYTL1", "COMP", "MATN1"}

ANGLE_DEG = -45
SPOT_SIZE = 0.85
MAX_CLIP_QUANTILE = 99.0
CMAP_NAME = "magma"
DARK_BG = "#101820"
ZERO_BG = "#26303A"
HIGHLIGHT_COLOR = "#20A39E"
DPI = 350

XMIN, YMIN = 2463.3256, 3502.2579
XMAX, YMAX = 10307.3373, 11346.2696

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.0,
    "axes.titlesize": 9.5,
    "axes.titleweight": "bold",
    "text.color": "#252525",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

print("✔ Cell 1 Status")

# =============================================================================
# CELL 2 | UMI MATRIX EXTRACTION & -45° ROTATION ALIGNMENT
# =============================================================================
print("1. Status...", flush=True)
adata = sc.read_h5ad(SPATIAL_H5AD, backed="r")
try:
    bc = pd.Index(adata.obs_names)
    xy = np.asarray(adata.obsm["spatial"][:, :2], dtype=np.float64)
finally:
    adata.file.close()

assert len(bc) == len(xy) == 69618, "Status"

print("2. Status 18 Status UMI Status...", flush=True)
raw = sc.read_10x_h5(RAW_H5)
raw.var_names_make_unique()

missing = [g for g in GENES if g not in raw.var_names]
if missing:
    raise ValueError("Status: " + ", ".join(missing))

idx = raw.obs_names.get_indexer(bc)
lib_size = np.asarray(raw.X[idx].sum(axis=1)).ravel().astype(np.float64)
cols = raw.var_names.get_indexer(GENES)
counts_18 = raw.X[idx, :][:, cols].toarray().astype(np.float64)

expr_all = np.log1p(counts_18 * (1e4 / np.maximum(lib_size[:, None], 1.0)))
del raw, counts_18, lib_size
gc.collect()

center = xy.mean(axis=0)
th = np.deg2rad(ANGLE_DEG)
rot_mat = np.array([
    [np.cos(th), -np.sin(th)],
    [np.sin(th), np.cos(th)],
])
xy_rot = (xy - center) @ rot_mat.T + center
x_all, y_all = xy_rot.T

inside = (x_all >= XMIN) & (x_all <= XMAX) & (y_all >= YMIN) & (y_all <= YMAX)
assert int(inside.sum()) == 61008, (
    f"ROI spot count mismatch: expected 61008, observed {inside.sum()}."
)

x = x_all[inside]
y = y_all[inside]
E = expr_all[inside]
del expr_all
gc.collect()

print(f"✔ Cell 2 data loaded; ROI contains {len(x):,} high-quality spatial spots.")

# =============================================================================
# CELL 3 | PLOTTING & MULTI-FORMAT EXPORT (ROUND 1 & ROUND 2)
# =============================================================================
def plot_gene_panel_3x3(gene_list, round_name):
    fig, axes = plt.subplots(
        3, 3, figsize=(7.08, 7.35), dpi=DPI, facecolor="white",
    )
    axes = axes.flatten()

    for i, gene in enumerate(gene_list):
        ax = axes[i]
        col_idx = GENES.index(gene)
        vals = E[:, col_idx]
        ax.set_facecolor(DARK_BG)

        zero_mask = vals <= 0
        ax.scatter(
            x[zero_mask], y[zero_mask],
            s=SPOT_SIZE * 0.7, c=ZERO_BG,
            marker="o", linewidths=0, rasterized=True,
        )

        pos_mask = vals > 0
        if np.any(pos_mask):
            v_max = np.percentile(vals[pos_mask], MAX_CLIP_QUANTILE)
            norm = Normalize(vmin=0, vmax=max(v_max, 1e-4), clip=True)
            sc_pos = ax.scatter(
                x[pos_mask], y[pos_mask],
                c=vals[pos_mask], cmap=CMAP_NAME, norm=norm,
                s=SPOT_SIZE, marker="o", linewidths=0, rasterized=True,
            )

            cax = ax.inset_axes([0.65, 0.08, 0.28, 0.04])
            cb = fig.colorbar(sc_pos, cax=cax, orientation="horizontal")
            cb.ax.tick_params(labelsize=5.2, length=1.5, pad=1, colors="#FFFFFF")
            cb.set_ticks([0, round(v_max, 1)])
            cb.outline.set_linewidth(0.4)
            cb.outline.set_color("#FFFFFF")

        ax.set_xlim(XMIN, XMAX)
        ax.set_ylim(YMAX, YMIN)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

        is_selected = gene in SELECTED_GENES
        title_box = dict(
            boxstyle="round,pad=0.25",
            facecolor=HIGHLIGHT_COLOR if is_selected else "#FFFFFF",
            edgecolor="#252525",
            linewidth=0.8 if is_selected else 0.5,
            alpha=0.95,
        )
        italic_gene = f"$\\it{{{gene}}}$"
        title_text = f"{italic_gene} (selected)" if is_selected else italic_gene
        ax.set_title(
            title_text, fontsize=8.5, fontweight="bold",
            color="#FFFFFF" if is_selected else "#252525",
            pad=4, bbox=title_box,
        )

    display_title = (
        "10x Visium marker screening: round 1 (literature candidates)"
        if "Round1" in round_name
        else "10x Visium marker screening: round 2 (matrix refinement)"
    )
    fig.suptitle(
        display_title, fontsize=12, fontweight="bold",
        y=0.985, color="#252525",
    )

    plt.subplots_adjust(
        left=0.04, right=0.96, top=0.93, bottom=0.03,
        wspace=0.06, hspace=0.10,
    )

    out_base = OUT_DIR / f"Visium_Marker_Screening_{round_name}"
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{out_base}.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(
        f"{out_base}.tiff", dpi=DPI, bbox_inches="tight",
        facecolor="white", pil_kwargs={"compression": "tiff_lzw"},
    )

    plt.show()
    plt.close(fig)
    print(f"✔ {round_name} saved to: {out_base.name} (PDF/SVG/TIFF)")


plot_gene_panel_3x3(ROUND_1, "Round1_LiteratureCandidates")
plot_gene_panel_3x3(ROUND_2, "Round2_MatrixRefinement")
