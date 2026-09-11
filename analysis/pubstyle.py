"""Shared publication-figure styling (Nature-portfolio conventions).

Import this module from any figure-producing code and call `pubstyle()`
once before drawing, so all figures in the project stay visually consistent.
"""

import matplotlib.pyplot as plt

# --- Figure sizing (inches), for print columns -----------------------------
COL1 = 3.50   # single column (89 mm)
COL2 = 7.20   # double column (183 mm)

# --- Colorblind-safe categorical palette (Wong, 2011, Nature Methods) ------
WONG = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "skyblue": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
    "grey": "#999999",
}
WONG_ORDER = ["blue", "orange", "green", "vermillion", "purple", "skyblue", "yellow", "grey"]
WONG_CYCLE = [WONG[k] for k in WONG_ORDER]

# Sequential colormap for heatmaps / confusion matrices with non-negative data.
SEQUENTIAL_CMAP = "Blues"
# Diverging colormap for correlation-style heatmaps (can be negative);
# RdBu_r is used instead of the default coolwarm for better colorblind safety
# while remaining a standard diverging Nature-style choice.
DIVERGING_CMAP = "RdBu_r"

PUB_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.pad": 1,
    "ytick.major.pad": 1,
    "axes.labelpad": 1,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.dpi": 150,
}


def pubstyle():
    """Apply the shared publication rcParams. Call once per figure/script."""
    plt.rcParams.update(PUB_RC)


def _despine(ax, keep=("left", "bottom")):
    """Remove all spines except those named in `keep`."""
    for side, spine in ax.spines.items():
        spine.set_visible(side in keep)
    ax.tick_params(top=False, right=False)


def suptitle(fig, text, fontsize=7.5):
    """Left-aligned figure suptitle, consistent across all figures."""
    fig.suptitle(text, x=0.01, ha="left", fontsize=fontsize, fontweight="bold")


def colorbar_kwargs():
    """Standard small, unobtrusive colorbar styling for heatmaps."""
    return {"fraction": 0.02, "pad": 0.01}


def style_colorbar(cbar, labelsize=6.5):
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=labelsize)


def annotation_color(value, vmin, vmax, dark="#333333", light="#FFFFFF"):
    """Pick white/dark-grey text color for a heatmap cell based on normalized value."""
    if vmax == vmin:
        norm = 0.5
    else:
        norm = (value - vmin) / (vmax - vmin)
    return light if norm > 0.6 or norm < 0.25 else dark
