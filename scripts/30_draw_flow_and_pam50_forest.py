#!/usr/bin/env python3
"""30_draw_flow_and_pam50_forest.py
Generate publication-quality study flowchart (Fig1_flow) and PAM50 subgroup forest plot (Fig_PAM50_forest).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, ArrowStyle, FancyArrowPatch
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------
# 1. REMARK-style Study Flowchart (Fig1_flow)
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

def draw_box(x, y, w, h, text, title="", bg="#F4F6F9", ec="#2C3E50", title_color="#1A252F", lw=1.2, fontsize=8.5):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.2",
                       fc=bg, ec=ec, lw=lw, zorder=2)
    ax.add_patch(p)
    if title:
        ax.text(x + w / 2, y + h - 2.8, title, ha="center", va="center",
                fontsize=fontsize + 0.8, fontweight="bold", color=title_color, zorder=3)
        ax.text(x + w / 2, y + (h - 2.8) / 2, text, ha="center", va="center",
                fontsize=fontsize, color="#2C3E50", zorder=3, linespacing=1.25)
    else:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fontsize, color="#2C3E50", zorder=3, linespacing=1.25)

def draw_arrow(x1, y1, x2, y2, ec="#34495E"):
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                            arrowstyle="-|>", mutation_scale=12,
                            color=ec, lw=1.5, zorder=1)
    ax.add_patch(arrow)

# Discovery Box
draw_box(26, 88, 48, 10,
         "Microarray (Affymetrix GPL570)\n104 breast tumors + 17 normal tissues (35 OS events)",
         title="Discovery Cohort: GSE42568", bg="#EBF5FB", ec="#2980B9")

# Step 1: Candidate Pool
draw_arrow(50, 88, 50, 77)
draw_box(26, 68, 48, 9,
         "Tumor vs. Normal DE (|log2FC| > 0.5, BH-FDR < 0.05)\nintersected with 11 KEGG immune & metabolic pathways\n-> 297 candidate genes",
         title="Feature Selection: Candidate Pool", bg="#FEF9E7", ec="#F39C12")

# Step 2: Univariable & LASSO
draw_arrow(50, 68, 50, 57)
draw_box(26, 48, 48, 9,
         "Full-cohort Cox screen (p < 0.01 -> 41 genes)\n5-fold CV tuned LASSO alpha in the screened set\n-> 14 genes; multivariable Cox refit",
         title="Model Construction: Penalized Cox", bg="#FEF9E7", ec="#F39C12")

# Step 3: Locked Signature
draw_arrow(50, 48, 50, 37)
draw_box(22, 28, 56, 9,
         "14 genes: 6 metabolic + 8 immune\nParameters recorded as locked; GSE20685 same-commit order unresolved\nEPV 2.5; planned reduction rule not followed",
         title="Locked 14-Gene Immune-Metabolic Score", bg="#EAFAF1", ec="#27AE60", fontsize=7.4)

# Branching arrows
draw_arrow(36, 28, 17, 19)
draw_arrow(50, 28, 50, 19)
draw_arrow(64, 28, 83, 19)

# Supportive same-platform cohort evaluation: GSE20685
draw_box(2, 3, 30, 16,
         "327 breast tumors (83 OS deaths)\nSupportive same-platform evaluation (GPL570)\nLock/analysis order is not established\nContinuous HR = 1.59 [1.29–1.98], p = 2.1e-5\nC-index = 0.656 [0.60–0.71]\nBinary cutoff HR = 1.95 [1.16–3.29]\nIncremental ΔC = +0.028",
         title="Supportive GSE20685 Evaluation", bg="#EBF5FB", ec="#2980B9", fontsize=7.2)

# Validation Cohort 2: SCAN-B GSE96058
draw_box(35, 3, 30, 16,
         "3,273 unique patients (336 OS deaths)\nLocked cross-platform test (RNA-seq)\nContinuous HR = 1.44 [1.31–1.59], p = 8.3e-14\nC-index = 0.595 [0.56–0.63]\nIncremental ΔC = +0.024\nPAM50 interaction p = 0.031 (exploratory)",
         title="Later Locked Application: SCAN-B", bg="#EBF5FB", ec="#2980B9", fontsize=7.5)

# Validation Cohort 3: METABRIC
draw_box(68, 3, 30, 16,
         "1,980 patients (1,143 OS deaths)\nIndependent Illumina microarray test\nContinuous HR = 1.13 [1.07–1.20], p = 2.2e-5\nC-index = 0.573 [0.55–0.59]\nLocked cutoff HR = 1.00 (non-transport)\nPAM50 interaction p = 0.42 (not replicated)",
         title="Later Locked Application: METABRIC", bg="#EBF5FB", ec="#2980B9", fontsize=7.5)

plt.tight_layout()
fig.savefig(FIG / "Fig1_flow.png", dpi=300)
fig.savefig(FIG / "Fig1_flow.pdf")
plt.close(fig)
print("Fig1_flow successfully generated.")

# -------------------------------------------------------------
# 2. PAM50 Subtype Forest Plot (Fig_PAM50_forest)
# -------------------------------------------------------------
subtypes = [
    "Basal-like",
    "HER2-enriched",
    "Normal-like",
    "Luminal A",
    "Luminal B",
    "Overall SCAN-B (Total)"
]
n_patients = [339, 327, 221, 1657, 729, 3273]
n_events = [67, 55, 16, 116, 82, 336]
hrs = [1.40, 1.45, 1.45, 0.99, 0.95, 1.44]
ci_lo = [1.10, 1.11, 0.92, 0.83, 0.76, 1.31]
ci_hi = [1.77, 1.91, 2.30, 1.19, 1.18, 1.59]
p_vals = ["0.005", "0.007", "0.11", "0.95", "0.62", "8.3e-14"]
q_bh = ["0.014", "0.014", "0.15", "0.95", "0.83", "—"]

y_pos = np.arange(len(subtypes))[::-1]

# Keep the estimates in a distinct panel. Placing them at fixed x coordinates on
# the forest axis caused the labels to collide with the confidence intervals when
# the figure was scaled into the journal-width manuscript PDF.
fig = plt.figure(figsize=(9.5, 4.3), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=(1.25, 1.9), wspace=0.08,
                      left=0.23, right=0.985, top=0.88, bottom=0.22)
ax = fig.add_subplot(gs[0, 0])
tab_ax = fig.add_subplot(gs[0, 1], sharey=ax)

ax.axvline(1.0, color="#7F8C8D", linestyle="--", linewidth=1.2, zorder=1)
colors = ["#C0392B", "#C0392B", "#7F8C8D", "#2980B9", "#2980B9", "#1B4F72"]
markers = ["s", "s", "o", "s", "s", "D"]
for i, y in enumerate(y_pos):
    ax.plot([ci_lo[i], ci_hi[i]], [y, y], color=colors[i], lw=2.2, zorder=2)
    ax.plot(hrs[i], y, marker=markers[i], markersize=8 if i < 5 else 9,
            color=colors[i], zorder=3)

ax.set_yticks(y_pos, subtypes, fontsize=9, fontweight="medium")
ax.set_xlabel("Hazard ratio per score SD (95% CI)", fontsize=9.5, labelpad=7)
ax.set_xlim(0.65, 2.4)
ax.set_ylim(-0.65, len(subtypes) - 0.15)
ax.grid(axis="x", linestyle=":", alpha=0.6)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

tab_ax.set_xlim(0, 1)
tab_ax.set_ylim(ax.get_ylim())
tab_ax.axis("off")
columns = [(0.02, "n / Events"), (0.31, "HR (95% CI)"), (0.76, "p-value (FDR)")]
for x, label in columns:
    tab_ax.text(x, len(subtypes) - 0.15, label, ha="left", va="bottom",
                fontsize=8, fontweight="bold")
for i, y in enumerate(y_pos):
    fdr_str = f" ({q_bh[i]})" if q_bh[i] != "—" else ""
    tab_ax.text(columns[0][0], y, f"{n_patients[i]} / {n_events[i]}",
                ha="left", va="center", fontsize=8)
    tab_ax.text(columns[1][0], y,
                f"{hrs[i]:.2f} [{ci_lo[i]:.2f}, {ci_hi[i]:.2f}]",
                ha="left", va="center", fontsize=8)
    tab_ax.text(columns[2][0], y, f"{p_vals[i]}{fdr_str}",
                ha="left", va="center", fontsize=8)

fig.text(0.52, 0.055,
         "Global risk × PAM50 interaction: LR χ² = 10.64, df = 4, p = 0.031",
         ha="center", fontsize=8.5, fontstyle="italic", color="#2C3E50")
fig.savefig(FIG / "Fig_PAM50_forest.png", dpi=300)
fig.savefig(FIG / "Fig_PAM50_forest.pdf")
plt.close(fig)
print("Fig_PAM50_forest successfully generated.")
