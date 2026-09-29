"""Rebuild Fig. S2 as the four-cohort locked-score distribution (FIGS2-001).

External review round 4: the published risk_distributions figure predates the
cross-platform extension — two untitled panels (GSE42568, GSE20685) while the
paper validates in four cohorts, and its caption said "discovery and validation
cohorts". The figure is also the direct visual evidence for the headline
transportability result: the frozen cutoff (-0.2880, derivation median) sits in
a different place in every cohort. All four scores are the locked linear
predictor with frozen coefficients and derivation means/SDs, read from the
published per-sample files (METABRIC/SCAN-B/GSE20685) or recomputed from the
z-scored training matrix with the locked model (GSE42568).

Outputs: figures/risk_distributions.png/.pdf (same filenames the renderers map)
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
DER = ROOT / "results/derived"
FIG = ROOT / "figures"
LOGS = ROOT / "logs"
logf = open(LOGS / "fig_risk_distributions4.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")


locked = json.load(open(DER / "locked_model.json"))
cutoff = locked["cutoff"]

# GSE42568 (derivation): recompute the locked score from the z matrix
zh = pd.read_csv(RES / "heatmap_train_14g_z.tsv", sep="\t", index_col=0)
genes = locked["genes"]
assert set(genes) <= set(zh.columns), set(genes) - set(zh.columns)
coefs = np.asarray(locked["coefs"], float)  # ordered as locked["genes"]
train_risk = zh[genes].values @ coefs

vr = pd.read_csv(RES / "validation_risk_GSE20685.tsv", sep="\t").risk.values
sc = pd.read_csv(RES / "rnaseq_risk_GSE96058.tsv", sep="\t").risk.values
mt = pd.read_csv(RES / "metabric_risk.tsv", sep="\t").risk.values

panels = [("GSE42568 (derivation)", train_risk), ("GSE20685 (supportive evaluation)", vr),
          ("SCAN-B (RNA-seq)", sc), ("METABRIC (Illumina)", mt)]
fig, ax = plt.subplots(1, 4, figsize=(15, 3.4))
for a, (nm, v) in zip(ax, panels):
    a.hist(v, bins=25, alpha=0.75, color="#1f77b4")
    a.axvline(cutoff, color="r", ls="--", lw=1.5)
    a.set_title(nm, fontsize=9)
    a.set_xlabel("Locked score")
    above = int((v > cutoff).sum())
    a.text(0.03, 0.93, f"n={len(v)}\nabove cutoff: {above}", transform=a.transAxes,
           va="top", fontsize=7)
    log(f"{nm}: n={len(v)} above locked cutoff {cutoff:.4f}: {above}")
ax[1].set_ylabel("samples")
plt.tight_layout()
plt.savefig(FIG / "risk_distributions.png", dpi=150)
plt.savefig(FIG / "risk_distributions.pdf")
log("WROTE figures/risk_distributions.png/.pdf (FIGS2-001)")
