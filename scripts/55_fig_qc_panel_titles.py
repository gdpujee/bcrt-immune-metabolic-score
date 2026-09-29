"""Rebuild Fig. S1 with explicit per-panel cohort titles (FIGS1-001).

External review round 4: the follow-up/expression QC figure (drawn by
scripts/03_qc_endpoints.py before the cross-platform extension) has two
untitled panels, and its supplement caption said "across cohorts" while it
shows only the two GPL570 cohorts. Same content, same inputs, labelled — the
caption is corrected in the supplement to name GSE42568 and GSE20685.

Outputs: figures/QC_followup_expr.png/.pdf (same filenames the renderers map)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata"
PROC = ROOT / "data/processed"
FIG = ROOT / "figures"
LOGS = ROOT / "logs"
logf = open(LOGS / "fig_qc_panel_titles.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")

# same cohorts as scripts/03_qc_endpoints.py
c425 = pd.read_csv(META / "GSE42568_clinical_curated.tsv", sep="\t")
train = c425[(c425["tissue"] == "breast cancer") & c425["OS_time_days"].notna()
             & c425["OS_event"].notna()].copy()
c206 = pd.read_csv(META / "GSE20685_clinical_curated.tsv", sep="\t")
valid = c206[c206["OS_years"].notna() & c206["OS_event"].notna()].copy()
log(f"train n={len(train)} (expect 104); valid n={len(valid)} (expect 327)")
assert (len(train), len(valid)) == (104, 327)

e425 = pd.read_csv(PROC / "GSE42568_expr_gene.tsv", sep="\t", index_col=0)
e206 = pd.read_csv(PROC / "GSE20685_expr_gene.tsv", sep="\t", index_col=0)
meds425 = np.nanmedian(e425[train.GSM.tolist()].values, axis=0)
meds206 = np.nanmedian(e206[valid.GSM.tolist()].values, axis=0)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].hist(train["OS_years"], bins=20, alpha=0.7, label="GSE42568 (derivation, n=104)")
ax[0].hist(valid["OS_years"], bins=20, alpha=0.5, label="GSE20685 (supportive evaluation, n=327)")
ax[0].set_xlabel("OS years"); ax[0].set_ylabel("n"); ax[0].legend(fontsize=8)
ax[0].set_title("Overall-survival follow-up (years)", fontsize=9)
ax[1].hist(meds425, bins=20, alpha=0.7, label="GSE42568")
ax[1].hist(meds206, bins=20, alpha=0.5, label="GSE20685")
ax[1].set_xlabel("per-sample median expression"); ax[1].legend(fontsize=8)
ax[1].set_title("Expression scale, GPL570 cohorts", fontsize=9)
plt.tight_layout()
plt.savefig(FIG / "QC_followup_expr.png", dpi=150)
plt.savefig(FIG / "QC_followup_expr.pdf")
log("WROTE figures/QC_followup_expr.png/.pdf (FIGS1-001)")
