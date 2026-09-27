#!/usr/bin/env python3
"""03_qc_endpoints.py — RUN-ID: QC-001
QC + endpoint harmonization + inclusion/exclusion freeze.
Inputs: metadata/*_clinical_raw.tsv, data/processed/*_expr_gene.tsv
Outputs: metadata/*_clinical_curated.tsv, results/raw/QC_summary.json, results/raw/QC_report.md, figures/QC_*.png
"""
import pandas as pd, numpy as np, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
META = ROOT/"metadata"; PROC = ROOT/"data/processed"; RES = ROOT/"results/raw"; FIG = ROOT/"figures"; LOGS = ROOT/"logs"
FIG.mkdir(exist_ok=True, parents=True); RES.mkdir(exist_ok=True, parents=True)

def load(gse):
    clin = pd.read_csv(META/f"{gse}_clinical_raw.tsv", sep="\t")
    expr = pd.read_csv(PROC/f"{gse}_expr_gene.tsv", sep="\t", index_col=0)
    return clin, expr

out = {}
report = ["# QC Report (QC-001, 2026-09-20)", ""]

# --- GSE42568 ---
c425, e425 = load("GSE42568")
# inclusion: breast cancer + OS non-missing
c425["OS_time_days"] = pd.to_numeric(c425["overall_survival_time_days"], errors="coerce")
c425["OS_event"] = pd.to_numeric(c425["overall_survival_event"], errors="coerce")
c425["RFS_time_days"] = pd.to_numeric(c425["relapse_free_survival_time_days"], errors="coerce")
c425["RFS_event"] = pd.to_numeric(c425["relapse_free_survival_event"], errors="coerce")
c425["OS_years"] = c425["OS_time_days"]/365.25
train = c425[(c425["tissue"]=="breast cancer") & c425["OS_time_days"].notna() & c425["OS_event"].notna()].copy()
excluded_425 = len(c425) - len(train) - (c425["tissue"]=="normal breast").sum()
normals_425 = (c425["tissue"]=="normal breast").sum()
out["GSE42568"] = {"n_total": int(len(c425)), "n_train_tumors": int(len(train)),
  "n_normals": int(normals_425), "n_excluded_other": int(excluded_425),
  "n_OS_events": int((train["OS_event"]==1).sum()), "n_OS_censored": int((train["OS_event"]==0).sum()),
  "median_OS_years": float(train["OS_years"].median()), "median_followup_censored_years": float(train.loc[train.OS_event==0,"OS_years"].median()),
  "age_median": float(pd.to_numeric(c425.loc[c425.tissue=='breast cancer',"age"], errors="coerce").median())}
report += [f"## GSE42568 (training): n_total={len(c425)}, tumors w/ OS={len(train)} (events={(train.OS_event==1).sum()}, censored={(train.OS_event==0).sum()}), normals={normals_425}",
 f"- OS years median {out['GSE42568']['median_OS_years']:.2f}; censored median {out['GSE42568']['median_followup_censored_years']:.2f}; max {train.OS_years.max():.2f}",
 f"- Expression: {e425.shape[0]} genes x {e425.shape[1]} samples; range [{np.nanmin(e425.values):.2f},{np.nanmax(e425.values):.2f}]"]
c425.to_csv(META/"GSE42568_clinical_curated.tsv", sep="\t", index=False)

# --- GSE20685 ---
c206, e206 = load("GSE20685")
c206["OS_years"] = pd.to_numeric(c206["follow_up_duration_years"], errors="coerce")
c206["OS_event"] = pd.to_numeric(c206["event_death"], errors="coerce")
c206["MET_event"] = pd.to_numeric(c206["event_metastasis"], errors="coerce")
valid = c206[c206["OS_years"].notna() & c206["OS_event"].notna()].copy()
out["GSE20685"] = {"n_total": int(len(c206)), "n_valid_OS": int(len(valid)),
  "n_OS_events": int((valid["OS_event"]==1).sum()), "n_OS_censored": int((valid["OS_event"]==0).sum()),
  "median_OS_years": float(valid["OS_years"].median()), "median_followup_censored": float(valid.loc[valid.OS_event==0,"OS_years"].median()),
  "n_metastasis_events": int((valid["MET_event"]==1).sum())}
report += [f"## GSE20685 (validation): n_total={len(c206)}, valid OS={len(valid)} (deaths={(valid.OS_event==1).sum()}, censored={(valid.OS_event==0).sum()})",
 f"- OS years median {out['GSE20685']['median_OS_years']:.2f}; censored median {out['GSE20685']['median_followup_censored']:.2f}",
 f"- Metastasis events {(valid.MET_event==1).sum()}; Expression {e206.shape[0]} genes x {e206.shape[1]} samples"]
c206.to_csv(META/"GSE20685_clinical_curated.tsv", sep="\t", index=False)

# --- GSE45827 ---
c458, e458 = load("GSE45827")
# exclude cell lines: diagnosis != Breast cancer and tumor_subtype NA? Use cell_line non-empty as cell-line flag
c458["is_cell_line"] = c458["cell_line"].fillna("").astype(str).str.len() > 0
c458_bio = c458[~c458["is_cell_line"]].copy()
n_tumor = (c458_bio["diagnosis"]=="Breast cancer").sum()
n_normal = ((c458_bio["diagnosis"]=="None (normal)") | (c458_bio["diagnosis"]=="None")).sum()
out["GSE45827"] = {"n_total": int(len(c458)), "n_excluded_cell_lines": int(c458["is_cell_line"].sum()),
  "n_bio": int(len(c458_bio)), "n_tumor": int(n_tumor), "n_normal": int(n_normal),
  "batches": sorted(c458_bio["batch"].dropna().unique().tolist())}
report += [f"## GSE45827 (biology-only): n_total={len(c458)}, excluded cell lines={out['GSE45827']['n_excluded_cell_lines']}, bio={len(c458_bio)} (tumor={n_tumor}, normal={n_normal})",
 f"- Batches {out['GSE45827']['batches']}; NO survival (by design); Expression {e458.shape[0]} genes x {e458.shape[1]} samples (pre-filtered, 14.5k genes)"]
c458.to_csv(META/"GSE45827_clinical_curated.tsv", sep="\t", index=False)

# harmonization note
report += ["", "## Endpoint harmonization",
 "- Primary OS: 42568 death/censor + days/365.25; 20685 event_death + follow_up_duration_years. Both all-cause OS from diagnosis/sample; units years for KM/Cox.",
 "- Secondary (separate, never pooled): 42568 RFS, 20685 metastasis. Definitions differ -> reported separately.",
 "- Inclusion frozen: 42568 train = 104 tumors w/ OS; 20685 validation = 327 w/ OS; 45827 biology = 141 (130+11) excl. 14 cell lines."]

with open(RES/"QC_summary.json","w") as f: json.dump(out, f, indent=2)
with open(RES/"QC_report.md","w") as f: f.write("\n".join(report)+"\n")

# figures: OS distribution + expression medians
fig, ax = plt.subplots(1,2, figsize=(10,4))
ax[0].hist(train["OS_years"], bins=20, alpha=0.7, label="42568 train")
ax[0].hist(valid["OS_years"], bins=20, alpha=0.5, label="20685 valid")
ax[0].set_xlabel("OS years"); ax[0].set_ylabel("n"); ax[0].legend()
meds425 = np.nanmedian(e425.values, axis=0); meds206 = np.nanmedian(e206.values, axis=0)
ax[1].hist(meds425, bins=20, alpha=0.7, label="42568")
ax[1].hist(meds206, bins=20, alpha=0.5, label="20685")
ax[1].set_xlabel("per-sample median expression"); ax[1].legend()
plt.tight_layout(); plt.savefig(FIG/"QC_followup_expr.png", dpi=150); plt.savefig(FIG/"QC_followup_expr.pdf")
print(json.dumps(out, indent=2))
print("WROTE QC_summary.json, QC_report.md, QC_followup_expr.png")
with open(LOGS/"qc.log", "w") as _lf:
    _lf.write(json.dumps(out, indent=2) + "\nWROTE QC_summary.json, QC_report.md, QC_followup_expr.png\n")
