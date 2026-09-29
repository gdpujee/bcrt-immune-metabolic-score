"""Materialise the documented METABRIC complete-case adjusted frame (MBADJ-FRAME-001).

External review round 4 item 10 asks for the adjusted-model cox.zph diagnostics
that GSE20685 and SCAN-B already report. The exact test runs in R (scripts/53),
so this script writes the analysis frame with the SAME construction as
scripts/33_metabric_adjusted.py (cBioPortal covariate caches, mm-derived T
category per D-039, per-SD risk) and asserts it reproduces the published
complete-case n/deaths before anything is fitted.

Output: metadata/METABRIC_adjusted_cc_frame.tsv, logs/metabric_adj_frame.log
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
logf = open(LOGS / "metabric_adj_frame.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")


risk = pd.read_csv(RES / "metabric_risk.tsv", sep="\t")
summ = json.load(open(RES / "metabric_summary.json"))
clin = pd.read_csv(META / "METABRIC_clinical_dl.tsv", sep="\t").drop_duplicates("patientId")
pc = pd.read_csv(META / "METABRIC_clinical_adj.tsv", sep="\t")
sc = pd.read_csv(META / "METABRIC_sample_clinical.tsv", sep="\t")
m = (risk.merge(clin, on="patientId", how="inner")
         .merge(pc, on="patientId", how="left", suffixes=("_dl", "_cbio"))
         .merge(sc, on="patientId", how="left"))
if {"AGE_AT_DIAGNOSIS_dl", "AGE_AT_DIAGNOSIS_cbio"} <= set(m.columns):
    a = pd.to_numeric(m.AGE_AT_DIAGNOSIS_dl, errors="coerce")
    b = pd.to_numeric(m.AGE_AT_DIAGNOSIS_cbio, errors="coerce")
    both = a.notna() & b.notna()
    dmax = float((a[both] - b[both]).abs().max()) if int(both.sum()) else 0.0
    log(f"AGE_AT_DIAGNOSIS cross-source check: {int(both.sum())} rows, max|diff|={dmax}")
    if dmax > 1e-9:
        raise SystemExit("age disagrees between the two cBioPortal downloads")
    m["AGE_AT_DIAGNOSIS"] = b.where(b.notna(), a)
m["OS_years"] = pd.to_numeric(m.OS_MONTHS, errors="coerce") / 12
m["OS_event"] = m.OS_STATUS.str.contains("DECEASED", na=False).astype(int)
m = m.dropna(subset=["OS_years", "OS_event", "risk"]).copy()

m["age"] = pd.to_numeric(m.AGE_AT_DIAGNOSIS, errors="coerce")
m["nodal_burden"] = pd.to_numeric(m.LYMPH_NODES_EXAMINED_POSITIVE, errors="coerce")
m["grade"] = pd.to_numeric(m.GRADE, errors="coerce")
m["er"] = m.ER_STATUS.map({"Positive": 1, "Negative": 0})
m["her2"] = m.HER2_STATUS.map({"Positive": 1, "Negative": 0})
_size_mm = pd.to_numeric(m.TUMOR_SIZE, errors="coerce")
m["tumor_size"] = pd.Series(
    np.select([_size_mm <= 20, (_size_mm > 20) & (_size_mm <= 50), _size_mm > 50],
              [1.0, 2.0, 3.0], default=np.nan), index=m.index)
m.loc[_size_mm.isna(), "tumor_size"] = np.nan
m["risk_sd"] = m.risk / m.risk.std()

COVS = ["age", "nodal_burden", "grade", "er", "her2", "tumor_size"]
cc = m.dropna(subset=COVS + ["risk_sd"]).copy()
log(f"complete case: n={len(cc)} deaths={int(cc.OS_event.sum())}")
adj = json.load(open(RES / "metabric_adjusted.json"))["complete_case"]
if (len(cc), int(cc.OS_event.sum())) != (adj["n"], adj["deaths"]):
    raise SystemExit(f"frame {len(cc)}/{int(cc.OS_event.sum())} does not match the "
                     f"published complete-case model {adj['n']}/{adj['deaths']}")
cols = ["OS_years", "OS_event", "risk_sd"] + COVS
cc[cols].to_csv(META / "METABRIC_adjusted_cc_frame.tsv", sep="\t", index=False)
log(f"WROTE metadata/METABRIC_adjusted_cc_frame.tsv ({len(cc)} rows) (MBADJ-FRAME-001)")
