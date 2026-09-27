"""Median-probe sensitivity for max-mean mapping caveat (MEDPROBE-001, exploratory).

Rebuilds microarray gene matrices with median-ranked probe per gene (probes
sorted by mean expression, middle rank; deterministic, no outcome use) instead
of max-mean, then applies the LOCKED genes/coefficients/scaling/cutoff
(no refitting) and recomputes GSE20685 validation metrics. Compares against
main locked results (continuous HR 1.59, binary HR 1.95, C 0.656).
SCAN-B (RNA-seq) is unaffected by probe mapping.
"""
import gzip
import json
import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
from pathlib import Path

ROOT = Path(".")
RAW = ROOT / "data/raw"
RES = ROOT / "results/raw"
DER = ROOT / "results/derived"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
SEED = 42
np.random.seed(SEED)


def log(m):
    print(m, flush=True)
    with open(LOGS / "medprobe.log", "a") as f:
        f.write(m + "\n")


def harrell_c(t, e, r):
    cc = tt = 0
    n = len(t)
    for i in range(n):
        if e[i] != 1:
            continue
        for j in range(n):
            if t[j] > t[i]:
                tt += 1
                cc += (r[j] < r[i]) + 0.5 * (r[j] == r[i])
    return cc / tt if tt else float("nan")


locked = json.load(open(DER / "locked_model.json"))
genes = locked["genes"]
coefs = np.array(locked["coefs"])
means = pd.Series(locked["scaling"]["means"])
sds = pd.Series(locked["scaling"]["sds"])
cutoff = locked["cutoff"]

# GPL570 annot probe->symbol (mirror of 02_parse_geo.py, incl. /// split)
probe2sym = {}
with gzip.open(RAW / "GPL570.annot.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("!") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3 and parts[2].strip():
            sym = parts[2].strip().strip('"').split("///")[0].strip()
            if sym:
                probe2sym[parts[0].strip().strip('"')] = sym

out = {}
for gse in ["GSE42568", "GSE20685"]:
    expr = pd.read_csv(RAW / f"{gse}_series_matrix.txt.gz", sep="\t",
                       compression="gzip", comment="!", header=0, index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr_mapped = expr[expr.index.isin(probe2sym)]
    probe_means = expr_mapped.mean(axis=1, skipna=True)
    tmp = pd.DataFrame({"probe": expr_mapped.index,
                        "symbol": [probe2sym[p] for p in expr_mapped.index],
                        "mean": probe_means.values})
    picks = {}
    for sym, grp in tmp.groupby("symbol"):
        grp = grp.sort_values("mean").reset_index(drop=True)
        picks[sym] = grp.loc[len(grp) // 2, "probe"]
    syms = sorted(picks)
    gene_expr = expr_mapped.loc[[picks[s] for s in syms]].copy()
    gene_expr.index = syms
    gene_expr.index.name = "symbol"
    log(f"{gse}: median-probe genes={len(picks)}")
    out[gse] = gene_expr

# GSE20685 locked evaluation on median-mapped data
clin = pd.read_csv(META / "GSE20685_clinical_curated.tsv", sep="\t")
X = out["GSE20685"].loc[genes, clin.GSM].T
Z = (X - means) / sds
risk = (Z.values @ coefs)
T = clin.follow_up_duration_years.values.astype(float)
E = clin.event_death.values.astype(int)
sd = risk.std()
rc = PHReg(T, (risk / sd).reshape(-1, 1), E).fit(disp=0)
hr = float(np.exp(rc.params[0]))
ci = [float(np.exp(rc.params[0] - 1.96 * rc.bse[0])),
      float(np.exp(rc.params[0] + 1.96 * rc.bse[0]))]
hi = (risk > cutoff).astype(int)
rb = PHReg(T, hi, E).fit(disp=0)
hrb = float(np.exp(rb.params[0]))
c = harrell_c(T, E, risk)
log(f"MEDIAN-PROBE GSE20685: cont HR={hr:.3f} [{ci[0]:.3f},{ci[1]:.3f}] "
    f"p={float(rc.pvalues[0]):.3g}; binary HR={hrb:.3f}; C={c:.3f}")
log(f"MAIN-LOCKED reference: cont HR=1.59 [1.29,1.98]; binary HR=1.95; C=0.656")
pd.DataFrame({"GSM": clin.GSM.values, "risk_medprobe": risk}).to_csv(
    RES / "medprobe_risk_GSE20685.tsv", sep="\t", index=False)
json.dump({"cont_HR": hr, "cont_CI": ci, "cont_p": float(rc.pvalues[0]),
           "bin_HR": hrb, "C": c, "seed": SEED},
          open(RES / "medprobe_summary.json", "w"), indent=2)
log("WROTE medprobe outputs (MEDPROBE-001)")
