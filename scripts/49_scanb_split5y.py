"""SCAN-B 5-year split-time sensitivity (SCANB-SPLIT-001).

The manuscript's SCAN-B time-stratified sentence (<=5y HR vs >5y HR) must be
readable from a raw output file like the GSE20685 and METABRIC splits are, so
this script computes both windows from results/raw/rnaseq_risk_GSE96058.tsv
with the same estimator as scripts/28_ph_diagnostics.py (statsmodels PHReg,
per-SD risk):
  early: administrative censoring at 5y on the FULL cohort
         (time=min(T,5), event=E&(T<=5)) — subsetting to T<=5 would drop
         at-risk person-time and bias the estimate;
  late:  delayed entry (conditional on survival past 5y, time=T-5).
Outputs: results/raw/scanb_split5y.json, logs/scanb_split5y.log
"""
import json
import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
LOGS = ROOT / "logs"

logf = open(LOGS / "scanb_split5y.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")


r = pd.read_csv(RES / "rnaseq_risk_GSE96058.tsv", sep="\t")
T = r.OS_years.values.astype(float)
E = r.OS_event.values.astype(int)
risk = (r.risk.values / r.risk.std()).reshape(-1, 1)
log(f"cohort n={len(T)} deaths={int(E.sum())}")

out = {"method": "statsmodels PHReg, per-SD locked risk; early = administrative "
        "censoring at 5y on the full cohort, late = delayed entry at 5y",
       "unit": "per SD of locked linear predictor, within-cohort SD"}

eT = np.minimum(T, 5)
eE = ((T <= 5) & (E == 1)).astype(int)
rr = PHReg(eT, risk, eE).fit(disp=0)
out["early"] = {"n": int(len(eT)), "deaths": int(eE.sum()),
                "HR": round(float(np.exp(rr.params[0])), 3),
                "CI": [round(float(np.exp(rr.params[0] - 1.96 * rr.bse[0])), 3),
                       round(float(np.exp(rr.params[0] + 1.96 * rr.bse[0])), 3)],
                "p": float(rr.pvalues[0])}

late = T > 5
rr2 = PHReg(T[late] - 5, risk[late], E[late]).fit(disp=0)
out["late"] = {"n": int(late.sum()), "deaths": int(E[late].sum()),
               "HR": round(float(np.exp(rr2.params[0])), 3),
               "CI": [round(float(np.exp(rr2.params[0] - 1.96 * rr2.bse[0])), 3),
                      round(float(np.exp(rr2.params[0] + 1.96 * rr2.bse[0])), 3)],
               "p": float(rr2.pvalues[0])}

log(f"early: n={out['early']['n']} deaths={out['early']['deaths']} "
    f"HR={out['early']['HR']} CI={out['early']['CI']}")
log(f"late:  n={out['late']['n']} deaths={out['late']['deaths']} "
    f"HR={out['late']['HR']} CI={out['late']['CI']}")
json.dump(out, open(RES / "scanb_split5y.json", "w"), indent=2)
log("WROTE scanb_split5y.json (SCANB-SPLIT-001)")
