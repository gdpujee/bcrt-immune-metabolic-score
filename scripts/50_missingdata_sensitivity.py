"""Complete-case sensitivity for the GSE20685 and SCAN-B adjusted models (MDSENS-001).

The missing-data policy paragraph (Methods, supplement, REMARK) states that the
GSE20685 and SCAN-B primary adjusted models are median-imputed full-cohort fits
with complete-case analyses "yielding consistent estimates".  That claim must be
readable from a raw output file like every other number, so this script fits the
the reported covariate sets on complete cases only and records the result. The
available materials do not establish when the covariate sets were selected.

Outputs: results/raw/missingdata_sensitivity.json, logs/missingdata_sensitivity.log
"""
import json
import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results/raw"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
logf = open(LOGS / "missingdata_sensitivity.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")


def fit_cc(df, covs, label, sd):
    cols = ["risk_unit"] + covs
    d = df.dropna(subset=cols).copy()
    X = np.column_stack([d[c].values.astype(float) for c in cols])
    r = PHReg(d.time.values, X, d.event.values.astype(int)).fit(disp=0)
    b, se = r.params[0], r.bse[0]
    k = float(sd)
    out = {"n": int(len(d)), "deaths": int(d.event.sum()),
           "excluded_missing": int(len(df) - len(d)),
           "per_SD_HR": round(float(np.exp(b * k)), 4),
           "per_SD_CI": [round(float(np.exp((b - 1.96 * se) * k)), 4),
                         round(float(np.exp((b + 1.96 * se) * k)), 4)],
           "per_SD_p": float(r.pvalues[0])}
    log(f"{label} complete case: n={out['n']} deaths={out['deaths']} "
        f"(excluded {out['excluded_missing']}) per-SD HR={out['per_SD_HR']:.3f} "
        f"[{out['per_SD_CI'][0]:.3f},{out['per_SD_CI'][1]:.3f}] p={out['per_SD_p']:.3g}")
    return out


out = {"method": "statsmodels PHReg, per-SD locked risk, reported adjusted "
        "covariate sets fitted on complete cases only (sensitivity; selection timing "
        "not established; the primary "
        "models are median-imputed full-cohort fits)"}

# ---- GSE20685: age + T + N ----
g = pd.read_csv(META / "GSE20685_clinical_curated.tsv", sep="\t")
gr = pd.read_csv(RAW / "validation_risk_GSE20685.tsv", sep="\t")
g = g.merge(gr[["GSM", "risk"]], on="GSM")
g = g.dropna(subset=["follow_up_duration_years", "event_death", "risk"])
g = g.rename(columns={"follow_up_duration_years": "time", "event_death": "event",
                      "risk": "risk_unit"})
for c in ("age_at_diagnosis", "t_stage", "n_stage"):
    g[c] = pd.to_numeric(g[c], errors="coerce")
out["GSE20685"] = fit_cc(g, ["age_at_diagnosis", "t_stage", "n_stage"],
                         "GSE20685 (age+T+N)", sd=g.risk_unit.std())

# ---- SCAN-B: age + ER + HER2 ----
rr = pd.read_csv(RAW / "rnaseq_risk_GSE96058.tsv", sep="\t")
cl = pd.read_csv(META / "GSE96058_clinical_raw.tsv", sep="\t")
merged = cl.merge(rr[["sample", "risk", "OS_years", "OS_event"]],
                  left_on="title", right_on="sample", how="inner")
sc = pd.DataFrame({
    "time": pd.to_numeric(merged.OS_years, errors="coerce"),
    "event": pd.to_numeric(merged.OS_event, errors="coerce"),
    "risk_unit": pd.to_numeric(merged.risk, errors="coerce"),
    "age_at_diagnosis": pd.to_numeric(merged.age_at_diagnosis, errors="coerce"),
    "er_status": pd.to_numeric(merged.er_status, errors="coerce"),
    "her2_status": pd.to_numeric(merged.her2_status, errors="coerce")})
sc = sc.dropna(subset=["time", "event", "risk_unit"])
out["SCANB"] = fit_cc(sc, ["age_at_diagnosis", "er_status", "her2_status"],
                      "SCAN-B (age+ER+HER2)", sd=sc.risk_unit.std())

json.dump(out, open(RAW / "missingdata_sensitivity.json", "w"), indent=2)
log("WROTE results/raw/missingdata_sensitivity.json (MDSENS-001)")
